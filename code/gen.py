"""
Générateur d'instances IRP calibrées par faits stylisés (Guinée).

Familles :
  BG  : Basse-Guinée, dépôt à Conakry (Kaloum)
  GF  : Guinée forestière, dépôt régional à N'Zérékoré
  UBG, UGF : familles témoins (mêmes paramètres, stations uniformes
             dans le rectangle englobant de la région)

Sources :
  - Populations : INS Guinée, RGPH-3 (2014), résultats définitifs.
  - Coordonnées des chefs-lieux : base GeoNames (paquet Python geonamescache).
"""
import hashlib, json, math, os, sys
import numpy as np

# --------------------------------------------------------------------------
# Localités : (nom, population RGPH 2014, latitude, longitude)  -- GeoNames
# Conakry : les 5 communes (Kaloum, Dixinn, Matam, Matoto, Ratoma) sont
# agrégées en une seule localité, centrée sur Dixinn, de rayon plus grand.
# --------------------------------------------------------------------------
REGIONS = {
    "BG": {
        "depot": ("Conakry (Kaloum)", 9.53795, -13.67729),
        "localites": [
            ("Conakry",     1660973,  9.55111, -13.67306),
            ("Coyah",        263861,  9.70643, -13.38465),
            ("Dubréka",      330548,  9.79111, -13.52333),
            ("Forécariah",   242942,  9.43056, -13.08806),
            ("Kindia",       439614, 10.05692, -12.86576),
            ("Télimélé",     284409, 10.90000, -13.03333),
            ("Boffa",        212583, 10.16667, -14.03333),
            ("Boké",         450278, 10.93217, -14.29055),
            ("Fria",          96700, 10.36676, -13.58253),
            ("Gaoual",       193612, 11.75000, -13.20000),
            ("Koundara",     129974, 12.48333, -13.30000),
        ],
    },
    "GF": {
        "depot": ("N'Zérékoré (dépôt régional)", 7.75624, -8.81790),
        "localites": [
            ("Beyla",        326082,  8.69011,  -8.64869),
            ("Guéckédou",    290611,  8.56744, -10.13360),
            ("Lola",         171561,  7.80000,  -8.53333),
            ("Macenta",      278456,  8.54351,  -9.47099),
            ("N'Zérékoré",   396949,  7.75624,  -8.81790),
            ("Yomou",        114371,  7.56748,  -9.26283),
        ],
    },
}
RAYON_CONAKRY_KM = 10.0   # agglomération étendue de la capitale

# Paramètres de base (scénario de référence)
BASE = dict(rho=5.0, phi=1.3, kappa=1.5, dbar=100.0)
ETA = {"L": 0.01, "H": 0.10}   # coût de détention relatif faible / élevé


def haversine(a, b):
    """Distance orthodromique (km) entre deux points (lat, lon)."""
    R = 6371.0
    la1, lo1, la2, lo2 = map(math.radians, (a[0], a[1], b[0], b[1]))
    h = math.sin((la2 - la1) / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin((lo2 - lo1) / 2) ** 2
    return 2 * R * math.asin(math.sqrt(h))


def offset(rng, lat, lon, radius_km):
    """Point uniforme dans un disque de rayon radius_km autour de (lat, lon)."""
    r = radius_km * math.sqrt(rng.random())
    th = 2 * math.pi * rng.random()
    dlat = (r * math.cos(th)) / 111.32
    dlon = (r * math.sin(th)) / (111.32 * math.cos(math.radians(lat)))
    return lat + dlat, lon + dlon


def plus_fort_reste(pops, N):
    """Répartit N stations proportionnellement aux populations."""
    tot = sum(pops)
    q = [N * p / tot for p in pops]
    n = [int(math.floor(v)) for v in q]
    reste = N - sum(n)
    ordre = sorted(range(len(pops)), key=lambda k: q[k] - n[k], reverse=True)
    for k in ordre[:reste]:
        n[k] += 1
    return n


def generer(famille, N, H, eta_code, rep, rho=None, phi=None, kappa=None, g_choices=(2, 3)):
    """Retourne (instance stylisée, instance témoin) sous forme de dictionnaires."""
    rho = BASE["rho"] if rho is None else rho
    phi = BASE["phi"] if phi is None else phi
    kappa = BASE["kappa"] if kappa is None else kappa
    reg = REGIONS[famille]
    # La graine ne dépend que de l'instance de base : les variantes de sensibilité
    # (kappa, phi, g) partagent ainsi exactement les mêmes stations et demandes.
    cle = f"{famille}-{N}-{H}-{eta_code}-{rep}-{BASE['rho']}-{BASE['phi']}-{BASE['kappa']}-{(2, 3)}"
    seed = int(hashlib.sha256(cle.encode()).hexdigest()[:8], 16)   # graine distincte par instance
    rng = np.random.default_rng(seed)

    locs = reg["localites"]
    pops = [l[1] for l in locs]
    n_par_loc = plus_fort_reste(pops, N)

    # 1. Positions et localité de rattachement
    pos, loc_of = [], []
    for k, (nom, P, lat, lon) in enumerate(locs):
        rad = RAYON_CONAKRY_KM if nom == "Conakry" else rho
        for _ in range(n_par_loc[k]):
            pos.append(offset(rng, lat, lon, rad))
            loc_of.append(k)

    # 3. Demande proportionnelle à la population par station
    xi = rng.uniform(0.8, 1.2, size=N)
    raw = np.array([pops[loc_of[i]] / n_par_loc[loc_of[i]] for i in range(N)]) * xi
    d = np.maximum(1, np.rint(BASE["dbar"] * raw / raw.mean())).astype(int)

    # 4. Stocks
    g = rng.choice(list(g_choices), size=N)
    C = (g * d).astype(int)
    I0 = (C - d).astype(int)

    # 5. Camion
    Q = int(max(round(kappa * d.sum()), C.max()))

    # Rectangle englobant (pour la famille témoin)
    lats = [l[2] for l in locs] + [reg["depot"][1]]
    lons = [l[3] for l in locs] + [reg["depot"][2]]
    bbox = (min(lats), max(lats), min(lons), max(lons))
    pos_u = [(rng.uniform(bbox[0], bbox[1]), rng.uniform(bbox[2], bbox[3])) for _ in range(N)]

    def assemble(positions, nom):
        pts = [(reg["depot"][1], reg["depot"][2])] + list(positions)
        n = len(pts)
        c = [[0 if i == j else int(round(phi * haversine(pts[i], pts[j]))) for j in range(n)] for i in range(n)]
        # garantir des coûts >= 1 entre sites distincts (inégalité triangulaire préservée à l'arrondi près)
        for i in range(n):
            for j in range(n):
                if i != j and c[i][j] == 0:
                    c[i][j] = 1
        cbar = float(np.mean([c[0][i] for i in range(1, n)]))
        h = ETA[eta_code] * cbar / BASE["dbar"]
        return dict(
            name=nom, family=nom.split("_")[0], N=N, H=H, eta=eta_code, rep=rep,
            rho=rho, phi=phi, kappa=kappa, g_choices=list(g_choices),
            coords=pts, c=c, d=[int(v) for v in d], C=[int(v) for v in C], I0=[int(v) for v in I0],
            Q=Q, r=int(d.sum()), B0=int(C.sum()), h0=0.0, h=[h] * N,
            loc=[locs[k][0] for k in loc_of],
        )

    tag = f"N{N}_H{H}_{eta_code}_r{rep}"
    extra = ""
    if (rho, phi, kappa, tuple(g_choices)) != (BASE["rho"], BASE["phi"], BASE["kappa"], (2, 3)):
        extra = f"_rho{rho}_phi{phi}_k{kappa}_g{'-'.join(map(str, g_choices))}"
    sty = assemble(pos, f"{famille}_{tag}{extra}")
    tem = assemble(pos_u, f"U{famille}_{tag}{extra}")
    return sty, tem


def main(outdir, reps=3):
    os.makedirs(outdir, exist_ok=True)
    count = 0
    for fam in ("BG", "GF"):
        for N in (5, 10, 15, 20, 25, 30):
            for H in (3, 6):
                for eta in ("L", "H"):
                    for rep in range(1, reps + 1):
                        for inst in generer(fam, N, H, eta, rep):
                            with open(os.path.join(outdir, inst["name"] + ".json"), "w") as f:
                                json.dump(inst, f)
                            count += 1
    print(count, "instances écrites dans", outdir)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "instances")
