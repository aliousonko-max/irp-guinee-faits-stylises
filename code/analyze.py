"""Analyse des résultats -> tableaux LaTeX (dossier ../tables) et figures (../figures)."""
import json, os, math, statistics as st
from collections import defaultdict
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB, FIG = os.path.join(ROOT, "tables"), os.path.join(ROOT, "figures")
os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)


def load(name):
    p = os.path.join(ROOT, name)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


def fr(x, nd=1):
    """Nombre au format français."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "--"
    s = f"{x:,.{nd}f}".replace(",", " ").replace(".", ",")
    return s.replace(" ", "\\,")


def gap_pct(r):
    if r.get("obj") is None or r.get("dual") is None:
        return None
    return 100 * max(0.0, r["obj"] - r["dual"]) / max(abs(r["obj"]), 1e-9)


def stations_par_tournee(r):
    lens = [len(v) - 1 for v in r["routes"].values() if len(v) > 1]
    return st.mean(lens) if lens else 0.0


def freq(r):
    """Fréquence de visite (visites / période) par station."""
    return [sum(v) / r["H"] for v in r["visits"].values()]


# ---------------------------------------------------------------- validation
val = load("results_validation.jsonl")
if val:
    rows = defaultdict(list)
    for r in val:
        rows[(r["N"], r["H"])].append(r)
    lines = []
    for (N, H) in sorted(rows, key=lambda k: (k[1], k[0])):
        g = rows[(N, H)]
        dmax = max(abs(r[k] - r["enum"]) for r in g for k in ("bc_1", "bc_0", "mtz_0"))
        lines.append(f"{H} & {N} & {len(g)} & {len(g)} & {dmax:.1e} \\\\".replace("e-0", "e-").replace(".", ","))
    with open(os.path.join(TAB, "validation.tex"), "w") as f:
        f.write("\n".join(lines))
    print("validation:", len(val), "instances, écart max",
          max(abs(r[k] - r["enum"]) for r in val for k in ("bc_1", "bc_0", "mtz_0")))

# ---------------------------------------------------------------- principal
main = [r for r in load("results_main.jsonl") if r.get("windows") and r.get("status") not in ("error",)]
by = defaultdict(list)
for r in main:
    by[(r["family"], r["H"], r["N"])].append(r)


def agg(g):
    opt = sum(r["status"] == "optimal" for r in g)
    t = st.mean(r["time"] for r in g)
    gaps = [gap_pct(r) for r in g if gap_pct(r) is not None]
    nodes = st.mean(r.get("nodes", 0) for r in g)
    cuts = st.mean(r.get("cuts", 0) for r in g)
    return opt, len(g), t, (st.mean(gaps) if gaps else None), nodes, cuts


for fam in ("BG", "GF"):
    lines = []
    for H in (3, 6):
        for N in (5, 10, 15, 20, 25, 30):
            a, b = by.get((fam, H, N)), by.get(("U" + fam, H, N))
            if not a or not b:
                continue
            A, Bt = agg(a), agg(b)
            lines.append(f"{H} & {N} & {A[0]}/{A[1]} & {fr(A[2])} & {fr(A[3],2)} & {fr(A[4],0)} & "
                         f"{Bt[0]}/{Bt[1]} & {fr(Bt[2])} & {fr(Bt[3],2)} & {fr(Bt[4],0)} \\\\")
        if H == 3:
            lines.append("\\midrule")
    with open(os.path.join(TAB, f"res_{fam}.tex"), "w") as f:
        f.write("\n".join(lines))

# comparaison appariée stylisée / témoin
pair = {}
for r in main:
    key = r["name"].lstrip("U") if r["family"].startswith("U") else r["name"]
    pair.setdefault(key, {})["tem" if r["family"].startswith("U") else "sty"] = r
summary = {}
for fam in ("BG", "GF"):
    P = [v for k, v in pair.items() if k.startswith(fam) and "sty" in v and "tem" in v]
    if not P:
        continue
    ts = np.array([v["sty"]["time"] for v in P]); tt = np.array([v["tem"]["time"] for v in P])
    ns = np.array([v["sty"].get("nodes", 0) for v in P]); nt = np.array([v["tem"].get("nodes", 0) for v in P])
    ratio = np.exp(np.mean(np.log((tt + 1) / (ts + 1))))
    try:
        p = stats.wilcoxon(tt, ts).pvalue
    except ValueError:
        p = float("nan")
    harder = int(np.sum(tt > ts))
    cost_ratio = np.mean([v["tem"]["obj"] / v["sty"]["obj"] for v in P if v["tem"]["obj"] and v["sty"]["obj"]])
    opt_s = sum(v["sty"]["status"] == "optimal" for v in P); opt_t = sum(v["tem"]["status"] == "optimal" for v in P)
    summary[fam] = dict(n=len(P), ratio=ratio, p=p, harder=harder, cost_ratio=cost_ratio, opt_s=opt_s, opt_t=opt_t,
                        nodes_ratio=float(np.exp(np.mean(np.log((nt + 1) / (ns + 1))))))

# structure des solutions (instances résolues à l'optimum)
struct = {}
for fam in ("BG", "UBG", "GF", "UGF"):
    R = [r for r in main if r["family"] == fam and r["status"] == "optimal" and r.get("obj")]
    if not R:
        continue
    share = st.mean(r["transport"] / r["obj"] for r in R)
    spt = st.mean(stations_par_tournee(r) for r in R)
    fq = st.mean(st.mean(freq(r)) for r in R)
    rho_d, rho_dist = [], []
    for r in R:
        f_ = freq(r)
        if len(set(f_)) > 1:
            rho_d.append(stats.spearmanr(f_, r["d"]).statistic)
            rho_dist.append(stats.spearmanr(f_, r["dist0"]).statistic)
    struct[fam] = dict(n=len(R), share=share, spt=spt, freq=fq,
                       rho_d=(st.mean(rho_d) if rho_d else None), rho_dist=(st.mean(rho_dist) if rho_dist else None),
                       share_L=st.mean([r["transport"] / r["obj"] for r in R if r["eta"] == "L"] or [float("nan")]),
                       share_H=st.mean([r["transport"] / r["obj"] for r in R if r["eta"] == "H"] or [float("nan")]))
if struct:
    names = {"BG": "BG", "UBG": "U(BG)", "GF": "GF", "UGF": "U(GF)"}
    lines = []
    for fam, s in struct.items():
        lines.append(f"{names[fam]} & {s['n']} & {fr(100*s['share_L'])} & {fr(100*s['share_H'])} & {fr(s['spt'],2)} & "
                     f"{fr(s['freq'],2)} & {fr(s['rho_d'],2)} & {fr(s['rho_dist'],2)} \\\\")
    with open(os.path.join(TAB, "structure.tex"), "w") as f:
        f.write("\n".join(lines))

# ---------------------------------------------------------------- fenêtres
win = load("results_windows.jsonl")
if win:
    ref = {r["name"]: r for r in main}
    rows = defaultdict(list)
    for r in win:
        w = ref.get(r["name"])
        if w is None:
            continue
        best = min(x for x in (w.get("obj"), r.get("obj")) if x is not None)
        rb = lambda z: z["dual"] if (z.get("root_dual") is None or z["root_dual"] > 1e19) else z["root_dual"]
        rg_w = max(0.0, 100 * (best - rb(w)) / best)
        rg_n = max(0.0, 100 * (best - rb(r)) / best)
        rows[(r["family"], r["H"])].append((rg_n, rg_w, r["time"], w["time"], r["nodes"], w["nodes"],
                                            r["status"] == "optimal", w["status"] == "optimal"))
    lines = []
    for (fam, H) in sorted(rows):
        g = rows[(fam, H)]
        m = lambda k: st.mean(x[k] for x in g)
        lines.append(f"{fam} & {H} & {len(g)} & {fr(m(0),2)} & {fr(m(1),2)} & {fr(m(2))} & {fr(m(3))} & "
                     f"{fr(m(4),0)} & {fr(m(5),0)} & {sum(x[6] for x in g)} & {sum(x[7] for x in g)} \\\\")
    with open(os.path.join(TAB, "windows.tex"), "w") as f:
        f.write("\n".join(lines))

# ---------------------------------------------------------------- sensibilité
sens = load("results_sens.jsonl")
if sens:
    base = [r for r in main if r["family"] == "BG" and r["N"] == 15 and r["H"] == 6 and r["eta"] == "L"]
    groups = defaultdict(list)
    for r in base:
        groups[("base", "")].append(r)
    for r in sens:
        if r["kappa"] != 1.5:
            groups[("kappa", r["kappa"])].append(r)
        elif r["g"] != [2, 3]:
            groups[("g", "-".join(map(str, r["g"])))].append(r)
        elif r["phi"] != 1.3:
            groups[("phi", r["phi"])].append(r)
    order = [("base", "")] + [("kappa", 1.0), ("kappa", 2.0), ("g", "3-4"), ("g", "4-6"), ("phi", 1.2), ("phi", 1.4)]
    label = {"base": "Référence ($\\kappa=1{,}5$, $g\\in\\{2,3\\}$, $\\phi=1{,}3$)", "kappa": "$\\kappa=$", "g": "$g_i\\in$", "phi": "$\\phi=$"}
    lines = []
    for k in order:
        g = groups.get(k)
        if not g:
            continue
        lab = label[k[0]] if k[0] == "base" else (label[k[0]] + (fr(k[1], 1) if k[0] != "g" else "\\{" + k[1].replace("-", ",") + "\\}"))
        R = [r for r in g if r.get("obj")]
        lines.append(f"{lab} & {sum(r['status']=='optimal' for r in g)}/{len(g)} & {fr(st.mean(r['obj'] for r in R),0)} & "
                     f"{fr(100*st.mean(r['transport']/r['obj'] for r in R))} & {fr(st.mean(stations_par_tournee(r) for r in R),2)} & "
                     f"{fr(st.mean(st.mean(freq(r)) for r in R),2)} & {fr(st.mean(r['time'] for r in g))} \\\\")
    with open(os.path.join(TAB, "sens.tex"), "w") as f:
        f.write("\n".join(lines))

# ---------------------------------------------------------------- figure
if main:
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#888888", "axes.labelcolor": "#333333", "xtick.color": "#555555",
                         "ytick.color": "#555555"})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
    col = {"BG": "#2a78d6", "GF": "#eb6834"}
    for ax, H in zip(axes, (3, 6)):
        for fam in ("BG", "GF"):
            for pref, ls, mk, lab in (("", "-", "o", fam), ("U", "--", "s", f"U({fam})")):
                Ns, Ts = [], []
                for N in (5, 10, 15, 20, 25, 30):
                    g = by.get((pref + fam, H, N))
                    if g:
                        Ns.append(N); Ts.append(np.exp(np.mean(np.log([r["time"] + 0.01 for r in g]))))
                if Ns:
                    ax.plot(Ns, Ts, ls=ls, marker=mk, ms=4.5, lw=1.6, color=col[fam], label=lab,
                            mfc=col[fam] if pref == "" else "white")
        ax.set_yscale("log"); ax.set_xlabel("Nombre de stations $N$"); ax.set_title(f"$H={H}$", fontsize=9)
        ax.grid(axis="y", color="#e6e6e6", lw=0.6); ax.set_xticks([5, 10, 15, 20, 25, 30])
        ax.axhline(600, color="#999999", lw=0.8, ls=":")
    axes[0].set_ylabel("Temps de résolution (s, moy. géom.)")
    axes[1].text(30, 640, "limite 600 s", ha="right", va="bottom", fontsize=7, color="#666666")
    axes[1].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "temps.pdf")); fig.savefig(os.path.join(FIG, "temps.png"), dpi=150)

json.dump(dict(summary=summary, struct=struct), open(os.path.join(ROOT, "summary.json"), "w"), indent=1, default=float)
print(json.dumps(dict(summary=summary, struct=struct), indent=1, default=float))
