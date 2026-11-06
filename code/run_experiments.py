"""
Lance les expérimentations (reprise possible : les instances déjà traitées
sont ignorées). Usage : python3 run_experiments.py <plan> <fichier_resultats> [n_workers]
Plans : main | windows | sens
"""
import json, os, sys
from multiprocessing import Pool
import gen
from irp import solve

TL = 600


def plan_main():
    for fam in ("BG", "GF"):
        for N in (5, 10, 15, 20, 25, 30):
            for H in (3, 6):
                for eta in ("L", "H"):
                    for rep in (1, 2, 3):
                        yield dict(fam=fam, N=N, H=H, eta=eta, rep=rep, which="both", windows=True)


def plan_windows():
    for fam in ("BG", "GF"):
        for N in (10, 15, 20):
            for H in (3, 6):
                for rep in (1, 2, 3):
                    yield dict(fam=fam, N=N, H=H, eta="L", rep=rep, which="sty", windows=False)


def plan_sens():
    base = dict(fam="BG", N=15, H=6, eta="L", which="sty", windows=True)
    for rep in (1, 2, 3):
        for kappa in (1.0, 2.0):
            yield dict(base, rep=rep, kappa=kappa)
        for g in ((3, 4), (4, 6)):
            yield dict(base, rep=rep, g=g)
        for phi in (1.2, 1.4):
            yield dict(base, rep=rep, phi=phi)


def instance_name(p):
    kw = {}
    if "kappa" in p: kw["kappa"] = p["kappa"]
    if "g" in p: kw["g_choices"] = tuple(p["g"])
    if "phi" in p: kw["phi"] = p["phi"]
    sty, tem = gen.generer(p["fam"], p["N"], p["H"], p["eta"], p["rep"], **kw)
    return tem["name"] if p["which"] == "tem" else sty["name"]


def job(p):
    """Résout une instance dans un sous-processus, avec arrêt forcé (garde-fou)."""
    import subprocess
    q = dict(p, TL=TL)
    try:
        out = subprocess.run([sys.executable, "solve_one.py", json.dumps(q)], capture_output=True, text=True,
                             timeout=TL + 180, cwd=os.path.dirname(os.path.abspath(__file__)))
        for line in out.stdout.splitlines():
            if line.startswith("RESULT "):
                return json.loads(line[7:])
        return dict(name=instance_name(p), windows=p["windows"], status="error", time=None, err=out.stderr[-500:])
    except subprocess.TimeoutExpired:
        return dict(name=instance_name(p), windows=p["windows"], status="killed", time=TL + 180, obj=None,
                    family=instance_name(p).split("_")[0], N=p["N"], H=p["H"], eta=p["eta"], rep=p["rep"])


if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor
    plan = {"main": plan_main, "windows": plan_windows, "sens": plan_sens}[sys.argv[1]]
    path = sys.argv[2]
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    done = set()
    if os.path.exists(path):
        for l in open(path):
            try:
                r = json.loads(l); done.add((r["name"], r["windows"]))
            except Exception:
                pass
    todo = []
    for p in plan():
        for w in (("sty", "tem") if p["which"] == "both" else ("sty",)):
            q = dict(p, which=w)
            if (instance_name(q), q["windows"]) not in done:
                todo.append(q)
    todo.sort(key=lambda p: -(p["N"] * p["H"]))
    print(len(todo), "instances", flush=True)
    from concurrent.futures import as_completed
    with ThreadPoolExecutor(nw) as ex, open(path, "a") as f:
        for fut in as_completed([ex.submit(job, p) for p in todo]):
            r = fut.result()
            f.write(json.dumps(r, default=str) + "\n"); f.flush()
            print(r["name"], r["status"], r.get("time"), flush=True)
