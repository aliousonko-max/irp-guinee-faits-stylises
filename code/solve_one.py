"""Résout une seule instance (appelé en sous-processus par run_experiments.py)."""
import json, sys
import gen
from irp import solve
p = json.loads(sys.argv[1])
kw = {}
if "kappa" in p: kw["kappa"] = p["kappa"]
if "g" in p: kw["g_choices"] = tuple(p["g"])
if "phi" in p: kw["phi"] = p["phi"]
sty, tem = gen.generer(p["fam"], p["N"], p["H"], p["eta"], p["rep"], **kw)
inst = tem if p["which"] == "tem" else sty
r = solve(inst, "bc", p["windows"], time_limit=p["TL"])
r.update(family=inst["name"].split("_")[0], N=inst["N"], H=inst["H"], eta=inst["eta"], rep=inst["rep"],
         kappa=inst["kappa"], phi=inst["phi"], g=inst["g_choices"], d=inst["d"],
         dist0=[inst["c"][0][i] for i in range(1, inst["N"] + 1)], loc=inst["loc"])
print("RESULT " + json.dumps(r, default=str))
