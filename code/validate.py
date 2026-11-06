"""Validation croisée : branch-and-cut (avec/sans fenêtres), MTZ et énumération."""
import json, sys, time
import gen
from irp import solve
from enum_check import optimum

out = open(sys.argv[1], "w")
cases = [(N, 3) for N in (3, 4, 5)] + [(N, 6) for N in (3, 4)]
for fam in ("BG", "GF"):
    for (N, H) in cases:
        for eta in ("L", "H"):
            for rep in (1, 2, 3):
                for inst in gen.generer(fam, N, H, eta, rep):
                    t = time.time(); e = optimum(inst); te = time.time() - t
                    row = dict(name=inst["name"], N=N, H=H, enum=e, t_enum=te)
                    for form, w in (("bc", True), ("bc", False), ("mtz", False)):
                        r = solve(inst, form, w, time_limit=300)
                        row[f"{form}_{int(w)}"] = r["obj"]; row[f"{form}_{int(w)}_st"] = r["status"]
                    out.write(json.dumps(row) + "\n"); out.flush()
