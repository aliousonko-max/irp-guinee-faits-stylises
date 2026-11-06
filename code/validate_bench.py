"""Validation sur les instances de référence d'Archetti et al. (2007), version multi-véhicule
(K=2, politique ML), contre les valeurs optimales du branch-and-cut publiées pour le DIMACS IRP Challenge."""
import json, sys, openpyxl
from irp_mv import read_dat, solve_mv
REPO = sys.argv[1]; out = open(sys.argv[2], "w")
wb = openpyxl.load_workbook(REPO + "/DIMACS_Results_IRP.xlsx", read_only=True, data_only=True)
ref = {}
for r in list(wb["Summary"].iter_rows(values_only=True))[1:]:
    if r[0] and r[0].strip() == "BC":
        ref[(r[1].strip(), r[2], r[3], r[4])] = (r[7], r[8], r[11])   # LB, UB, Optimal
folders = {"High": "Istanze0105", "Low": "Istanze001005"}
for H, Ns in ((3, (5, 10, 15)), (6, (5, 10))):
    for cls in ("High", "Low"):
        for N in Ns:
            for a in range(1, 6):
                path = f"{REPO}/Data/Small/{folders[cls]}h{H}/abs{a}n{N}_1.dat"
                inst = read_dat(path, 2)
                res = solve_mv(inst, "ML", "same", 600)
                lb, ub, opt = ref.get((f"{cls}T{H}", a, N, 2), (None, None, None))
                res.update(cls=cls, H=H, N=N, abs=a, ref_lb=lb, ref_ub=ub, ref_opt=opt,
                           obj_dimacs=(res["obj"] + res["const0"]) if res["obj"] is not None else None)
                out.write(json.dumps(res) + "\n"); out.flush()
                print(res["name"], cls, H, res["status"], round(res["obj_dimacs"], 2), "réf", ub, opt, round(res["time"], 1), flush=True)
