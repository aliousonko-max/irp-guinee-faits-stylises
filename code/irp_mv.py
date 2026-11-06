"""
Variante multi-véhicule (K camions homogènes) du modèle, utilisée uniquement
pour la validation sur les instances de référence d'Archetti et al. (2007),
dans leur version multi-véhicule (Coelho et Laporte, 2013), pour lesquelles
des valeurs optimales sont publiées. Politique 'ML' (maximum level) ou 'OU'.
"""
import math, time
import networkx as nx
from pyscipopt import Model, Conshdlr, SCIP_RESULT, quicksum

EPS = 1e-6


def read_dat(path, K):
    """Lit un fichier .dat (format Archetti) ; la capacité lue est celle d'un camion."""
    L = [l.split() for l in open(path) if l.strip()]
    n1, H, Q = int(L[0][0]), int(L[0][1]), int(float(L[0][2]))
    dep = L[1]
    pts = [(float(dep[1]), float(dep[2]))]
    B0, r0, h0 = float(dep[3]), float(dep[4]), float(dep[5])
    I0, U, Lo, d, h = [], [], [], [], []
    for row in L[2:2 + n1 - 1]:
        pts.append((float(row[1]), float(row[2])))
        I0.append(float(row[3])); U.append(float(row[4])); Lo.append(float(row[5]))
        d.append(float(row[6])); h.append(float(row[7]))
    n = len(pts)
    c = [[int(math.sqrt((pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2) + 0.5) for j in range(n)] for i in range(n)]
    return dict(name=path.split("/")[-1], N=n - 1, H=H, K=K, Q=Q, c=c, d=d, C=U, L=Lo, I0=I0,
                r=r0, B0=B0, h0=h0, h=h)


class GSECmv(Conshdlr):
    def __init__(self, data):
        self.d = data; self.ncuts = 0

    def _viol(self, sol, integral):
        x, y, V, Vc, T, Ks = (self.d[k] for k in ("x", "y", "V", "Vc", "T", "Ks"))
        m = self.model; cuts = []
        for t in T:
            for k in Ks:
                G = nx.DiGraph(); G.add_nodes_from(V)
                for i in V:
                    for j in V:
                        if i != j:
                            v = m.getSolVal(sol, x[i, j, k, t])
                            if v > EPS:
                                G.add_edge(i, j, capacity=v)
                if integral:
                    for comp in nx.connected_components(G.to_undirected()):
                        if 0 in comp:
                            continue
                        mm = max(comp, key=lambda a: m.getSolVal(sol, y[a, k, t]))
                        if m.getSolVal(sol, y[mm, k, t]) > 0.5:
                            cuts.append((frozenset(comp), mm, k, t))
                else:
                    done = set()
                    for mm in Vc:
                        ym = m.getSolVal(sol, y[mm, k, t])
                        if ym < 0.01 or mm in done:
                            continue
                        if not nx.has_path(G, 0, mm):
                            S = (set(nx.ancestors(G, mm)) | {mm}) - {0}
                            cuts.append((frozenset(S), mm, k, t)); done |= S; continue
                        val, (A, B) = nx.minimum_cut(G, 0, mm)
                        if val < ym - 1e-3:
                            cuts.append((frozenset(B), mm, k, t)); done |= set(B)
        return cuts

    def _add(self, cuts):
        x, y, V = self.d["x"], self.d["y"], self.d["V"]
        for S, mm, k, t in set(cuts):
            self.model.addCons(quicksum(x[i, j, k, t] for i in V if i not in S for j in S) >= y[mm, k, t], removable=True)
            self.ncuts += 1

    def conscheck(self, constraints, solution, checkintegrality, checklprows, printreason, completely, **kw):
        return {"result": SCIP_RESULT.INFEASIBLE if self._viol(solution, True) else SCIP_RESULT.FEASIBLE}

    def consenfolp(self, constraints, nusefulconss, solinfeasible):
        c = self._viol(None, True)
        if c:
            self._add(c); return {"result": SCIP_RESULT.CONSADDED}
        return {"result": SCIP_RESULT.FEASIBLE}

    def consenfops(self, constraints, nusefulconss, solinfeasible, objinfeasible):
        return {"result": SCIP_RESULT.INFEASIBLE if self._viol(None, True) else SCIP_RESULT.FEASIBLE}

    def conssepalp(self, constraints, nusefulconss):
        c = self._viol(None, False)
        if c:
            self._add(c); return {"result": SCIP_RESULT.CONSADDED}
        return {"result": SCIP_RESULT.DIDNOTFIND}

    def conslock(self, constraint, locktype, nlockspos, nlocksneg):
        for v in list(self.d["x"].values()) + list(self.d["y"].values()):
            self.model.addVarLocks(v, nlockspos + nlocksneg, nlockspos + nlocksneg)


def solve_mv(inst, policy="ML", supplier="same", time_limit=600):
    """supplier='same' : la production de la période t est disponible en t ;
       supplier='before' : les livraisons de t sont limitées au stock de début de période."""
    t0 = time.time()
    N, H, K = inst["N"], inst["H"], inst["K"]
    Vc = list(range(1, N + 1)); V = [0] + Vc; T = list(range(1, H + 1)); Ks = list(range(1, K + 1))
    c, d, U, Lo, I0, Q = inst["c"], inst["d"], inst["C"], inst["L"], inst["I0"], inst["Q"]
    m = Model(); m.hideOutput(); m.setParam("limits/time", time_limit); m.setParam("parallel/maxnthreads", 1)
    x = {(i, j, k, t): m.addVar(vtype="B") for i in V for j in V if i != j for k in Ks for t in T}
    y = {(i, k, t): m.addVar(vtype="B") for i in V for k in Ks for t in T}
    q = {(i, k, t): m.addVar(lb=0) for i in Vc for k in Ks for t in T}
    I = {(i, t): m.addVar(lb=Lo[i - 1]) for i in Vc for t in T}
    B = {t: m.addVar(lb=0) for t in T}
    Ip = lambda i, t: I0[i - 1] if t == 1 else I[i, t - 1]
    Bp = lambda t: inst["B0"] if t == 1 else B[t - 1]
    m.setObjective(quicksum(c[i][j] * v for (i, j, k, t), v in x.items())
                   + quicksum(inst["h0"] * B[t] for t in T)
                   + quicksum(inst["h"][i - 1] * I[i, t] for i in Vc for t in T), "minimize")
    for t in T:
        tot = quicksum(q[i, k, t] for i in Vc for k in Ks)
        m.addCons(B[t] == Bp(t) + inst["r"] - tot)
        if supplier == "before":
            m.addCons(tot <= Bp(t))
        for i in Vc:
            qi = quicksum(q[i, k, t] for k in Ks); yi = quicksum(y[i, k, t] for k in Ks)
            m.addCons(I[i, t] == Ip(i, t) + qi - d[i - 1])
            m.addCons(Ip(i, t) + qi <= U[i - 1])
            m.addCons(yi <= 1)
            for k in Ks:
                m.addCons(q[i, k, t] <= min(U[i - 1], Q) * y[i, k, t])
                m.addCons(y[i, k, t] <= y[0, k, t])
            if policy == "OU":
                m.addCons(qi >= U[i - 1] * yi - Ip(i, t))
        for k in Ks:
            m.addCons(quicksum(q[i, k, t] for i in Vc) <= Q * y[0, k, t])
            for i in V:
                m.addCons(quicksum(x[i, j, k, t] for j in V if j != i) == y[i, k, t])
                m.addCons(quicksum(x[j, i, k, t] for j in V if j != i) == y[i, k, t])
            if k < K:
                m.addCons(y[0, k, t] >= y[0, k + 1, t])
    hdl = GSECmv(dict(x=x, y=y, V=V, Vc=Vc, T=T, Ks=Ks))
    m.includeConshdlr(hdl, "GSECmv", "subtours", sepapriority=1000000, enfopriority=-1000000,
                      chckpriority=-1000000, sepafreq=1, needscons=False)
    m.optimize()
    const0 = sum(inst["h"][i] * I0[i] for i in range(N)) + inst["h0"] * inst["B0"]
    return dict(name=inst["name"], K=K, policy=policy, supplier=supplier, status=m.getStatus(),
                obj=(m.getObjVal() if m.getNSols() else None), dual=m.getDualbound(),
                time=time.time() - t0, const0=const0)
