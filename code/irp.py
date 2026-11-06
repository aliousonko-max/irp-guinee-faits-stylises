"""
IRP mono-véhicule, politique order-up-to-level (OU).
Résolution par branch-and-cut avec SCIP (PySCIPOpt) :
  - contraintes d'élimination des sous-tours (GSEC) ajoutées dynamiquement
    par un gestionnaire de contraintes (séparation entière : composantes
    connexes ; séparation fractionnaire : flot maximal / coupe minimale) ;
  - inégalités de fenêtre (proposition 4) ajoutées a priori (option).
Variante compacte MTZ pour la validation croisée.
"""
import json, math, sys, time
import networkx as nx
from pyscipopt import Model, Conshdlr, SCIP_RESULT, quicksum

EPS = 1e-6


class GSEC(Conshdlr):
    """Gestionnaire des contraintes  sum_{i notin S, j in S} x_ijt >= y_mt."""

    def __init__(self, data):
        self.d = data          # dict avec x, y, V, Vc, T
        self.ncuts = 0

    # -- séparation ---------------------------------------------------------
    def _violations(self, sol, integral):
        x, y, V, Vc, T = self.d["x"], self.d["y"], self.d["V"], self.d["Vc"], self.d["T"]
        m = self.model
        cuts = []
        for t in T:
            G = nx.DiGraph()
            G.add_nodes_from(V)
            for (i, j, tt), var in x.items():
                if tt != t:
                    continue
                v = m.getSolVal(sol, var)
                if v > EPS:
                    G.add_edge(i, j, capacity=v)
            if integral:
                und = G.to_undirected()
                for comp in nx.connected_components(und):
                    if 0 in comp:
                        continue
                    S = set(comp)
                    mm = max(S, key=lambda k: m.getSolVal(sol, y[k, t]))
                    if m.getSolVal(sol, y[mm, t]) > 0.5:
                        cuts.append((frozenset(S), mm, t))
            else:
                done = set()
                for mm in Vc:
                    ym = m.getSolVal(sol, y[mm, t])
                    if ym < 0.01 or mm in done:
                        continue
                    if mm not in G or not nx.has_path(G, 0, mm):
                        S = set(nx.ancestors(G, mm)) | {mm} if mm in G else {mm}
                        S.discard(0)
                        cuts.append((frozenset(S), mm, t)); done |= S
                        continue
                    val, (A, B) = nx.minimum_cut(G, 0, mm)
                    if val < ym - 1e-3:
                        S = set(B)
                        cuts.append((frozenset(S), mm, t)); done |= S
        return cuts

    def _add(self, cuts):
        x, y, V = self.d["x"], self.d["y"], self.d["V"]
        m = self.model
        seen = set()
        for S, mm, t in cuts:
            if (S, mm, t) in seen:
                continue
            seen.add((S, mm, t))
            expr = quicksum(x[i, j, t] for i in V if i not in S for j in S if i != j)
            m.addCons(expr >= y[mm, t], removable=True, name=f"gsec_{t}_{self.ncuts}")
            self.ncuts += 1

    def _is_integral(self, sol):
        m = self.model
        return all(abs(m.getSolVal(sol, v) - round(m.getSolVal(sol, v))) < 1e-5 for v in self.d["x"].values())

    # -- callbacks SCIP -----------------------------------------------------
    def conscheck(self, constraints, solution, checkintegrality, checklprows, printreason, completely, **kw):
        cuts = self._violations(solution, integral=True)
        return {"result": SCIP_RESULT.INFEASIBLE if cuts else SCIP_RESULT.FEASIBLE}

    def consenfolp(self, constraints, nusefulconss, solinfeasible):
        cuts = self._violations(None, integral=True)
        if cuts:
            self._add(cuts)
            return {"result": SCIP_RESULT.CONSADDED}
        return {"result": SCIP_RESULT.FEASIBLE}

    def consenfops(self, constraints, nusefulconss, solinfeasible, objinfeasible):
        cuts = self._violations(None, integral=True)
        return {"result": SCIP_RESULT.INFEASIBLE if cuts else SCIP_RESULT.FEASIBLE}

    def conssepalp(self, constraints, nusefulconss):
        cuts = self._violations(None, integral=False)
        if cuts:
            self._add(cuts)
            return {"result": SCIP_RESULT.CONSADDED}
        return {"result": SCIP_RESULT.DIDNOTFIND}

    def conslock(self, constraint, locktype, nlockspos, nlocksneg):
        for v in list(self.d["x"].values()) + list(self.d["y"].values()):
            self.model.addVarLocks(v, nlockspos + nlocksneg, nlockspos + nlocksneg)


def build(inst, form="bc", windows=True, logical=True, time_limit=600, verbose=False):
    N, H = inst["N"], inst["H"]
    Vc = list(range(1, N + 1)); V = [0] + Vc; T = list(range(1, H + 1))
    c, d, C, I0, Q = inst["c"], inst["d"], inst["C"], inst["I0"], inst["Q"]
    h, h0, r, B0 = inst["h"], inst["h0"], inst["r"], inst["B0"]
    m = Model(f"IRP_{inst['name']}")
    if not verbose:
        m.hideOutput()
    m.setParam("limits/time", time_limit)
    m.setParam("parallel/maxnthreads", 1)

    x = {(i, j, t): m.addVar(vtype="B", name=f"x_{i}_{j}_{t}") for i in V for j in V if i != j for t in T}
    y = {(i, t): m.addVar(vtype="B", name=f"y_{i}_{t}") for i in V for t in T}
    q = {(i, t): m.addVar(lb=0, name=f"q_{i}_{t}") for i in Vc for t in T}
    I = {(i, t): m.addVar(lb=0, name=f"I_{i}_{t}") for i in Vc for t in T}
    B = {t: m.addVar(lb=0, name=f"B_{t}") for t in T}
    Iprev = lambda i, t: I0[i - 1] if t == 1 else I[i, t - 1]
    Bprev = lambda t: B0 if t == 1 else B[t - 1]

    m.setObjective(
        quicksum(c[i][j] * x[i, j, t] for (i, j, t) in x)
        + quicksum(h0 * B[t] for t in T)
        + quicksum(h[i - 1] * I[i, t] for i in Vc for t in T), "minimize")

    for t in T:
        m.addCons(B[t] == Bprev(t) + r - quicksum(q[i, t] for i in Vc))                       # (2)
        for i in Vc:
            m.addCons(I[i, t] == Iprev(i, t) + q[i, t] - d[i - 1])                             # (3)
            m.addCons(Iprev(i, t) + q[i, t] <= C[i - 1])                                       # (4)
            m.addCons(q[i, t] >= C[i - 1] * y[i, t] - Iprev(i, t))                             # (5)
            m.addCons(q[i, t] <= C[i - 1] * y[i, t])                                           # (6)
        m.addCons(quicksum(q[i, t] for i in Vc) <= Q * y[0, t])                                # (7)
        for i in V:
            m.addCons(quicksum(x[i, j, t] for j in V if j != i) == y[i, t])                    # (8)
            m.addCons(quicksum(x[j, i, t] for j in V if j != i) == y[i, t])
        if logical:
            for i in Vc:
                m.addCons(y[i, t] <= y[0, t])
                for j in V:
                    if j != i:
                        m.addCons(x[i, j, t] <= y[i, t]); m.addCons(x[j, i, t] <= y[i, t])

    nwin = 0
    if windows:                                                                                 # (11)-(12)
        for i in Vc:
            Ci, di = C[i - 1], d[i - 1]
            for t in T:
                rhs = math.ceil((t * di - I0[i - 1]) / Ci - 1e-9)
                if rhs >= 1:
                    m.addCons(quicksum(y[i, u] for u in range(1, t + 1)) >= rhs); nwin += 1
                for s in range(1, t):
                    rhs = math.ceil(((t - s + 1) * di - Ci) / Ci - 1e-9)
                    if rhs >= 1:
                        m.addCons(quicksum(y[i, u] for u in range(s + 1, t + 1)) >= rhs); nwin += 1

    data = dict(x=x, y=y, V=V, Vc=Vc, T=T)
    hdl = None
    if form == "bc":
        hdl = GSEC(data)
        m.includeConshdlr(hdl, "GSEC", "subtour elimination", sepapriority=1000000,
                          enfopriority=-1000000, chckpriority=-1000000, sepafreq=1, needscons=False)
    elif form == "mtz":
        u = {(i, t): m.addVar(lb=1, ub=N, name=f"u_{i}_{t}") for i in Vc for t in T}
        for t in T:
            for i in Vc:
                for j in Vc:
                    if i != j:
                        m.addCons(u[i, t] - u[j, t] + N * x[i, j, t] <= N - 1)
    return m, data, hdl, nwin


def solve(inst, form="bc", windows=True, time_limit=600, verbose=False):
    t0 = time.time()
    m, data, hdl, nwin = build(inst, form, windows, time_limit=time_limit, verbose=verbose)
    m.optimize()
    wall = time.time() - t0
    status = m.getStatus()
    res = dict(name=inst["name"], form=form, windows=windows, status=status, time=wall,
               nodes=m.getNNodes(), cuts=(hdl.ncuts if hdl else 0), nwin=nwin,
               root_dual=None)
    try:
        res["root_dual"] = m.getDualboundRoot()
    except Exception:
        pass
    try:
        res["dual"] = m.getDualbound()
    except Exception:
        res["dual"] = None
    if m.getNSols() > 0:
        sol = m.getBestSol()
        res["obj"] = m.getSolObjVal(sol)
        res["gap"] = m.getGap()
        V, Vc, T = data["V"], data["Vc"], data["T"]
        x, y = data["x"], data["y"]
        c = inst["c"]
        tr = sum(c[i][j] for (i, j, t), v in x.items() if m.getSolVal(sol, v) > 0.5)
        res["transport"] = tr
        res["holding"] = res["obj"] - tr
        res["visits"] = {i: [int(m.getSolVal(sol, y[i, t]) > 0.5) for t in T] for i in Vc}
        res["routes"] = {}
        for t in T:
            succ = {i: j for (i, j, tt), v in x.items() if tt == t and m.getSolVal(sol, v) > 0.5}
            route, cur = [0], succ.get(0)
            while cur not in (None, 0) and len(route) <= len(V):
                route.append(cur); cur = succ.get(cur)
            res["routes"][t] = route
    else:
        res["obj"] = None; res["gap"] = None
    try:
        res["root_lp"] = m.getFirstLpTime()
    except Exception:
        pass
    return res


if __name__ == "__main__":
    inst = json.load(open(sys.argv[1]))
    print(json.dumps(solve(inst, verbose=True), default=str, indent=1))
