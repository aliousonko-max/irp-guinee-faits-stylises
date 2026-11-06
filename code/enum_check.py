"""
Validation exacte par énumération (fondée sur la proposition 3) :
pour chaque calendrier de visites réalisable, les quantités et les stocks
sont déterminés ; le coût de routage de chaque période est un TSP sur les
stations visitées, résolu exactement par programmation dynamique.
Utilisable pour N <= 6 environ.
"""
import itertools, functools, json, sys


def held_karp(c, nodes):
    nodes = list(nodes)
    if not nodes:
        return 0
    n = len(nodes)
    idx = {v: k for k, v in enumerate(nodes)}
    INF = float("inf")
    dp = {}
    for k, v in enumerate(nodes):
        dp[(1 << k, k)] = c[0][v]
    for mask in range(1, 1 << n):
        for k in range(n):
            if not (mask >> k) & 1 or (mask, k) not in dp:
                continue
            base = dp[(mask, k)]
            for l in range(n):
                if (mask >> l) & 1:
                    continue
                key = (mask | (1 << l), l)
                val = base + c[nodes[k]][nodes[l]]
                if val < dp.get(key, INF):
                    dp[key] = val
    full = (1 << n) - 1
    return min(dp[(full, k)] + c[nodes[k]][0] for k in range(n))


def optimum(inst):
    N, H = inst["N"], inst["H"]
    c, d, C, I0, Q, h = inst["c"], inst["d"], inst["C"], inst["I0"], inst["Q"], inst["h"]
    r, B0, h0 = inst["r"], inst["B0"], inst["h0"]
    # motifs réalisables par station (stock >= 0)
    pats = []
    for i in range(N):
        ok = []
        for pat in itertools.product((0, 1), repeat=H):
            I, qs, hold, feas = I0[i], [], 0.0, True
            for t in range(H):
                q = (C[i] - I) if pat[t] else 0
                I = I + q - d[i]
                if I < 0:
                    feas = False; break
                qs.append(q); hold += h[i] * I
            if feas:
                ok.append((pat, qs, hold))
        pats.append(ok)
    tsp = functools.lru_cache(maxsize=None)(lambda S: held_karp(c, S))
    best, arg = float("inf"), None
    for combo in itertools.product(*pats):
        load = [sum(p[1][t] for p in combo) for t in range(H)]
        if any(l > Q for l in load):
            continue
        B, hold0, feas = B0, 0.0, True
        for t in range(H):
            B = B + r - load[t]
            if B < 0:
                feas = False; break
            hold0 += h0 * B
        if not feas:
            continue
        cost = sum(p[2] for p in combo) + hold0
        if cost >= best:
            continue
        for t in range(H):
            S = frozenset(i + 1 for i in range(N) if combo[i][0][t])
            cost += tsp(S)
            if cost >= best:
                break
        if cost < best:
            best, arg = cost, combo
    return best


if __name__ == "__main__":
    print(optimum(json.load(open(sys.argv[1]))))
