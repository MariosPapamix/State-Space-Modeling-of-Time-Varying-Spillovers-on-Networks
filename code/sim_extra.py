"""
Additional simulation studies:
  S1d  Correctly specified Gaussian run (state path drawn from the random walk,
       filter uses the true Q, sigma): verifies Theorem 2(i) exact coverage and
       Theorem 2(v) steady-state level in their own regime (the S1 jump design
       deliberately violates them, so it cannot).
  S4x  Three further design points for the Poisson-Laplace calibration (dense
       graph; large q; small N,T), at further design points.
  S7   Size and power of the randomization test when unit heterogeneity is
       ALIGNED WITH THE GRAPH: unrestricted
       label test vs the stratified conditional test, and power under a
       genuine spillover.  Gaussian NSSM, q fixed and disclosed.
  S1b  Table fragment for the five-topology study (previously narrative only).
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, math, sys, time
import numpy as np
sys.path.insert(0, _HERE)
from nssm import row_normalise, poisson_laplace_filter

t0 = time.time(); OUT = _os.path.join(_HERE, "out"); TAB = _os.path.join(_ROOT, "paper", "tables")
rng0 = np.random.default_rng(0)

def er_graph(N, d, rng):
    A = (rng.random((N, N)) < d / (N - 1)).astype(float)
    A = np.triu(A, 1); A = A + A.T
    iso = np.where(A.sum(1) == 0)[0]
    for i in iso:
        j = rng.integers(N - 1); j = j + (j >= i)
        A[i, j] = A[j, i] = 1
    return row_normalise(A)

def kalman_rw(Y, X_list, Q, s2, m0, P0):
    T1 = len(X_list); K = len(m0)
    m, P = m0.copy(), P0.copy()
    ms = np.zeros((T1, K)); Ps = np.zeros((T1, K, K)); ll = 0.0
    for t in range(T1):
        P = P + Q
        X = X_list[t]; y = Y[t]
        Sm = X @ P @ X.T + s2 * np.eye(len(y))
        Kg = np.linalg.solve(Sm, X @ P).T
        r = y - X @ m
        ll += float(-0.5 * (len(y) * math.log(2 * math.pi) + np.linalg.slogdet(Sm)[1]
                            + r @ np.linalg.solve(Sm, r)))
        m = m + Kg @ r; P = P - Kg @ Sm @ Kg.T
        ms[t] = m; Ps[t] = P
    return ms, Ps, ll

# ---------------- S1d: correctly specified ----------------
def s1d(R=200, N=200, T=60, q=(1e-4, 1e-3, 1e-3), sig=1.0, seed=1):
    rng = np.random.default_rng(seed)
    cov = []; ratio = []
    Q = np.diag(q)
    for rep in range(R):
        W = er_graph(N, 4, rng)
        th = np.array([0.0, 0.40, 0.30]); y = rng.normal(0, 1, N)
        Ys, Xs, ths = [], [], []
        for t in range(T):
            X = np.column_stack([np.ones(N), W @ y, y])
            y = X @ th + rng.normal(0, sig, N)
            Ys.append(y); Xs.append(X); ths.append(th.copy())
            th = th + rng.normal(0, np.sqrt(np.array(q)))
        ms, Ps, _ = kalman_rw(np.array(Ys), Xs, Q, sig ** 2, np.zeros(3), np.eye(3))
        tr = np.array(ths)
        z = (tr[5:, 1] - ms[5:, 1]) / np.sqrt(Ps[5:, 1, 1])
        cov.append(np.mean(np.abs(z) <= 1.96))
        Ibar = np.mean([np.linalg.inv(np.cov(np.array([Xs[t][:, 1]]))) if False else
                        (Xs[t][:, 1] @ Xs[t][:, 1] - (Xs[t][:, 1] @ np.column_stack(
                            [np.ones(N), Xs[t][:, 2]])) @ np.linalg.lstsq(
                            np.column_stack([np.ones(N), Xs[t][:, 2]]),
                            Xs[t][:, 1], rcond=None)[0]) / sig ** 2 for t in range(T // 2, T)])
        Pst = (math.sqrt(q[1] ** 2 * Ibar ** 2 + 4 * q[1] * Ibar) - q[1] * Ibar) / (2 * Ibar)
        ratio.append(Ps[-1, 1, 1] / Pst)
    return float(np.mean(cov)), float(np.std(cov) / math.sqrt(R)), float(np.mean(ratio))
c95, c95se, vr = s1d()
print(f"[S1d] coverage {c95:.3f} (MC-SE {c95se:.3f}); v_T/P* = {vr:.2f} ({time.time()-t0:.0f}s)", flush=True)

# ---------------- S4x: Poisson-Laplace calibration sweep ----------------
def s4_point(N, T, dexp, q, R=150, seed=2):
    rng = np.random.default_rng(seed)
    cov = []; rmse = []
    for rep in range(R):
        W = er_graph(N, dexp, rng)
        th = np.array([-0.3, 0.30, 0.30]); y = rng.poisson(2.0, N).astype(float)
        Ys, Xs, ths = [], [], []
        for t in range(T):
            X = np.column_stack([np.ones(N), W @ np.log1p(y), np.log1p(y)])
            lam = np.exp(np.clip(X @ th, -20, 20))
            y = rng.poisson(lam).astype(float)
            Ys.append(y); Xs.append(X); ths.append(th.copy())
            th = th + rng.normal(0, math.sqrt(q), 3)
        fit = poisson_laplace_filter(np.array(Ys, dtype=int), Xs, q * np.eye(3),
                                     np.zeros(3), np.eye(3))
        tr = np.array(ths); ms = fit["m_filt"]; Ps = fit["P_filt"]
        z = (tr[5:, 1] - ms[5:, 1]) / np.sqrt(Ps[5:, 1, 1])
        cov.append(np.mean(np.abs(z) <= 1.96))
        rmse.append(np.sqrt(np.mean((tr[5:, 1] - ms[5:, 1]) ** 2)))
    return float(np.mean(cov)), float(np.std(cov) / math.sqrt(R)), float(np.mean(rmse))
pts = [("sparse, $N=200$, $q=10^{-3}$ (S4)", 200, 60, 4, 1e-3),
       ("dense ($d=N/4$), $N=200$", 200, 60, 50, 1e-3),
       ("sparse, $N=200$, $q=10^{-2}$", 200, 60, 4, 1e-2),
       ("sparse, $N=100$, $T=40$", 100, 40, 4, 1e-3)]
s4rows = []
for lab, N_, T_, d_, q_ in pts:
    c, cse, rm = s4_point(N_, T_, d_, q_)
    s4rows.append((lab, c, cse, rm))
    print(f"[S4x] {lab}: cov {c:.3f} ({cse:.3f}) rmse {rm:.3f} ({time.time()-t0:.0f}s)", flush=True)
open(f"{TAB}/tab_s4ext.tex", "w").write("\n".join(
    f"{lab} & {c:.3f} & ({cse:.3f}) & {rm:.3f} \\\\" for lab, c, cse, rm in s4rows) + "\n")

# ---------------- S7: randomization size/power ----------------
def s7(M=150, N=60, T=60, q=1e-3, seed=3):
    rng = np.random.default_rng(seed)
    out = {}
    for scen, beta1 in [("size", 0.0), ("power", 0.25)]:
        rej_u = 0; rej_s = 0
        for rep in range(M):
            W = er_graph(N, 4, rng)
            a = 2.0 * np.linalg.solve(np.eye(N) - 0.8 * W, rng.normal(0, 1, N))  # graph-smooth effects
            y = rng.normal(0, 1, N)
            Ys, Xs = [], []
            for t in range(T):
                X = np.column_stack([np.ones(N), W @ y, y])
                y = a + beta1 * (W @ y) + 0.30 * y + rng.normal(0, 1, N)
                Ys.append(y); Xs.append(X)
            Ya = np.array(Ys)
            Xo = [np.column_stack([np.ones(N), Ya[t - 1]]) for t in range(1, T)]
            mo, _, _ = kalman_rw(Ya[1:], Xo, q * np.eye(2), 1.0, np.zeros(2), np.eye(2))
            base = np.array([0.5 * np.sum((Ya[t + 1] - Xo[t] @ mo[t]) ** 2)
                             for t in range(T - 13, T - 1)]).sum()
            def plugin_gain(Wm):
                Xl = [np.column_stack([np.ones(N), Wm @ Ya[t - 1], Ya[t - 1]]) for t in range(1, T)]
                mn, _, _ = kalman_rw(Ya[1:], Xl, q * np.eye(3), 1.0, np.zeros(3), np.eye(3))
                g = base
                for t in range(T - 13, T - 1):
                    g += -0.5 * np.sum((Ya[t + 1] - Xl[t] @ mn[t]) ** 2)
                return g
            obs = plugin_gain(W)
            trmean = Ya[:T - 12].mean(0)
            strata = np.digitize(trmean, np.quantile(trmean, [.2, .4, .6, .8]))
            vu, vs = [], []
            for _ in range(99):
                p = rng.permutation(N); vu.append(plugin_gain(W[np.ix_(p, p)]))
                p2 = np.arange(N)
                for sgrp in range(5):
                    idx = np.where(strata == sgrp)[0]; p2[idx] = rng.permutation(idx)
                vs.append(plugin_gain(W[np.ix_(p2, p2)]))
            if (1 + sum(v >= obs for v in vu)) / 100.0 <= 0.05: rej_u += 1
            if (1 + sum(v >= obs for v in vs)) / 100.0 <= 0.05: rej_s += 1
        out[scen] = dict(unrestricted=rej_u / M, stratified=rej_s / M)
        print(f"[S7 {scen}] unrestricted {rej_u/M:.3f} stratified {rej_s/M:.3f} ({time.time()-t0:.0f}s)", flush=True)
    return out
s7res = s7()
open(f"{TAB}/tab_s7.tex", "w").write(
    f"size ($\\beta_1=0$, graph-aligned effects) & {s7res['size']['unrestricted']:.3f} & {s7res['size']['stratified']:.3f} \\\\\n"
    f"power ($\\beta_1=0.25$) & {s7res['power']['unrestricted']:.3f} & {s7res['power']['stratified']:.3f} \\\\\n")

# ---------------- S1b table from stored results ----------------
try:
    ad = json.load(open(f"{_os.path.join(_HERE,'out')}/addendum_sim.json"))
    key = [k for k in ad if "1b" in k.lower() or "topol" in k.lower()]
    print("[S1b] keys:", list(ad.keys())[:12], "->", key, flush=True)
    if key:
        rows = ad[key[0]]
        frag = []
        if isinstance(rows, list):
            for r in rows:
                nm = r.get("family", r.get("name", "?"))
                frag.append(f"{nm} & {r.get('tau', r.get('tauW', float('nan'))):.1f} & "
                            f"{r.get('info', float('nan')):.1f} & {r.get('mse', float('nan')):.4f} \\\\")
        elif isinstance(rows, dict):
            for nm, r in rows.items():
                frag.append(f"{nm} & {r.get('tau', r.get('tauW', float('nan'))):.1f} & "
                            f"{r.get('info', float('nan')):.1f} & {r.get('mse', float('nan')):.4f} \\\\")
        if frag:
            open(f"{TAB}/tab_s1b.tex", "w").write("\n".join(frag) + "\n")
except Exception as e:
    print("[S1b] skipped:", e, flush=True)

json.dump(dict(s1d=dict(cov=c95, se=c95se, vratio=vr),
               s4x=[dict(label=l, cov=c, se=s, rmse=r) for l, c, s, r in s4rows],
               s7=s7res), open(f"{OUT}/sim_extra.json", "w"))
print("DONE", round(time.time() - t0), "s", flush=True)
