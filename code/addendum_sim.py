"""
Simulation addendum:
S4  Validation of the Poisson Laplace filter under the true model (coverage/RMSE).
S1b Realised information vs tau_W across graph families (collapse check).
S1c Spatially smooth lagged fields destroy information on a fixed sparse graph.
S5  Log-scale S-scan under a persistent Gaussian coefficient draw: divergence at
    h>=3 visible when the posterior is diffuse, invisible when concentrated.
S6  Cap sensitivity of the raw-count S-scan at h=3.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, time, os, sys, math
import numpy as np
sys.path.insert(0, _HERE)
from nssm import (row_normalise, er_graph, complete_graph, tau_W, kalman_rw,
                  poisson_laplace_filter, profile_information)
t0 = time.time(); OUT = _os.path.join(_HERE, "out"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(5)

# ---------------- S4: Laplace filter validation under the true Poisson NSSM ----
N, T = 200, 60
A = er_graph(N, 4, rng); W = row_normalise(A)
Q = np.diag([1e-4, 2.5e-4, 2.5e-4]); Lq = np.linalg.cholesky(Q)
R = 200
cov_hits, err2, infos = [], [], []
n_skip = 0
for rep in range(R):
    th = np.array([0.3, 0.25, 0.25]) + 0.05 * rng.standard_normal(3)
    y = rng.poisson(2.0, N).astype(float)
    ths, Xs, Ys = [], [], []
    for t in range(T):
        th = th + rng.standard_normal(3) @ Lq.T
        l = np.log1p(y)
        X = np.column_stack([np.ones(N), W @ l, l])
        lam = np.exp(np.clip(X @ th, -25, 25))
        y = rng.poisson(np.minimum(lam, 1e9)).astype(float)
        ths.append(th.copy()); Xs.append(X); Ys.append(y)
    try:
        fit = poisson_laplace_filter(np.array(Ys), Xs, Q, m0=np.zeros(3), P0=np.eye(3))
    except np.linalg.LinAlgError:
        n_skip += 1; continue
    B = np.array(ths)
    Pd = fit["P_filt"][:, 1, 1]
    if not np.all(np.isfinite(Pd)) or np.any(Pd <= 0):
        n_skip += 1; continue
    e1 = fit["m_filt"][:, 1] - B[:, 1]; sd1 = np.sqrt(Pd)
    keep = np.arange(T) >= 10
    cov_hits.append(np.mean((np.abs(e1) <= 1.959963985 * sd1)[keep]))
    err2.append(np.mean(e1[keep] ** 2))
    Xk = Xs
    infos.append(float(np.mean([profile_information(Xk[t], np.exp(np.clip(Xk[t] @ fit["m_filt"][t], -25, 25)), 1)[0]
                                for t in range(10, T, 10)])))
S4 = dict(coverage=float(np.mean(cov_hits)), cov_se=float(np.std(cov_hits) / math.sqrt(len(cov_hits))),
          rmse=float(math.sqrt(np.mean(err2))), mean_info=float(np.mean(infos)), n_skip=n_skip, n_used=len(cov_hits))
print(f"[S4] Laplace coverage={S4['coverage']:.3f}+-{S4['cov_se']:.3f} rmse={S4['rmse']:.4f} "
      f"info={S4['mean_info']:.0f} ({time.time()-t0:.0f}s)", flush=True)

# ---------------- S1b: graph families, realised info vs tau_W -------------------
def ring_lattice(N, k=2):
    A = np.zeros((N, N))
    for i in range(N):
        for d in range(1, k + 1):
            A[i, (i + d) % N] = A[i, (i - d) % N] = 1
    return A
def watts_strogatz(N, k, p, r):
    A = ring_lattice(N, k)
    edges = np.argwhere(np.triu(A, 1) > 0)
    for e in edges:
        if r.random() < p:
            i = e[0]
            cand = [j for j in range(N) if j != i and A[i, j] == 0]
            if cand:
                j = cand[r.integers(len(cand))]
                A[e[0], e[1]] = A[e[1], e[0]] = 0
                A[i, j] = A[j, i] = 1
    iso = np.where(A.sum(1) == 0)[0]
    for i in iso:
        j = (i + 1) % N; A[i, j] = A[j, i] = 1
    return A
def barabasi_albert(N, m, r):
    A = np.zeros((N, N))
    for i in range(m + 1):
        for j in range(i + 1, m + 1): A[i, j] = A[j, i] = 1
    deg = A.sum(1)
    for v in range(m + 1, N):
        probs = deg[:v] / deg[:v].sum()
        tgt = r.choice(v, size=m, replace=False, p=probs)
        for u in tgt: A[v, u] = A[u, v] = 1
        deg = A.sum(1)
    return A

def sim_gauss_info(W, Tt=60, R=50, sigma=1.0, seed=1):
    r = np.random.default_rng(seed)
    Nn = W.shape[0]; Q = np.diag([1e-4, 1e-3, 1e-3])
    infos, mses = [], []
    for rep in range(R):
        b0 = np.cumsum(r.normal(0, 0.02, Tt)); b1 = np.full(Tt, 0.30); b2 = np.full(Tt, 0.30)
        Yl = np.zeros((Tt + 1, Nn)); Yl[0] = r.normal(0, 1, Nn)
        for t in range(1, Tt + 1):
            Yl[t] = b0[t-1] + b1[t-1] * (W @ Yl[t-1]) + b2[t-1] * Yl[t-1] + r.normal(0, sigma, Nn)
        X_list = [np.column_stack([np.ones(Nn), W @ Yl[t-1], Yl[t-1]]) for t in range(1, Tt + 1)]
        fit = kalman_rw(Yl[1:], X_list, sigma**2, Q, np.zeros(3), np.eye(3))
        e1 = fit["m_filt"][:, 1] - b1
        keep = np.arange(Tt) >= 10
        mses.append(np.mean(e1[keep] ** 2))
        infos.append(np.mean([profile_information(X_list[t], np.full(Nn, 1/sigma**2), 1)[0] for t in range(10, Tt)]))
    return float(np.mean(infos)), float(np.mean(mses))

Ng = 200; fams = {}
fams["ring lattice ($d=4$)"] = ring_lattice(Ng, 2)
fams["Watts--Strogatz"] = watts_strogatz(Ng, 2, 0.1, rng)
fams["Barab\\'asi--Albert"] = barabasi_albert(Ng, 2, rng)
fams["Erd\\H os--R\\'enyi ($d=4$)"] = er_graph(Ng, 4, rng)
fams["dense ER ($d=N/4$)"] = er_graph(Ng, Ng / 4, rng)
fams["complete"] = complete_graph(Ng)
S1b = []
for nm, Ad in fams.items():
    Wf = row_normalise(Ad); tw = tau_W(Wf)
    info, mse = sim_gauss_info(Wf, R=50, seed=3)
    S1b.append(dict(family=nm, tau=tw, info=info, mse=mse))
    print(f"[S1b] {nm:28s} tau={tw:8.2f} info={info:8.2f} ratio={info/max(tw,1e-12):6.2f} mse={mse:.5f}", flush=True)

# ---------------- S1c: smooth fields on a fixed sparse graph --------------------
Wf = row_normalise(er_graph(Ng, 4, np.random.default_rng(9))); tw = tau_W(Wf)
S1c = []
for phi in [0.0, 0.5, 0.9, 0.99]:
    Minv = np.linalg.inv(np.eye(Ng) - phi * Wf)
    r = np.random.default_rng(13); infos = []
    for rep in range(40):
        x = Minv @ r.normal(0, 1, Ng)
        X = np.column_stack([np.ones(Ng), Wf @ x, x])
        infos.append(profile_information(X, col=1)[0])
    S1c.append(dict(phi=phi, info=float(np.mean(infos))))
    print(f"[S1c] phi={phi}: realised info={np.mean(infos):8.2f} (tau={tw:.2f})", flush=True)

# ---------------- S5: log-scale S-scan under persistent Gaussian coefficient ----
N5 = 10; A5 = er_graph(N5, 3, np.random.default_rng(3)); W5 = row_normalise(A5)
y0 = np.full(N5, 3.0); mtheta = np.array([0.3, 0.3, 0.4])
S5 = {}
for sdev in [0.3, 1.0]:
    r = np.random.default_rng(17)
    res = {}
    for S in [10**3, 10**4, 10**5, 10**6]:
        th = mtheta[None, :] + sdev * r.standard_normal((S, 3))
        y = np.repeat(y0[None, :], S, 0)
        means = {}
        for h in range(1, 5):
            l = np.log1p(y)
            eta = np.clip(th[:, 0:1] + th[:, 1:2] * (l @ W5.T) + th[:, 2:3] * l, -700, 700)
            lam = np.exp(eta)
            means[h] = float(lam.mean())
            if h < 4: y = r.poisson(np.minimum(lam, 1e12)).astype(float)
        res[S] = means
        print(f"[S5] sd={sdev} S={S:7d} " + " ".join(f"h{h}:{means[h]:.4g}" for h in means), flush=True)
    S5[str(sdev)] = {str(k): v for k, v in res.items()}

# ---------------- S6: cap sensitivity of the raw-count S-scan at h=3 ------------
b0, b1, b2 = 0.3, 0.3, 0.4
S6 = {}
for clip in [200.0, 700.0]:
    r = np.random.default_rng(23); res = {}
    for S in [10**4, 10**5]:
        y = np.repeat(y0[None, :], S, 0)
        means = {}
        for h in range(1, 4):
            eta = np.clip(b0 + b1 * (y @ W5.T) + b2 * y, -clip, clip)
            lam = np.exp(eta)
            means[h] = float(lam.mean())
            if h < 3: y = r.poisson(np.minimum(lam, 1e12)).astype(float)
        res[S] = means
        print(f"[S6] clip={clip} S={S}: " + " ".join(f"h{h}:{means[h]:.4g}" for h in means), flush=True)
    S6[str(clip)] = {str(k): v for k, v in res.items()}

json.dump(dict(S4=S4, S1b=S1b, S1c=S1c, S5=S5, S6=S6, runtime=time.time() - t0),
          open(f"{OUT}/addendum_sim.json", "w"))
print("DONE", round(time.time() - t0), "s", flush=True)
