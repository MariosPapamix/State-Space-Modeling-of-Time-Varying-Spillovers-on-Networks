"""
Additional Chicago analyses:
A. Point forecasts that exist: predictive MEDIAN forecasts (all horizons) and
   frozen-coefficient (plug-in) Monte Carlo means (finite by Theorem 4(c)),
   alongside the sampled-coefficient means (population-finite only at h<=2).
B. Exact randomization test: 99 node-label permutations (+60 degree-preserving
   rewirings as a second reference set), statistic = test-window one-step
   plug-in log-score gain of the network model over the no-network model.
C. Negative-binomial observation layer at h=1: mode filter, r tuned on the
   training window, PIT and log-score comparison with the Poisson layer.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, time, os, sys, math
import numpy as np, pandas as pd, scipy.io as sio
from scipy.special import gammaln, logsumexp
from scipy.stats import chisquare
sys.path.insert(0, _HERE)
from nssm import (row_normalise, poisson_laplace_filter, poisson_logpmf, profile_information)

t0 = time.time()
OUT = _os.path.join(_HERE, "out"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20260823)

df = pd.read_csv(_os.path.join(_ROOT, "data", "crime.csv"), index_col=0)
Y = df.values.T.astype(int); T, N = Y.shape
A = sio.mmread(_os.path.join(_ROOT, "data", "neighborhood.mtx")).tocsr()
A = ((A + A.T) > 0).astype(float).toarray(); np.fill_diagonal(A, 0)
W = row_normalise(A); ell = np.log1p(Y)
TRAIN_END = 60; HORIZONS = [1, 2, 4, 6]; S_MAIN = 400
prev = json.load(open(f"{OUT}/chicago_results.json"))
tuned = {k: float(v) for k, v in prev["tuned_q"].items()}

def make_builder(Wmat, kind):
    def build(l_prev):
        if kind == "net":   return np.column_stack([np.ones(N), Wmat @ l_prev, l_prev])
        if kind == "nonet": return np.column_stack([np.ones(N), l_prev])
        raise ValueError(kind)
    return build

MODELS = {"net": 3, "nonet": 2}

def run_filter(builder, q, K, Ymat=Y):
    X_list = [None] + [builder(np.log1p(Ymat[t - 1])) for t in range(1, T)]
    fit = poisson_laplace_filter(Ymat[1:], X_list[1:], Q=q * np.eye(K), m0=np.zeros(K), P0=np.eye(K))
    return fit, X_list

# ------------------------------------------------------------------ A. forecasts
def predict_paths(m, P, Q, h, S, rng_l, builder, y_last, frozen=False):
    K = len(m)
    if frozen:
        th = np.repeat(m[None, :], S, 0)
    else:
        th = m[None, :] + rng_l.standard_normal((S, K)) @ np.linalg.cholesky(P + 1e-12 * np.eye(K)).T
    Lq = np.linalg.cholesky(Q + 1e-14 * np.eye(K))
    y_prev = np.repeat(y_last[None, :], S, 0).astype(float)
    for k in range(h):
        if k > 0 and not frozen:
            th = th + rng_l.standard_normal((S, K)) @ Lq.T
        lam = np.zeros_like(y_prev)
        for s in range(S):
            lam[s] = np.exp(np.clip(builder(np.log1p(y_prev[s])) @ th[s], -30, 30))
        if k < h - 1:
            y_prev = rng_l.poisson(np.minimum(lam, 1e12)).astype(float)
    return lam

def mixture_median(lam_paths):
    S, Nn = lam_paths.shape
    ymax = int(min(max(50, np.quantile(lam_paths, 0.999).max() * 3 + 20), 5000))
    grid = np.arange(0, ymax + 1)
    pmf = np.zeros((Nn, len(grid)))
    for s in range(S):
        pmf += np.exp(poisson_logpmf(grid[None, :], np.minimum(lam_paths[s], 1e6)[:, None]))
    cdf = np.cumsum(pmf / S, axis=1)
    return (cdf >= 0.5).argmax(axis=1).astype(float)

res = {}
for name, K in MODELS.items():
    b = make_builder(W, name); q = tuned[name]
    fit, _ = run_filter(b, q, K)
    rows = {h: [] for h in HORIZONS}
    rng_l = np.random.default_rng(101)
    for h in HORIZONS:
        for t_or in range(TRAIN_END - 1, T - h):
            idx = t_or - 1
            m, P = fit["m_filt"][idx], fit["P_filt"][idx]
            lam_full = predict_paths(m, P, q * np.eye(K), h, S_MAIN, rng_l, b, Y[t_or], frozen=False)
            lam_frozen = predict_paths(m, P, q * np.eye(K), h, S_MAIN, rng_l, b, Y[t_or], frozen=True)
            med = mixture_median(lam_full)
            ytrue = Y[t_or + h]
            rows[h].append(dict(
                mae_sampledmean=float(np.mean(np.abs(ytrue - lam_full.mean(0)))),
                mae_median=float(np.mean(np.abs(ytrue - med))),
                mae_frozenmean=float(np.mean(np.abs(ytrue - lam_frozen.mean(0)))),
            ))
    res[name] = rows
    print(f"[A] {name} done ({time.time()-t0:.0f}s)", flush=True)

pairsA = {}
for h in HORIZONS:
    d = {}
    for key in ["mae_sampledmean", "mae_median", "mae_frozenmean"]:
        a = np.array([r[key] for r in res["net"][h]]); bb = np.array([r[key] for r in res["nonet"][h]])
        d[key] = dict(net=float(a.mean()), nonet=float(bb.mean()), dmae=float((a - bb).mean()))
    pairsA[h] = d
    print(f"[A] h={h}: " + " | ".join(f"{k}: net={v['net']:.3f} nonet={v['nonet']:.3f} d={v['dmae']:+.4f}" for k, v in d.items()), flush=True)

# ------------------------------------------------------------------ B. randomization test
def plugin_logscore_test(builder, q, K):
    fit, X_list = run_filter(builder, q, K)
    ls = 0.0
    for t in range(TRAIN_END, T):
        eta = np.clip(X_list[t] @ fit["m_pred"][t - 1], -30, 30)
        ls += poisson_logpmf(Y[t], np.exp(eta)).sum()
    return float(ls)

ls_nonet = plugin_logscore_test(make_builder(W, "nonet"), tuned["nonet"], 2)
ls_obs = plugin_logscore_test(make_builder(W, "net"), tuned["net"], 3)
stat_obs = ls_obs - ls_nonet
r2 = np.random.default_rng(7)
perm_stats = []
for p in range(99):
    perm = r2.permutation(N)
    Wp = W[np.ix_(perm, perm)]
    perm_stats.append(plugin_logscore_test(make_builder(Wp, "net"), tuned["net"], 3) - ls_nonet)
    if p % 20 == 0: print(f"[B] perm {p} ({time.time()-t0:.0f}s)", flush=True)
perm_stats = np.array(perm_stats)
p_perm = (1 + int(np.sum(perm_stats >= stat_obs))) / (1 + len(perm_stats))

def rewire(Amat, n_swaps, r):
    A2 = Amat.copy(); edges = np.argwhere(np.triu(A2, 1) > 0)
    E = len(edges); done = 0; tries = 0
    while done < n_swaps and tries < 50 * n_swaps:
        tries += 1
        i1, i2 = r.choice(E, 2, replace=False)
        a_, b_ = edges[i1]; c_, d_ = edges[i2]
        if len({a_, b_, c_, d_}) < 4 or A2[a_, d_] or A2[c_, b_]: continue
        A2[a_, b_] = A2[b_, a_] = 0; A2[c_, d_] = A2[d_, c_] = 0
        A2[a_, d_] = A2[d_, a_] = 1; A2[c_, b_] = A2[b_, c_] = 1
        edges[i1] = [a_, d_]; edges[i2] = [c_, b_]; done += 1
    return A2

rw_stats = []
for p in range(60):
    Wr = row_normalise(rewire(A, int(A.sum()), r2))
    rw_stats.append(plugin_logscore_test(make_builder(Wr, "net"), tuned["net"], 3) - ls_nonet)
rw_stats = np.array(rw_stats)
p_rw = (1 + int(np.sum(rw_stats >= stat_obs))) / (1 + len(rw_stats))
print(f"[B] stat_obs={stat_obs:+.1f}; perms: max={perm_stats.max():+.1f} p={p_perm:.3f}; "
      f"rewire: max={rw_stats.max():+.1f} p={p_rw:.3f} ({time.time()-t0:.0f}s)", flush=True)

# ------------------------------------------------------------------ C. negative binomial at h=1
def nb_logpmf(y, mu, r):
    mu = np.maximum(mu, 1e-300)
    return (gammaln(y + r) - gammaln(r) - gammaln(y + 1.0)
            + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))

def nb_mode_filter(Ymat, X_list, Q, m0, P0, r, iters=30):
    Tt, Nn = Ymat.shape; K = len(m0)
    m_pred = np.zeros((Tt, K)); m_filt = np.zeros((Tt, K)); P_filt = np.zeros((Tt, K, K))
    m, P = np.asarray(m0, float).copy(), np.asarray(P0, float).copy()
    for t in range(Tt):
        mp, Pp = m, P + Q
        X = X_list[t]; y = Ymat[t]
        Ppi = np.linalg.inv(Pp); th = mp.copy()
        for it in range(iters):
            eta = np.clip(X @ th, -30, 30); mu = np.exp(eta)
            g = X.T @ (y - (y + r) * mu / (r + mu)) - Ppi @ (th - mp)
            wts = r * mu / (r + mu)                       # expected information weights
            H = X.T @ (X * wts[:, None]) + Ppi
            step = np.linalg.solve(H, g)
            s = 1.0
            def obj(tt):
                e = np.clip(X @ tt, -30, 30); m_ = np.exp(e)
                return float(np.sum(y * e - (y + r) * np.log(r + m_)) - 0.5 * (tt - mp) @ Ppi @ (tt - mp))
            f_old = obj(th); th_new = th + step
            while s > 1e-4:
                th_new = th + s * step
                if obj(th_new) >= f_old - 1e-10: break
                s *= 0.5
            th = th_new
            if np.max(np.abs(s * step)) < 1e-9: break
        eta = np.clip(X @ th, -30, 30); mu = np.exp(eta)
        wts = r * mu / (r + mu)
        Pf = np.linalg.inv(X.T @ (X * wts[:, None]) + Ppi); Pf = 0.5 * (Pf + Pf.T)
        m_pred[t], m_filt[t], P_filt[t] = mp, th, Pf
        m, P = th, Pf
    return dict(m_pred=m_pred, m_filt=m_filt, P_filt=P_filt)

nbres = {}
for name, K in MODELS.items():
    b = make_builder(W, name)
    X_list = [None] + [b(ell[t - 1]) for t in range(1, T)]
    best = None
    for r_ in [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 1e8]:
        fit = nb_mode_filter(Y[1:], X_list[1:], tuned[name] * np.eye(K), np.zeros(K), np.eye(K), r_)
        ls = sum(nb_logpmf(Y[t], np.exp(np.clip(X_list[t] @ fit["m_pred"][t - 1], -30, 30)), r_).sum()
                 for t in range(12, TRAIN_END))
        if best is None or ls > best[1]: best = (r_, ls, fit)
    r_star, _, fit = best
    rng_l = np.random.default_rng(11)
    pits, dls_rows, cov_rows = [], [], []
    for t_or in range(TRAIN_END - 1, T - 1):
        idx = t_or - 1
        m, P = fit["m_filt"][idx], fit["P_filt"][idx]
        th = m[None, :] + rng_l.standard_normal((S_MAIN, K)) @ np.linalg.cholesky(P + tuned[name] * np.eye(K)).T
        mu = np.exp(np.clip(th @ b(ell[t_or]).T, -30, 30))          # S x N
        y = Y[t_or + 1]
        lp = nb_logpmf(y[None, :], mu, r_star)
        ls_sum = float((logsumexp(lp, axis=0) - math.log(S_MAIN)).sum())
        ymax = int(min(max(50, np.quantile(mu, .999).max() * 4 + 30), 4000))
        grid = np.arange(ymax + 1)
        pmf = np.zeros((N, ymax + 1))
        for s in range(S_MAIN): pmf += np.exp(nb_logpmf(grid[None, :], mu[s][:, None], r_star))
        cdf = np.cumsum(pmf / S_MAIN, 1)
        iy = np.minimum(y, ymax)
        Fy = cdf[np.arange(N), iy]; Fm = np.where(y > 0, cdf[np.arange(N), np.maximum(iy - 1, 0)], 0.0)
        pits.append(Fm + rng_l.random(N) * (Fy - Fm))
        lo = (cdf >= 0.05).argmax(1); hi = (cdf >= 0.95).argmax(1)
        cov_rows.append(float(np.mean((y >= lo) & (y <= hi))))
        dls_rows.append(ls_sum)
    pits = np.concatenate(pits)
    cnt, _ = np.histogram(pits, bins=10, range=(0, 1))
    chi2s, pval = chisquare(cnt)
    nbres[name] = dict(r=r_star, ls_rows=dls_rows, pit_counts=cnt.tolist(),
                       chi2=float(chi2s), cov90=float(np.mean(cov_rows)))
    print(f"[C] NB {name}: r={r_star:g} chi2={chi2s:.1f} cov90={np.mean(cov_rows):.3f} ({time.time()-t0:.0f}s)", flush=True)
d_nb = float(np.mean(np.array(nbres["net"]["ls_rows"]) - np.array(nbres["nonet"]["ls_rows"])))
print(f"[C] NB net-vs-nonet dLS h=1: {d_nb:+.2f}", flush=True)

json.dump(dict(pointforecasts={str(h): v for h, v in pairsA.items()},
               randtest=dict(stat_obs=stat_obs, perm=perm_stats.tolist(), p_perm=p_perm,
                             rewire=rw_stats.tolist(), p_rw=p_rw, ls_nonet=ls_nonet),
               nb=nbres, nb_dls_h1=d_nb, runtime=time.time() - t0),
          open(f"{OUT}/addendum_chicago.json", "w"))
print("DONE", round(time.time() - t0), "s", flush=True)
