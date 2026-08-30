"""
Benchmarks from the literature on the Chicago data, one protocol for all models:
expanding-window refits at each origin t0 in {59..70}, one-step (and two-step)
forecasts of month t0+h, log score = sum over the 552 nodes of the log predictive
pmf, MAE over nodes of the model's point forecast.

Models:
  NSSM net / nonet (Poisson, full mixture predictive, S=400)     [ours]
  NSSM net (NB layer, r tuned)                                    [ours, from addendum]
  NSSM net + seasonal covariates (NB layer)                       [ours, new]
  PNAR(1) linear   (Armillotta & Fokianos 2024)                   [benchmark]
  PNAR(1) log-linear                                              [benchmark]
  GNAR-type log-linear, two neighbourhood orders                  [benchmark]
  hhh4-type NB endemic-epidemic: unit intercepts + seasonality    [benchmark]
  kernel tvNAR (one-sided local-constant Poisson)                 [benchmark]
  EWMA-NB per node                                                [benchmark]

Self-tests: PNAR-linear MLE and the hhh4 coordinate ascent are checked on
simulated data before use.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, time, os, sys, math
import numpy as np, pandas as pd, scipy.io as sio
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp
sys.path.insert(0, _HERE)
from nssm import row_normalise, poisson_laplace_filter, poisson_logpmf

t0c = time.time()
OUT = _os.path.join(_HERE, "out"); TAB = _os.path.join(_ROOT, "paper", "tables")
rng = np.random.default_rng(31)

df = pd.read_csv(_os.path.join(_ROOT, "data", "crime.csv"), index_col=0)
Y = df.values.T.astype(int); T, N = Y.shape
A = sio.mmread(_os.path.join(_ROOT, "data", "neighborhood.mtx")).tocsr()
A = ((A + A.T) > 0).astype(float).toarray(); np.fill_diagonal(A, 0)
W = row_normalise(A); ell = np.log1p(Y)
TRAIN_END = 60; ORIG = list(range(TRAIN_END - 1, T - 1))          # 12 origins, h=1
prev = json.load(open(f"{OUT}/chicago_results.json"))
q_net, q_non = float(prev["tuned_q"]["net"]), float(prev["tuned_q"]["nonet"])
month = np.arange(T) % 12
SIN = np.sin(2 * np.pi * month / 12); COS = np.cos(2 * np.pi * month / 12)

from nb_filter import nb_mode_filter, nb_logpmf

# ---------------------------------------------------------------- self-tests
def pnar_lin_fit(Ys, Wm, t_hi):
    """MLE of linear PNAR(1) on data up to t_hi (exclusive)."""
    Xr = np.stack([np.ones((t_hi - 1, Wm.shape[0])), (Ys[:t_hi - 1] @ Wm.T), Ys[:t_hi - 1]], -1)
    yv = Ys[1:t_hi]
    def nll(b):
        lam = np.maximum(Xr @ b, 1e-8)
        return -(np.sum(yv * np.log(lam) - lam))
    res = minimize(nll, np.array([0.5, 0.2, 0.3]), method="L-BFGS-B",
                   bounds=[(1e-8, None)] * 3)
    return res.x

from nssm import er_graph
Wsim = row_normalise(er_graph(60, 4, rng))
y = rng.poisson(2.0, 60).astype(float); rows = [y]
for t in range(299):
    lam = 0.6 + 0.25 * (Wsim @ rows[-1]) + 0.35 * rows[-1]
    rows.append(rng.poisson(lam).astype(float))
Ysim = np.array(rows)
bhat = pnar_lin_fit(Ysim, Wsim, 300)
print(f"[selftest] PNAR-lin MLE on sim: {np.round(bhat,3)} vs (0.6,0.25,0.35)", flush=True)
assert np.max(np.abs(bhat - np.array([0.6, 0.25, 0.35]))) < 0.06

def hhh4_fit(Ys, Wm, t_hi, sin_v, cos_v, outer=6):
    """NB endemic-epidemic: mu_it = exp(a_i + g*sin_t + d*cos_t) + lam*y_{i,t-1} + phi*(W y_{t-1})_i."""
    Nn = Ys.shape[1]
    ylag = Ys[:t_hi - 1]; ycur = Ys[1:t_hi]
    wy = ylag @ Wm.T
    a = np.log(ycur.mean(0) + 0.05)
    lam, phi, g, d, logr = 0.1, 0.1, 0.0, 0.0, math.log(2.0)
    s_v, c_v = sin_v[1:t_hi], cos_v[1:t_hi]
    def mu_of(a_, lam_, phi_, g_, d_):
        endm = np.exp(np.clip(a_[None, :] + g_ * s_v[:, None] + d_ * c_v[:, None], -30, 30))
        return endm + lam_ * ylag + phi_ * wy, endm
    def ll(a_, lam_, phi_, g_, d_, logr_):
        mu, _ = mu_of(a_, lam_, phi_, g_, d_)
        return float(nb_logpmf(ycur, mu, math.exp(logr_)).sum())
    for it in range(3):                                       # warm start: coordinate ascent
        r_ = math.exp(logr)
        for _ in range(3):
            mu, endm = mu_of(a, lam, phi, g, d)
            w1 = ycur / mu - (ycur + r_) / (mu + r_)
            w2 = -ycur / mu**2 + (ycur + r_) / (mu + r_)**2
            gvec = (w1 * endm).sum(0)
            hvec = (w2 * endm**2 + w1 * endm).sum(0)
            step = np.where(hvec < -1e-12, -gvec / hvec, 0.0)
            a = a + np.clip(step, -1.0, 1.0)
        def nll5(p):
            return -ll(a, p[0], p[1], p[2], p[3], p[4])
        res = minimize(nll5, np.array([lam, phi, g, d, logr]), method="L-BFGS-B",
                       bounds=[(0, None), (0, None), (-3, 3), (-3, 3), (math.log(0.2), math.log(500))])
        lam, phi, g, d, logr = res.x
    # joint refinement over all N+5 parameters with analytic gradient
    from scipy.special import digamma
    def pack(a_, lam_, phi_, g_, d_, logr_): return np.concatenate([a_, [lam_, phi_, g_, d_, logr_]])
    def nll_grad(p):
        a_ = p[:Nn]; lam_, phi_, g_, d_, logr_ = p[Nn:]
        r_ = math.exp(logr_)
        mu, endm = mu_of(a_, lam_, phi_, g_, d_)
        llv = nb_logpmf(ycur, mu, r_).sum()
        w1 = ycur / mu - (ycur + r_) / (mu + r_)
        ga = (w1 * endm).sum(0)
        gl = float((w1 * ylag).sum()); gp = float((w1 * wy).sum())
        gg = float((w1 * endm * s_v[:, None]).sum()); gd = float((w1 * endm * c_v[:, None]).sum())
        glr = float((r_ * (digamma(ycur + r_) - digamma(r_) + np.log(r_ / (r_ + mu))
                           + 1.0 - (ycur + r_) / (r_ + mu))).sum())
        return -llv, -np.concatenate([ga, [gl, gp, gg, gd, glr]])
    bounds = [(None, None)] * Nn + [(0, None), (0, None), (-3, 3), (-3, 3), (math.log(0.2), math.log(500))]
    res = minimize(nll_grad, pack(a, lam, phi, g, d, logr), jac=True, method="L-BFGS-B",
                   bounds=bounds, options=dict(maxiter=400))
    a = res.x[:Nn]; lam, phi, g, d, logr = res.x[Nn:]
    return dict(a=a, lam=lam, phi=phi, g=g, d=d, r=math.exp(logr), ll=-res.fun)

y0 = rng.poisson(3.0, 60).astype(float); rows = [y0]
a_true = np.log(rng.uniform(0.5, 2.0, 60))
for t in range(299):
    endm = np.exp(a_true + 0.4 * SIN[(t + 1) % 12 if False else 0] * 0)  # constant endemic in sim
    lam = endm + 0.3 * rows[-1] + 0.2 * (Wsim @ rows[-1])
    rows.append(rng.poisson(lam).astype(float))
Ysim2 = np.array(rows)
fit = hhh4_fit(Ysim2, Wsim, 300, np.zeros(300), np.zeros(300))
print(f"[selftest] hhh4 on sim: lam={fit['lam']:.3f} (0.3) phi={fit['phi']:.3f} (0.2) "
      f"r={fit['r']:.0f} (Poisson-limit)", flush=True)
assert abs(fit["lam"] - 0.3) < 0.05 and abs(fit["phi"] - 0.2) < 0.06

# ---------------------------------------------------------------- shared scoring
def score_plugin(lam_or_mu, ytrue, kind="pois", r=None):
    if kind == "pois":
        ls = float(poisson_logpmf(ytrue, np.maximum(lam_or_mu, 1e-10)).sum())
    else:
        ls = float(nb_logpmf(ytrue, lam_or_mu, r).sum())
    return ls, float(np.mean(np.abs(ytrue - lam_or_mu)))

def mc_score_fixedtheta(build_lam1, step_lam, ytrue, S=400, kind="pois", r=None, seed=0):
    """Two-step predictive under fixed parameters: simulate y_{t+1}, score mixture at t+2."""
    rl = np.random.default_rng(seed)
    lam1 = build_lam1()
    y1 = (rl.poisson(np.tile(np.minimum(lam1, 1e9), (S, 1))) if kind == "pois"
          else rl.negative_binomial(r, r / (r + np.tile(lam1, (S, 1)))))
    lam2 = np.stack([step_lam(y1[s].astype(float)) for s in range(S)])
    lp = (poisson_logpmf(ytrue[None, :], np.maximum(lam2, 1e-10)) if kind == "pois"
          else nb_logpmf(ytrue[None, :], lam2, r))
    ls = float((logsumexp(lp, axis=0) - math.log(S)).sum())
    return ls, float(np.mean(np.abs(ytrue - lam2.mean(0))))

RES = {}  # model -> {"h1_ls": [...], "h1_mae": [...], "h2_ls": [...]}
def add(model, h, ls, mae=None):
    RES.setdefault(model, {}).setdefault(f"h{h}_ls", []).append(ls)
    if mae is not None: RES[model].setdefault(f"h{h}_mae", []).append(mae)

# ---------------------------------------------------------------- our NSSM (Poisson mixture)
OFF = np.log(Y[:TRAIN_END].mean(0) + 0.05)          # unit offsets from the training window only

def nssm_eval(kind, q, seasonal=False, nb_r=None, label=None, horizons=(1, 2), S=400, offsets=False):
    K = (3 if kind == "net" else 2) + (2 if seasonal else 0) + (1 if offsets else 0)
    def build(l_prev, t):
        base = [np.ones(N)]
        if seasonal: base += [np.full(N, SIN[t]), np.full(N, COS[t])]
        if offsets: base += [OFF]
        if kind == "net": base += [W @ l_prev, l_prev]
        else: base += [l_prev]
        return np.column_stack(base)
    X_list = [None] + [build(np.log1p(Y[t - 1]), t) for t in range(1, T)]
    if nb_r is None:
        fit = poisson_laplace_filter(Y[1:], X_list[1:], q * np.eye(K), np.zeros(K), np.eye(K))
    else:
        sys.path.insert(0, _HERE)
        from nb_filter import nb_mode_filter  # noqa
        fit = nb_mode_filter(Y[1:], X_list[1:], q * np.eye(K), np.zeros(K), np.eye(K), nb_r)
    rl = np.random.default_rng(5)
    for h in horizons:
        for t_or in range(TRAIN_END - 1, T - h):
            idx = t_or - 1
            m, P = fit["m_filt"][idx], fit["P_filt"][idx]
            th = m[None, :] + rl.standard_normal((S, K)) @ np.linalg.cholesky(P + q * np.eye(K) + 1e-12 * np.eye(K)).T
            Lq = np.linalg.cholesky(q * np.eye(K) + 1e-14 * np.eye(K))
            y_prev = np.repeat(Y[t_or][None, :].astype(float), S, 0)
            for k in range(h):
                if k > 0: th = th + rl.standard_normal((S, K)) @ Lq.T
                lam = np.zeros_like(y_prev)
                for s in range(S):
                    lam[s] = np.exp(np.clip(build(np.log1p(y_prev[s]), t_or + k + 1) @ th[s], -30, 30))
                if k < h - 1:
                    y_prev = (rl.poisson(np.minimum(lam, 1e12)).astype(float) if nb_r is None
                              else rl.negative_binomial(nb_r, nb_r / (nb_r + np.minimum(lam, 1e12))).astype(float))
            ytrue = Y[t_or + h]
            lp = (poisson_logpmf(ytrue[None, :], np.minimum(lam, 1e8)) if nb_r is None
                  else nb_logpmf(ytrue[None, :], np.minimum(lam, 1e8), nb_r))
            ls = float((logsumexp(lp, axis=0) - math.log(S)).sum())
            add(label, h, ls, float(np.mean(np.abs(ytrue - np.median(lam, 0)))))
    print(f"[{label}] done ({time.time()-t0c:.0f}s)", flush=True)

nssm_eval("net", q_net, label="NSSM net (Poisson)")
nssm_eval("nonet", q_non, label="NSSM no-network (Poisson)")

# NB rows from the addendum (same origins, h=1, absolute LS sums)
addc = json.load(open(f"{OUT}/addendum_chicago.json"))
RES["NSSM net (NB)"] = {"h1_ls": addc["nb"]["net"]["ls_rows"]}
RES["NSSM no-network (NB)"] = {"h1_ls": addc["nb"]["nonet"]["ls_rows"]}

# NSSM net + seasonal, NB layer: small (q, r) tuning on the training window
def train_plugin_nb(kind, q, r, seasonal, offsets=False):
    K = (3 if kind == "net" else 2) + (2 if seasonal else 0) + (1 if offsets else 0)
    def build(l_prev, t):
        base = [np.ones(N)]
        if seasonal: base += [np.full(N, SIN[t]), np.full(N, COS[t])]
        if offsets: base += [OFF]
        if kind == "net": base += [W @ l_prev, l_prev]
        else: base += [l_prev]
        return np.column_stack(base)
    X_list = [None] + [build(np.log1p(Y[t - 1]), t) for t in range(1, T)]
    from nb_filter import nb_mode_filter
    fit = nb_mode_filter(Y[1:TRAIN_END], X_list[1:TRAIN_END], q * np.eye(K), np.zeros(K), np.eye(K), r)
    ls = sum(nb_logpmf(Y[t], np.exp(np.clip(X_list[t] @ fit["m_pred"][t - 1], -30, 30)), r).sum()
             for t in range(12, TRAIN_END))
    return float(ls)

best = None
for q in [1e-4, 3e-4, 1e-3]:
    for r in [2.0, 5.0, 10.0]:
        v = train_plugin_nb("net", q, r, True)
        if best is None or v > best[0]: best = (v, q, r)
print(f"[seasonal NSSM] tuned q={best[1]:g} r={best[2]:g} ({time.time()-t0c:.0f}s)", flush=True)
nssm_eval("net", best[1], seasonal=True, nb_r=best[2], label="NSSM net + seasonal (NB)", horizons=(1,))

for kind, lab in [("net", "NSSM net + offsets, seasonal (NB)"),
                  ("nonet", "NSSM no-network + offsets, seasonal (NB)")]:
    bo = None
    for q in [1e-4, 3e-4, 1e-3]:
        for r in [2.0, 5.0, 10.0]:
            v = train_plugin_nb(kind, q, r, True, offsets=True)
            if bo is None or v > bo[0]: bo = (v, q, r)
    print(f"[offset NSSM {kind}] tuned q={bo[1]:g} r={bo[2]:g} ({time.time()-t0c:.0f}s)", flush=True)
    nssm_eval(kind, bo[1], seasonal=True, nb_r=bo[2], label=lab, horizons=(1,), offsets=True)

# ---------------------------------------------------------------- PNAR linear / log-linear / GNAR2
def glm_pois_fit(Xr, yv, b0):
    def nll(b):
        eta = np.clip(Xr @ b, -30, 30); lam = np.exp(eta)
        return float(np.sum(lam - yv * eta))
    return minimize(nll, b0, method="L-BFGS-B").x

A2 = ((A @ A) > 0).astype(float); np.fill_diagonal(A2, 0); A2 = np.maximum(A2 - A, 0)
rs = A2.sum(1); W2 = A2 / np.where(rs > 0, rs, 1.0)[:, None]

for t_or in ORIG:
    yl = Y[:t_or].astype(float)
    # PNAR linear
    b = pnar_lin_fit(Y.astype(float), W, t_or + 1)
    lam1 = np.maximum(b[0] + b[1] * (W @ Y[t_or]) + b[2] * Y[t_or], 1e-8)
    ls, mae = score_plugin(lam1, Y[t_or + 1]); add("PNAR(1) linear", 1, ls, mae)
    if t_or <= T - 3:
        ls2, _ = mc_score_fixedtheta(
            lambda: lam1, lambda y1: np.maximum(b[0] + b[1] * (W @ y1) + b[2] * y1, 1e-8),
            Y[t_or + 2], seed=t_or); add("PNAR(1) linear", 2, ls2)
    # PNAR log-linear
    Xr = np.column_stack([np.ones(((t_or) * N)),
                          (np.log1p(Y[:t_or]) @ W.T).ravel(), np.log1p(Y[:t_or]).ravel()])
    bl = glm_pois_fit(Xr, Y[1:t_or + 1].ravel(), np.array([-0.3, 0.4, 0.5]))
    lam1 = np.exp(np.clip(bl[0] + bl[1] * (W @ np.log1p(Y[t_or])) + bl[2] * np.log1p(Y[t_or]), -30, 30))
    ls, mae = score_plugin(lam1, Y[t_or + 1]); add("PNAR(1) log-linear", 1, ls, mae)
    if t_or <= T - 3:
        ls2, _ = mc_score_fixedtheta(
            lambda: lam1,
            lambda y1: np.exp(np.clip(bl[0] + bl[1] * (W @ np.log1p(y1)) + bl[2] * np.log1p(y1), -30, 30)),
            Y[t_or + 2], seed=t_or); add("PNAR(1) log-linear", 2, ls2)
    # GNAR-type two orders
    Xr = np.column_stack([np.ones(((t_or) * N)), (np.log1p(Y[:t_or]) @ W.T).ravel(),
                          (np.log1p(Y[:t_or]) @ W2.T).ravel(), np.log1p(Y[:t_or]).ravel()])
    bg = glm_pois_fit(Xr, Y[1:t_or + 1].ravel(), np.array([-0.3, 0.3, 0.1, 0.5]))
    lam1 = np.exp(np.clip(bg[0] + bg[1] * (W @ np.log1p(Y[t_or])) + bg[2] * (W2 @ np.log1p(Y[t_or]))
                          + bg[3] * np.log1p(Y[t_or]), -30, 30))
    ls, mae = score_plugin(lam1, Y[t_or + 1]); add("GNAR-type, two orders", 1, ls, mae)
print(f"[PNAR/GNAR] done ({time.time()-t0c:.0f}s)", flush=True)

# ---------------------------------------------------------------- hhh4-type
for t_or in ORIG:
    f = hhh4_fit(Y.astype(float), W, t_or + 1, SIN, COS)
    endm = np.exp(np.clip(f["a"] + f["g"] * SIN[t_or + 1] + f["d"] * COS[t_or + 1], -30, 30))
    mu1 = endm + f["lam"] * Y[t_or] + f["phi"] * (W @ Y[t_or])
    ls, mae = score_plugin(mu1, Y[t_or + 1], kind="nb", r=f["r"]); add("hhh4-type NB", 1, ls, mae)
    if t_or <= T - 3:
        endm2 = np.exp(np.clip(f["a"] + f["g"] * SIN[t_or + 2] + f["d"] * COS[t_or + 2], -30, 30))
        ls2, _ = mc_score_fixedtheta(
            lambda: mu1, lambda y1: endm2 + f["lam"] * y1 + f["phi"] * (W @ y1),
            Y[t_or + 2], kind="nb", r=f["r"], seed=t_or); add("hhh4-type NB", 2, ls2)
RES["hhh4-type NB"]["fitted"] = dict(lam=float(f["lam"]), phi=float(f["phi"]), r=float(f["r"]))
print(f"[hhh4] done, last fit lam={f['lam']:.3f} phi={f['phi']:.3f} r={f['r']:.1f} ({time.time()-t0c:.0f}s)", flush=True)

# ---------------------------------------------------------------- kernel tvNAR
def ktv_fit(t_hi, bw):
    wts = np.exp(-0.5 * ((t_hi - 1 - np.arange(1, t_hi)) / bw) ** 2)
    Xr = np.column_stack([np.ones(((t_hi - 1) * N)),
                          (np.log1p(Y[:t_hi - 1]) @ W.T).ravel(), np.log1p(Y[:t_hi - 1]).ravel()])
    yv = Y[1:t_hi].ravel(); wv = np.repeat(wts, N)
    def nll(b):
        eta = np.clip(Xr @ b, -30, 30)
        return float(np.sum(wv * (np.exp(eta) - yv * eta)))
    return minimize(nll, np.array([-0.3, 0.4, 0.5]), method="L-BFGS-B").x

bw_best, sc_best = None, -np.inf
for bw in [4.0, 8.0, 12.0, 16.0]:
    sc = 0.0
    for s in range(30, TRAIN_END - 1):
        bb = ktv_fit(s + 1, bw)
        lam = np.exp(np.clip(bb[0] + bb[1] * (W @ np.log1p(Y[s])) + bb[2] * np.log1p(Y[s]), -30, 30))
        sc += poisson_logpmf(Y[s + 1], lam).sum()
    if sc > sc_best: bw_best, sc_best = bw, sc
print(f"[ktv] bandwidth {bw_best} months ({time.time()-t0c:.0f}s)", flush=True)
for t_or in ORIG:
    bb = ktv_fit(t_or + 1, bw_best)
    lam1 = np.exp(np.clip(bb[0] + bb[1] * (W @ np.log1p(Y[t_or])) + bb[2] * np.log1p(Y[t_or]), -30, 30))
    ls, mae = score_plugin(lam1, Y[t_or + 1]); add("kernel tvNAR", 1, ls, mae)
    if t_or <= T - 3:
        ls2, _ = mc_score_fixedtheta(
            lambda: lam1,
            lambda y1: np.exp(np.clip(bb[0] + bb[1] * (W @ np.log1p(y1)) + bb[2] * np.log1p(y1), -30, 30)),
            Y[t_or + 2], seed=t_or); add("kernel tvNAR", 2, ls2)

# ---------------------------------------------------------------- EWMA-NB
best = None
for aa in [0.1, 0.2, 0.3, 0.5, 0.7]:
    for r in [1.0, 2.0, 5.0, 10.0, 20.0]:
        m = Y[0].astype(float) + 0.1; sc = 0.0
        for t in range(1, TRAIN_END):
            if t >= 12: sc += nb_logpmf(Y[t], np.maximum(m, 1e-3), r).sum()
            m = aa * Y[t] + (1 - aa) * m
        if best is None or sc > best[0]: best = (sc, aa, r)
_, aa, r_ew = best
for t_or in ORIG:
    mm = Y[0].astype(float) + 0.1
    for t in range(1, t_or + 1): mm = aa * Y[t] + (1 - aa) * mm
    ls, mae = score_plugin(np.maximum(mm, 1e-3), Y[t_or + 1], kind="nb", r=r_ew)
    add("EWMA-NB", 1, ls, mae)
print(f"[EWMA] a={aa} r={r_ew} ({time.time()-t0c:.0f}s)", flush=True)

# ---------------------------------------------------------------- table + stats
REF1, REF2 = "NSSM net + offsets, seasonal (NB)", "NSSM net (Poisson)"
def block_boot_ci(diff, B=4000, L=3, seed=2):
    rl = np.random.default_rng(seed); n = len(diff); out = []
    for _ in range(B):
        idx = []
        while len(idx) < n:
            s0 = rl.integers(0, n); idx += [(s0 + k) % n for k in range(L)]
        out.append(np.mean(np.array(diff)[idx[:n]]))
    return float(np.quantile(out, .025)), float(np.quantile(out, .975))

order = ["NSSM net + offsets, seasonal (NB)", "NSSM no-network + offsets, seasonal (NB)",
         "NSSM net + seasonal (NB)", "NSSM net (NB)", "NSSM no-network (NB)",
         "NSSM net (Poisson)", "NSSM no-network (Poisson)", "hhh4-type NB", "EWMA-NB",
         "kernel tvNAR", "GNAR-type, two orders", "PNAR(1) linear", "PNAR(1) log-linear"]
rows = []
for mdl in order:
    d = RES[mdl]
    ls1 = np.mean(d["h1_ls"]); mae1 = np.mean(d.get("h1_mae", [np.nan]))
    if mdl == REF1:
        dtxt = "---"
    else:
        n = min(len(d["h1_ls"]), len(RES[REF1]["h1_ls"]))
        diff = np.array(RES[REF1]["h1_ls"][:n]) - np.array(d["h1_ls"][:n])
        lo, hi = block_boot_ci(diff)
        star = "$^{*}$" if lo > 0 or hi < 0 else ""
        dtxt = f"${np.mean(diff):+.1f}$ [{lo:+.1f}, {hi:+.1f}]{star}"
    ls2 = np.mean(d["h2_ls"]) if "h2_ls" in d else np.nan
    ls2txt = f"{ls2:.1f}" if np.isfinite(ls2) else "--"
    maetxt = f"{mae1:.3f}" if np.isfinite(mae1) else "--"
    rows.append(f"{mdl} & {ls1:.1f} & {dtxt} & {ls2txt} & {maetxt} \\\\")
open(f"{TAB}/tab_chicago_bench.tex", "w").write("\n".join(rows) + "\n")

json.dump({k: {kk: (list(map(float, vv)) if isinstance(vv, list) else vv) for kk, vv in v.items()} for k, v in RES.items()},
          open(f"{OUT}/benchmarks.json", "w"))
for mdl in order:
    d = RES[mdl]
    print(f"{mdl:32s} LS1={np.mean(d['h1_ls']):8.1f}  LS2={np.mean(d.get('h2_ls',[np.nan])):8.1f}  "
          f"MAE1={np.mean(d.get('h1_mae',[np.nan])):.3f}", flush=True)
print("DONE", round(time.time() - t0c), "s", flush=True)
