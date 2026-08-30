import re
"""Time-varying competitors fitted to the identical design.

Adds two competitors that let coefficients move, so the comparison isolates the
coefficient dynamics and nothing else:
  (a) rolling-window negative-binomial regression with the same regressors as
      the final state-space model (window width tuned on the training period),
      run on both the Texas and the Chicago data;
  (b) a penalised-spline hhh4-type model on Texas, with the endemic component
      carrying unit intercepts and seasonality and the epidemic own- and
      network-coefficients varying over time through cubic B-splines
      (number of interior knots tuned on the training period, light ridge on
      the spline deviations, refitted at every forecast origin).
Also computes the real-time coefficient paths of all three approaches and
their lead-lag relation to statewide epidemic growth, and draws the figure.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, math, sys, time
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln, digamma, expit
from scipy.interpolate import BSpline
from scipy.stats import nbinom
sys.path.insert(0, _HERE)
from nssm import row_normalise
from nb_filter import nb_mode_filter, nb_logpmf

t0 = time.time(); OUT = _os.path.join(_HERE, "out"); TAB = _os.path.join(_ROOT, "paper", "tables")
FIG = _os.path.join(_ROOT, "paper", "figs")

# ---------------- Texas data (identical to the main pipeline) ----------------
raw = pd.read_csv(_os.path.join(_ROOT, "data", "covid_us.csv"))
lk = pd.read_csv(_os.path.join(_ROOT, "data", "uid_lookup.csv"))
tx = raw[(raw.Province_State == "Texas") & raw.Admin2.notna() & (raw.Lat != 0)].copy()
tx = tx.merge(lk[["FIPS", "Population"]], on="FIPS", how="left")
tx = tx[tx.Population > 0].reset_index(drop=True)
dc = [c for c in tx.columns if c.count("/") == 2]
cum = tx[dc].values.astype(float)
dates = pd.to_datetime(dc, format="%m/%d/%y")
start = np.searchsorted(dates, pd.Timestamp("2020-03-01"))
new = np.clip(np.diff(cum[:, start - 1:], axis=1), 0, None)
nw = new.shape[1] // 7
wk = new[:, :nw * 7].reshape(len(tx), nw, 7).sum(2)
T = min(nw, 92); Yw = wk[:, :T].T.astype(int); N = Yw.shape[1]
pop = tx.Population.values.astype(float); lat, lon = tx.Lat.values, tx.Long_.values
D = np.sqrt((lat[:, None] - lat[None, :]) ** 2 + (lon[:, None] - lon[None, :]) ** 2)
np.fill_diagonal(D, np.inf)
A = np.zeros((N, N))
for i in range(N): A[i, np.argsort(D[i])[:5]] = 1
W5 = row_normalise(((A + A.T) > 0).astype(float))
ell = np.log1p(Yw); TRAIN = 72; TARGETS = list(range(TRAIN, T)); p_off = np.log(pop / 1e5)

def design_tx(t):
    l = ell[t - 1]
    return np.column_stack([np.ones(N), W5 @ l, l, p_off,
                            np.full(N, math.sin(2 * math.pi * t / 52)),
                            np.full(N, math.cos(2 * math.pi * t / 52))])

# ---------------- generic rolling-window NB fit (analytic gradient) ----------------
def nb_nll_grad(par, X, y):
    b, lr = par[:-1], par[-1]; r = math.exp(lr)
    mu = np.exp(np.clip(X @ b, -30, 30))
    ll = np.sum(gammaln(y + r) - gammaln(r) - gammaln(y + 1)
                + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))
    dmu = y / np.clip(mu, 1e-12, None) - (y + r) / (mu + r)
    gb = X.T @ (dmu * mu)
    glr = r * np.sum(digamma(y + r) - digamma(r) + np.log(r / (r + mu)) + (mu - y) / (mu + r))
    return -ll, -np.concatenate([gb, [glr]])

def rolling_fit(designs, Y, t_or, w, par0):
    lo = max(1, t_or - w + 1)
    X = np.vstack([designs(t) for t in range(lo, t_or + 1)])
    y = np.concatenate([Y[t] for t in range(lo, t_or + 1)]).astype(float)
    res = minimize(nb_nll_grad, par0, args=(X, y), jac=True, method="L-BFGS-B",
                   options=dict(maxiter=300))
    return res.x

def rolling_eval(designs, Y, K, train_hi, T_all, targets, windows, tag):
    """Tune w on the training window, then score one-step plug-in on the targets,
    and record the full real-time spillover path."""
    par0 = np.zeros(K + 1); par0[-1] = math.log(2.0)
    s0 = max(max(windows), 10)          # common origin set: identical for every w
    best = None
    for w in windows:
        p = par0.copy(); tr = []
        for t_or in range(s0, train_hi - 1):
            p = rolling_fit(designs, Y, t_or, w, p)
            mu = np.exp(np.clip(designs(t_or + 1) @ p[:-1], -30, 30))
            tr.append(float(nb_logpmf(Y[t_or + 1], mu, math.exp(p[-1])).sum()))
        m = float(np.mean(tr))
        if best is None or m > best[0]: best = (m, w)
        print(f"[{tag} tune] w={w} mean per-origin train LS {m:.1f} over {len(tr)} common origins "
              f"({time.time()-t0:.0f}s)", flush=True)
    w = best[1]
    return (w,) + rolling_run(designs, Y, K, T_all, targets, w)

def rolling_run(designs, Y, K, T_all, targets, w):
    p = np.zeros(K + 1); p[-1] = math.log(2.0)
    ls = []; mm = []; ma = []; path = {}
    for t_or in range(min(max(w, 10), 12), T_all - 1):
        p = rolling_fit(designs, Y, t_or, w, p)
        path[t_or] = float(p[1])
        if t_or + 1 in targets:
            mu = np.exp(np.clip(designs(t_or + 1) @ p[:-1], -30, 30)); r = math.exp(p[-1])
            ls.append(float(nb_logpmf(Y[t_or + 1], mu, r).sum()))
            mm.append(float(np.mean(np.abs(Y[t_or + 1] - nbinom.ppf(0.5, r, r / (r + mu))))))
            ma.append(float(np.mean(np.abs(Y[t_or + 1] - mu))))
    return float(np.mean(ls)), float(np.mean(mm)), float(np.mean(ma)), path, list(map(float, ls))

roll_w, roll_ls, roll_mm, roll_ma, roll_path, roll_ls_list = rolling_eval(
    design_tx, Yw, 6, TRAIN, T, set(TARGETS), [8, 12, 16, 24, 36, 48], "TX-roll")
exp_ls, exp_mm, exp_ma, exp_path, exp_ls_list = rolling_run(design_tx, Yw, 6, T, set(TARGETS), 10**6)
print(f"[TX expanding] LS/wk {exp_ls:.1f} medMAE {exp_mm:.1f} ({time.time()-t0:.0f}s)", flush=True)
print(f"[TX rolling] w={roll_w} LS/wk {roll_ls:.1f} medMAE {roll_mm:.1f} ({time.time()-t0:.0f}s)", flush=True)

# ---------------- penalised-spline hhh4 on Texas ----------------
def spline_basis(Kk):
    kn = np.linspace(1, T, Kk + 2)[1:-1]
    knots = np.concatenate([[1, 1, 1, 1], kn, [T, T, T, T]])
    def B(t):
        return BSpline.design_matrix(np.atleast_1d(float(t)), knots, 3).toarray()[0]
    m = Kk + 4
    return B, m

def hhh4tv_fit(Kk, t_hi, par0=None, ridge=1e-3):
    B, m = spline_basis(Kk)
    Bt = np.array([B(t) for t in range(1, t_hi)])          # (t_hi-1, m)
    Yl = Yw[:t_hi - 1]; Yt = Yw[1:t_hi]; WY = Yl @ W5.T
    sw = np.array([math.sin(2 * math.pi * t / 52) for t in range(1, t_hi)])
    cw = np.array([math.cos(2 * math.pi * t / 52) for t in range(1, t_hi)])
    npar = N + 2 + 2 * m + 1
    if par0 is None:
        par0 = np.zeros(npar); par0[:N] = np.log(Yw[:min(20, t_hi)].mean(0) + 0.1)
        par0[N + 2] = -1.0; par0[N + 2 + m] = -2.0; par0[-1] = math.log(2.0)
    sl = dict(a=slice(0, N), gs=N, gc=N + 1, cl=slice(N + 2, N + 2 + m),
              cf=slice(N + 2 + m, N + 2 + 2 * m), lr=npar - 1)
    def nllg(p):
        a = p[sl["a"]]; gs = p[sl["gs"]]; gc = p[sl["gc"]]
        cl = p[sl["cl"]]; cf = p[sl["cf"]]; r = math.exp(p[sl["lr"]])
        eta_l = Bt @ cl; eta_f = Bt @ cf
        lam = np.log1p(np.exp(np.clip(eta_l, -30, 30)))      # softplus
        phi = np.log1p(np.exp(np.clip(eta_f, -30, 30)))
        nu = np.exp(np.clip(a[None, :] + gs * sw[:, None] + gc * cw[:, None], -30, 30))
        mu = np.clip(nu + lam[:, None] * Yl + phi[:, None] * WY, 1e-8, None)
        ll = np.sum(gammaln(Yt + r) - gammaln(r) - gammaln(Yt + 1)
                    + r * np.log(r / (r + mu)) + Yt * np.log(mu / (r + mu)))
        dmu = Yt / mu - (Yt + r) / (mu + r)
        g = np.zeros(npar)
        g[sl["a"]] = (dmu * nu).sum(0)
        g[sl["gs"]] = np.sum(dmu * nu * sw[:, None]); g[sl["gc"]] = np.sum(dmu * nu * cw[:, None])
        sig_l = expit(np.clip(eta_l, -30, 30)); sig_f = expit(np.clip(eta_f, -30, 30))
        g[sl["cl"]] = Bt.T @ ((dmu * Yl).sum(1) * sig_l)
        g[sl["cf"]] = Bt.T @ ((dmu * WY).sum(1) * sig_f)
        g[sl["lr"]] = r * np.sum(digamma(Yt + r) - digamma(r) + np.log(r / (r + mu)) + (mu - Yt) / (mu + r))
        pen = ridge * (np.sum(np.diff(cl) ** 2) + np.sum(np.diff(cf) ** 2))
        gp = np.zeros(npar)
        gp[sl["cl"]][:-1] += -2 * ridge * np.diff(cl); gp[sl["cl"]][1:] += 2 * ridge * np.diff(cl)
        gp[sl["cf"]][:-1] += -2 * ridge * np.diff(cf); gp[sl["cf"]][1:] += 2 * ridge * np.diff(cf)
        return -(ll - pen), -(g - gp)
    res = minimize(nllg, par0, jac=True, method="L-BFGS-B", options=dict(maxiter=400))
    return res.x, sl, B, m

def hhh4tv_predict(p, sl, B, tgt):
    a = p[sl["a"]]; gs = p[sl["gs"]]; gc = p[sl["gc"]]
    lam = math.log1p(math.exp(min(30, float(B(tgt) @ p[sl["cl"]]))))
    phi = math.log1p(math.exp(min(30, float(B(tgt) @ p[sl["cf"]]))))
    nu = np.exp(np.clip(a + gs * math.sin(2 * math.pi * tgt / 52) + gc * math.cos(2 * math.pi * tgt / 52), -30, 30))
    mu = np.clip(nu + lam * Yw[tgt - 1] + phi * (W5 @ Yw[tgt - 1]), 1e-8, None)
    return mu, math.exp(p[sl["lr"]]), phi

# tune K on the training window (fit to 64, score 65..71)
bestK = None
for Kk in [3, 6]:
    p, sl, B, m = hhh4tv_fit(Kk, 65)
    v = 0.0
    for tgt in range(65, TRAIN):
        mu, r, _ = hhh4tv_predict(p, sl, B, tgt)
        v += float(nb_logpmf(Yw[tgt], mu, r).sum())
    print(f"[TX hhh4tv tune] K={Kk} val LS {v:.0f} ({time.time()-t0:.0f}s)", flush=True)
    if bestK is None or v > bestK[0]: bestK = (v, Kk)
Kk = bestK[1]
p, sl, B, m = hhh4tv_fit(Kk, TRAIN)
tv_ls = []; tv_mm = []; tv_cov = []; tv_path = {}
for t_or in range(20, T - 1):
    p, sl, B, m = hhh4tv_fit(Kk, t_or + 1, par0=p)
    _, _, phi_now = hhh4tv_predict(p, sl, B, t_or)
    tv_path[t_or] = phi_now
    tgt = t_or + 1
    if tgt in set(TARGETS):
        mu, r, _ = hhh4tv_predict(p, sl, B, tgt)
        tv_ls.append(float(nb_logpmf(Yw[tgt], mu, r).sum()))
        tv_mm.append(float(np.mean(np.abs(Yw[tgt] - nbinom.ppf(0.5, r, r / (r + mu))))))
        lo, hi = nbinom.ppf(0.025, r, r / (r + mu)), nbinom.ppf(0.975, r, r / (r + mu))
        tv_cov.append(float(np.mean((Yw[tgt] >= lo) & (Yw[tgt] <= hi))))
print(f"[TX hhh4tv] K={Kk} LS/wk {np.mean(tv_ls):.1f} medMAE {np.mean(tv_mm):.1f} ({time.time()-t0:.0f}s)", flush=True)

# ---------------- filter real-time path (final model, as in the main pipeline) ----------------
Q_, R_ = 1e-3, 2.0
XL = [design_tx(t) for t in range(1, T)]
fit = nb_mode_filter(Yw[1:T], XL, Q_ * np.eye(6), np.zeros(6), np.eye(6), R_)
filt_path = {t: float(fit["m_filt"][t - 1, 1]) for t in range(1, T)}
filt_ls = [float(nb_logpmf(Yw[t], np.exp(np.clip(XL[t - 1] @ fit["m_pred"][t - 1], -30, 30)), R_).sum())
           for t in TARGETS]

# block-bootstrap CI for filter minus rolling, and filter minus spline (weekly LS)
def mbb_ci(d, Bn=4000, blk=4, seed=7):
    d = np.array(d); n = len(d); rng = np.random.default_rng(seed); out = []
    for _ in range(Bn):
        idx = []
        while len(idx) < n:
            s0 = rng.integers(0, n - blk + 1); idx += list(range(s0, s0 + blk))
        out.append(d[idx[:n]].mean())
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))
d_roll = [a - b for a, b in zip(filt_ls, roll_ls_list)]
d_tv = [a - b for a, b in zip(filt_ls, tv_ls)]
ci_roll = mbb_ci(d_roll); ci_tv = mbb_ci(d_tv)
print(f"[TX deltas] filter-rolling {np.mean(d_roll):+.1f} CI {ci_roll}; "
      f"filter-spline {np.mean(d_tv):+.1f} CI {ci_tv}", flush=True)

# ---------------- lead-lag against statewide growth ----------------
tot = Yw.sum(1).astype(float)
G = np.diff(np.log1p(tot))            # growth, index t-1 -> value at week t
ts = sorted(set(filt_path) & set(roll_path) & set(tv_path))
Gv = np.array([G[t - 1] for t in ts])
def xcorr(path):
    v = np.array([path[t] for t in ts])
    v = (v - v.mean()) / (v.std() + 1e-12); g = (Gv - Gv.mean()) / (Gv.std() + 1e-12)
    cs = {}
    for k in range(-4, 5):
        if k >= 0: cs[k] = float(np.corrcoef(v[:len(v) - k or None], g[k:])[0, 1]) if len(v) - k > 3 else np.nan
        else: cs[k] = float(np.corrcoef(v[-k:], g[:k])[0, 1])
    kbest = max(cs, key=lambda k: cs[k])
    return cs, kbest
cs_f, k_f = xcorr(filt_path); cs_r, k_r = xcorr(roll_path); cs_t, k_t = xcorr(tv_path)
print(f"[lead-lag] filter lag0 {cs_f[0]:+.2f} best k={k_f}; rolling lag0 {cs_r[0]:+.2f} best k={k_r}; "
      f"spline lag0 {cs_t[0]:+.2f} best k={k_t}", flush=True)

# ---------------- figure ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True,
                       gridspec_kw=dict(height_ratios=[1, 1.6], hspace=0.12))
wkx = np.arange(T)
ax[0].fill_between(wkx, tot, color="0.82")
ax[0].set_yscale("log"); ax[0].set_ylabel("statewide cases")
ax[0].axvspan(TRAIN, T - 1, color="0.93")
tt = np.array(ts)
ax[1].plot(tt, [filt_path[t] for t in ts], "-", color="black", lw=1.6, label="filter (state space)")
ax[1].plot(tt, [roll_path[t] for t in ts], "--", color="0.35", lw=1.4, label=f"rolling window ($w={roll_w}$)")
ax[1].plot(tt, [tv_path[t] for t in ts], ":", color="0.35", lw=1.8, label="spline hhh4 ($\\phi_t$)")
ax[1].axhline(0, color="0.7", lw=0.8)
ax[1].axvspan(TRAIN, T - 1, color="0.93")
ax[1].set_xlabel("week"); ax[1].set_ylabel("network coefficient (real time)")
ax[1].legend(frameon=False, fontsize=8, loc="upper right")
fig.savefig(_os.path.join(FIG, "fig_tvpaths.pdf"), bbox_inches="tight")
print("[figure] written", flush=True)

# ---------------- Chicago rolling row ----------------
dfc = pd.read_csv(_os.path.join(_ROOT, "data", "crime.csv"), index_col=0)
Yc = dfc.values.T.astype(int); Tc, Nc = Yc.shape
import scipy.io as sio
Ac = sio.mmread(_os.path.join(_ROOT, "data", "neighborhood.mtx")).tocsr()
Ac = ((Ac + Ac.T) > 0).astype(float).toarray(); np.fill_diagonal(Ac, 0)
Wc = row_normalise(Ac); ellc = np.log1p(Yc)
TRAINC = 60; oc = np.log(Yc[:TRAINC].mean(0) + 0.05)
def design_ch(t):
    l = ellc[t - 1]
    return np.column_stack([np.ones(Nc), Wc @ l, l, oc,
                            np.full(Nc, math.sin(2 * math.pi * (t % 12) / 12)),
                            np.full(Nc, math.cos(2 * math.pi * (t % 12) / 12))])
cw, cls, cmm, cma, cpath, _ = rolling_eval(design_ch, Yc, 6, TRAINC, Tc,
                                      set(range(TRAINC, Tc)), [12, 24, 36, 48], "CH-roll")
cels, cemm, cema, _, _ = rolling_run(design_ch, Yc, 6, Tc, set(range(TRAINC, Tc)), 10**6)
print(f"[CH expanding] LS/mo {cels:.1f} meanMAE {cema:.3f} ({time.time()-t0:.0f}s)", flush=True)
print(f"[CH rolling] w={cw} LS/mo {cls:.1f} medMAE {cmm:.2f} ({time.time()-t0:.0f}s)", flush=True)

# ---------------- fragments + json ----------------
bt = open(_os.path.join(TAB, "tab_covid_bench.tex")).read()
if "Rolling-window" not in bt:
    lines = bt.rstrip().split("\n")
    ins = [f"Rolling-window NB, same design ($w={roll_w}$) & {roll_ls:.1f} & {roll_mm:.1f} & -- \\\\",
           f"Spline hhh4-type NB (time-varying $\\lambda_t,\\phi_t$) & {np.mean(tv_ls):.1f} & {np.mean(tv_mm):.1f} & -- \\\\"]
    lines = lines[:3] + ins + lines[3:]
    open(_os.path.join(TAB, "tab_covid_bench.tex"), "w").write("\n".join(lines) + "\n")
    print("[frag] covid bench rows inserted", flush=True)
bc = open(_os.path.join(TAB, "tab_chicago_bench.tex")).read()
if "Rolling-window" not in bc:
    lines = bc.rstrip().split("\n")
    lines.insert(2, f"Rolling-window NB, same design ($w={cw}$) & {cls:.1f} & -- & -- & {cmm:.2f} \\\\")
    open(_os.path.join(TAB, "tab_chicago_bench.tex"), "w").write("\n".join(lines) + "\n")
    print("[frag] chicago bench row inserted", flush=True)

json.dump(dict(
    texas=dict(rolling=dict(w=int(roll_w), ls=roll_ls, medmae=roll_mm),
               spline=dict(K=int(Kk), ls=float(np.mean(tv_ls)), medmae=float(np.mean(tv_mm))),
               filter_ls=float(np.mean(filt_ls)),
               delta_roll=dict(mean=float(np.mean(d_roll)), ci=ci_roll),
               delta_spline=dict(mean=float(np.mean(d_tv)), ci=ci_tv),
               leadlag=dict(filter=dict(lag0=cs_f[0], best=int(k_f)),
                            rolling=dict(lag0=cs_r[0], best=int(k_r)),
                            spline=dict(lag0=cs_t[0], best=int(k_t))),
               amp=dict(filter=[float(min(filt_path[t] for t in ts)), float(max(filt_path[t] for t in ts))],
                        rolling=[float(min(roll_path[t] for t in ts)), float(max(roll_path[t] for t in ts))],
                        spline=[float(min(tv_path[t] for t in ts)), float(max(tv_path[t] for t in ts))])),
    chicago=dict(rolling=dict(w=int(cw), ls=cls, medmae=cmm))),
    open(_os.path.join(OUT, "tvcompare.json"), "w"))
print("DONE", round(time.time() - t0), "s", flush=True)

# ---------------- score-driven (GAS) time-varying spillover, Texas ----------------
def gas_fit_predict():
    def unpack(p):
        b0, b2, b3, bs, bc, om, a, zb, lr = p
        return b0, b2, b3, bs, bc, om, a, 0.999 / (1 + math.exp(-zb)), math.exp(lr)
    def run(p, t_hi, score_from):
        b0, b2, b3, bs, bc, om, a, bpers, r = unpack(p)
        f = om / max(1 - bpers, 1e-3); tot = 0.0; ls_t = []; fpath = {}
        for t in range(1, t_hi):
            X = design_tx(t); xnet = X[:, 1]
            eta = b0 + f * xnet + b2 * X[:, 2] + b3 * X[:, 3] + bs * X[:, 4] + bc * X[:, 5]
            mu = np.exp(np.clip(eta, -30, 30))
            lp = float(nb_logpmf(Yw[t], mu, r).sum())
            if t >= score_from: tot += lp; ls_t.append(lp)
            fpath[t] = f
            wgt = r * mu / (r + mu)
            sc = float(np.sum((Yw[t] - mu) * (r / (r + mu)) * xnet))
            info = float(np.sum(xnet ** 2 * wgt)) + 1e-6
            f = om + bpers * f + a * sc / info
        return tot, ls_t, fpath
    def nll(p):
        v, _, _ = run(p, TRAIN, 2)
        return -v
    p0 = np.array([-1.0, 0.5, 0.5, 0.0, 0.0, 0.02, 0.05, 2.0, math.log(2.0)])
    res = minimize(nll, p0, method="L-BFGS-B", options=dict(maxiter=400))
    _, ls_test, fpath = run(res.x, T, TRAIN)
    b0, b2, b3, bs, bc, om, a, bpers, r = unpack(res.x)
    mm = []
    f = fpath  # recompute medians along the same recursion
    _, _, fp2 = run(res.x, T, TRAIN)
    for j, t in enumerate(TARGETS):
        X = design_tx(t)
        eta = b0 + fp2[t] * X[:, 1] + b2 * X[:, 2] + b3 * X[:, 3] + bs * X[:, 4] + bc * X[:, 5]
        mu = np.exp(np.clip(eta, -30, 30))
        mm.append(float(np.mean(np.abs(Yw[t] - nbinom.ppf(0.5, r, r / (r + mu))))))
    return float(np.mean(ls_test)), float(np.mean(mm)), fp2, ls_test, dict(persistence=bpers, alpha=a)
gas_ls, gas_mm, gas_path, gas_ls_list, gas_par = gas_fit_predict()
print(f"[TX GAS] LS/wk {gas_ls:.1f} medMAE {gas_mm:.1f} persistence {gas_par['persistence']:.2f} "
      f"({time.time()-t0:.0f}s)", flush=True)
cs_g, k_g = xcorr({t: gas_path[t] for t in ts if t in gas_path})
print(f"[GAS lead-lag] lag0 {cs_g[0]:+.2f} best {k_g}", flush=True)

# ---------------- native hhh4 per-week scores (for the CI on the incumbent margin) ----------------
import pandas as _pd
_logs = _pd.read_csv(_os.path.join(_ROOT, "data", "tx_hhh4_logs.csv")).values
_pred = _pd.read_csv(_os.path.join(_ROOT, "data", "tx_hhh4_pred.csv")).values
_psi = _pd.read_csv(_os.path.join(_ROOT, "data", "tx_hhh4_psi.csv")).iloc[:, -1].values
nat_ls_list = []
for j, t in enumerate(TARGETS):
    row = -_logs[j]; bad = ~np.isfinite(row)
    if bad.any():
        r_ = max(float(np.exp(_psi[j])), 1e-6)
        row[bad] = nb_logpmf(Yw[t][bad], np.clip(_pred[j][bad], 1e-300, None), r_)
    nat_ls_list.append(float(np.sum(row)))
d_nat = [a - b for a, b in zip(filt_ls, nat_ls_list)]
d_exp = [a - b for a, b in zip(filt_ls, exp_ls_list)]
d_gas = [a - b for a, b in zip(filt_ls, gas_ls_list)]
ci_nat = mbb_ci(d_nat, seed=11); ci_exp = mbb_ci(d_exp, seed=12); ci_gas = mbb_ci(d_gas, seed=13)
print(f"[TX deltas 2] vs native hhh4 {np.mean(d_nat):+.1f} CI {ci_nat}; "
      f"vs expanding {np.mean(d_exp):+.1f} CI {ci_exp}; vs GAS {np.mean(d_gas):+.1f} CI {ci_gas}", flush=True)

# ---------------- calibration (goodness of absolute fit) ----------------
def coverage_pearson(Y, XLl, fitl, r, targets):
    c80 = []; c95 = []; pe = []
    for t in targets:
        mu = np.exp(np.clip(XLl[t - 1] @ fitl["m_pred"][t - 1], -30, 30))
        lo80, hi80 = nbinom.ppf(0.10, r, r / (r + mu)), nbinom.ppf(0.90, r, r / (r + mu))
        lo95, hi95 = nbinom.ppf(0.025, r, r / (r + mu)), nbinom.ppf(0.975, r, r / (r + mu))
        c80.append(np.mean((Y[t] >= lo80) & (Y[t] <= hi80)))
        c95.append(np.mean((Y[t] >= lo95) & (Y[t] <= hi95)))
        pe.append(np.mean((Y[t] - mu) ** 2 / (mu * (1 + mu / r) + 1e-9)))
    return float(np.mean(c80)), float(np.mean(c95)), float(np.mean(pe))
tx_cov = coverage_pearson(Yw, XL, fit, R_, TARGETS)
XLc = [design_ch(t) for t in range(1, Tc)]
fitc = nb_mode_filter(Yc[1:Tc], XLc, 3e-4 * np.eye(6), np.zeros(6), np.eye(6), 5.0)
ch_cov = coverage_pearson(Yc, XLc, fitc, 5.0, range(TRAINC, Tc))
print(f"[calibration] TX cov80/95 {tx_cov[0]:.2f}/{tx_cov[1]:.2f} pearson {tx_cov[2]:.2f}; "
      f"CH cov80/95 {ch_cov[0]:.2f}/{ch_cov[1]:.2f} pearson {ch_cov[2]:.2f}; "
      f"spline cov95 {np.mean(tv_cov):.2f}", flush=True)

# ---------------- rewrite both benchmark fragments from one source ----------------
c2 = json.load(open(_os.path.join(OUT, "covid2.json")))
rows_tx = [
 ("NSSM net (NB, pop.\\ offsets, seasonal)", c2["ls"]["net_h1"], c2["medmae"]["net"], c2["ls"]["net_h2"]),
 ("NSSM net (NB, plug-in at filtered mode)", float(np.mean(filt_ls)), None, None),
 ("NSSM no-network (NB, pop.\\ offsets, seasonal)", c2["ls"]["nonet_h1"], c2["medmae"]["nonet"], c2["ls"]["nonet_h2"]),
 ("hhh4, native \\texttt{surveillance} fit", float(np.mean(nat_ls_list)), 98.1, None),
 ("Spline hhh4-type NB (time-varying $\\lambda_t,\\phi_t$)", float(np.mean(tv_ls)), float(np.mean(tv_mm)), None),
 ("Score-driven spillover (GAS $\\phi_t$), same design", gas_ls, gas_mm, None),
 (f"Rolling-window NB, same design ($w={roll_w}$)", roll_ls, roll_mm, None),
 ("Expanding-window NB, same design", exp_ls, exp_mm, None),
 ("EWMA--NB (per node)", c2["bench"]["ewma"][0], c2["bench"]["ewma"][1], None),
 ("PNAR(1) log-linear ($+\\log$ pop.)", c2["bench"]["pnar_log"][0], c2["bench"]["pnar_log"][1], None),
 ("PNAR(1) linear", c2["bench"]["pnar_lin"][0], c2["bench"]["pnar_lin"][1], None)]
open(_os.path.join(TAB, "tab_covid_bench.tex"), "w").write("\n".join(
    f"{a} & {b:.1f} & " + (f"{c:.1f}" if c is not None else "--") + " & " +
    (f"{d:.1f}" if d is not None else "--") + " \\\\" for a, b, c, d in rows_tx) + "\n")
bc = open(_os.path.join(TAB, "tab_chicago_bench.tex")).read()
bc = re.sub(r"Rolling-window NB, same design \(\$w=\d+\$\) & -?[\d.]+ & -- & -- & [\d.]+ \\\\",
            f"Rolling-window NB, same design ($w={cw}$) & {cls:.1f} & -- & -- & {cma:.3f} \\\\\n"
            f"Expanding-window NB, same design & {cels:.1f} & -- & -- & {cema:.3f} \\\\", bc)
open(_os.path.join(TAB, "tab_chicago_bench.tex"), "w").write(bc)
print("[frags] rewritten", flush=True)

j = json.load(open(_os.path.join(OUT, "tvcompare.json")))
j.update(dict(v2=dict(
    tx_rolling=dict(w=int(roll_w), ls=roll_ls, medmae=roll_mm),
    tx_expanding=dict(ls=exp_ls, medmae=exp_mm),
    tx_gas=dict(ls=gas_ls, medmae=gas_mm, persistence=gas_par["persistence"], lag0=cs_g[0], best=int(k_g)),
    tx_native=dict(ls=float(np.mean(nat_ls_list))),
    deltas=dict(native=dict(mean=float(np.mean(d_nat)), ci=ci_nat),
                expanding=dict(mean=float(np.mean(d_exp)), ci=ci_exp),
                gas=dict(mean=float(np.mean(d_gas)), ci=ci_gas)),
    calib=dict(tx=tx_cov, ch=ch_cov, spline_cov95=float(np.mean(tv_cov))),
    ch_rolling=dict(w=int(cw), ls=cls, meanmae=cma),
    ch_expanding=dict(ls=cels, meanmae=cema))))
json.dump(j, open(_os.path.join(OUT, "tvcompare.json"), "w"))
print("DONE-V2", round(time.time() - t0), flush=True)
