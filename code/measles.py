"""Third dataset and implementation cross-check.

(1) Weekly measles counts in the 17 districts of the Weser-Ems region, the
    showcase data of the surveillance package: the incumbent hhh4 model is fitted
    NATIVELY in R with rolling one-step predictions scored by the package's own
    scoring rules, and raced against the network state-space model, a
    rolling-window fit of the identical design, and an EWMA baseline, all under
    one protocol.
(2) Cross-check of the Python hhh4-type implementation against the R package on
    both the measles and the Texas data: coefficients, dispersion and training
    log-likelihoods side by side, plus the package's native one-step test score
    on Texas next to the Python implementation's.
(3) Calibration of the state-space model on the measles test window: empirical
    coverage of central 80 and 95 percent one-step predictive intervals and the
    mean squared Pearson residual.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, math, sys, time
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln, digamma
from scipy.stats import nbinom
sys.path.insert(0, _HERE)
from nssm import row_normalise
from nb_filter import nb_mode_filter, nb_logpmf
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

t0 = time.time(); OUT = _os.path.join(_HERE, "out"); TAB = _os.path.join(_ROOT, "paper", "tables")
DATA = _os.path.join(_ROOT, "data")

# ---------------- load measles ----------------
Ym = pd.read_csv(_os.path.join(DATA, "measles_counts.csv")).values.astype(int)   # T x N
NBo = pd.read_csv(_os.path.join(DATA, "measles_nborder.csv")).values
popm = pd.read_csv(_os.path.join(DATA, "measles_pop.csv"))["pop"].values.astype(float)
T, N = Ym.shape; W = row_normalise((NBo == 1).astype(float))
ell = np.log1p(Ym); TRAIN = 52; TARGETS = list(range(TRAIN, T))
p_off = np.log(popm / popm.mean())
print(f"[measles] T={T} N={N} mean count {Ym.mean():.2f} zeros {np.mean(Ym==0):.2f}", flush=True)

def design_m(t):
    l = ell[t - 1]
    return np.column_stack([np.ones(N), W @ l, l, p_off,
                            np.full(N, math.sin(2 * math.pi * t / 52)),
                            np.full(N, math.cos(2 * math.pi * t / 52))])

def fit_nssm(q, r, net=True, T_hi=T):
    K = 6 if net else 5
    XL = [design_m(t)[:, [0, 2, 3, 4, 5]] if not net else design_m(t) for t in range(1, T_hi)]
    return nb_mode_filter(Ym[1:T_hi], XL, q * np.eye(K), np.zeros(K), np.eye(K), r), XL
def pls(fit, XL, lo, hi, r):
    return sum(nb_logpmf(Ym[t], np.exp(np.clip(XL[t - 1] @ fit["m_pred"][t - 1], -30, 30)), r).sum()
               for t in range(lo, hi))
def tune(net=True):
    best = None
    for q in [1e-4, 3e-4, 1e-3, 3e-3]:
        for r in [0.5, 1.0, 2.0, 5.0]:
            f, XL = fit_nssm(q, r, net)
            v = pls(f, XL, 8, TRAIN, r)
            if best is None or v > best[0]: best = (v, q, r)
    return best[1], best[2]
qn, rn = tune(True); qo, ro = tune(False)
fN, XLn = fit_nssm(qn, rn, True); fO, XLo = fit_nssm(qo, ro, False)
def test_stats(f, XL, r):
    ls = []; mm = []; cov80 = []; cov95 = []; pear = []
    for t in TARGETS:
        mu = np.exp(np.clip(XL[t - 1] @ f["m_pred"][t - 1], -30, 30))
        ls.append(float(nb_logpmf(Ym[t], mu, r).sum()))
        mm.append(float(np.mean(np.abs(Ym[t] - nbinom.ppf(0.5, r, r / (r + mu))))))
        lo80, hi80 = nbinom.ppf(0.10, r, r / (r + mu)), nbinom.ppf(0.90, r, r / (r + mu))
        lo95, hi95 = nbinom.ppf(0.025, r, r / (r + mu)), nbinom.ppf(0.975, r, r / (r + mu))
        cov80.append(np.mean((Ym[t] >= lo80) & (Ym[t] <= hi80)))
        cov95.append(np.mean((Ym[t] >= lo95) & (Ym[t] <= hi95)))
        pear.append(np.mean((Ym[t] - mu) ** 2 / (mu * (1 + mu / r) + 1e-9)))
    return (float(np.mean(ls)), float(np.mean(mm)), float(np.mean(cov80)),
            float(np.mean(cov95)), float(np.mean(pear)))
net_stats = test_stats(fN, XLn, rn); non_stats = test_stats(fO, XLo, ro)
print(f"[NSSM] net (q={qn},r={rn}) LS {net_stats[0]:.1f} medMAE {net_stats[1]:.2f} "
      f"cov80/95 {net_stats[2]:.2f}/{net_stats[3]:.2f} pearson {net_stats[4]:.2f}; "
      f"nonet LS {non_stats[0]:.1f} ({time.time()-t0:.0f}s)", flush=True)

# rolling same design (incl. expanding = w >= T)
best = None
for w in [16, 26, 52, 200]:
    p = np.zeros(7); p[-1] = 0.0; tr = 0.0
    for t_or in range(12, TRAIN - 1):
        p = rolling_fit(design_m, Ym, t_or, w, p)
        mu = np.exp(np.clip(design_m(t_or + 1) @ p[:-1], -30, 30))
        tr += float(nb_logpmf(Ym[t_or + 1], mu, math.exp(p[-1])).sum())
    if best is None or tr > best[0]: best = (tr, w)
wsel = best[1]
p = np.zeros(7); p[-1] = 0.0; rls = []; rmm = []
for t_or in range(12, T - 1):
    p = rolling_fit(design_m, Ym, t_or, wsel, p)
    if t_or + 1 in set(TARGETS):
        mu = np.exp(np.clip(design_m(t_or + 1) @ p[:-1], -30, 30)); r = math.exp(p[-1])
        rls.append(float(nb_logpmf(Ym[t_or + 1], mu, r).sum()))
        rmm.append(float(np.mean(np.abs(Ym[t_or + 1] - nbinom.ppf(0.5, r, r / (r + mu))))))
roll_row = (wsel, float(np.mean(rls)), float(np.mean(rmm)))
print(f"[rolling] w={wsel} LS {roll_row[1]:.1f} medMAE {roll_row[2]:.2f} ({time.time()-t0:.0f}s)", flush=True)

# EWMA-NB
bestE = None
for a in [0.1, 0.2, 0.4, 0.6]:
    for r in [0.5, 1.0, 2.0]:
        m = Ym[:4].mean(0) + 0.2; tr = 0.0
        for t in range(4, TRAIN):
            tr += float(nb_logpmf(Ym[t], np.clip(m, 1e-8, None), r).sum()); m = (1 - a) * m + a * Ym[t]
        if bestE is None or tr > bestE[0]: bestE = (tr, a, r)
_, a, rE = bestE
m = Ym[:4].mean(0) + 0.2
for t in range(4, TRAIN): m = (1 - a) * m + a * Ym[t]
els = []; emm = []
for t in TARGETS:
    mu = np.clip(m, 1e-8, None)
    els.append(float(nb_logpmf(Ym[t], mu, rE).sum()))
    emm.append(float(np.mean(np.abs(Ym[t] - nbinom.ppf(0.5, rE, rE / (rE + mu))))))
    m = (1 - a) * m + a * Ym[t]
ew_row = (float(np.mean(els)), float(np.mean(emm)))

# ---------------- native hhh4 rows via the package's own scores ----------------
def native_row(prefix, Y, targets):
    logs = pd.read_csv(_os.path.join(DATA, f"{prefix}_hhh4_logs.csv")).values
    pred = pd.read_csv(_os.path.join(DATA, f"{prefix}_hhh4_pred.csv")).values
    psi = pd.read_csv(_os.path.join(DATA, f"{prefix}_hhh4_psi.csv")).iloc[:, -1].values
    # detect the size parameterisation by matching the package's own log scores
    size_vec = np.exp(psi)
    agree = []
    for j, t in enumerate(targets):
        lp = nb_logpmf(Y[t], np.clip(pred[j], 1e-9, None), max(size_vec[j], 1e-6))
        fin = np.isfinite(logs[j])
        if fin.any(): agree.append(float(np.max(np.abs(-lp[fin] - logs[j][fin]))))
    print(f"[{prefix} native] scoring agreement with package on finite cells: "
          f"max |diff| {max(agree):.2e}; pred NA {int(np.isnan(pred).sum())}", flush=True)
    ls = []; mm = []
    for j, t in enumerate(targets):
        mu = np.clip(pred[j], 1e-300, None); r = max(size_vec[j], 1e-6)
        row = -logs[j]
        bad = ~np.isfinite(row)
        if bad.any(): row[bad] = nb_logpmf(Y[t][bad], mu[bad], r)
        ls.append(float(np.sum(row)))
        mm.append(float(np.mean(np.abs(Y[t] - nbinom.ppf(0.5, r, r / (r + mu))))))
    return float(np.mean(ls)), float(np.mean(mm)), "exp(psi)"
mn_ls, mn_mm, _ = native_row("measles", Ym, TARGETS)
print(f"[measles native hhh4] LS {mn_ls:.1f} medMAE {mn_mm:.2f}", flush=True)
Yt = pd.read_csv(_os.path.join(DATA, "tx_counts.csv")).values.astype(int)
tx_ls, tx_mm, _ = native_row("tx", Yt, list(range(72, Yt.shape[0])))
print(f"[texas native hhh4] LS {tx_ls:.1f} medMAE {tx_mm:.1f}", flush=True)

# ---------------- cross-check: python hhh4 on both datasets ----------------
def py_hhh4(Y, Wm, TRAIN_, period=52):
    Tn, Nn = Y.shape
    a0 = np.log(Y[:TRAIN_].mean(0) + 0.1)
    sw = np.array([math.sin(2 * math.pi * t / period) for t in range(1, TRAIN_)])
    cw = np.array([math.cos(2 * math.pi * t / period) for t in range(1, TRAIN_)])
    Yl = Y[:TRAIN_ - 1]; Yt_ = Y[1:TRAIN_]; WY = Yl @ Wm.T
    def unpack(p): return p[:Nn], p[Nn], p[Nn + 1], math.exp(p[Nn + 2]), math.exp(p[Nn + 3]), math.exp(p[Nn + 4])
    def nll(p):
        a, gs, gc, lam, phi, r = unpack(p)
        nu = np.exp(np.clip(a[None, :] + gs * sw[:, None] + gc * cw[:, None], -30, 30))
        mu = np.clip(nu + lam * Yl + phi * WY, 1e-9, None)
        return float(-np.sum(gammaln(Yt_ + r) - gammaln(r) - gammaln(Yt_ + 1)
                             + r * np.log(r / (r + mu)) + Yt_ * np.log(mu / (r + mu))))
    p0 = np.concatenate([a0, [0, 0], [math.log(0.5), math.log(0.02), math.log(1.5)]])
    res = minimize(nll, p0, method="L-BFGS-B", options=dict(maxiter=500))
    a, gs, gc, lam, phi, r = unpack(res.x)
    return dict(loglik=-float(res.fun), lam=float(lam), phi=float(phi), size=float(r))
pyM = py_hhh4(Ym, W, 52)
Wt = pd.read_csv(_os.path.join(DATA, "tx_W.csv")).values
pyT = py_hhh4(Yt, Wt, 72)
Rco_m = pd.read_csv(_os.path.join(DATA, "measles_hhh4_coef.csv"))
Rco_t = pd.read_csv(_os.path.join(DATA, "tx_hhh4_coef.csv"))
def rget(df, nm): return float(df[df.name == nm].value.iloc[0])
cross = dict(
  measles=dict(R=dict(loglik=rget(Rco_m, "loglik"), lam=math.exp(rget(Rco_m, "ar.1")), phi=math.exp(rget(Rco_m, "ne.1")),
                      size=1.0 / rget(Rco_m, "overdisp")),
               PY=pyM),
  texas=dict(R=dict(loglik=rget(Rco_t, "loglik"), lam=math.exp(rget(Rco_t, "ar.1")), phi=math.exp(rget(Rco_t, "ne.1")),
                    size=1.0 / rget(Rco_t, "overdisp")),
             PY=pyT))
print("[crosscheck]", json.dumps(cross), flush=True)

# ---------------- fragment + json ----------------
rows = [("NSSM net (NB, offsets, seasonal)", net_stats[0], net_stats[1], net_stats[3]),
        ("NSSM no-network", non_stats[0], non_stats[1], non_stats[3]),
        ("hhh4, native \\texttt{surveillance} fit", mn_ls, mn_mm, None),
        (f"Rolling-window NB, same design ($w={roll_row[0]}$)", roll_row[1], roll_row[2], None),
        ("EWMA--NB (per node)", ew_row[0], ew_row[1], None)]
open(_os.path.join(TAB, "tab_measles.tex"), "w").write("\n".join(
    f"{a} & {b:.1f} & {c:.2f} & " + (f"{d:.2f}" if d is not None else "--") + " \\\\"
    for a, b, c, d in rows) + "\n")
json.dump(dict(nssm=dict(q=qn, r=rn, ls=net_stats[0], medmae=net_stats[1],
                         cov80=net_stats[2], cov95=net_stats[3], pearson=net_stats[4]),
               nonet=dict(ls=non_stats[0]), native_hhh4=dict(ls=mn_ls, medmae=mn_mm),
               rolling=dict(w=roll_row[0], ls=roll_row[1], medmae=roll_row[2]),
               ewma=dict(ls=ew_row[0], medmae=ew_row[1]),
               texas_native_hhh4=dict(ls=tx_ls, medmae=tx_mm),
               crosscheck=cross, dims=dict(T=T, N=N)),
          open(_os.path.join(OUT, "measles.json"), "w"))
print("DONE", round(time.time() - t0), flush=True)
