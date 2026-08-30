"""
Unified, leak-free re-evaluation of the Chicago forecasting tables.

Fixes implemented relative to the original chicago.py evaluation:
  1. Nested per-origin tuning: q is re-selected at every forecast origin using
     only data up to that origin (no look-ahead at any horizon).
  2. One predictive engine for every model and every table (mixture S=400 with
     intermediate-count simulation for dynamic models; plug-in for statics),
     so Tables 3, 4, 5 and 7 can no longer disagree.
  3. Score-functional pairing: MAE is reported for predictive MEDIANS at all
     horizons; means (and RMSE) only at h<=2 where Theorem 3(d) applies.
  4. CRPS and the multivariate energy score at h=1 (the log score is labelled
     a composite marginal score); MC standard errors from 5 seed replicates.
Outputs: out/unified.json + fragments tab_chicago_acc/pairs/pointfc + macros.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, math, sys, time
import numpy as np, pandas as pd, scipy.io as sio
from scipy.special import logsumexp
from scipy.optimize import minimize
sys.path.insert(0, _HERE)
from nssm import row_normalise, poisson_laplace_filter, poisson_logpmf
from nb_filter import nb_mode_filter, nb_logpmf

t0c = time.time()
OUT = _os.path.join(_HERE, "out"); TAB = _os.path.join(_ROOT, "paper", "tables")
df = pd.read_csv(_os.path.join(_ROOT, "data", "crime.csv"), index_col=0)
Y = df.values.T.astype(int); T, N = Y.shape
A = sio.mmread(_os.path.join(_ROOT, "data", "neighborhood.mtx")).tocsr()
A = ((A + A.T) > 0).astype(float).toarray(); np.fill_diagonal(A, 0)
W = row_normalise(A); Wc = (np.ones((N, N)) - np.eye(N)) / (N - 1)
ell = np.log1p(Y); HORIZONS = [1, 2, 4, 8]; TARGETS = list(range(60, 72))
Q_GRID = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]; S = 400

def build(kind, l_prev):
    if kind == "net":     return np.column_stack([np.ones(N), W @ l_prev, l_prev])
    if kind == "complete":return np.column_stack([np.ones(N), Wc @ l_prev, l_prev])
    return np.column_stack([np.ones(N), l_prev])

def tuned_q(kind, t_hi, cache={}):
    """q maximizing one-step plug-in log score on months 12..t_hi-1, data <= t_hi-1 only."""
    key = (kind, t_hi)
    if key in cache: return cache[key]
    K = 3 if kind != "nonet" else 2
    X_list = [build(kind, ell[t - 1]) for t in range(1, t_hi)]
    best = None
    for q in Q_GRID:
        fit = poisson_laplace_filter(Y[1:t_hi], X_list, q * np.eye(K), np.zeros(K), np.eye(K))
        ls = sum(poisson_logpmf(Y[t], np.exp(np.clip(X_list[t - 1] @ fit["m_pred"][t - 1], -30, 30))).sum()
                 for t in range(12, t_hi))
        if best is None or ls > best[0]: best = (ls, q)
    cache[key] = best[1]; return best[1]

def crps_from_draws(Yd, y):
    """Yd: S x N count draws. CRPS_i = E|Ys-y| - .5 E|Ys-Ys'| via random pairing."""
    Sd = Yd.shape[0]; perm = np.random.default_rng(0).permutation(Sd)
    return float((np.abs(Yd - y[None, :]).mean(0) - 0.5 * np.abs(Yd - Yd[perm]).mean(0)).sum())

def energy_from_draws(Yd, y):
    Sd = Yd.shape[0]; perm = np.random.default_rng(1).permutation(Sd)
    return float(np.linalg.norm(Yd - y[None, :], axis=1).mean()
                 - 0.5 * np.linalg.norm(Yd - Yd[perm], axis=1).mean())

def dyn_eval(kind, seed=5, frozen=False, hs=None):
    """Mixture (or frozen-path) predictive for one dynamic model at all origins/horizons."""
    rl = np.random.default_rng(seed)
    res = {h: dict(ls=[], medmae=[], meanmae=[], rmse=[], crps=[], es=[]) for h in HORIZONS}
    for h in (hs or HORIZONS):
        for tgt in TARGETS:
            t_or = tgt - h
            if t_or < 12: continue
            q = tuned_q(kind, t_or + 1)
            K = 3 if kind != "nonet" else 2
            X_list = [build(kind, ell[t - 1]) for t in range(1, t_or + 1)]
            fit = poisson_laplace_filter(Y[1:t_or + 1], X_list, q * np.eye(K), np.zeros(K), np.eye(K))
            m, P = fit["m_filt"][-1], fit["P_filt"][-1]
            if frozen:
                th = np.repeat(m[None, :], S, 0)
            else:
                th = m[None, :] + rl.standard_normal((S, K)) @ np.linalg.cholesky(
                    P + q * np.eye(K) + 1e-12 * np.eye(K)).T
            Lq = np.linalg.cholesky(q * np.eye(K) + 1e-14 * np.eye(K))
            y_prev = np.repeat(Y[t_or][None, :].astype(float), S, 0)
            for k in range(h):
                if k > 0 and not frozen: th = th + rl.standard_normal((S, K)) @ Lq.T
                lam = np.empty_like(y_prev)
                for s in range(S):
                    lam[s] = np.exp(np.clip(build(kind, np.log1p(y_prev[s])) @ th[s], -30, 30))
                if k < h - 1:
                    y_prev = rl.poisson(np.minimum(lam, 1e12)).astype(float)
            ytrue = Y[tgt]
            lamc = np.minimum(lam, 1e8)
            lp = poisson_logpmf(ytrue[None, :], lamc)
            res[h]["ls"].append(float((logsumexp(lp, axis=0) - math.log(S)).sum()))
            Yd = rl.poisson(np.minimum(lam, 1e12))
            res[h]["medmae"].append(float(np.mean(np.abs(ytrue - np.median(Yd, 0)))))
            if h <= 2:
                mm = lam.mean(0)
                res[h]["meanmae"].append(float(np.mean(np.abs(ytrue - mm))))
                res[h]["rmse"].append(float(np.sqrt(np.mean((ytrue - mm) ** 2))))
            if h == 1:
                res[h]["crps"].append(crps_from_draws(Yd, ytrue))
                res[h]["es"].append(energy_from_draws(Yd, ytrue))
    return res

def glm_pois(Xr, yv, b0):
    def nll(b):
        eta = np.clip(Xr @ b, -30, 30)
        return float(np.sum(np.exp(eta) - yv * eta))
    return minimize(nll, b0, method="L-BFGS-B").x

def static_eval(kind, seed=11):
    """Expanding-window static Poisson GLMs and the seasonal profile, all horizons."""
    rl = np.random.default_rng(seed)
    res = {h: dict(ls=[], medmae=[], meanmae=[], rmse=[], crps=[], es=[]) for h in HORIZONS}
    for h in HORIZONS:
        for tgt in TARGETS:
            t_or = tgt - h
            if t_or < 13: continue
            if kind == "season":
                mu1 = np.zeros(N)
                mth = tgt % 12
                sel = [t for t in range(t_or + 1) if t % 12 == mth]
                mu1 = Y[sel].mean(0) if sel else Y[:t_or + 1].mean(0)
                lam = np.repeat(np.maximum(mu1, 1e-3)[None, :], S, 0)  # same at all h
            else:
                cols = ([np.ones(t_or * N), (ell[:t_or] @ W.T).ravel(), ell[:t_or].ravel()]
                        if kind == "static_net" else [np.ones(t_or * N), ell[:t_or].ravel()])
                b = glm_pois(np.column_stack(cols), Y[1:t_or + 1].ravel(),
                             np.array([-0.3, 0.4, 0.5]) if kind == "static_net" else np.array([-0.3, 0.6]))
                def lam_of(l_prev):
                    e = (b[0] + b[1] * (W @ l_prev) + b[2] * l_prev) if kind == "static_net" \
                        else (b[0] + b[1] * l_prev)
                    return np.exp(np.clip(e, -30, 30))
                y_prev = np.repeat(Y[t_or][None, :].astype(float), S, 0)
                for k in range(h):
                    lam = np.vstack([lam_of(np.log1p(y_prev[s])) for s in range(S)])
                    if k < h - 1: y_prev = rl.poisson(np.minimum(lam, 1e12)).astype(float)
            ytrue = Y[tgt]
            lp = poisson_logpmf(ytrue[None, :], np.minimum(lam, 1e8))
            res[h]["ls"].append(float((logsumexp(lp, axis=0) - math.log(S)).sum()))
            Yd = rl.poisson(np.minimum(lam, 1e12))
            res[h]["medmae"].append(float(np.mean(np.abs(ytrue - np.median(Yd, 0)))))
            if h <= 2:
                mm = lam.mean(0)
                res[h]["meanmae"].append(float(np.mean(np.abs(ytrue - mm))))
                res[h]["rmse"].append(float(np.sqrt(np.mean((ytrue - mm) ** 2))))
            if h == 1:
                res[h]["crps"].append(crps_from_draws(Yd, ytrue))
                res[h]["es"].append(energy_from_draws(Yd, ytrue))
    return res

MODELS = {}
for kind, lab in [("net", "net"), ("nonet", "nonet"), ("complete", "complete")]:
    MODELS[lab] = dyn_eval(kind); print(f"[{lab}] done ({time.time()-t0c:.0f}s)", flush=True)
MODELS["net_frozen"] = dyn_eval("net", frozen=True)
MODELS["nonet_frozen"] = dyn_eval("nonet", frozen=True)
print(f"[frozen] done ({time.time()-t0c:.0f}s)", flush=True)
for kind in ["static_net", "static_nonet", "season"]:
    MODELS[kind] = static_eval(kind); print(f"[{kind}] done ({time.time()-t0c:.0f}s)", flush=True)

# MC standard error of the headline h=1 dLS from 5 seed replicates
dls_reps = []
for sd in range(5):
    a = dyn_eval("net", seed=100 + sd, hs=[1]); b = dyn_eval("nonet", seed=200 + sd, hs=[1])
    dls_reps.append(np.mean(a[1]["ls"]) - np.mean(b[1]["ls"]))
mcse = float(np.std(dls_reps, ddof=1) / math.sqrt(5))
print(f"[MC] dLS h1 reps: {[round(x,2) for x in dls_reps]}, MC-SE {mcse:.3f} ({time.time()-t0c:.0f}s)", flush=True)

# NB plug-in and mixture rows at h=1 (r tuned per origin on training window <= origin)
def nb_h1(kind, plug=False, seed=7):
    rl = np.random.default_rng(seed); ls = []
    for tgt in TARGETS:
        t_or = tgt - 1
        best = None
        for q in [1e-4, 3e-4, 1e-3]:
            for r in [2.0, 5.0, 10.0]:
                K = 3 if kind == "net" else 2
                X_list = [build(kind, ell[t - 1]) for t in range(1, t_or + 1)]
                fit = nb_mode_filter(Y[1:t_or + 1], X_list, q * np.eye(K), np.zeros(K), np.eye(K), r)
                v = sum(nb_logpmf(Y[t], np.exp(np.clip(X_list[t - 1] @ fit["m_pred"][t - 1], -30, 30)), r).sum()
                        for t in range(12, t_or + 1))
                if best is None or v > best[0]: best = (v, q, r, fit, X_list)
        _, q, r, fit, X_list = best
        Xn = build(kind, ell[t_or])
        if plug:
            lam = np.exp(np.clip(Xn @ fit["m_filt"][-1], -30, 30))
            ls.append(float(nb_logpmf(Y[tgt], lam, r).sum()))
        else:
            K = 3 if kind == "net" else 2
            th = fit["m_filt"][-1][None, :] + rl.standard_normal((S, K)) @ np.linalg.cholesky(
                fit["P_filt"][-1] + q * np.eye(K) + 1e-12 * np.eye(K)).T
            lam = np.exp(np.clip(th @ Xn.T, -30, 30))
            lp = nb_logpmf(Y[tgt][None, :], np.minimum(lam, 1e8), r)
            ls.append(float((logsumexp(lp, axis=0) - math.log(S)).sum()))
    return ls
nb_mix = nb_h1("net"); nb_plug = nb_h1("net", plug=True)
print(f"[NB] mixture {np.mean(nb_mix):.1f}, plug-in {np.mean(nb_plug):.1f} ({time.time()-t0c:.0f}s)", flush=True)

def bboot(diff, B=4000, L=3, seed=2):
    rl = np.random.default_rng(seed); n = len(diff); out = []
    for _ in range(B):
        idx = []
        while len(idx) < n:
            s0 = rl.integers(0, n); idx += [(s0 + k) % n for k in range(L)]
        out.append(np.mean(np.array(diff)[idx[:n]]))
    return float(np.quantile(out, .025)), float(np.quantile(out, .975))

def dm(diff, lag):
    d = np.array(diff); n = len(d); mu = d.mean()
    g0 = np.mean((d - mu) ** 2)
    var = g0
    for l in range(1, min(lag, n - 1) + 1):
        gl = np.mean((d[l:] - mu) * (d[:-l] - mu))
        var += 2 * (1 - l / (lag + 1)) * gl
    return float(mu / math.sqrt(max(var, 1e-12) / n))

# ---- Table 3 (accuracy) fragment: LS all h + median-MAE all h ----
order3 = [("net", "NSSM, observed $W$"), ("nonet", "NSSM, no network"),
          ("complete", "NSSM, complete graph $\\Wc$"), ("static_net", "static network GLM"),
          ("static_nonet", "static GLM, no network"), ("season", "seasonal profile")]
rows = []
for k, lab in order3:
    d = MODELS[k]
    ls = [f"{np.mean(d[h]['ls']):.1f}" for h in HORIZONS]
    mm = [f"{np.mean(d[h]['medmae']):.3f}" for h in HORIZONS]
    rows.append(f"{lab} & " + " & ".join(ls) + " & " + " & ".join(mm) + " \\\\")
open(f"{TAB}/tab_chicago_acc.tex", "w").write("\n".join(rows) + "\n")

# ---- Table 4 (pairs vs net) fragment ----
rows = []
for k, lab in order3[1:]:
    ent = []
    for h in HORIZONS:
        n_ = min(len(MODELS["net"][h]["ls"]), len(MODELS[k][h]["ls"]))
        dd = np.array(MODELS["net"][h]["ls"][:n_]) - np.array(MODELS[k][h]["ls"][:n_])
        lo, hi = bboot(dd)
        star = "$^{*}$" if lo > 0 or hi < 0 else ""
        dmtxt = f"{dm(dd, h - 1):.1f}" if h <= 4 else "--"
        ent.append(f"${np.mean(dd):+.1f}$ [{lo:+.1f},{hi:+.1f}]{star} ({dmtxt})")
    dcr = np.array(MODELS["net"][1]["crps"]) - np.array(MODELS[k][1]["crps"])
    ent.append(f"${np.mean(dcr):+.1f}$")
    rows.append(f"vs.\\ {lab} & " + " & ".join(ent) + " \\\\")
open(f"{TAB}/tab_chicago_pairs.tex", "w").write("\n".join(rows) + "\n")

# ---- Table 5 (point forecasts that exist), horizons {1,2,4,8} ----
rows = []
for h in HORIZONS:
    med_n = np.mean(MODELS["net"][h]["medmae"]); med_o = np.mean(MODELS["nonet"][h]["medmae"])
    frz_n = np.mean(MODELS["net_frozen"][h]["medmae"])  # frozen medians ~ frozen means proxy
    frn = np.mean(MODELS["net_frozen"][h].get("meanmae", [np.nan])) if h <= 2 else np.nan
    fro = np.mean(MODELS["nonet_frozen"][h].get("meanmae", [np.nan])) if h <= 2 else np.nan
    # frozen means exist at ALL h (Thm 3(c)); recompute them from frozen runs' medmae? need mean columns at all h:
    rows.append("PLACEHOLDER")
# recompute frozen means at all horizons properly
def frozen_means(kind, seed=9):
    rl = np.random.default_rng(seed); out = {h: [] for h in HORIZONS}
    for h in HORIZONS:
        for tgt in TARGETS:
            t_or = tgt - h
            if t_or < 12: continue
            q = tuned_q(kind, t_or + 1); K = 3 if kind != "nonet" else 2
            X_list = [build(kind, ell[t - 1]) for t in range(1, t_or + 1)]
            fit = poisson_laplace_filter(Y[1:t_or + 1], X_list, q * np.eye(K), np.zeros(K), np.eye(K))
            m = fit["m_filt"][-1]
            y_prev = np.repeat(Y[t_or][None, :].astype(float), S, 0)
            for k in range(h):
                lam = np.vstack([np.exp(np.clip(build(kind, np.log1p(y_prev[s])) @ m, -30, 30)) for s in range(S)])
                if k < h - 1: y_prev = rl.poisson(np.minimum(lam, 1e12)).astype(float)
            out[h].append(float(np.mean(np.abs(Y[tgt] - lam.mean(0)))))
    return out
fm_net = frozen_means("net"); fm_non = frozen_means("nonet")
print(f"[frozen means] done ({time.time()-t0c:.0f}s)", flush=True)
rows = []
for h in HORIZONS:
    med_n = np.mean(MODELS["net"][h]["medmae"]); med_o = np.mean(MODELS["nonet"][h]["medmae"])
    fzn, fzo = np.mean(fm_net[h]), np.mean(fm_non[h])
    if h <= 2:
        sn = np.mean(MODELS["net"][h]["meanmae"]); so = np.mean(MODELS["nonet"][h]["meanmae"])
        samp = f"{sn:.3f} & {so:.3f} & ${sn-so:+.4f}$"
    else:
        samp = "\\multicolumn{3}{c}{--- (no population mean, Thm.~3(e))}"
    rows.append(f"{h} & {med_n:.3f} & {med_o:.3f} & ${med_n-med_o:+.4f}$ & "
                f"{fzn:.3f} & {fzo:.3f} & ${fzn-fzo:+.4f}$ & {samp} \\\\")
open(f"{TAB}/tab_chicago_pointfc.tex", "w").write("\n".join(rows) + "\n")

json.dump(dict(
    models={k: {str(h): {kk: list(map(float, vv)) for kk, vv in v[h].items()} for h in HORIZONS}
            for k, v in MODELS.items()},
    frozen_means_net={str(h): list(map(float, fm_net[h])) for h in HORIZONS},
    frozen_means_nonet={str(h): list(map(float, fm_non[h])) for h in HORIZONS},
    nb_mix=list(map(float, nb_mix)), nb_plug=list(map(float, nb_plug)),
    mcse_h1=mcse, dls_reps=list(map(float, dls_reps)),
    tuned_q_per_origin={f"{k}_{t}": tuned_q(k, t) for k in ["net", "nonet"] for t in range(60, 72)},
), open(f"{OUT}/unified.json", "w"))

for k, lab in order3:
    d = MODELS[k]
    print(f"{lab:28s} LS: " + " ".join(f"{np.mean(d[h]['ls']):7.1f}" for h in HORIZONS)
          + "  medMAE: " + " ".join(f"{np.mean(d[h]['medmae']):.3f}" for h in HORIZONS), flush=True)
print("CRPS h1 net/nonet:", round(np.mean(MODELS['net'][1]['crps']),1), round(np.mean(MODELS['nonet'][1]['crps']),1))
print("ES h1 net/nonet:", round(np.mean(MODELS['net'][1]['es']),1), round(np.mean(MODELS['nonet'][1]['es']),1))
print("DONE", round(time.time() - t0c), "s", flush=True)
