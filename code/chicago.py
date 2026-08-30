"""
Chicago burglary application (Poisson NSSM).  Reproducible end-to-end:
data from github.com/nick3703/Chicago-Data (552 regions x 72 months, 2010-2015).
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, time, sys, os
import numpy as np, pandas as pd, scipy.io as sio
from scipy.stats import chisquare
from scipy.special import logsumexp
sys.path.insert(0, _HERE)
from nssm import (row_normalise, poisson_laplace_filter, poisson_predict, mixture_scores,
                  poisson_glm, two_way_loglinear, profile_information, graph_roughness,
                  poisson_logpmf, op_norm, tau_W)

t_start = time.time()
OUT = _os.path.join(_HERE, "out"); os.makedirs(OUT, exist_ok=True)

# ---------------- data ----------------
df = pd.read_csv(_os.path.join(_ROOT, "data", "crime.csv"), index_col=0)
Y = df.values.T.astype(int)                        # T x N
T, N = Y.shape
A = sio.mmread(_os.path.join(_ROOT, "data", "neighborhood.mtx")).tocsr()
A = ((A + A.T) > 0).astype(float).toarray(); np.fill_diagonal(A, 0)
W = row_normalise(A)
Wc = np.ones((N, N)) / (N - 1); np.fill_diagonal(Wc, 0)      # complete graph
ell = np.log1p(Y)
TARGETS = list(range(60, 72))                        # last 12 months are the targets
HORIZONS = [1, 2, 4, 8]
S_MAIN, S_STRESS = 400, 200
Q_GRID = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]
deg = A.sum(1)
data_summary = dict(N=int(N), T=int(T), n_edges=int(A.sum() / 2), deg_min=float(deg.min()),
                    deg_mean=float(deg.mean()), deg_max=float(deg.max()),
                    mean_count=float(Y.mean()), zero_frac=float((Y == 0).mean()),
                    max_count=int(Y.max()), monthly_total_first=int(Y[0].sum()), monthly_total_last=int(Y[-1].sum()),
                    tau_W=tau_W(W), tau_Wc=tau_W(Wc), tau_regular_bound=float(np.sum(1 / deg)))
print("data:", data_summary)

def make_builder(Wmat, kind):
    def build(l_prev):
        if kind == "net":
            return np.column_stack([np.ones(N), Wmat @ l_prev, l_prev])
        if kind == "nonet":
            return np.column_stack([np.ones(N), l_prev])
        if kind == "complete":
            return np.column_stack([np.ones(N), Wc @ l_prev, l_prev])
        raise ValueError(kind)
    return build

MODELS = {"net": 3, "nonet": 2, "complete": 3}

def designs(builder, transform=np.log1p):
    Z = transform(Y)
    return [None] + [builder(Z[t - 1]) for t in range(1, T)]

def run_filter(builder, q, K, transform=np.log1p):
    X_list = designs(builder, transform)
    fit = poisson_laplace_filter(Y[1:], X_list[1:], Q=q * np.eye(K), m0=np.zeros(K), P0=np.eye(K))
    return fit, X_list

def plugin_logscore_train(fit, X_list, t_lo=12, t_hi=60):
    ls = 0.0
    for t in range(t_lo, t_hi):
        eta = np.clip(X_list[t] @ fit["m_pred"][t - 1], -30, 30)
        ls += poisson_logpmf(Y[t], np.exp(eta)).sum()
    return ls

def tune_q(builder, K, transform=np.log1p):
    best, table = None, {}
    for q in Q_GRID:
        fit, X_list = run_filter(builder, q, K, transform)
        ls = plugin_logscore_train(fit, X_list)
        table[q] = ls
        if best is None or ls > best[1]:
            best = (q, ls)
    return best[0], table

fits, tuned, tune_tables = {}, {}, {}
for name, K in MODELS.items():
    b = make_builder(W, name)
    q, tab = tune_q(b, K)
    tuned[name] = q; tune_tables[name] = tab
    fits[name], _ = run_filter(b, q, K)
    print(f"[tune] {name}: q={q:g}  train plug-in LS by q: " + ", ".join(f"{k:g}:{v:.1f}" for k, v in tab.items()))

# ---------------- rolling evaluation (fixed targets, origin = target - h) ----------------
def evaluate_dynamic(name, Wmat, q, S, horizons, seed, pit_h1=False, transform=np.log1p, light=False):
    builder = make_builder(Wmat, name)
    K = MODELS[name]
    fit, X_list = run_filter(builder, q, K, transform)
    rng_l = np.random.default_rng(seed)
    res = {h: [] for h in horizons}; pits = []
    for h in horizons:
        for tau in TARGETS:
            t0 = tau - h
            m, P = fit["m_filt"][t0 - 1], fit["P_filt"][t0 - 1]
            lam, coef = poisson_predict(m, P, q * np.eye(K), h, S, rng_l, builder, Y[t0], transform=transform)
            explosive = float(np.mean(np.any(lam > 1e4, axis=1)))
            rho = coef[:, :, 1] + coef[:, :, -1] if name in ("net", "complete") else coef[:, :, 1]
            row = dict(target=int(tau), origin=int(t0), explosive=explosive,
                       max_lam=float(lam.max()), rho_ge1=float(np.mean(np.any(rho >= 1.0, axis=0))))
            if light:
                lp = poisson_logpmf(Y[tau][None, :], np.minimum(lam, 1e12))
                mean = lam.mean(0); err = np.abs(Y[tau] - mean)
                row.update(mae=float(err.mean()), medae=float(np.median(err)), mse=float((err ** 2).mean()),
                           ls=float((logsumexp(lp, axis=0) - np.log(S)).sum()), cov=float("nan"))
            else:
                sc = mixture_scores(Y[tau], lam, rng_l)
                err = np.abs(Y[tau] - sc["mean"])
                row.update(mae=float(err.mean()), medae=float(np.median(err)), mse=float((err ** 2).mean()),
                           ls=sc["logscore_sum"], cov=sc["coverage"], err_by_region=err.tolist())
                if h == 1 and pit_h1:
                    pits.append(sc["pit"])
            res[h].append(row)
    return fit, res, (np.concatenate(pits) if pits else None)

results, betas, pit_store = {}, {}, {}
for name in MODELS:
    fit, res, pit = evaluate_dynamic(name, W, tuned[name], S_MAIN, HORIZONS, seed=1, pit_h1=True)
    results[name] = res; pit_store[name] = pit
    betas[name] = dict(m=fit["m_filt"].tolist(), sd=np.sqrt(np.einsum("tkk->tk", fit["P_filt"])).tolist(),
                       info=fit["info"].tolist())
    print(f"[eval] {name} done ({time.time()-t_start:.0f}s)")

# ---------------- raw-count specification (Theorem 4 diagnostic) ----------------
ident = lambda x: x
q_raw, tab_raw = tune_q(make_builder(W, "net"), 3, transform=ident)
fit_raw, res_raw, _ = evaluate_dynamic("net", W, q_raw, S_MAIN, [1, 2, 4], seed=5, transform=ident, light=True)
results["raw_net"] = res_raw
betas["raw_net"] = dict(m=fit_raw["m_filt"].tolist(), sd=np.sqrt(np.einsum("tkk->tk", fit_raw["P_filt"])).tolist())
print(f"[raw] q={q_raw:g} done ({time.time()-t_start:.0f}s)")

# S-scan: Monte Carlo predictive mean of the city total at origin 59 (target 59+h) as S grows
def s_scan(name, q, transform, t0, hs, S_list, seed=9):
    builder = make_builder(W, name); K = MODELS[name]
    fit, _ = run_filter(builder, q, K, transform)
    m, P = fit["m_filt"][t0 - 1], fit["P_filt"][t0 - 1]
    out = {}
    for h in hs:
        out[h] = {}
        for S in S_list:
            r = np.random.default_rng(seed)
            lam, _ = poisson_predict(m, P, q * np.eye(K), h, S, r, builder, Y[t0], transform=transform)
            tot = lam.sum(1)
            out[h][S] = dict(mean_total=float(tot.mean()), median_total=float(np.median(tot)),
                             max_total=float(tot.max()), frac_gt_1e4=float(np.mean(tot > 1e4)))
    return out
sscan = dict(raw=s_scan("net", q_raw, ident, 59, [1, 2, 3, 4], [100, 1000, 5000]),
             log1p=s_scan("net", tuned["net"], np.log1p, 59, [1, 2, 3, 4], [100, 1000, 5000]))
print("[sscan] done", f"({time.time()-t_start:.0f}s)")

# ---------------- static baselines (expanding-window refit at each origin) ----------------
def evaluate_static(kind, S, horizons, seed):
    rng_l = np.random.default_rng(seed)
    res = {h: [] for h in horizons}
    b = make_builder(W, kind) if kind != "season" else None
    Xs = designs(b) if kind != "season" else None
    for h in horizons:
        for tau in TARGETS:
            t0 = tau - h
            if kind == "season":
                a, bm = two_way_loglinear(Y[:t0 + 1])
            else:
                Xtr = np.vstack(Xs[1:t0 + 1]); ytr = Y[1:t0 + 1].ravel()
                beta = poisson_glm(Xtr, ytr)
            y_prev = np.repeat(Y[t0][None, :], S, axis=0).astype(float)
            for k in range(h):
                if kind == "season":
                    lam = np.exp(a + bm[(t0 + k + 1) % 12])[None, :].repeat(S, 0)
                else:
                    lam = np.vstack([np.exp(np.clip(b(np.log1p(y_prev[s])) @ beta, -30, 30)) for s in range(S)])
                if k < h - 1:
                    y_prev = rng_l.poisson(lam).astype(float)
            sc = mixture_scores(Y[tau], lam, rng_l)
            err = np.abs(Y[tau] - sc["mean"])
            res[h].append(dict(target=int(tau), origin=int(t0), mae=float(err.mean()), medae=float(np.median(err)),
                               mse=float((err ** 2).mean()), ls=sc["logscore_sum"], cov=sc["coverage"],
                               explosive=0.0, rho_ge1=float("nan"), err_by_region=err.tolist()))
    return res

for kind in ["net", "nonet", "season"]:
    results["static_" + kind] = evaluate_static(kind, S_MAIN, HORIZONS, seed=2)
    print(f"[static] {kind} done ({time.time()-t_start:.0f}s)")

# ---------------- paired comparisons: block bootstrap over targets ----------------
def block_boot_ci(d, B=4000, block=3, seed=0):
    d = np.asarray(d, float); n = len(d); r = np.random.default_rng(seed)
    nb = int(np.ceil(n / block)); starts = np.arange(0, n - block + 1)
    means = np.empty(B)
    for b in range(B):
        idx = np.concatenate([np.arange(s, s + block) for s in r.choice(starts, nb)])[:n]
        means[b] = d[idx].mean()
    return float(d.mean()), float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))

def dm_stat(d, lag):
    d = np.asarray(d, float); n = len(d); x = d - d.mean()
    v = (x @ x) / n
    for k in range(1, min(lag, n - 1) + 1):
        v += 2 * (1 - k / (lag + 1)) * (x[k:] @ x[:-k]) / n
    se = np.sqrt(max(v, 1e-300) / n)
    return float(d.mean() / se)

def summarise(res):
    out = {}
    for h, rows in res.items():
        out[h] = {k: float(np.nanmean([r[k] for r in rows])) for k in ["mae", "medae", "mse", "ls", "cov", "explosive", "rho_ge1", "max_lam"] if k in rows[0]}
        out[h]["n"] = len(rows)
    return out

summary = {k: summarise(v) for k, v in results.items()}
pairs = {}
for other in ["nonet", "complete", "static_net", "static_nonet", "static_season"]:
    key = "net_vs_" + other; pairs[key] = {}
    for h in HORIZONS:
        a = results["net"][h]; b = results[other][h]
        d_mae = [ra["mae"] - rb["mae"] for ra, rb in zip(a, b)]
        d_ls = [ra["ls"] - rb["ls"] for ra, rb in zip(a, b)]
        m, lo, hi = block_boot_ci(d_mae); ml, llo, lhi = block_boot_ci(d_ls, seed=1)
        pairs[key][h] = dict(dmae=m, dmae_ci=[lo, hi], dmae_dm=dm_stat(d_mae, h - 1),
                             dls=ml, dls_ci=[llo, lhi], dls_dm=dm_stat(d_ls, h - 1),
                             prop_targets_better_ls=float(np.mean(np.array(d_ls) > 0)))

E_net = np.array([r["err_by_region"] for r in results["net"][1]]).mean(0)
E_non = np.array([r["err_by_region"] for r in results["nonet"][1]]).mean(0)
dreg = E_net - E_non
regionwise = dict(q10=float(np.quantile(dreg, .1)), q25=float(np.quantile(dreg, .25)), median=float(np.median(dreg)),
                  q75=float(np.quantile(dreg, .75)), q90=float(np.quantile(dreg, .9)), prop_better=float(np.mean(dreg < 0)))

pit_tests = {}
for name in MODELS:
    p = pit_store[name]
    counts, _ = np.histogram(p, bins=10, range=(0, 1))
    stat, pval = chisquare(counts)
    pit_tests[name] = dict(counts=counts.tolist(), chi2=float(stat), pval=float(pval), n=int(len(p)))

# ---------------- identification diagnostics (Theorem 1) ----------------
identd = dict(roughness=[], info_net=[], info_complete=[], lam_min_net=[], month=[], noise_info=[])
fit_net = fits["net"]; bnet = make_builder(W, "net"); bc = make_builder(W, "complete")
for t in range(1, T):
    X = bnet(ell[t - 1]); Xc = bc(ell[t - 1])
    lam = np.exp(np.clip(X @ fit_net["m_filt"][t - 1], -30, 30))
    i_net, lmin = profile_information(X, weights=lam, col=1)
    i_c, _ = profile_information(Xc, weights=lam, col=1)
    identd["roughness"].append(graph_roughness(W, ell[t - 1]) / N)
    identd["info_net"].append(i_net); identd["info_complete"].append(i_c); identd["lam_min_net"].append(lmin)
    identd["month"].append(t)

# ---------------- network stress / placebo tests (h = 1, 2) ----------------
def rewire_degree_preserving(Amat, n_swaps, r):
    A2 = Amat.copy(); edges = np.argwhere(np.triu(A2, 1) > 0)
    E = len(edges); done = 0; tries = 0
    while done < n_swaps and tries < 50 * n_swaps:
        tries += 1
        i1, i2 = r.choice(E, 2, replace=False)
        a, b = edges[i1]; c, d = edges[i2]
        if len({a, b, c, d}) < 4 or A2[a, d] or A2[c, b]:
            continue
        A2[a, b] = A2[b, a] = 0; A2[c, d] = A2[d, c] = 0
        A2[a, d] = A2[d, a] = 1; A2[c, b] = A2[b, c] = 1
        edges[i1] = [a, d]; edges[i2] = [c, b]; done += 1
    return A2

stress, perturbations = {}, {}
r = np.random.default_rng(7)
perturbations["original"] = W
for a_ in [0.25, 0.5, 0.75, 1.0]:
    perturbations[f"mix_uniform_{a_}"] = (1 - a_) * W + a_ * Wc
A2 = A.copy(); edges = np.argwhere(np.triu(A2, 1) > 0)
drop = r.choice(len(edges), int(0.2 * len(edges)), replace=False)
for e in edges[drop]:
    A2[e[0], e[1]] = A2[e[1], e[0]] = 0
perturbations["edge_delete_0.2"] = row_normalise(A2)
perturbations["rewire_degseq"] = row_normalise(rewire_degree_preserving(A, int(2 * A.sum() / 2), r))
for p in range(3):
    perm = r.permutation(N)
    perturbations[f"permute_labels_{p}"] = W[np.ix_(perm, perm)]

for pname, Wp in perturbations.items():
    q = tuned["net"]
    fitp, resp, _ = evaluate_dynamic("net", Wp, q, S_STRESS, [1, 2], seed=11)
    row = dict(opnorm_diff=float(op_norm(Wp - W)), tau=tau_W(Wp))
    for h in [1, 2]:
        a = resp[h]; b = results["nonet"][h]
        d_ls = [ra["ls"] - rb["ls"] for ra, rb in zip(a, b)]; d_mae = [ra["mae"] - rb["mae"] for ra, rb in zip(a, b)]
        row[f"dmae_h{h}"] = float(np.mean(d_mae)); row[f"dls_h{h}"] = float(np.mean(d_ls))
        row[f"dls_h{h}_ci"] = list(block_boot_ci(d_ls, seed=3)[1:])
    bp = make_builder(Wp, "net"); infos = []
    for t in range(1, T):
        X = bp(ell[t - 1]); lam = np.exp(np.clip(X @ fitp["m_filt"][t - 1], -30, 30))
        infos.append(profile_information(X, weights=lam, col=1)[0])
    row["mean_profile_info"] = float(np.mean(infos))
    row["mean_beta1"] = float(fitp["m_filt"][59:, 1].mean()); row["mean_sd_beta1"] = float(np.sqrt(fitp["P_filt"][59:, 1, 1]).mean())
    stress[pname] = row
    print(f"[stress] {pname}: " + ", ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in row.items()))

out = dict(data=data_summary, targets=TARGETS, horizons=HORIZONS, S_main=S_MAIN, S_stress=S_STRESS,
           tuned_q=tuned, q_raw=q_raw, tune_tables={k: {str(q): v for q, v in tab.items()} for k, tab in tune_tables.items()},
           tune_table_raw={str(q): v for q, v in tab_raw.items()},
           summary=summary, pairs=pairs, regionwise=regionwise, pit_tests=pit_tests, ident=identd, stress=stress,
           betas=betas, sscan={k: {str(h): {str(S): v for S, v in d.items()} for h, d in dd.items()} for k, dd in sscan.items()},
           per_target={k: {str(h): [{kk: vv for kk, vv in r_.items() if kk != "err_by_region"} for r_ in rows] for h, rows in v.items()} for k, v in results.items()},
           runtime_sec=time.time() - t_start)
with open(f"{OUT}/chicago_results.json", "w") as f:
    json.dump(out, f)
print("TOTAL RUNTIME", round(time.time() - t_start), "s")
