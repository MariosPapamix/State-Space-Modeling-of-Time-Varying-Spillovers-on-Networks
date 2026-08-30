"""
Second application: weekly county-level COVID-19 case counts, Texas (JHU CSSE),
k-nearest-neighbour contiguity graph built from the county coordinates shipped
with the data. Replicates the paper's pipeline in a different domain and count
regime: identification diagnostics, NB-NSSM net vs no-network forecasting,
exact randomization test, and the equivalence-inflation check.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, time, os, sys, math, urllib.request
import numpy as np, pandas as pd
from scipy.special import gammaln, logsumexp
sys.path.insert(0, _HERE)
from nssm import row_normalise, profile_information
from nb_filter import nb_mode_filter, nb_logpmf

t0c = time.time()
OUT = _os.path.join(_HERE, "out"); TAB = _os.path.join(_ROOT, "paper", "tables")
DATA = _os.path.join(_ROOT, "data", "covid_us.csv")
if not os.path.exists(DATA):
    u = ("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/"
         "csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
    urllib.request.urlretrieve(u, DATA)
    print(f"[covid] downloaded {os.path.getsize(DATA)/1e6:.1f} MB ({time.time()-t0c:.0f}s)", flush=True)

df = pd.read_csv(DATA)
tx = df[(df.Province_State == "Texas") & (df.Lat != 0) & df.Admin2.notna()
        & (~df.Admin2.isin(["Unassigned", "Out of TX"]))].reset_index(drop=True)
datecols = [c for c in tx.columns if "/" in c]
cum = tx[datecols].values.astype(float)
daily = np.maximum(np.diff(cum, axis=1), 0.0)
dates = pd.to_datetime(datecols[1:], format="%m/%d/%y")
start = np.searchsorted(dates, pd.Timestamp("2020-03-02"))
daily = daily[:, start:]; dates = dates[start:]
n_weeks = daily.shape[1] // 7
Yw = daily[:, :n_weeks * 7].reshape(len(tx), n_weeks, 7).sum(2)
T_use = min(92, n_weeks)
Y = np.round(Yw[:, :T_use]).astype(int).T                      # T x N
T, N = Y.shape
lat, lon = tx.Lat.values, tx.Long_.values
print(f"[covid] N={N} counties, T={T} weeks, mean count={Y.mean():.1f}, max={Y.max()}", flush=True)

# k-nearest-neighbour graph (k=5), symmetrised
K_NN = 5
D2 = (lat[:, None] - lat[None, :])**2 + (np.cos(np.radians(31)) * (lon[:, None] - lon[None, :]))**2
np.fill_diagonal(D2, np.inf)
A = np.zeros((N, N))
for i in range(N):
    A[i, np.argsort(D2[i])[:K_NN]] = 1
A = np.maximum(A, A.T); W = row_normalise(A)
tauW = float(np.trace(W.T @ (np.eye(N) - np.ones((N, N)) / N) @ W))
TRAIN = 72; ORIG = list(range(TRAIN - 1, T - 1)); S = 300
ell = np.log1p(Y)
print(f"[covid] tau_W={tauW:.1f} (complete graph: {N/(N-1)-1:.4f}; sum 1/n_i={np.sum(1/A.sum(1)):.1f})", flush=True)

def build(kind, l_prev):
    return (np.column_stack([np.ones(N), W @ l_prev, l_prev]) if kind == "net"
            else np.column_stack([np.ones(N), l_prev]))

def run_nb(kind, q, r, t_hi=None):
    t_hi = t_hi or T
    K = 3 if kind == "net" else 2
    X_list = [build(kind, ell[t - 1]) for t in range(1, t_hi)]
    fit = nb_mode_filter(Y[1:t_hi], X_list, q * np.eye(K), np.zeros(K), np.eye(K), r)
    return fit, X_list

def train_score(kind, q, r):
    fit, X_list = run_nb(kind, q, r, TRAIN)
    return float(sum(nb_logpmf(Y[t], np.exp(np.clip(X_list[t - 1] @ fit["m_pred"][t - 1], -30, 30)), r).sum()
                     for t in range(10, TRAIN)))

tuned = {}
for kind in ["net", "nonet"]:
    best = None
    for q in [1e-4, 3e-4, 1e-3, 3e-3]:
        for r in [2.0, 5.0, 10.0, 20.0]:
            v = train_score(kind, q, r)
            if best is None or v > best[0]: best = (v, q, r)
    tuned[kind] = (best[1], best[2])
    print(f"[covid] {kind}: q={best[1]:g} r={best[2]:g} ({time.time()-t0c:.0f}s)", flush=True)

res = {}
for kind in ["net", "nonet"]:
    q, r = tuned[kind]; K = 3 if kind == "net" else 2
    fit, X_list = run_nb(kind, q, r)
    rl = np.random.default_rng(9)
    ls_rows, mae_rows, cov_rows = [], [], []
    for t_or in ORIG:
        idx = t_or - 1
        m, P = fit["m_filt"][idx], fit["P_filt"][idx]
        th = m[None, :] + rl.standard_normal((S, K)) @ np.linalg.cholesky(P + q * np.eye(K) + 1e-12 * np.eye(K)).T
        mu = np.exp(np.clip(th @ build(kind, ell[t_or]).T, -30, 30))
        y = Y[t_or + 1]
        lp = nb_logpmf(y[None, :], np.minimum(mu, 1e8), r)
        ls_rows.append(float((logsumexp(lp, axis=0) - math.log(S)).sum()))
        mae_rows.append(float(np.mean(np.abs(y - np.median(mu, 0)))))
    if kind == "net":
        b1 = fit["m_filt"][:, 1]; sd1 = np.sqrt(fit["P_filt"][:, 1, 1])
        wts = None
        infos = []
        for t in range(9, T - 1):
            mu_t = np.exp(np.clip(X_list[t] @ fit["m_filt"][t], -30, 30))
            lamw = r * mu_t / (r + mu_t)
            infos.append(profile_information(X_list[t], lamw, 1)[0])
        res["beta1"] = dict(first=float(np.mean(b1[9:19])), last=float(np.mean(b1[-10:])),
                            minv=float(b1[9:].min()), maxv=float(b1[9:].max()),
                            sd=float(np.mean(sd1[9:])))
        res["info_mean"] = float(np.mean(infos)); res["info_min"] = float(np.min(infos)); res["info_max"] = float(np.max(infos))
    res[kind] = dict(q=q, r=r, ls=ls_rows, mae=mae_rows)
    print(f"[covid] {kind}: LS1={np.mean(ls_rows):.1f} medMAE={np.mean(mae_rows):.1f} ({time.time()-t0c:.0f}s)", flush=True)

dls = np.array(res["net"]["ls"]) - np.array(res["nonet"]["ls"])
def block_boot_ci(diff, B=4000, L=3, seed=2):
    rl = np.random.default_rng(seed); n = len(diff); out = []
    for _ in range(B):
        idx = []
        while len(idx) < n:
            s0 = rl.integers(0, n); idx += [(s0 + k) % n for k in range(L)]
        out.append(np.mean(np.array(diff)[idx[:n]]))
    return float(np.quantile(out, .025)), float(np.quantile(out, .975))
lo, hi = block_boot_ci(dls)
dm = float(np.mean(dls) / (np.std(dls, ddof=1) / math.sqrt(len(dls))))
res["dls"] = dict(mean=float(np.mean(dls)), lo=lo, hi=hi, dm=dm, n=len(dls))
print(f"[covid] dLS net-nonet: {np.mean(dls):+.1f} [{lo:+.1f},{hi:+.1f}] DM={dm:.2f}", flush=True)

# exact randomization test (99 label permutations; plug-in NB test-window statistic)
def stat_of(Wm):
    K = 3
    X_list = [np.column_stack([np.ones(N), Wm @ ell[t - 1], ell[t - 1]]) for t in range(1, T)]
    q, r = tuned["net"]
    fit = nb_mode_filter(Y[1:], X_list, q * np.eye(K), np.zeros(K), np.eye(K), r)
    return float(sum(nb_logpmf(Y[t], np.exp(np.clip(X_list[t - 1] @ fit["m_pred"][t - 1], -30, 30)), r).sum()
                     for t in range(TRAIN, T)))
q, r = tuned["nonet"]
Xn = [np.column_stack([np.ones(N), ell[t - 1]]) for t in range(1, T)]
fitn = nb_mode_filter(Y[1:], Xn, q * np.eye(2), np.zeros(2), np.eye(2), r)
base = float(sum(nb_logpmf(Y[t], np.exp(np.clip(Xn[t - 1] @ fitn["m_pred"][t - 1], -30, 30)), r).sum()
                 for t in range(TRAIN, T)))
s_obs = stat_of(W) - base
r2 = np.random.default_rng(7); perm_stats = []
for p in range(99):
    pm = r2.permutation(N)
    perm_stats.append(stat_of(W[np.ix_(pm, pm)]) - base)
    if p % 25 == 0: print(f"[covid] perm {p} ({time.time()-t0c:.0f}s)", flush=True)
perm_stats = np.array(perm_stats)
p_perm = (1 + int(np.sum(perm_stats >= s_obs))) / 100
res["randtest"] = dict(stat=s_obs, perm_max=float(perm_stats.max()), p=p_perm)
print(f"[covid] randtest: obs={s_obs:+.1f} perm_max={perm_stats.max():+.1f} p={p_perm:.3f}", flush=True)

# equivalence inflation at alpha=0.5
Wc = (np.ones((N, N)) - np.eye(N)) / (N - 1)
Wa = 0.5 * W + 0.5 * Wc
q, r = tuned["net"]
Xa = [np.column_stack([np.ones(N), Wa @ ell[t - 1], ell[t - 1]]) for t in range(1, T)]
fita = nb_mode_filter(Y[1:], Xa, q * np.eye(3), np.zeros(3), np.eye(3), r)
fit0, _ = run_nb("net", q, r)
b_obs = float(np.mean(fit0["m_filt"][9:, 1])); b_mix = float(np.mean(fita["m_filt"][9:, 1]))
res["equiv"] = dict(b_obs=b_obs, b_mix=b_mix, ratio=b_mix / b_obs)
print(f"[covid] equivalence: beta1 {b_obs:.3f} -> {b_mix:.3f} (ratio {b_mix/b_obs:.2f}, predicted 2.00)", flush=True)

json.dump(res, open(f"{OUT}/covid_results.json", "w"))
rows = [
    f"nodes $N$ / weeks $T$ / test targets & {N} / {T} / {len(ORIG)} \\\\",
    f"$\\tau_W$ (KNN graph) / $\\sum_i n_i^{{-1}}$ & {tauW:.1f} / {np.sum(1/A.sum(1)):.1f} \\\\",
    f"weighted information $\\bar{{\\calI}}_t$ (min--max) & {res['info_mean']:.0f} ({res['info_min']:.0f}--{res['info_max']:.0f}) \\\\",
    f"filtered $\\hat\\beta_{{1,t}}$: range (mean s.d.) & {res['beta1']['minv']:.2f}--{res['beta1']['maxv']:.2f} ({res['beta1']['sd']:.3f}) \\\\",
    f"$\\Delta$LS ($h=1$, net $-$ no-net), 95\\% CI, DM & ${res['dls']['mean']:+.1f}$ [{lo:+.1f}, {hi:+.1f}], {dm:.1f} \\\\",
    f"randomization $p$ (99 label permutations) & {p_perm:.3f} \\\\",
    f"$\\hat\\beta_1$ under $0.5W+0.5W_c$: observed / predicted & {b_mix:.2f} / {2*b_obs:.2f} \\\\",
]
open(f"{TAB}/tab_covid.tex", "w").write("\n".join(rows) + "\n")
print("DONE", round(time.time() - t0c), "s", flush=True)
