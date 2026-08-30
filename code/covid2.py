"""
Texas COVID replication, extended to a full parallel pipeline:
  - county POPULATION offsets (JHU UID lookup table, downloaded at run time);
  - external benchmarks on this second dataset too (PNAR linear/log-linear,
    hhh4-type NB with unit intercepts [which absorb population exposure],
    EWMA-NB), one protocol, h=1 and h=2 where the predictive exists;
  - median-MAE alongside the composite log score; PIT (orientation);
  - conditional randomization test stratified by population quintiles with
    (q, r) retuned inside every permutation; unrestricted-retuned comparison;
  - k=8 nearest-neighbour graph sensitivity.
"""
import os as _os, urllib.request
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, math, sys, time
import numpy as np, pandas as pd
from scipy.special import logsumexp, gammaln
from scipy.optimize import minimize
sys.path.insert(0, _HERE)
from nssm import row_normalise, poisson_logpmf
from nb_filter import nb_mode_filter, nb_logpmf

t0 = time.time(); OUT = _os.path.join(_HERE, "out"); TAB = _os.path.join(_ROOT, "paper", "tables")
DATA = _os.path.join(_ROOT, "data")
CSV = _os.path.join(DATA, "covid_us.csv")
if not _os.path.exists(CSV):
    urllib.request.urlretrieve("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/"
        "csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv", CSV)
LOOK = _os.path.join(DATA, "uid_lookup.csv")
if not _os.path.exists(LOOK):
    urllib.request.urlretrieve("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/"
        "csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv", LOOK)
raw = pd.read_csv(CSV); lk = pd.read_csv(LOOK)
tx = raw[(raw.Province_State == "Texas") & raw.Admin2.notna() & (raw.Lat != 0)].copy()
tx = tx.merge(lk[["FIPS", "Population"]], on="FIPS", how="left")
tx = tx[tx.Population > 0].reset_index(drop=True)
datecols = [c for c in tx.columns if c.count("/") == 2]
cum = tx[datecols].values.astype(float)
dates = pd.to_datetime(datecols, format="%m/%d/%y")
start = np.searchsorted(dates, pd.Timestamp("2020-03-01"))
new = np.clip(np.diff(cum[:, start - 1:], axis=1), 0, None)
nweeks = new.shape[1] // 7
wk = new[:, :nweeks * 7].reshape(len(tx), nweeks, 7).sum(2)
Tfull = min(nweeks, 92); Yw = wk[:, :Tfull].T.astype(int)  # T x N
N = Yw.shape[1]; T = Yw.shape[0]
pop = tx.Population.values.astype(float)
lat, lon = tx.Lat.values, tx.Long_.values
def knn_graph(k):
    D = np.sqrt((lat[:, None] - lat[None, :]) ** 2 + (lon[:, None] - lon[None, :]) ** 2)
    np.fill_diagonal(D, np.inf)
    A = np.zeros((N, N))
    for i in range(N):
        A[i, np.argsort(D[i])[:k]] = 1
    A = ((A + A.T) > 0).astype(float)
    return row_normalise(A)
W5, W8 = knn_graph(5), knn_graph(8)
ell = np.log1p(Yw); TRAIN = 72; TARGETS = list(range(TRAIN, T))
p_off = np.log(pop / 1e5); S = 300
print(f"[data] N={N} T={T} targets={len(TARGETS)} ({time.time()-t0:.0f}s)", flush=True)

def design(t, Wm, net=True):
    l_prev = ell[t - 1]; cols = [np.ones(N)]
    if net: cols.append(Wm @ l_prev)
    cols += [l_prev, p_off, np.full(N, math.sin(2 * math.pi * t / 52)),
             np.full(N, math.cos(2 * math.pi * t / 52))]
    return np.column_stack(cols)

QG, RG = [1e-4, 3e-4, 1e-3], [2.0, 5.0, 10.0]
def fit_nssm(Wm, q, r, net=True, T_hi=T):
    K = 6 if net else 5
    XL = [design(t, Wm, net) for t in range(1, T_hi)]
    return nb_mode_filter(Yw[1:T_hi], XL, q * np.eye(K), np.zeros(K), np.eye(K), r), XL
def plug_ls(fit, XL, lo, hi, r):
    return sum(nb_logpmf(Yw[t], np.exp(np.clip(XL[t - 1] @ fit["m_pred"][t - 1], -30, 30)), r).sum()
               for t in range(lo, hi))
def tune(Wm, net=True):
    best = None
    for q in QG:
        for r in RG:
            fit, XL = fit_nssm(Wm, q, r, net)
            v = plug_ls(fit, XL, 8, TRAIN, r)
            if best is None or v > best[0]: best = (v, q, r)
    return best[1], best[2]

qn, rn = tune(W5, True); qo, ro = tune(W5, False)
fit_net, XLn = fit_nssm(W5, qn, rn, True); fit_non, XLo = fit_nssm(W5, qo, ro, False)
print(f"[tuned] net (q={qn}, r={rn}) nonet (q={qo}, r={ro}) ({time.time()-t0:.0f}s)", flush=True)

def nssm_scores(fit, XL_unused, Wm, net, q, r, seed=3):
    rl = np.random.default_rng(seed)
    out = {1: dict(ls=[], medmae=[]), 2: dict(ls=[])}
    K = 6 if net else 5
    pit_u = []
    for h in (1, 2):
        for tgt in TARGETS:
            t_or = tgt - h
            if t_or < 8: continue
            fitc, XLc = fit_nssm(Wm, q, r, net, T_hi=t_or + 1)
            m, P = fitc["m_filt"][-1], fitc["P_filt"][-1]
            th = m[None, :] + rl.standard_normal((S, K)) @ np.linalg.cholesky(
                P + q * np.eye(K) + 1e-12 * np.eye(K)).T
            y_prev = np.repeat(Yw[t_or][None, :].astype(float), S, 0)
            Lq = np.linalg.cholesky(q * np.eye(K) + 1e-14 * np.eye(K))
            for kstep in range(h):
                if kstep > 0: th = th + rl.standard_normal((S, K)) @ Lq.T
                lam = np.empty_like(y_prev)
                for s in range(S):
                    lp2 = np.log1p(y_prev[s]); cols = [np.ones(N)]
                    if net: cols.append(Wm @ lp2)
                    cols += [lp2, p_off, np.full(N, math.sin(2 * math.pi * (t_or + kstep + 1) / 52)),
                             np.full(N, math.cos(2 * math.pi * (t_or + kstep + 1) / 52))]
                    lam[s] = np.exp(np.clip(np.column_stack(cols) @ th[s], -30, 30))
                if kstep < h - 1:
                    g = rl.gamma(r, 1 / r, size=lam.shape)
                    y_prev = rl.poisson(np.minimum(lam * g, 1e12)).astype(float)
            lp = nb_logpmf(Yw[tgt][None, :], np.minimum(lam, 1e8), r)
            out[h]["ls"].append(float((logsumexp(lp, axis=0) - math.log(S)).sum()))
            if h == 1:
                g = rl.gamma(r, 1 / r, size=lam.shape)
                Yd = rl.poisson(np.minimum(lam * g, 1e12))
                out[1]["medmae"].append(float(np.mean(np.abs(Yw[tgt] - np.median(Yd, 0)))))
    return out

from scipy.stats import nbinom
def nbcdf(y, mu, r):
    return nbinom.cdf(y, r, r / (r + mu))

res_net = nssm_scores(fit_net, XLn, W5, True, qn, rn)
res_non = nssm_scores(fit_non, XLo, W5, False, qo, ro)
print(f"[NSSM] net h1 {np.mean(res_net[1]['ls']):.1f} h2 {np.mean(res_net[2]['ls']):.1f}; "
      f"nonet h1 {np.mean(res_non[1]['ls']):.1f} h2 {np.mean(res_non[2]['ls']):.1f} ({time.time()-t0:.0f}s)", flush=True)

# PIT (orientation only) for net h=1 via randomized PIT on plug-in one-step
u_all = []
rlp = np.random.default_rng(8)
for tgt in TARGETS:
    mu = np.exp(np.clip(XLn[tgt - 1] @ fit_net["m_pred"][tgt - 1], -30, 30))
    lo = nbcdf(Yw[tgt] - 1, mu, rn); hi = nbcdf(Yw[tgt], mu, rn)
    u_all.append(lo + rlp.uniform(size=N) * (hi - lo))
u_all = np.concatenate(u_all)
cnt, _ = np.histogram(u_all, bins=10, range=(0, 1))
pit_chi = float(((cnt - len(u_all) / 10) ** 2 / (len(u_all) / 10)).sum())

# ---------------- external benchmarks on the same protocol ----------------
def glm_fit(nll, b0):
    return minimize(nll, b0, method="L-BFGS-B", options=dict(maxiter=400)).x
def pnar_lin_eval():
    ls, mm = [], []
    rl = np.random.default_rng(21)
    for tgt in TARGETS:
        t_or = tgt - 1
        Xr = np.column_stack([np.ones((t_or) * N), (Yw[:t_or] @ W5.T).ravel(), Yw[:t_or].ravel()])
        yv = Yw[1:t_or + 1].ravel()
        def nll(b):
            lam = np.clip(Xr @ b, 1e-8, None)
            return float(np.sum(lam - yv * np.log(lam)))
        b = minimize(nll, np.array([.5, .3, .5]), method="L-BFGS-B",
                     bounds=[(1e-8, None)] * 3).x
        lam1 = np.clip(b[0] + b[1] * (W5 @ Yw[t_or]) + b[2] * Yw[t_or], 1e-8, None)
        ls.append(float(poisson_logpmf(Yw[tgt], lam1).sum()))
        mm.append(float(np.mean(np.abs(Yw[tgt] - np.floor(lam1 + 1/3 - 0.02/np.maximum(lam1,1e-6))))))
    return np.mean(ls), np.mean(mm)
def pnar_log_eval():
    ls, mm = [], []
    for tgt in TARGETS:
        t_or = tgt - 1
        Xr = np.column_stack([np.ones(t_or * N), (ell[:t_or] @ W5.T).ravel(), ell[:t_or].ravel(),
                              np.repeat(p_off[None, :], t_or, 0).ravel()])
        yv = Yw[1:t_or + 1].ravel()
        def nll(b):
            eta = np.clip(Xr @ b, -30, 30)
            return float(np.sum(np.exp(eta) - yv * eta))
        b = glm_fit(nll, np.array([-.5, .5, .5, .5]))
        lam1 = np.exp(np.clip(b[0] + b[1] * (W5 @ ell[t_or]) + b[2] * ell[t_or] + b[3] * p_off, -30, 30))
        ls.append(float(poisson_logpmf(Yw[tgt], lam1).sum()))
        mm.append(float(np.mean(np.abs(Yw[tgt] - np.floor(lam1 + 1/3)))))
    return np.mean(ls), np.mean(mm)
def hhh4_eval():
    # unit intercepts (absorb population), annual seasonality, own + network epidemic, NB r
    a0 = np.log(Yw[:TRAIN].mean(0) + 0.1)
    par0 = np.concatenate([a0, [0.0, 0.0], [0.3, 0.2], [math.log(5.0)]])
    sw = np.array([math.sin(2 * math.pi * t / 52) for t in range(1, TRAIN)])
    cw = np.array([math.cos(2 * math.pi * t / 52) for t in range(1, TRAIN)])
    Yl = Yw[:TRAIN - 1]; Yt = Yw[1:TRAIN]; WY = Yl @ W5.T
    def unpack(p): return p[:N], p[N], p[N + 1], p[N + 2], p[N + 3], math.exp(p[N + 4])
    def nll(p):
        a, gs, gc, lam_e, phi_e, r = unpack(p)
        nu = np.exp(np.clip(a[None, :] + gs * sw[:, None] + gc * cw[:, None], -30, 30))
        mu = nu + max(lam_e, 0) * Yl + max(phi_e, 0) * WY
        mu = np.clip(mu, 1e-8, None)
        return float(-np.sum(gammaln(Yt + r) - gammaln(r) - gammaln(Yt + 1)
                             + r * np.log(r / (r + mu)) + Yt * np.log(mu / (r + mu))))
    p = minimize(nll, par0, method="L-BFGS-B", options=dict(maxiter=300)).x
    a, gs, gc, lam_e, phi_e, r = unpack(p)
    ls, mm = [], []
    for tgt in TARGETS:
        swt, cwt = math.sin(2 * math.pi * tgt / 52), math.cos(2 * math.pi * tgt / 52)
        mu = np.clip(np.exp(np.clip(a + gs * swt + gc * cwt, -30, 30))
                     + max(lam_e, 0) * Yw[tgt - 1] + max(phi_e, 0) * (W5 @ Yw[tgt - 1]), 1e-8, None)
        ls.append(float(nb_logpmf(Yw[tgt], mu, r).sum()))
        med = nbinom.ppf(0.5, r, r / (r + mu))
        mm.append(float(np.mean(np.abs(Yw[tgt] - med))))
    return np.mean(ls), np.mean(mm), float(lam_e), float(phi_e), float(r)
def ewma_eval():
    best = None
    for a in [0.2, 0.4, 0.6, 0.8]:
        for r in [1.0, 2.0, 5.0]:
            m = Yw[:4].mean(0) + 0.5; v = 0.0
            tr_ls = 0.0
            for t in range(4, TRAIN):
                tr_ls += float(nb_logpmf(Yw[t], np.clip(m, 1e-8, None), r).sum())
                m = (1 - a) * m + a * Yw[t]
            if best is None or tr_ls > best[0]: best = (tr_ls, a, r)
    _, a, r = best
    m = Yw[:4].mean(0) + 0.5
    for t in range(4, TRAIN): m = (1 - a) * m + a * Yw[t]
    ls, mm = [], []
    for tgt in TARGETS:
        mu = np.clip(m, 1e-8, None)
        ls.append(float(nb_logpmf(Yw[tgt], mu, r).sum()))
        mm.append(float(np.mean(np.abs(Yw[tgt] - nbinom.ppf(0.5, r, r / (r + mu))))))
        m = (1 - a) * m + a * Yw[tgt]
    return np.mean(ls), np.mean(mm), a, r
b_lin = pnar_lin_eval(); print(f"[PNAR lin] {b_lin} ({time.time()-t0:.0f}s)", flush=True)
b_log = pnar_log_eval(); print(f"[PNAR log] {b_log} ({time.time()-t0:.0f}s)", flush=True)
b_h4 = hhh4_eval(); print(f"[hhh4] {b_h4[:2]} lam/phi/r {b_h4[2:]} ({time.time()-t0:.0f}s)", flush=True)
b_ew = ewma_eval(); print(f"[EWMA] {b_ew[:2]} a={b_ew[2]} r={b_ew[3]} ({time.time()-t0:.0f}s)", flush=True)

# NSSM plug-in row + medians
plug_net = float(np.mean([nb_logpmf(Yw[t], np.exp(np.clip(XLn[t - 1] @ fit_net["m_pred"][t - 1], -30, 30)), rn).sum() for t in TARGETS]))
med_net = float(np.mean(res_net[1]["medmae"])); med_non = float(np.mean(res_non[1]["medmae"]))

# ---------------- randomization: stratified by population + retuned ----------------
test_non_pl = plug_ls(fit_non, XLo, TRAIN, T, ro)
def stat_of(Wm):
    # q re-tuned inside every permutation; r fixed at its tuned value so the
    # statistic is one fixed function of (W, data) at a third of the cost
    best = None
    for q in QG:
        fit, XL = fit_nssm(Wm, q, rn, True)
        v = plug_ls(fit, XL, 8, TRAIN, rn)
        if best is None or v > best[0]: best = (v, fit, XL)
    _, fit, XL = best
    return plug_ls(fit, XL, TRAIN, T, rn) - test_non_pl
obs = stat_of(W5)
strata = np.digitize(p_off, np.quantile(p_off, [.2, .4, .6, .8]))
res_rt = {}
for label, stratified in [("stratified", True), ("unrestricted", False)]:
    rng = np.random.default_rng(51 if stratified else 52)
    vals = []
    for _ in range(99):
        p = np.arange(N)
        if stratified:
            for sgrp in range(5):
                idx = np.where(strata == sgrp)[0]; p[idx] = rng.permutation(idx)
        else:
            p = rng.permutation(N)
        vals.append(stat_of(W5[np.ix_(p, p)]))
    pval = (1 + sum(v >= obs for v in vals)) / 100.0
    res_rt[label] = dict(obs=float(obs), max_perm=float(max(vals)), p=float(pval))
    print(f"[randtest {label}] obs {obs:.1f} max {max(vals):.1f} p={pval:.3f} ({time.time()-t0:.0f}s)", flush=True)

# ---------------- k=8 sensitivity + final path ----------------
q8, r8 = tune(W8, True)
fit8, XL8 = fit_nssm(W8, q8, r8, True)
d8 = plug_ls(fit8, XL8, TRAIN, T, r8) - test_non_pl
d5 = plug_ls(fit_net, XLn, TRAIN, T, rn) - test_non_pl
b1 = fit_net["m_filt"][:, 1]; sd1 = np.sqrt(fit_net["P_filt"][:, 1, 1])
path = dict(first=float(b1[8:16].mean()), last=float(b1[-8:].mean()),
            lo=float(b1[8:].min()), hi=float(b1[8:].max()), sd_mean=float(sd1[8:].mean()))
print(f"[k-sens] plug-in gain k=5 {d5:.1f} vs k=8 {d8:.1f}; path {path} ({time.time()-t0:.0f}s)", flush=True)

rows_bench = [
    ("NSSM net (NB, offsets, seasonal)", np.mean(res_net[1]["ls"]), med_net, np.mean(res_net[2]["ls"])),
    ("NSSM net (NB, plug-in)", plug_net, None, None),
    ("NSSM no-network (NB, offsets, seasonal)", np.mean(res_non[1]["ls"]), med_non, np.mean(res_non[2]["ls"])),
    ("hhh4-type NB (unit intercepts, seasonal)", b_h4[0], b_h4[1], None),
    ("PNAR(1) log-linear (+ log-pop)", b_log[0], b_log[1], None),
    ("PNAR(1) linear", b_lin[0], b_lin[1], None),
    ("EWMA--NB (per node)", b_ew[0], b_ew[1], None),
]
frag = []
for lab, l1, m1, l2 in rows_bench:
    frag.append(f"{lab} & {l1:.1f} & " + (f"{m1:.2f}" if m1 is not None else "--")
                + " & " + (f"{l2:.1f}" if l2 is not None else "--") + " \\\\")
open(f"{TAB}/tab_covid_bench.tex", "w").write("\n".join(frag) + "\n")

json.dump(dict(N=int(N), T=int(T), targets=len(TARGETS), tuned=dict(q=qn, r=rn),
    ls=dict(net_h1=float(np.mean(res_net[1]["ls"])), nonet_h1=float(np.mean(res_non[1]["ls"])),
            net_h2=float(np.mean(res_net[2]["ls"])), nonet_h2=float(np.mean(res_non[2]["ls"])),
            plug_net=plug_net),
    medmae=dict(net=med_net, nonet=med_non),
    bench=dict(pnar_lin=b_lin, pnar_log=b_log, hhh4=b_h4[:2], ewma=b_ew[:2]),
    pit_chi=pit_chi, randtest=res_rt, ksens=dict(d5=float(d5), d8=float(d8)),
    path=path, pop=dict(lo=float(pop.min()), hi=float(pop.max()))),
    open(f"{OUT}/covid2.json", "w"))
print("DONE", round(time.time() - t0), "s", flush=True)
