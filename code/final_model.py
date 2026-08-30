"""
Diagnostics under the FINAL model (NB + unit offsets + seasonality), plus the
conditional randomization tests. Contents:
  - central beta_{1,t} path, information and test recomputed under the final
    model rather than the pooled Poisson model;
  - profile information with the ENLARGED nuisance design [1, offset, own lag]
    (seasonal columns are cross-sectionally degenerate: verified by rank check);
  - Q with q_{beta1}=0 ablation: is spillover drift data-preferred?
  - conditional randomization test: labels permuted WITHIN offset-quintile
    strata, and the (q, r) tuning rerun inside every permutation (exactness);
  - unrestricted-but-retuned test for comparison;
  - tau_W of the 60 degree-preserving rewirings (non-invariance, item B9).
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, math, sys, time
import numpy as np, pandas as pd, scipy.io as sio
sys.path.insert(0, _HERE)
from nssm import row_normalise
from nb_filter import nb_mode_filter, nb_logpmf

t0 = time.time()
OUT = _os.path.join(_HERE, "out")
df = pd.read_csv(_os.path.join(_ROOT, "data", "crime.csv"), index_col=0)
Y = df.values.T.astype(int); T, N = Y.shape
A0 = sio.mmread(_os.path.join(_ROOT, "data", "neighborhood.mtx")).tocsr()
A0 = ((A0 + A0.T) > 0).astype(float).toarray(); np.fill_diagonal(A0, 0)
W = row_normalise(A0); ell = np.log1p(Y)
TRAIN = 60; o = np.log(Y[:TRAIN].mean(0) + 0.05)
mth = lambda t: 2 * math.pi * (t % 12) / 12.0

def design(t, Wm, net=True):
    l_prev = ell[t - 1]
    cols = [np.ones(N)]
    if net: cols.append(Wm @ l_prev)
    cols += [l_prev, o, np.full(N, math.sin(mth(t))), np.full(N, math.cos(mth(t)))]
    return np.column_stack(cols)

# seasonal degeneracy: at fixed t the seasonal columns are scalar multiples of 1
assert np.linalg.matrix_rank(np.column_stack([np.ones(N), np.full(N, .3), np.full(N, .9)])) == 1

def fit_model(Wm, q, r, net=True, qzero_spill=False, T_hi=T):
    K = 6 if net else 5
    Q = q * np.eye(K)
    if qzero_spill and net: Q[1, 1] = 0.0
    X_list = [design(t, Wm, net) for t in range(1, T_hi)]
    return nb_mode_filter(Y[1:T_hi], X_list, Q, np.zeros(K), np.eye(K), r), X_list

def plug_ls(fit, X_list, lo, hi, r):
    return sum(nb_logpmf(Y[t], np.exp(np.clip(X_list[t - 1] @ fit["m_pred"][t - 1], -30, 30)), r).sum()
               for t in range(lo, hi))

QG, RG = [1e-4, 3e-4, 1e-3], [2.0, 5.0, 10.0]
def tune(Wm, net=True, qzero=False):
    best = None
    for q in QG:
        for r in RG:
            fit, XL = fit_model(Wm, q, r, net, qzero)
            v = plug_ls(fit, XL, 12, TRAIN, r)
            if best is None or v > best[0]: best = (v, q, r, fit, XL)
    return best

# ---------- 1. final-model fits, path, ablation ----------
v_net = tune(W, net=True); v_non = tune(W, net=False); v_qz = tune(W, net=True, qzero=True)
_, qn, rn, fit_net, XLn = v_net
_, qo, ro, fit_non, XLo = v_non
_, qz, rz, fit_qz, XLz = v_qz
test_net = plug_ls(fit_net, XLn, TRAIN, T, rn)
test_non = plug_ls(fit_non, XLo, TRAIN, T, ro)
test_qz = plug_ls(fit_qz, XLz, TRAIN, T, rz)
b1 = fit_net["m_filt"][:, 1]; sd1 = np.sqrt(fit_net["P_filt"][:, 1, 1])
path = dict(first=float(b1[9:19].mean()), last=float(b1[-10:].mean()),
            lo=float(b1[9:].min()), hi=float(b1[9:].max()),
            sd_mean=float(sd1[9:].mean()), sd_max=float(sd1[9:].max()))
print(f"[final] q={qn} r={rn}; test LS net {test_net:.1f} nonet {test_non:.1f} "
      f"qzero {test_qz:.1f}; b1 {path['first']:.2f}->{path['last']:.2f} "
      f"[{path['lo']:.2f},{path['hi']:.2f}] sd {path['sd_mean']:.3f} ({time.time()-t0:.0f}s)", flush=True)

# ---------- 2. weighted information, enlarged nuisance [1, o, ell] ----------
def winfo(fit, Wm, months, spill_idx=1):
    out = []
    for t in months:
        X = np.column_stack([np.ones(N), Wm @ ell[t - 1], ell[t - 1], o])
        mu = np.exp(np.clip(design(t, Wm) @ fit["m_pred"][t - 1], -30, 30))
        lam_w = rn * mu / (rn + mu)
        G = X.T @ (X * lam_w[:, None])
        try: out.append(1.0 / np.linalg.inv(G)[spill_idx, spill_idx])
        except np.linalg.LinAlgError: out.append(0.0)
    return np.array(out)
months = range(12, T)
I_enl = winfo(fit_net, W, months)
# base three-column weighted info under the same fitted weights, for comparison
def winfo_base(fit, Wm, months):
    out = []
    for t in months:
        X = np.column_stack([np.ones(N), Wm @ ell[t - 1], ell[t - 1]])
        mu = np.exp(np.clip(design(t, Wm) @ fit["m_pred"][t - 1], -30, 30))
        lam_w = rn * mu / (rn + mu)
        G = X.T @ (X * lam_w[:, None])
        out.append(1.0 / np.linalg.inv(G)[1, 1])
    return np.array(out)
I_base = winfo_base(fit_net, W, months)
print(f"[info] enlarged mean {I_enl.mean():.1f} vs base {I_base.mean():.1f}; "
      f"enlarged<=base at all t: {bool(np.all(I_enl <= I_base + 1e-8))} ({time.time()-t0:.0f}s)", flush=True)

# ---------- 3. randomization: stratified + retuned, and unrestricted + retuned ----------
strata = np.digitize(o, np.quantile(o, [.2, .4, .6, .8]))
def stat_of(Wm):
    _, q, r, fit, XL = tune(Wm, net=True)
    return plug_ls(fit, XL, TRAIN, T, r) - test_non
obs = stat_of(W)
def perm_W(rng, stratified):
    p = np.arange(N)
    if stratified:
        for s in range(5):
            idx = np.where(strata == s)[0]; p[idx] = rng.permutation(idx)
    else:
        p = rng.permutation(N)
    return W[np.ix_(p, p)]
res_rt = {}
for label, stratified in [("stratified", True), ("unrestricted", False)]:
    rng = np.random.default_rng(41 if stratified else 42)
    vals = [stat_of(perm_W(rng, stratified)) for _ in range(99)]
    p = (1 + sum(v >= obs for v in vals)) / 100.0
    res_rt[label] = dict(obs=float(obs), max_perm=float(max(vals)), p=float(p),
                         vals=list(map(float, vals)))
    print(f"[randtest {label}] obs {obs:.1f} max perm {max(vals):.1f} p={p:.3f} "
          f"({time.time()-t0:.0f}s)", flush=True)

# ---------- 4. rewiring tau_W non-invariance ----------
def tau_of(Am):
    Wm = row_normalise(Am)
    return float(np.sum(Wm * Wm) - np.sum(Wm.sum(0) ** 2) / N)
def rewire(Am, n_swaps, rng):
    A = Am.copy(); edges = np.argwhere(np.triu(A) > 0)
    for _ in range(n_swaps):
        i = rng.integers(len(edges)); j = rng.integers(len(edges))
        (a, b), (c, d) = edges[i], edges[j]
        if len({a, b, c, d}) < 4: continue
        if A[a, d] or A[c, b]: continue
        A[a, b] = A[b, a] = A[c, d] = A[d, c] = 0
        A[a, d] = A[d, a] = A[c, b] = A[b, c] = 1
        edges[i] = [min(a, d), max(a, d)]; edges[j] = [min(c, b), max(c, b)]
    return A
taus = [tau_of(rewire(A0, 4 * int(A0.sum() // 2), np.random.default_rng(s))) for s in range(60)]
tau_obs = tau_of(A0)
print(f"[rewire] tau_W obs {tau_obs:.2f}; rewired range [{min(taus):.2f}, {max(taus):.2f}]", flush=True)

json.dump(dict(tuned=dict(q=qn, r=rn), test_ls=dict(net=float(test_net), nonet=float(test_non),
        qzero_spill=float(test_qz), delta_net=float(test_net - test_non),
        delta_drift=float(test_net - test_qz)), path=path,
    info=dict(enlarged_mean=float(I_enl.mean()), base_mean=float(I_base.mean()),
              enlarged_min=float(I_enl.min()), enlarged_max=float(I_enl.max())),
    randtest=res_rt, rewire_tau=dict(obs=tau_obs, lo=float(min(taus)), hi=float(max(taus)))),
    open(f"{OUT}/final_model.json", "w"))
print("DONE", round(time.time() - t0), "s", flush=True)
