"""
Independent numerical verification of every mathematical claim in the paper:
Lemmas A.1-A.5, Theorem 1 (+Corollary), Theorem 2, Proposition 3, Theorem 4,
and all explicit constants. Prints PASS/FAIL per check.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import numpy as np, sys, math
from scipy.stats import chi2
from scipy.special import gammaln, logsumexp, comb
sys.path.insert(0, _HERE)
from nssm import row_normalise, er_graph, complete_graph, op_norm, tau_W, profile_information

rng = np.random.default_rng(0)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)

def MZ(Z):
    return np.eye(Z.shape[0]) - Z @ np.linalg.pinv(Z)

# ---------------------------------------------------------------- Lemma A.1 (Schur)
ok = True
for trial in range(200):
    K = rng.integers(3, 6)
    Amat = rng.normal(size=(K, K)); M = Amat @ Amat.T + 0.1 * np.eye(K)
    s_block = M[1, 1] - M[1, np.r_[0, 2:K]] @ np.linalg.inv(M[np.ix_(np.r_[0, 2:K], np.r_[0, 2:K])]) @ M[np.r_[0, 2:K], 1]
    ok &= abs(1.0 / np.linalg.inv(M)[1, 1] - s_block) < 1e-8 * abs(s_block)
    # variational form
    E = np.zeros((K, K - 1)); idx = [j for j in range(K) if j != 1]
    for c, j in enumerate(idx): E[j, c] = 1
    a_star = -np.linalg.solve(M[np.ix_(idx, idx)], M[idx, 1])
    v = np.zeros(K); v[1] = 1; v += E @ a_star
    ok &= abs(v @ M @ v - s_block) < 1e-8 * abs(s_block)
    # random a gives larger value
    for _ in range(5):
        vr = np.zeros(K); vr[1] = 1; vr += E @ rng.normal(size=K - 1)
        ok &= vr @ M @ vr >= s_block - 1e-10
check("Lemma A.1(a): s(M)=1/(M^{-1})_22 = variational min", ok)

ok = True
for trial in range(200):
    K = 4
    A1 = rng.normal(size=(K, K)); A1 = A1 @ A1.T
    A2 = rng.normal(size=(K, K)); A2 = A2 @ A2.T
    def svar(M):
        idx = [j for j in range(K) if j != 1]
        Mi = M[np.ix_(idx, idx)] + 1e-12 * np.eye(K - 1)
        return M[1, 1] - M[1, idx] @ np.linalg.solve(Mi, M[idx, 1])
    ok &= svar(A1 + A2) >= svar(A1) + svar(A2) - 1e-8
check("Lemma A.1(b): superadditivity s(A+G)>=s(A)+s(G)", ok)

ok = True
for trial in range(200):
    N = rng.integers(5, 30)
    x = rng.normal(size=N); w = rng.normal(size=N); Z = np.column_stack([np.ones(N), x])
    info = np.min([np.sum((w + Z @ a) ** 2) for a in (rng.normal(size=(400, 2)) * 3)])
    exact = np.sum((MZ(Z) @ w) ** 2)
    a_star = -np.linalg.pinv(Z) @ w
    ok &= np.sum((w + Z @ a_star) ** 2) <= exact + 1e-9 and info >= exact - 1e-9
    X = np.column_stack([np.ones(N), w, x])
    pi, _ = profile_information(X, col=1)
    ok &= abs(pi - exact) < 1e-8 * max(1, exact)
check("Lemma A.1(c): profile info = ||M_Z w||^2, attained at -Z^+ w", ok)

# ---------------------------------------------------------------- Lemma A.2 (tau_W)
ok = True
for trial in range(60):
    N = int(rng.integers(6, 60)); A = er_graph(N, 4, rng); W = row_normalise(A)
    n_i = A.sum(1)
    tau = tau_W(W)
    ok &= abs(tau - (np.sum(W ** 2) - np.sum(W.T @ np.ones(N)) ** 0 * np.linalg.norm(W.T @ np.ones(N)) ** 2 / N)) < 1e-8
    ok &= abs(np.sum(W ** 2) - np.sum(1 / n_i)) < 1e-8
    ok &= tau <= np.sum(1 / n_i) + 1e-8 and tau <= N + 1e-8
# d-regular: cycle graph (d=2) and complete
for N in [5, 12, 40]:
    Ac = np.zeros((N, N))
    for i in range(N): Ac[i, (i + 1) % N] = Ac[i, (i - 1) % N] = 1
    ok &= abs(tau_W(row_normalise(Ac)) - (N / 2 - 1)) < 1e-8
    ok &= abs(tau_W(complete_graph(N) / (N - 1)) - 1 / (N - 1)) < 1e-8
# dense: all degrees > cN => tau <= 1/c
N = 60; A = er_graph(N, N * 0.5, rng); A[A.sum(1) < 0.3 * N] = 1; A = np.maximum(A, A.T); np.fill_diagonal(A, 0)
c = A.sum(1).min() / N
ok &= tau_W(row_normalise(A)) <= 1 / c + 1e-8
check("Lemma A.2: tau identity, =sum 1/n_i, regular N/d-1, complete 1/(N-1), <=1/c", ok)

# ---------------------------------------------------------------- Theorem 1(i)
ok = True
for trial in range(300):
    N = int(rng.integers(5, 40)); W = row_normalise(er_graph(N, 3, rng))
    x = rng.normal(size=N) * rng.uniform(0.5, 3); sig = rng.uniform(0.5, 2)
    Z = np.column_stack([np.ones(N), x]); Mz = MZ(Z)
    I1 = np.sum((Mz @ (W @ x)) ** 2) / sig ** 2
    I2 = np.sum((Mz @ ((np.eye(N) - W) @ x)) ** 2) / sig ** 2
    ok &= abs(I1 - I2) < 1e-8 * max(1, I1)
    ok &= I1 <= np.sum(((np.eye(N) - W) @ x) ** 2) / sig ** 2 + 1e-8
    xt = x - x.mean(); M1 = np.eye(N) - np.ones((N, N)) / N
    I3 = (np.sum((M1 @ W @ x) ** 2) - (xt @ (W @ x)) ** 2 / np.sum(xt ** 2)) / sig ** 2
    ok &= abs(I1 - I3) < 1e-7 * max(1, I1)
check("Thm 1(i): roughness identity, projection bound, two-term formula", ok)

# ---------------------------------------------------------------- Theorem 1(ii)
N = 25; Wc = complete_graph(N) / (N - 1)
ok = max(np.sum((MZ(np.column_stack([np.ones(N), x])) @ (Wc @ x)) ** 2)
         for x in rng.normal(size=(200, N))) < 1e-18
# only-if: any random row-stochastic zero-diagonal W != Wc has positive info for some basis vector
worst = np.inf
for trial in range(300):
    N = int(rng.integers(4, 20))
    W = rng.random((N, N)); np.fill_diagonal(W, 0); W = W / W.sum(1, keepdims=True)
    Wc = complete_graph(N) / (N - 1)
    if np.max(np.abs(W - Wc)) < 1e-6: continue
    m = max(np.sum((MZ(np.column_stack([np.ones(N), e])) @ (W @ e)) ** 2)
            for e in np.eye(N))
    worst = min(worst, m / np.max(np.abs(W - Wc)) ** 2)
ok &= worst > 1e-6
# structural step: sum_{j != i} a_j = 1 for all i forces a_j = 1/(N-1)
for N in [4, 9, 17]:
    Amat = np.ones((N, N)) - np.eye(N)
    a = np.linalg.solve(Amat, np.ones(N))
    ok &= np.max(np.abs(a - 1 / (N - 1))) < 1e-10
check("Thm 1(ii): complete graph <=> zero information (both directions)", ok, f"min info ratio off-Wc = {worst:.2e}")

# ---------------------------------------------------------------- Theorem 1(iii)
N = 50; W = row_normalise(er_graph(N, 4, rng)); tau = tau_W(W); sig = 1.3; c0 = 0.7
M1 = np.eye(N) - np.ones((N, N)) / N
R = 120000
E_all = rng.normal(0, sig, size=(R, N))
X = c0 + E_all                              # mu = c 1
MW = (M1 @ W)
num = np.einsum("ri,ri->r", E_all @ MW.T, E_all @ MW.T)        # ||M1 W eps||^2
quad = np.einsum("ri,ri->r", E_all, E_all @ MW.T)               # eps' M1 W eps  (since (M1W)' eps dot eps)
den = np.einsum("ri,ri->r", E_all @ M1.T, E_all @ M1.T)
I_mc = (num - quad ** 2 / den) / sig ** 2
cN = tau * (2 / (N - 1) + math.sqrt(3) * math.exp(-(N - 1) / 32)) + 2 * N / (N - 1)
m_hat, se = I_mc.mean(), I_mc.std() / math.sqrt(R)
ok = (m_hat <= tau + 4 * se) and (m_hat >= tau - cN - 4 * se)
check("Thm 1(iii): tau - c_N <= E[I] <= tau  (mu = c*1, Monte Carlo)", ok,
      f"E[I]={m_hat:.3f}+-{se:.3f}, tau={tau:.3f}, tau-c_N={tau-cN:.3f}")
# general mu upper bound
mu = rng.normal(0, 2, N)
Xg = mu + E_all
Zb = np.stack([np.column_stack([np.ones(N), x]) for x in Xg[:20000]])
I_g = np.array([np.sum((MZ(Zb[r]) @ (W @ Xg[r])) ** 2) for r in range(20000)]) / sig ** 2
ub = tau + np.sum((M1 @ W @ mu) ** 2) / sig ** 2
ok = I_g.mean() <= ub + 4 * I_g.std() / math.sqrt(20000)
check("Thm 1(iii): E[I] <= tau + ||M1 W mu||^2/sig^2 (general mu, MC)", ok,
      f"E[I]={I_g.mean():.3f}, bound={ub:.3f}")
# moment identities used in the proof
Bq = W.T @ M1 @ W
lhs = np.mean(np.einsum("ri,ri->r", E_all @ Bq.T, E_all) ** 2)
rhs = sig ** 4 * (np.trace(Bq) ** 2 + 2 * np.trace(Bq @ Bq))
ok = abs(lhs - rhs) < 5 * np.std(np.einsum("ri,ri->r", E_all @ Bq.T, E_all) ** 2) / math.sqrt(R)
Cs = (MW + MW.T) / 2
lhs2 = np.mean(quad ** 2)
rhs2 = sig ** 4 * (np.trace(W @ W) + tau)
ok &= abs(lhs2 - rhs2) < 5 * np.std(quad ** 2) / math.sqrt(R)
ok &= abs(np.trace(Cs) + 1) < 1e-9 and abs(2 * np.trace(Cs @ Cs) - (np.trace(W @ W) - 1 + tau)) < 1e-8
ok &= np.trace(W @ W) <= N + 1e-9
check("Thm 1(iii) proof moments: E(e'Be)^2, tr identities, tr(W^2)<=N", ok)
# Laurent-Massart specialisation
ok = all(chi2.cdf(n / 2, n) <= math.exp(-n / 16) + 1e-15 for n in range(1, 1001))
check("Laurent-Massart: P(chi2_n <= n/2) <= exp(-n/16), n=1..1000", ok)
# constants c_N
Ns = np.arange(3, 5000); cb = 4 + 4 / (Ns - 1) + math.sqrt(3) * Ns * np.exp(-(Ns - 1) / 32)
check("Thm 1(iii): c_N bound <= 26, max at N=32, -> 4", cb.max() <= 26 and abs(cb[-1] - 4) < 0.01,
      f"max={cb.max():.2f} at N={Ns[cb.argmax()]}")

# ---------------------------------------------------------------- Corollary (aggregation)
N = 12; W = row_normalise(er_graph(N, 3, rng)); th = np.array([0.3, 0.4, 0.2])
y0 = rng.normal(size=N); eps = rng.normal(0, 0.5, N)
y1 = th[0] + th[1] * W @ y0 + th[2] * y0 + eps
cvec = W.T @ np.ones(N)
ok = abs(y1.mean() - (th[0] + th[1] * (cvec @ y0) / N + th[2] * y0.mean() + eps.mean())) < 1e-10
Ac = np.zeros((N, N))
for i in range(N): Ac[i, (i + 1) % N] = Ac[i, (i - 1) % N] = 1
Wr = row_normalise(Ac); y1r = th[0] + th[1] * Wr @ y0 + th[2] * y0 + eps
ok &= abs(y1r.mean() - (th[0] + (th[1] + th[2]) * y0.mean() + eps.mean())) < 1e-10
check("Corollary 1: aggregation identities (general and doubly stochastic)", ok)

# ---------------------------------------------------------------- Theorem 2
# (i) exact coverage by simulation, batched Kalman filter
R, T, N = 40000, 6, 8
W = row_normalise(er_graph(N, 3, rng)); sig = 1.0; Q = np.diag([2e-3, 5e-3, 5e-3]); m0 = np.zeros(3); P0 = np.eye(3) * 0.4
th = m0 + rng.normal(size=(R, 3)) @ np.linalg.cholesky(P0).T
y = rng.normal(0, 1, size=(R, N))
m = np.repeat(m0[None, :], R, 0); P = np.repeat(P0[None, :, :], R, 0)
cov_hits = []; joint_hits = []; sandwich_ok = True; consist_ok = True
Lq = np.linalg.cholesky(Q)
for t in range(T):
    th = th + rng.normal(size=(R, 3)) @ Lq.T
    Xb = np.stack([np.ones((R, N)), y @ W.T, y], axis=2)          # R x N x 3
    y = np.einsum("rnk,rk->rn", Xb, th) + rng.normal(0, sig, size=(R, N))
    Pp = P + Q
    G = np.einsum("rnk,rnl->rkl", Xb, Xb) / sig ** 2
    Ppinv = np.linalg.inv(Pp)
    Pf = np.linalg.inv(Ppinv + G)
    m = np.einsum("rkl,rl->rk", Pf, np.einsum("rkl,rl->rk", Ppinv, m) + np.einsum("rnk,rn->rk", Xb, y) / sig ** 2)
    P = Pf
    e1 = th[:, 1] - m[:, 1]; sd1 = np.sqrt(P[:, 1, 1])
    cov_hits.append(np.mean(np.abs(e1) <= 1.959963985 * sd1))
    d = th - m
    q95 = chi2.ppf(0.95, 3)
    joint_hits.append(np.mean(np.einsum("rk,rkl,rl->r", d, np.linalg.inv(P), d) <= q95))
    # (ii) sandwich for a subsample
    for r in range(0, R, R // 50):
        Zr = np.column_stack([np.ones(N), Xb[r, :, 2]]); wr = Xb[r, :, 1]
        It = np.sum((MZ(Zr) @ wr) ** 2) / sig ** 2
        ustar = np.linalg.pinv(Zr) @ wr
        lo = 1.0 / (It + (1 + ustar @ ustar) / np.linalg.eigvalsh(Q).min())
        up = np.inf if It < 1e-14 else 1.0 / It
        sandwich_ok &= (P[r, 1, 1] <= up + 1e-9) and (P[r, 1, 1] >= lo - 1e-9)
        # (iii)
        lmin = np.linalg.eigvalsh(Xb[r].T @ Xb[r]).min()
        if lmin > 1e-10:
            consist_ok &= np.trace(P[r]) <= 3 * sig ** 2 / lmin + 1e-9
cov_hits = np.array(cov_hits); joint_hits = np.array(joint_hits)
se = math.sqrt(0.95 * 0.05 / R)
check("Thm 2(i): exact 95% coverage of beta1 credible intervals (all t)",
      np.all(np.abs(cov_hits - 0.95) < 4 * se), f"coverage by t: {np.round(cov_hits,4)}")
check("Thm 2(i): exact 95% coverage of joint credible ellipsoid (all t)",
      np.all(np.abs(joint_hits - 0.95) < 4 * se), f"{np.round(joint_hits,4)}")
check("Thm 2(ii): 1/(I+(1+||u*||^2)/lmin(Q)) <= Var(beta1|F_t) <= 1/I", sandwich_ok)
check("Thm 2(iii): tr P_{t|t} <= 3 sig^2 / lmin(X'X)", consist_ok)
# (iv) complete graph null direction
N = 15; Wc = complete_graph(N) / (N - 1); yv = rng.normal(size=N)
Xc = np.column_stack([np.ones(N), Wc @ yv, yv])
v = np.array([-N * yv.mean() / (N - 1), 1.0, 1.0 / (N - 1)])
ok = np.max(np.abs(Xc @ v)) < 1e-10
Qc = np.diag([1e-3, 4e-3, 2e-3]); Pp = rng.normal(size=(3, 3)); Pp = Pp @ Pp.T + 0.05 * np.eye(3); Pp = Pp + Qc
Pf = np.linalg.inv(np.linalg.inv(Pp) + Xc.T @ Xc)
ok &= np.max(np.abs(np.linalg.inv(Pf) @ v - np.linalg.inv(Pp) @ v)) < 1e-8
ok &= v @ Pf @ v >= np.linalg.eigvalsh(Qc).min() * (v @ v) - 1e-10 and v @ v >= 1
check("Thm 2(iv): X_t v_t = 0, precision invariance, variance floor lmin(Q)||v||^2", ok)

# ---------------------------------------------------------------- Proposition 3
ok = True
for trial in range(100):
    N = int(rng.integers(6, 30)); W = row_normalise(er_graph(N, 3, rng))
    Wt = row_normalise(er_graph(N, 4, rng))
    b0, b1, b2 = rng.normal(0, 0.3, 3)
    B = b1 * W + b2 * np.eye(N); Bt = b1 * Wt + b2 * np.eye(N)
    y = rng.normal(0, 2, N); M1 = np.eye(N) - np.ones((N, N)) / N
    for h in [1, 2, 3, 5, 8]:
        mu = y.copy(); mut = y.copy()
        for k in range(h):
            mu = b0 + B @ mu; mut = b0 + Bt @ mut
        d = mu - mut
        d2 = (np.linalg.matrix_power(B, h) - np.linalg.matrix_power(Bt, h)) @ (M1 @ y)
        ok &= np.max(np.abs(d - d2)) < 1e-8 * max(1, np.max(np.abs(d)))
        rb = max(op_norm(B), op_norm(Bt))
        ok &= np.linalg.norm(d) <= abs(b1) * op_norm(W - Wt) * h * rb ** (h - 1) * np.linalg.norm(M1 @ y) + 1e-8
        if rb < 1:
            ok &= np.linalg.norm(d) <= abs(b1) * op_norm(W - Wt) * np.linalg.norm(M1 @ y) / (1 - rb) ** 2 + 1e-8
    ok &= np.max(np.abs(np.linalg.eigvals(B))) >= abs(b1 + b2) - 1e-9
check("Prop 3(i)-(iii): exact identity, telescoping bound, sup bound, rho(B)>=|b1+b2|", ok)
# (iii) sharpness on a common eigenvector (cycle graph + Wc mixture)
ok = True
for N in [8, 20]:
    Ac = np.zeros((N, N))
    for i in range(N): Ac[i, (i + 1) % N] = Ac[i, (i - 1) % N] = 1
    W = row_normalise(Ac); Wc = complete_graph(N) / (N - 1)
    alpha = 0.3; Wt = (1 - alpha) * W + alpha * Wc
    evals, evecs = np.linalg.eigh(W)                      # symmetric here
    j = 1                                                 # a non-Perron eigenvector
    v = evecs[:, j]; lam = evals[j]
    lamt = (1 - alpha) * lam - alpha / (N - 1)
    ok &= np.max(np.abs(Wt @ v - lamt * v)) < 1e-8       # common eigenvector
    b0, b1, b2 = 0.1, 0.6, 0.55
    r = b1 * lam + b2; rt = b1 * lamt + b2
    B = b1 * W + b2 * np.eye(N); Bt = b1 * Wt + b2 * np.eye(N)
    cc = 1.7; y = 2.0 + cc * v
    for h in [1, 2, 4, 8]:
        mu = y.copy(); mut = y.copy()
        for k in range(h): mu = b0 + B @ mu; mut = b0 + Bt @ mut
        ok &= abs(np.linalg.norm(mu - mut) - abs(cc) * abs(r ** h - rt ** h)) < 1e-7
        if r * rt > 0:
            ok &= abs(r ** h - rt ** h) >= h * min(abs(r), abs(rt)) ** (h - 1) * abs(b1) * abs(lam - lamt) - 1e-9
    ok &= op_norm(W - Wt) >= abs(lam - lamt) - 1e-9
check("Prop 3(iii): common-eigenvector equality and geometric lower bound", ok)

# ---------------------------------------------------------------- Theorem 4
# (a) two-step formula by Monte Carlo
N = 6; A = er_graph(N, 3, rng); W = row_normalise(A)
b0, b1, b2 = 0.15, 0.30, 0.25
S = np.array([[b1 * W[i, j] + (b2 if i == j else 0.0) for j in range(N)] for i in range(N)])
yt = np.array([1., 0., 2., 0., 1., 3.])
lam1 = np.exp(b0 + S @ yt)
i0 = 0
formula = math.exp(b0) * math.exp(np.sum(lam1 * (np.exp(S[i0]) - 1)))
Rmc = 4_000_000
y1 = rng.poisson(lam1, size=(Rmc, N))
lam2_i = np.exp(b0 + y1 @ S[i0])
mc, mcse = lam2_i.mean(), lam2_i.std() / math.sqrt(Rmc)
check("Thm 4(a): two-step formula (eq. 8) vs Monte Carlo",
      abs(mc - formula) < 4 * mcse and abs(mc / formula - 1) < 0.02,
      f"MC={mc:.4f}+-{mcse:.4f}, formula={formula:.4f}")
# (a) h=3 divergence: exact lower-bound series terms explode (log scale)
def log_terms_h3(b0, b1, b2, W, yt, i0, m_max=250):
    N = W.shape[0]
    S = b1 * W + b2 * np.eye(N)
    lam1 = np.exp(b0 + S @ yt); Lam = lam1.sum()
    logs = []
    for m in range(m_max):
        # E_m: y_{i0,t+1}=m, others 0 -> lam_{j,t+2} = exp(b0 + S[j,i0] m)
        lam2 = np.exp(np.minimum(b0 + S[:, i0] * m, 700))
        Phi = np.sum(np.minimum(lam2, 1e290) * (np.exp(S[i0]) - 1))
        logs.append(-Lam + m * math.log(lam1[i0]) - gammaln(m + 1) + b0 + min(Phi, 1e290))
    return np.array(logs)
lt1 = log_terms_h3(0.1, 0.3, 0.25, W, yt, 0)                     # beta2>0
lt2 = log_terms_h3(0.1, 0.5, -0.4, W, yt, 0)                     # beta1>0, beta2<0 (neighbour case)
check("Thm 4(a): h=3 lower-bound series diverges (beta2>0)", lt1.max() > 1e6, f"max log term={lt1.max():.3g}")
check("Thm 4(a): h=3 lower-bound series diverges (beta1>0, beta2<0)", lt2.max() > 1e6, f"max log term={lt2.max():.3g}")
# (a) non-positive case bound
b0n, b1n, b2n = 0.4, -0.3, -0.2
Sn = b1n * W + b2n * np.eye(N)
yv = yt.copy(); ok = True
for rep in range(2000):
    y = yt.copy()
    for h in range(6):
        lam = np.exp(b0n + Sn @ y)
        ok &= np.all(lam <= math.exp(b0n) + 1e-12)
        y = rng.poisson(lam).astype(float)
check("Thm 4(a): beta1,beta2<=0 => intensities (hence means) <= e^{beta0}", ok)
# (b) divergence with random-walk coefficients: lower-bound series in log scale
N3 = 3
A3 = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0.]]); W3 = row_normalise(A3)
y3 = np.ones(3); i0 = 0
z_all = np.column_stack([np.ones(3), W3 @ y3, y3])               # z_{j,t}
mtheta = np.array([0.1, 0.2, 0.15]); Ptheta = 0.01 * np.eye(3)
Qb = 0.4 * np.eye(3); qmin = np.linalg.eigvalsh(Qb).min()
draws = mtheta + rng.normal(size=(300000, 3)) @ np.linalg.cholesky(Ptheta).T
log_lam_all = draws @ z_all.T                                     # log lambda_j(theta)
Lam_all = np.exp(log_lam_all).sum(1)
lt = []
for m in range(0, 120):
    expo = draws[:, 0] + m * draws[:, 2] + m * log_lam_all[:, i0] - Lam_all
    log_E = logsumexp(expo) - math.log(len(draws))
    lt.append(0.5 * qmin * m * m + log_E - gammaln(m + 1))
lt = np.array(lt)
check("Thm 4(b): h=2 lower-bound series diverges under RW coefficients (Q>0)",
      lt.max() > 1e3 and lt[-1] > lt[20], f"log term at m=119: {lt[-1]:.1f}")
# (b) h=1 finite under Gaussian coefficients: lognormal mean formula vs MC
z0 = z_all[i0]
exact1 = math.exp(z0 @ mtheta + 0.5 * z0 @ Ptheta @ z0)
mc1 = np.exp(draws @ z0).mean()
check("Thm 4(b): one-step mean finite and equals lognormal formula",
      abs(mc1 / exact1 - 1) < 0.01, f"MC={mc1:.5f}, exact={exact1:.5f}")
# (c) log1p: intensity bound and stable moments
ok = True
b0c, b1c, b2c = 0.2, 0.45, 0.35; rho = abs(b1c) + abs(b2c)
maxs = []
for rep in range(4000):
    y = np.array([3., 1., 6., 0., 2., 4.])
    for h in range(6):
        l = np.log1p(y)
        lam = np.exp(b0c + b1c * (W @ l) + b2c * l)
        ok &= np.all(lam <= math.exp(abs(b0c)) * (1 + y.max()) ** rho + 1e-9)
        y = rng.poisson(lam).astype(float)
    maxs.append(y.max())
maxs = np.array(maxs)
ok &= np.isfinite(np.mean((1 + maxs) ** 4.0))
check("Thm 4(c): log1p intensity bound e^{|b0|}(1+max y)^rho; moments finite", ok,
      f"E(1+M_6)^4 = {np.mean((1+maxs)**4):.1f}")
# (d) full inequality chain numerically on a random instance
N5 = 5; W5 = row_normalise(er_graph(N5, 2, rng))
y5 = np.array([2., 0., 3., 1., 1.]); L = math.log(1 + y5.max())
qbar = 0.05; Q5 = qbar * np.eye(3); P5 = 0.02 * np.eye(3); m5 = np.array([0.1, 0.2, 0.15])
ct = math.sqrt(3) * (1 + L); dN = 2 + math.log(1 + N5)
kap = max(6.0, math.sqrt(3) * ct + 2 * qbar * ct * ct)
assert 2 * kap * np.linalg.eigvalsh(P5).max() < 1
Rd = 2_000_000
thd = m5 + rng.normal(size=(Rd, 3)) @ np.linalg.cholesky(P5).T
l5 = np.log1p(y5); Z5 = np.column_stack([np.ones(N5), W5 @ l5, l5])
lam1d = np.exp(thd @ Z5.T)
y1d = rng.poisson(lam1d)
th2 = thd + rng.normal(size=(Rd, 3)) @ np.linalg.cholesky(Q5).T
l1d = np.log1p(y1d)
i0 = 0
z2 = np.stack([np.ones(Rd), l1d @ W5[i0], l1d[:, i0]], axis=1)
lam2d = np.exp(np.einsum("rk,rk->r", z2, th2))
lhs, lhs_se = lam2d.mean(), lam2d.std() / math.sqrt(Rd)
nrm = np.linalg.norm(thd, axis=1)
rhs = math.exp(qbar / 2) * (math.exp(2 * qbar * dN * dN) + 1.55) * np.mean(np.exp(kap * nrm ** 2 + math.sqrt(3) * (1 + dN) * nrm))
check("Thm 4(d): E[y_{t+2}] (MC) <= analytic bound with kappa_t, d_N", lhs + 4 * lhs_se <= rhs,
      f"LHS={lhs:.3f}+-{lhs_se:.3f}, RHS={rhs:.1f}")

# ---------------------------------------------------------------- Lemma A.3 (Poisson moments)
Bell = [1, 1, 2, 5, 15, 52]
ok = True
for lam in [0.3, 1.0, 5.0, 20.0]:
    kmax = int(lam + 30 * math.sqrt(lam) + 200)
    ks = np.arange(kmax); logp = ks * math.log(lam) - lam - gammaln(ks + 1)
    for p in [1.0, 2.0, 3.0, 2.5, 4.7]:
        lhs = np.sum(np.exp(logp) * (1 + ks) ** p)
        pc = math.ceil(p)
        Cp = sum(comb(pc, k) * Bell[k] for k in range(pc + 1))
        rhs = Cp ** (p / pc) * (1 + lam) ** p
        ok &= lhs <= rhs * (1 + 1e-10)
check("Lemma A.3: E(1+Y)^p <= C_p (1+lambda)^p (Bell-number bound)", ok)

# ---------------------------------------------------------------- Lemma A.4 (Poisson exp of log^2)
ok = True; worst_ratio = 0
for Lam in [0.5, 2.0, 10.0, 50.0]:
    kmax = int(Lam + 40 * math.sqrt(Lam) + 3000)
    ks = np.arange(kmax); logp = ks * math.log(Lam) - Lam - gammaln(ks + 1)
    lg = np.log1p(ks)
    for a in [0.0, 0.5, 1.0, 2.0, 3.0]:
        for q in [0.0, 0.1, 0.25, 0.375]:
            lhs = np.sum(np.exp(logp + a * lg + q * lg ** 2))
            ell = 2 + math.log(1 + Lam)
            rhs = math.exp(a * ell + q * ell ** 2) + 1.55 * math.exp(2 * a * a)
            ok &= lhs <= rhs * (1 + 1e-10)
            worst_ratio = max(worst_ratio, lhs / rhs)
check("Lemma A.4: E exp{a log(1+S)+q log^2(1+S)} <= exp(a l + q l^2)+1.55 e^{2a^2}", ok,
      f"worst LHS/RHS = {worst_ratio:.3f}")
# inequality steps inside the proof
ks = np.arange(0, 200000)
ok = np.all(np.log1p(ks) <= np.sqrt(ks) + 1e-12)
ok &= 1 / (math.exp(0.5) - 1) <= 1.55
for Lam in [0.3, 1, 7, 40]:
    K = math.ceil(math.e ** 2 * Lam)
    ok &= 1 + K <= math.e ** 2 * (1 + Lam) + 1e-9
    for k in range(K, K + 50):
        ok &= (math.e * Lam / k) ** k <= math.exp(-k) * (1 + 1e-12) if k > 0 else True
check("Lemma A.4 proof steps: log(1+k)<=sqrt(k), 1/(e^{1/2}-1)<=1.55, tail bounds", ok)

# ---------------------------------------------------------------- Lemma A.5 (Gaussian)
kap, b = 3.0, 2.0
P = np.diag([0.12, 0.05, 0.02]); m = np.array([0.5, -0.3, 0.2])
assert 2 * kap * np.linalg.eigvalsh(P).max() < 1
# closed form of the proof's dominating integral:
#   E exp{(k+e)||xi||^2 + 2k m'xi} = det(I-2(k+e)P)^{-1/2} exp(2 k^2 m'A^{-1}m), A=P^{-1}-2(k+e)I
def dom_integral(eps):
    Aq = np.linalg.inv(P) - 2 * (kap + eps) * np.eye(3)
    if np.linalg.eigvalsh(Aq).min() <= 0: return np.inf
    return np.linalg.det(np.eye(3) - 2 * (kap + eps) * P) ** -0.5 * math.exp(2 * kap ** 2 * m @ np.linalg.solve(Aq, m))
# (a) integral formula itself vs stable MC at a light-tailed exponent (identity is exponent-free)
xi = rng.normal(size=(2_000_000, 3)) @ np.linalg.cholesky(P).T
kap_s = 0.8
Aq_s = np.linalg.inv(P) - 2 * kap_s * np.eye(3)
cf = np.linalg.det(np.eye(3) - 2 * kap_s * P) ** -0.5 * math.exp(2 * kap ** 2 * m @ np.linalg.solve(Aq_s, m))
mc = np.exp(kap_s * np.sum(xi ** 2, 1) + 2 * kap * xi @ m).mean()
ok = abs(mc / cf - 1) < 0.02
# (b) pointwise inequality chain of the proof: k||th||^2+b||th|| <= c_e + (k+e)||xi||^2 + 2k m'xi
th = m + xi[:1_000_000]
lhs = kap * np.sum(th ** 2, 1) + b * np.linalg.norm(th, axis=1)
for eps in [0.05, 0.2, 0.5]:
    c_e = kap * m @ m + b * np.linalg.norm(m) + b * b / (4 * eps)
    rhs = c_e + (kap + eps) * np.sum(xi[:1_000_000] ** 2, 1) + 2 * kap * xi[:1_000_000] @ m
    ok &= np.all(lhs <= rhs + 1e-9)
# (c) hence E exp{k||th||^2+b||th||} <= min_e exp(c_e) * dom_integral(e) < infinity; MC below it
bound = min(math.exp(kap * m @ m + b * np.linalg.norm(m) + b * b / (4 * e)) * dom_integral(e)
            for e in [0.02, 0.05, 0.1, 0.2, 0.35, 0.49])
vals = np.exp(kap * np.sum((m + xi) ** 2, 1) + b * np.linalg.norm(m + xi, axis=1))
ok &= np.isfinite(bound) and vals.mean() <= bound
# (d) sharpness of the condition: 2 k lmax(P) >= 1 makes the dominating integral infinite
ok &= dom_integral((1 / (2 * P.max()) - kap) + 1e-9) == np.inf
check("Lemma A.5: closed-form bound verifies finiteness iff 2k lmax(P)<1",
      ok, f"MC/closed-form={mc/cf:.4f}; E<=bound: {vals.mean():.1f} <= {bound:.1f}")



# ================================================================================
# ADDITIONS: verification of the five results added in revision
# (Theorem 2(v)+Lemma A.6, Corollary floor, Prop. weighted, Prop. equiv, Thm 3(e))
# ================================================================================
from nssm import kalman_rw
from scipy.special import gammaln as _lg

def Pstar(I, q):
    return (math.sqrt(q*q*I*I + 4*q*I) - q*I) / (2*I)
def phi(v, I, q):
    return 1.0 / (1.0/(v+q) + I)

# ---- 32. Lemma A.6: fixed point, monotone convergence, bounds ----
ok = True
for _ in range(300):
    q = 10 ** rng.uniform(-5, 0); I = 10 ** rng.uniform(-2, 3)
    P = Pstar(I, q)
    ok &= abs(phi(P, I, q) - P) < 1e-12 * P
    for v0 in [0.0, P/3, 3*P + 1]:
        v = v0; prev = None; mono = True
        for t in range(400000):
            vn = phi(v, I, q)
            if prev is not None: mono &= (vn - v) * prev >= -1e-15
            step = vn - v; prev = step; v = vn
            if abs(step) < 1e-15 * max(P, 1e-300): break
        ok &= mono and abs(v - P) < 1e-7 * max(P, 1e-12)
    r = P / math.sqrt(q / I)
    ok &= r <= 1 + 1e-12
    if q * I <= 1: ok &= r >= 2/(1+math.sqrt(5)) - 1e-12
check("Lemma A.6 (Riccati): fixed point exact, monotone convergence, P*/sqrt(q/I) in [0.618,1]", ok)

# ---- 33. Thm 2(v): per-step envelope + envelope sequence on a simulated Kalman filter ----
ok = True
for rep in range(4):
    N = 100; A = er_graph(N, 4, rng); W = row_normalise(A)
    Q = np.diag([1e-4, 1e-3, 1e-3]); Lq = np.linalg.cholesky(Q); T = 80
    th = np.array([0.1, 0.3, 0.3]); y = rng.normal(0, 1, N)
    Xs, Ys, infos = [], [], []
    for t in range(T):
        th = th + rng.standard_normal(3) @ Lq.T
        X = np.column_stack([np.ones(N), W @ y, y])
        y = X @ th + rng.normal(0, 1, N)
        Xs.append(X); Ys.append(y); infos.append(profile_information(X, col=1)[0])
    fit = kalman_rw(np.array(Ys), Xs, 1.0, Q, np.zeros(3), np.eye(3))
    v_prev = 1.0; env = 1.0; Imin = min(infos)
    for t in range(T):
        v = fit["P_filt"][t, 1, 1]
        ok &= v <= phi(v_prev, infos[t], 1e-3) + 1e-10          # per-step inequality
        env = phi(env, Imin, 1e-3)
        ok &= v <= env + 1e-10                                   # dominated by envelope sequence
        v_prev = v
    ok &= fit["P_filt"][-1, 1, 1] <= Pstar(Imin, 1e-3) * 1.02 + 1e-9
check("Thm 2(v): filter variance obeys v_t <= phi_{I_t}(v_{t-1}) and the P* envelope", ok)

# ---- 34. Thm 2(v): exactness for orthogonal designs ----
ok = True
for I in [5.0, 40.0, 300.0]:
    N = 90
    v1 = np.zeros(N); v1[0] = 1; v1 -= v1.mean(); v1 /= np.linalg.norm(v1)
    v2 = np.zeros(N); v2[1] = 1
    v2 -= v2 @ np.ones(N)/N * np.ones(N); v2 -= (v2 @ v1) * v1; v2 /= np.linalg.norm(v2)
    X = np.column_stack([np.ones(N), math.sqrt(I) * v1, v2])
    G = X.T @ X
    ok &= np.max(np.abs(G - np.diag(np.diag(G)))) < 1e-10        # orthogonal columns
    ok &= abs(profile_information(X, col=1)[0] - I) < 1e-8 * I   # s(G)=I
    Q = np.diag([1e-4, 1e-3, 1e-3])
    fit = kalman_rw(np.zeros((250, N)), [X]*250, 1.0, Q, np.zeros(3), np.eye(3))
    v = 1.0; exact = True
    for t in range(250):
        v = phi(v, I, 1e-3)
        exact &= abs(fit["P_filt"][t, 1, 1] - v) < 1e-11
    ok &= exact and abs(fit["P_filt"][-1, 1, 1] - Pstar(I, 1e-3)) < 1e-10
check("Thm 2(v): orthogonal designs attain the scalar recursion and P* exactly", ok)

# ---- 35. Simulation claim: S1 sparse rows follow the sqrt(q/I) law (slope, levels) ----
import json
S1 = json.load(open(_os.path.join(_HERE, "out", "sim_results.json")))["S1"]
sp = [r for r in S1 if "sparse" in str(r.get("family", "")).lower()]
info = np.array([r["info"] for r in sp], float)
mse = np.array([r["mse_beta1"] for r in sp], float)
sd2 = np.array([r["post_sd_beta1"] for r in sp], float) ** 2
Ps = np.array([Pstar(i, 1e-3) for i in info])
sl_v = np.polyfit(np.log(info), np.log(sd2), 1)[0]
sl_m = np.polyfit(np.log(info), np.log(mse), 1)[0]
sl_P = np.polyfit(np.log(info), np.log(Ps), 1)[0]
ok = abs(sl_v - sl_P) < 0.03 and np.all(sd2 <= Ps * 1.02) and np.all(sd2 >= 0.80 * Ps) \
     and np.all(mse <= sd2) and (-0.66 <= sl_m <= -0.54)
check("S1 vs Thm 2(v): posterior variance tracks P* (slopes within 0.03), inside envelope at 82-90%; MSE below it",
      ok, f"sd2 slope={sl_v:.3f} vs P* slope={sl_P:.3f}; sd2/P* {min(sd2/Ps):.2f}-{max(sd2/Ps):.2f}; mse slope={sl_m:.3f}")

# ---- 36. Prop equiv (i)-(ii): span equality and exact coefficient map ----
ok = True
for _ in range(60):
    N = 40; W = row_normalise(er_graph(N, 4, rng)); Wc = row_normalise(complete_graph(N))
    x = rng.normal(0, 1, N); al = rng.uniform(0.05, 0.9)
    Wa = (1-al)*W + al*Wc
    X1 = np.column_stack([np.ones(N), W @ x, x]); X2 = np.column_stack([np.ones(N), Wa @ x, x])
    P1 = X1 @ np.linalg.pinv(X1); P2 = X2 @ np.linalg.pinv(X2)
    ok &= np.max(np.abs(P1 - P2)) < 1e-8                          # same column space
    b0, b1, b2 = rng.normal(0, 1, 3)
    ca = al*b1/((1-al)*(N-1))
    tb1, tb2, tb0 = b1/(1-al), b2 + ca, b0 - ca*N*x.mean()
    m1 = b0 + b1*(W @ x) + b2*x; m2 = tb0 + tb1*(Wa @ x) + tb2*x
    ok &= np.max(np.abs(m1 - m2)) < 1e-9 * (1 + np.max(np.abs(m1)))
check("Prop equiv (i)-(ii): span{1,W_a x,x}=span{1,Wx,x}; coefficient map reproduces the mean exactly", ok)

# ---- 37. Prop equiv (iii): B_a = B + c_a J; mu1 equal; mu2 shift formula ----
ok = True
for _ in range(60):
    N = 35; W = row_normalise(er_graph(N, 4, rng)); Wc = row_normalise(complete_graph(N))
    y = np.abs(rng.normal(1, 1, N)); al = rng.uniform(0.05, 0.9)
    b0, b1, b2 = 0.2, 0.5, 0.3
    ca = al*b1/((1-al)*(N-1)); J = np.ones((N, N))
    B = b1*W + b2*np.eye(N)
    Wa = (1-al)*W + al*Wc
    Ba = (b1/(1-al))*Wa + (b2+ca)*np.eye(N)
    ok &= np.max(np.abs(Ba - (B + ca*J))) < 1e-10
    tb0 = b0 - ca*N*y.mean()
    mu1 = b0 + B @ y; tmu1 = tb0 + Ba @ y
    ok &= np.max(np.abs(tmu1 - mu1)) < 1e-9
    mu2 = b0 + B @ mu1; tmu2 = tb0 + Ba @ tmu1
    pred = ca*N*(mu1.mean() - y.mean())
    ok &= np.max(np.abs((tmu2 - mu2) - pred)) < 1e-8
check("Prop equiv (iii): B_a=B+c_aJ, one-step means equal, two-step shift = c_a N (mu1bar-ybar) 1", ok)

# ---- 38. Prop weighted: identity, roughness, weight-free iff ----
ok = True
for _ in range(80):
    N = 30; W = row_normalise(er_graph(N, 3, rng))
    lam = np.exp(rng.normal(0, 0.7, N)); x = rng.normal(0, 1, N)
    X = np.column_stack([np.ones(N), W @ x, x]); Z = np.column_stack([np.ones(N), x])
    sl = np.sqrt(lam)
    MzWx = sl*(W @ x) - (sl[:, None]*Z) @ np.linalg.lstsq(sl[:, None]*Z, sl*(W @ x), rcond=None)[0]
    Iw = MzWx @ MzWx
    ok &= abs(Iw - 1.0/np.linalg.inv(X.T @ (X*lam[:, None]))[1, 1]) < 1e-7 * max(Iw, 1)
    MzRx = sl*((np.eye(N)-W) @ x) - (sl[:, None]*Z) @ np.linalg.lstsq(sl[:, None]*Z, sl*((np.eye(N)-W) @ x), rcond=None)[0]
    ok &= abs(Iw - MzRx @ MzRx) < 1e-7 * max(Iw, 1)
    Wc = row_normalise(complete_graph(N))
    Xc = np.column_stack([np.ones(N), Wc @ x, x])
    resid = sl*(Wc @ x) - (sl[:, None]*Z) @ np.linalg.lstsq(sl[:, None]*Z, sl*(Wc @ x), rcond=None)[0]
    ok &= resid @ resid < 1e-18 * (1 + x @ x)                    # complete graph: zero for ANY weights
    ok &= Iw > 1e-6                                              # non-complete: positive
check("Prop weighted: weighted identity + roughness form; complete-graph iff is weight-free", ok)

# ---- 39. Prop weighted (iii): Laplace-precision sandwich ----
ok = True
for _ in range(100):
    N = 30; W = row_normalise(er_graph(N, 3, rng))
    lam = np.exp(rng.normal(0, 0.7, N)); x = rng.normal(0, 1, N)
    X = np.column_stack([np.ones(N), W @ x, x]); Z = np.column_stack([np.ones(N), x])
    Amat = rng.normal(size=(3, 3)); A = Amat @ Amat.T + 0.3*np.eye(3)
    Pf = np.linalg.inv(A + X.T @ (X*lam[:, None]))
    sl = np.sqrt(lam)
    ustar = np.linalg.lstsq(sl[:, None]*Z, sl*(W @ x), rcond=None)[0]
    r = sl*(W @ x) - (sl[:, None]*Z) @ ustar; Iw = r @ r
    ok &= Pf[1, 1] <= 1.0/Iw + 1e-10
    ok &= Pf[1, 1] >= 1.0/(Iw + np.linalg.eigvalsh(A)[-1]*(1 + ustar @ ustar)) - 1e-10
check("Prop weighted (iii): (P_Laplace)_22 sandwiched by weighted information", ok)

# ---- 40. Thm 3(e) part 1: divergence with the actual constants + the chain on simulated paths ----
N = 6; W = row_normalise(er_graph(N, 3, rng)); i0 = 0
ell = np.zeros(N); ell[i0] = math.log(2.0)                       # y_{i,t}=1, others 0
gam = ell[i0] + (W @ ell)[i0]
s2 = 0.09; C0 = 1.0
from scipy.stats import norm as _norm
logf = []
for b in np.arange(5.0, 202.0, 4.0):
    # stable log lower bound on the box probability: min density on [b,b+1] times length 1
    lp_box = 2*(-(b+1)**2/(2*s2) - 0.5*math.log(2*math.pi*s2)) \
             + np.log(_norm.cdf(C0/math.sqrt(s2)) - _norm.cdf(-C0/math.sqrt(s2)))
    logf.append(lp_box + math.log(0.25) + gam*b**3/4)
logf = np.array(logf)
ok = np.all(np.diff(logf[25:]) > 0) and logf[-1] > 1e5
detail = f"log lower bound at b=201: {logf[-1]:.0f} (crossover where cubic beats Gaussian cost ~ b=64)"
for b in [6.0, 8.0]:
    th = np.array([0.0, b, b]); R = 4000
    lam1 = np.exp(th[0] + th[1]*(W @ ell) + th[2]*ell)
    y1 = rng.poisson(np.tile(lam1, (R, 1)))
    D1 = y1[:, i0] >= lam1[i0]/2
    ok &= D1.mean() >= 0.5
    l1 = np.log1p(y1[D1].astype(float))
    lam2 = np.exp(th[0] + th[1]*(l1 @ W.T) + th[2]*l1)
    # exact Poisson where numpy can sample it; normal approximation (error o(1) for the
    # half-intensity event) where lambda exceeds the sampler's range
    y2 = np.where(lam2 < 1e12, rng.poisson(np.minimum(lam2, 1e12)).astype(float),
                  np.maximum(lam2 + np.sqrt(lam2)*rng.standard_normal(lam2.shape), 0.0))
    D2 = y2[:, i0] >= lam2[:, i0]/2
    ok &= D2.mean() >= 0.5
    l2 = np.log1p(y2[D2].astype(float))
    loglam3 = th[0] + th[1]*(l2 @ W.T)[:, i0] + th[2]*l2[:, i0]
    ok &= np.min(loglam3) >= gam*b**3/4
check("Thm 3(e): lower-bound diverges with actual constants; chain events hold at rate >= 1/2 on paths",
      ok, detail)

# ---- 41. Thm 3(e) sharpness: series terms grow iff 2 s^2 L0 > 1; (d)-regime finite numerically ----
L0 = math.log(2.0)
def log_term(m, s2):
    bs = np.linspace(math.log(m)/L0, (math.log(m) + 1.0/m)/L0, 60)
    integ = -bs**2/(2*s2) - 0.5*math.log(2*math.pi*s2) - np.exp(bs*L0) + m*bs*L0 - _lg(m+1) + bs*np.log1p(m)
    return logsumexp(integ) + math.log(bs[1]-bs[0])
s2A = 1.0/L0                                                      # 2 s^2 L0 = 2 > 1
msA = np.array([10**k for k in range(2, 7)], dtype=float)
ltA = np.array([log_term(int(m), s2A) for m in msA])
coefA = np.polyfit(np.log(msA)**2, ltA, 1)[0]
target = (2*s2A*L0 - 1)/(2*s2A*L0**2)
okA = np.all(np.diff(ltA) > 0) and ltA[-1] > 50 and abs(coefA - target) < 0.15*target
s2B = 0.02                                                        # 2 kappa s^2 = 0.24 < 1: (d) finite
bg = np.linspace(-10*math.sqrt(s2B), 10*math.sqrt(s2B), 4001)
lamb = np.exp(bg*L0); inner = np.zeros_like(bg)
for j, b in enumerate(bg):
    ys = np.arange(0, 200)
    inner[j] = np.exp(logsumexp(-lamb[j] + ys*math.log(max(lamb[j], 1e-300)) - _lg(ys+1) + b*np.log1p(ys)))
EyB = np.trapezoid(inner*np.exp(-bg**2/(2*s2B))/math.sqrt(2*math.pi*s2B), bg)
okB = np.isfinite(EyB) and 0.5 < EyB < 5
check("Thm 3(e) sharpness: log term_m ~ c log^2 m with c=(2s^2L0-1)/(2s^2L0^2) when >1; (d)-regime integral finite",
      okA and okB, f"fitted c={coefA:.3f} vs {target:.3f}; E[y_(t+2)] under (d) = {EyB:.3f}")

n_fail = sum(1 for _, o in results if not o)
print(f"\n{len(results)} checks, {n_fail} failures")

# ---- additional checks (42-49) ----
def check_42_hub_tau():
    N=50; W=np.zeros((N,N)); W[0,1]=1.0
    for i in range(1,N): W[i,0]=1.0
    tau=np.sum(W*W)-np.sum(W.sum(0)**2)/N
    ok=abs(tau-(2-2/N))<1e-12
    return ok, f"directed hub tau_W={tau:.6f} vs 2-2/N={2-2/N:.6f}"
def check_43_riccati_large_qI():
    q,I=1.0,1e4
    P=(np.sqrt(q*q*I*I+4*q*I)-q*I)/(2*I)
    ok=abs(P*I-1)<2e-2
    return ok, f"large-qI regime: P*·I={P*I:.4f} (→1)"
def check_44_enlarged_leq_base():
    rng=np.random.default_rng(0); ok=True
    for _ in range(30):
        N=80; x=rng.normal(size=N); w=rng.normal(size=N); o=rng.normal(size=N); lam=rng.uniform(.2,3,N)
        def info(cols):
            X=np.column_stack(cols); G=X.T@(X*lam[:,None])
            return 1/np.linalg.inv(G)[1,1]
        base=info([np.ones(N),w,x]); enl=info([np.ones(N),w,x,o])
        ok&= enl<=base+1e-9
    return ok, "weighted profile info: enlarged nuisance <= base (30 random designs)"
def check_45_nb_concentration():
    from scipy.stats import nbinom, gamma
    ok=True; msg=[]
    for r in [0.5,2.0,10.0]:
        pr=0.5*(gamma.cdf(2,r,scale=1/r)-gamma.cdf(0.75,r,scale=1/r))
        for lam in [64,200,1000]:
            p=1-nbinom.cdf(int(lam/2)-1,r,r/(r+lam))
            ok&= p>=pr
            msg.append(f"r={r},lam={lam}: P(Y>=lam/2)={p:.3f}>=p_r={pr:.3f}")
    return ok, "; ".join(msg[:3])+" ..."
def check_46_nb_two_step_bound():
    from scipy.stats import nbinom, gamma
    ok=True
    for B,L0,r in [(4,np.log(4),2.0),(6,np.log(4),5.0)]:
        mu=np.exp(B*L0); pr=0.5*(gamma.cdf(2,r,scale=1/r)-gamma.cdf(0.75,r,scale=1/r))
        ks=np.arange(0,int(mu*8)); pmf=nbinom.pmf(ks,r,r/(r+mu))
        lhs=np.sum((1.0+ks.astype(float))**B*pmf)
        rhs=pr*(1+mu/2)**B
        ok&= lhs>=rhs
    return ok, "NB payoff bound E[(1+Y)^B] >= p_r (1+mu/2)^B at test points"
def check_47_induction_constants():
    ok=True
    for h,gam,C0 in [(4,1.0,1.0),(5,0.5,2.0)]:
        b=max(4*h,2**(h+3)*(C0+h+1)/gam)
        L=[ -(C0+2)+(b-2)*gam ]
        ok&= L[0]>=gam*b/2
        for k in range(1,h):
            nxt=-(C0+h)+(b-h)*(L[-1]-1 if k>1 else L[-1]-1)
            nxt=-(C0+h)+(b-h)*(gam*b**k/2**k-1)
            ok&= nxt>=gam*b**(k+1)/2**(k+1)
            L.append(nxt)
    return ok, "induction display verified numerically at (h,gamma,C0)=(4,1,1),(5,.5,2)"
def check_48_rewire_tau_varies():
    import json,os
    f=os.path.join(os.path.dirname(os.path.abspath(__file__)),'out','final_model.json')
    d=json.load(open(f))['rewire_tau']
    ok=(d['hi']-d['lo'])>1e-3 and abs(d['obs']-126.47)<0.05
    return ok, f"rewired tau_W spans [{d['lo']:.2f},{d['hi']:.2f}] vs obs {d['obs']:.2f}: not invariant"
def check_49_nb_h2_integrand():
    L0=np.log(4); s=np.sqrt(1.2/(2*L0))
    def part(a,b):
        B=np.linspace(a,b,4000)
        _tz=getattr(np,'trapezoid',getattr(np,'trapz',None))
        return _tz(np.exp(L0*B*B-B*np.log(2)-B*B/(2*s*s)),B)
    p1,p2,p3=part(2,6),part(6,10),part(10,14)
    ok=p2>5*p1 and p3>5*p2
    return ok, f"2s^2L0=1.2: partial integrals grow {p1:.2e},{p2:.2e},{p3:.2e} (divergent)"
_extra=[check_42_hub_tau,check_43_riccati_large_qI,check_44_enlarged_leq_base,
        check_45_nb_concentration,check_46_nb_two_step_bound,check_47_induction_constants,
        check_48_rewire_tau_varies,check_49_nb_h2_integrand]
_fails=0
for _f in _extra:
    _ok,_msg=_f()
    print(("[PASS] " if _ok else "[FAIL] ")+_f.__name__+": "+_msg)
    _fails+=(not _ok)
print(f"8 additional checks, {_fails} failures (49 total)")

def check_50_nb_h2_finiteness():
    L0=np.log(4); s=np.sqrt(0.8/(2*L0)); r=2.0
    def part(a,b):
        B=np.linspace(a,b,4000)
        _tz=getattr(np,'trapezoid',getattr(np,'trapz',None))
        return _tz(np.exp(L0*B*B+B*np.log(4*B/r)-B-B*B/(2*s*s)),B)
    p1,p2,p3=part(2,8),part(8,14),part(14,20)
    ok=p2<p1/3 and p3<p2/3
    return ok, f"2s^2L0=0.8: bound partial integrals shrink {p1:.2e},{p2:.2e},{p3:.2e} (convergent)"
_extra.append(check_50_nb_h2_finiteness)
_ok,_msg=check_50_nb_h2_finiteness()
print(("[PASS] " if _ok else "[FAIL] ")+"check_50_nb_h2_finiteness: "+_msg)
print(f"{'9' if _ok else '?'} additional checks total across suites (50 checks overall)")
