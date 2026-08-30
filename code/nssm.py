"""
nssm.py -- Network state-space models (NSSM): filters and diagnostics.

Gaussian random-walk Kalman filter, Poisson Laplace (mode/curvature) filter,
the spillover-identification diagnostic (profile information for beta_1),
and predictive-distribution utilities.  Pure numpy/scipy.
"""
import numpy as np
from scipy.special import gammaln, logsumexp


# ----------------------------------------------------------------------------
# Network utilities
# ----------------------------------------------------------------------------
def row_normalise(A):
    """Row-stochastic weight matrix W = D^{-1} A; zero rows stay zero."""
    A = np.asarray(A, dtype=float).copy()
    np.fill_diagonal(A, 0.0)
    d = A.sum(1)
    W = np.zeros_like(A)
    nz = d > 0
    W[nz] = A[nz] / d[nz, None]
    return W


def er_graph(N, expected_degree, rng, directed=False):
    p = min(1.0, expected_degree / (N - 1))
    A = (rng.random((N, N)) < p).astype(float)
    if not directed:
        A = np.triu(A, 1)
        A = A + A.T
    np.fill_diagonal(A, 0)
    # guarantee no isolated node (connect to a random other node)
    iso = np.where(A.sum(1) == 0)[0]
    for i in iso:
        j = rng.integers(N)
        while j == i:
            j = rng.integers(N)
        A[i, j] = A[j, i] = 1
    return A


def complete_graph(N):
    A = np.ones((N, N))
    np.fill_diagonal(A, 0)
    return A


def op_norm(M):
    return np.linalg.norm(M, 2)


def tau_W(W):
    """tau_W = tr(W' M_1 W), M_1 = I - 11'/N  (noise-driven information about beta_1)."""
    N = W.shape[0]
    M1 = np.eye(N) - np.ones((N, N)) / N
    return float(np.trace(W.T @ M1 @ W))


# ----------------------------------------------------------------------------
# Identification diagnostic (Theorem 1): profile information for beta_1
# ----------------------------------------------------------------------------
def profile_information(X, weights=None, col=1):
    """
    Conditional Fisher information for the coefficient in column `col` of the
    design X (N x K) after profiling out all other coefficients, for a
    (generalised) linear model with observation weights `weights`
    (Gaussian: 1/sigma^2; Poisson log-link: lambda_i).
    Equals  || M_{-col} sqrt(Omega) x_col ||^2 = 1 / (G^{-1})_{col,col},
    G = X' Omega X.  Returns (profile information, lambda_min(G)).
    """
    X = np.asarray(X, float)
    N, K = X.shape
    w = np.ones(N) if weights is None else np.asarray(weights, float)
    Xw = X * np.sqrt(w)[:, None]
    G = Xw.T @ Xw
    others = [j for j in range(K) if j != col]
    Xo = Xw[:, others]
    xc = Xw[:, col]
    coef, *_ = np.linalg.lstsq(Xo, xc, rcond=None)
    resid = xc - Xo @ coef
    info = float(resid @ resid)
    lam_min = float(np.linalg.eigvalsh(G).min())
    return info, lam_min


def graph_roughness(W, y):
    """|| (I - W) y ||^2 : roughness of the lagged field along the graph."""
    r = y - W @ y
    return float(r @ r)


# ----------------------------------------------------------------------------
# Gaussian random-walk Kalman filter
# ----------------------------------------------------------------------------
def kalman_rw(Y, X_list, sigma2, Q, m0, P0):
    """
    Y: T x N, X_list: list of T arrays (N x K), theta_t = theta_{t-1} + u_t,
    u_t ~ N(0, Q), Y_t = X_t theta_t + eps_t, eps_t ~ N(0, sigma2 I).
    Returns dict with m_pred, P_pred, m_filt, P_filt, yhat, loglik (per t),
    and G_t = X_t' X_t / sigma2 (information matrices).
    """
    T, N = Y.shape
    K = len(m0)
    m_pred = np.zeros((T, K)); P_pred = np.zeros((T, K, K))
    m_filt = np.zeros((T, K)); P_filt = np.zeros((T, K, K))
    yhat = np.zeros((T, N)); loglik = np.zeros(T)
    G = np.zeros((T, K, K))
    m, P = np.asarray(m0, float).copy(), np.asarray(P0, float).copy()
    for t in range(T):
        mp, Pp = m, P + Q
        X = X_list[t]
        yh = X @ mp
        resid = Y[t] - yh
        Gt = X.T @ X / sigma2
        Pinv = np.linalg.inv(Pp) + Gt
        Pf = np.linalg.inv(Pinv)
        mf = Pf @ (np.linalg.solve(Pp, mp) + X.T @ Y[t] / sigma2)
        logdetS = N * np.log(sigma2) + np.linalg.slogdet(np.eye(K) + Pp @ Gt)[1]
        Xr = X.T @ resid / sigma2
        quad = (resid @ resid) / sigma2 - Xr @ Pf @ Xr
        loglik[t] = -0.5 * (N * np.log(2 * np.pi) + logdetS + quad)
        m_pred[t], P_pred[t], m_filt[t], P_filt[t] = mp, Pp, mf, Pf
        yhat[t] = yh; G[t] = Gt
        m, P = mf, Pf
    return dict(m_pred=m_pred, P_pred=P_pred, m_filt=m_filt, P_filt=P_filt,
                yhat=yhat, loglik=loglik, G=G)


# ----------------------------------------------------------------------------
# Poisson NSSM: Laplace (mode + curvature) filter with random-walk state
# ----------------------------------------------------------------------------
def poisson_laplace_filter(Y, X_list, Q, m0, P0, newton_iters=25, eta_clip=30.0):
    """
    Y_it | theta_t ~ Poisson(exp(X_t theta_t)_i), theta_t = theta_{t-1} + u_t.
    At each t: Gaussian prior N(m_{t|t-1}, P_{t|t-1}); posterior approximated
    by N(mode, inverse negative Hessian) (Laplace / iterated EKF).
    """
    T, N = Y.shape
    K = len(m0)
    m_pred = np.zeros((T, K)); P_pred = np.zeros((T, K, K))
    m_filt = np.zeros((T, K)); P_filt = np.zeros((T, K, K))
    info = np.zeros(T)      # profile information for beta_1 at the mode
    m, P = np.asarray(m0, float).copy(), np.asarray(P0, float).copy()
    for t in range(T):
        mp, Pp = m, P + Q
        X = X_list[t]; y = Y[t]
        Pp_inv = np.linalg.inv(Pp)
        th = mp.copy()
        for it in range(newton_iters):
            eta = np.clip(X @ th, -eta_clip, eta_clip)
            lam = np.exp(eta)
            g = X.T @ (y - lam) - Pp_inv @ (th - mp)
            H = X.T @ (X * lam[:, None]) + Pp_inv
            step = np.linalg.solve(H, g)
            s = 1.0
            f_old = (y @ eta - lam.sum()) - 0.5 * (th - mp) @ Pp_inv @ (th - mp)
            th_new = th + step
            while s > 1e-4:
                th_new = th + s * step
                eta_n = np.clip(X @ th_new, -eta_clip, eta_clip)
                f_new = (y @ eta_n - np.exp(eta_n).sum()) - 0.5 * (th_new - mp) @ Pp_inv @ (th_new - mp)
                if f_new >= f_old - 1e-12:
                    break
                s *= 0.5
            th = th_new
            if np.max(np.abs(s * step)) < 1e-9:
                break
        eta = np.clip(X @ th, -eta_clip, eta_clip); lam = np.exp(eta)
        H = X.T @ (X * lam[:, None]) + Pp_inv
        Pf = np.linalg.inv(H); Pf = 0.5 * (Pf + Pf.T)
        m_pred[t], P_pred[t], m_filt[t], P_filt[t] = mp, Pp, th, Pf
        if K >= 2:
            info[t] = profile_information(X, weights=lam, col=1)[0]
        m, P = th, Pf
    return dict(m_pred=m_pred, P_pred=P_pred, m_filt=m_filt, P_filt=P_filt, info=info)


def poisson_logpmf(y, lam):
    lam = np.maximum(lam, 1e-300)
    return y * np.log(lam) - lam - gammaln(y + 1.0)


# ----------------------------------------------------------------------------
# Predictive distributions by Monte Carlo (Poisson NSSM)
# ----------------------------------------------------------------------------
def poisson_predict(m, P, Q, h, S, rng, build_X, y_last, transform=np.log1p, eta_clip=30.0):
    """
    Simulate S coefficient paths h steps ahead from the filtered state (m, P)
    (Gaussian random walk with innovation covariance Q) and future counts.
    build_X(transformed_prev_counts) -> N x K design.  y_last: counts at the
    origin.  Returns lam_paths (S x N) at horizon h and the coefficient paths.
    """
    K = len(m)
    # theta_{t+1} | F_t ~ N(m_{t|t}, P_{t|t} + Q); each further step adds an independent N(0, Q) innovation
    L = np.linalg.cholesky(P + Q + 1e-12 * np.eye(K))
    th = m[None, :] + rng.standard_normal((S, K)) @ L.T
    Lq = np.linalg.cholesky(Q + 1e-14 * np.eye(K))
    y_prev = np.repeat(y_last[None, :], S, axis=0).astype(float)
    lam = None
    coef_paths = np.zeros((h, S, K))
    for k in range(h):
        if k > 0:
            th = th + rng.standard_normal((S, K)) @ Lq.T
        coef_paths[k] = th
        lam = np.zeros_like(y_prev)
        for s in range(S):
            X = build_X(transform(y_prev[s]))
            lam[s] = np.exp(np.clip(X @ th[s], -eta_clip, eta_clip))
        if k < h - 1:
            y_prev = rng.poisson(np.minimum(lam, 1e12)).astype(float)
    return lam, coef_paths


def mixture_scores(y, lam_paths, rng, alpha=0.10):
    """
    Given counts y (N,) and S x N simulated intensities at the target horizon,
    compute: predictive mean, sum over nodes of marginal log scores
    log( mean_s Pois(y_i | lam_si) ), randomized PIT values, and coverage of
    central (1-alpha) intervals computed from the Poisson mixture.
    """
    S, N = lam_paths.shape
    lp = poisson_logpmf(y[None, :], lam_paths)             # S x N
    marg_ls = logsumexp(lp, axis=0) - np.log(S)             # N
    mean = lam_paths.mean(0)
    ymax = int(min(max(y.max(), np.quantile(lam_paths, 0.999).max() * 3 + 20), 5000))
    grid = np.arange(0, ymax + 1)
    pmf = np.zeros((N, len(grid)))
    for s in range(S):
        pmf += np.exp(poisson_logpmf(grid[None, :], np.minimum(lam_paths[s], 1e6)[:, None]))
    pmf /= S
    cdf = np.cumsum(pmf, axis=1)
    idx = np.minimum(y, ymax)
    F_y = cdf[np.arange(N), idx]
    F_ym1 = np.where(y > 0, cdf[np.arange(N), np.maximum(idx - 1, 0)], 0.0)
    u = rng.random(N)
    pit = F_ym1 + u * (F_y - F_ym1)
    lo = (cdf >= alpha / 2).argmax(axis=1)
    hi = (cdf >= 1 - alpha / 2).argmax(axis=1)
    covered = (y >= lo) & (y <= hi)
    return dict(mean=mean, logscore_sum=float(marg_ls.sum()), pit=pit,
                coverage=float(covered.mean()), lo=lo, hi=hi)


# ----------------------------------------------------------------------------
# Static Poisson GLM (IRLS) for baselines
# ----------------------------------------------------------------------------
def poisson_glm(X, y, iters=50, ridge=1e-8):
    X = np.asarray(X, float); y = np.asarray(y, float)
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.mean(), 1e-3))
    for _ in range(iters):
        eta = np.clip(X @ beta, -30, 30); lam = np.exp(eta)
        g = X.T @ (y - lam) - ridge * beta
        H = X.T @ (X * lam[:, None]) + ridge * np.eye(X.shape[1])
        step = np.linalg.solve(H, g)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-8:
            break
    return beta


def two_way_loglinear(Ymat, iters=200):
    """
    Poisson model log lambda_{t i} = a_i + b_{month(t)} fitted by iterative
    proportional fitting (exact MLE for the additive log-linear model).
    Ymat: T x N counts (training window).  Returns (a, b) with b of length 12.
    """
    T, N = Ymat.shape
    months = np.arange(T) % 12
    a = np.log(np.maximum(Ymat.mean(0), 1e-3)); b = np.zeros(12)
    for _ in range(iters):
        for mth in range(12):
            rows = months == mth
            denom = np.exp(a).sum() * rows.sum()
            b[mth] = np.log(max(Ymat[rows].sum(), 1e-6) / denom)
        a = np.log(np.maximum(Ymat.sum(0), 1e-6) / np.exp(b[months]).sum())
    return a, b
