"""Negative-binomial mode filter (expected-information weights) and NB log-pmf."""
import numpy as np, math
from scipy.special import gammaln

def nb_logpmf(y, mu, r):
    mu = np.maximum(mu, 1e-300)
    return (gammaln(y + r) - gammaln(r) - gammaln(y + 1.0)
            + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))

def nb_mode_filter(Ymat, X_list, Q, m0, P0, r, iters=30):
    Tt, Nn = Ymat.shape; K = len(m0)
    m_pred = np.zeros((Tt, K)); m_filt = np.zeros((Tt, K)); P_filt = np.zeros((Tt, K, K))
    m, P = np.asarray(m0, float).copy(), np.asarray(P0, float).copy()
    for t in range(Tt):
        mp, Pp = m, P + Q
        X = X_list[t]; y = Ymat[t]
        Ppi = np.linalg.inv(Pp); th = mp.copy()
        for it in range(iters):
            eta = np.clip(X @ th, -30, 30); mu = np.exp(eta)
            g = X.T @ (y - (y + r) * mu / (r + mu)) - Ppi @ (th - mp)
            wts = r * mu / (r + mu)                       # expected information weights
            H = X.T @ (X * wts[:, None]) + Ppi
            step = np.linalg.solve(H, g)
            s = 1.0
            def obj(tt):
                e = np.clip(X @ tt, -30, 30); m_ = np.exp(e)
                return float(np.sum(y * e - (y + r) * np.log(r + m_)) - 0.5 * (tt - mp) @ Ppi @ (tt - mp))
            f_old = obj(th); th_new = th + step
            while s > 1e-4:
                th_new = th + s * step
                if obj(th_new) >= f_old - 1e-10: break
                s *= 0.5
            th = th_new
            if np.max(np.abs(s * step)) < 1e-9: break
        eta = np.clip(X @ th, -30, 30); mu = np.exp(eta)
        wts = r * mu / (r + mu)
        Pf = np.linalg.inv(X.T @ (X * wts[:, None]) + Ppi); Pf = 0.5 * (Pf + Pf.T)
        m_pred[t], m_filt[t], P_filt[t] = mp, th, Pf
        m, P = th, Pf
    return dict(m_pred=m_pred, m_filt=m_filt, P_filt=P_filt)

