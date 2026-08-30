"""
Simulation study for the rewritten paper.
S1: Theorems 1-2 -- profile information for beta_1 tracks tau_W = tr(W'M_1 W);
    filter MSE ~ 1/information; coverage of Kalman 95% intervals; three graph
    families (sparse, dense, complete).
S2: Proposition 3 -- plug-in h-step mean error is linear in ||W - W~||_op and is
    amplified geometrically in h when the propagation operator is not a contraction.
S3: Theorem 4 -- Poisson NSSM with raw-count regressors: Monte Carlo estimates of
    the h-step predictive mean diverge with the number of draws S for h>=3 (fixed
    coefficients) while the median path stays bounded; log(1+y) regressors do not.
Usage: python3 sim_study.py [S1] [S2] [S3]   (default: all three; results merged into out/sim_results.json)
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, time, sys, os
import numpy as np
sys.path.insert(0, _HERE)
from nssm import row_normalise, er_graph, complete_graph, kalman_rw, profile_information, graph_roughness, op_norm, tau_W
OUT = _os.path.join(_HERE, "out"); os.makedirs(OUT, exist_ok=True)
t0 = time.time()

def simulate_gaussian(N, W, T, rng, sigma=1.0, break_t=30):
    b0 = np.cumsum(rng.normal(0, 0.02, T)); b1 = np.where(np.arange(T) < break_t, 0.30, 0.60); b2 = np.full(T, 0.30)
    Y = np.zeros((T + 1, N)); Y[0] = rng.normal(0, 1, N)
    for t in range(1, T + 1):
        Y[t] = b0[t - 1] + b1[t - 1] * (W @ Y[t - 1]) + b2[t - 1] * Y[t - 1] + rng.normal(0, sigma, N)
    return Y, np.column_stack([b0, b1, b2])

def run_S1(R=100, T=60, sigma=1.0):
    rng = np.random.default_rng(1)
    Ns = [25, 50, 100, 200, 400, 800]
    families = ["sparse (d=4)", "dense (d=N/4)", "complete"]
    Q = np.diag([1e-4, 1e-3, 1e-3])
    rows = []
    for fam in families:
        for N in Ns:
            mse, cov, cov_ex, info, lam_min, rough, tau, psd = [], [], [], [], [], [], [], []
            for r in range(R):
                if fam.startswith("sparse"):
                    A = er_graph(N, 4, rng)
                elif fam.startswith("dense"):
                    A = er_graph(N, N / 4, rng)
                else:
                    A = complete_graph(N)
                W = row_normalise(A); tau.append(tau_W(W))
                Y, B = simulate_gaussian(N, W, T, rng, sigma)
                X_list = [np.column_stack([np.ones(N), W @ Y[t - 1], Y[t - 1]]) for t in range(1, T + 1)]
                fit = kalman_rw(Y[1:], X_list, sigma ** 2, Q, m0=np.zeros(3), P0=np.eye(3))
                e1 = fit["m_filt"][:, 1] - B[:, 1]
                sd1 = np.sqrt(fit["P_filt"][:, 1, 1])
                inside = np.abs(e1) <= 1.96 * sd1
                keep = np.arange(T) >= 10
                keep_ex = keep & ~((np.arange(T) >= 30) & (np.arange(T) < 36))
                mse.append(np.mean(e1[keep_ex] ** 2)); cov.append(np.mean(inside[keep])); cov_ex.append(np.mean(inside[keep_ex]))
                psd.append(np.mean(sd1[keep_ex]))
                infos = [profile_information(X_list[t], weights=np.full(N, 1 / sigma ** 2), col=1) for t in range(10, T)]
                info.append(np.mean([i[0] for i in infos])); lam_min.append(np.mean([i[1] for i in infos]))
                rough.append(np.mean([graph_roughness(W, Y[t - 1]) / sigma ** 2 for t in range(10, T)]))
            rows.append(dict(family=fam, N=N, tau_W=float(np.mean(tau)), info=float(np.mean(info)),
                             lam_min=float(np.mean(lam_min)), roughness=float(np.mean(rough)),
                             mse_beta1=float(np.mean(mse)), mse_beta1_se=float(np.std(mse) / np.sqrt(R)),
                             post_sd_beta1=float(np.mean(psd)),
                             cov95=float(np.mean(cov)), cov95_ex_break=float(np.mean(cov_ex))))
            print(f"[S1] {fam:14s} N={N:4d} tau={rows[-1]['tau_W']:7.2f} info={rows[-1]['info']:8.2f} "
                  f"lam_min={rows[-1]['lam_min']:8.2f} MSE={rows[-1]['mse_beta1']:.5f} psd={rows[-1]['post_sd_beta1']:.4f} cov={rows[-1]['cov95']:.3f} "
                  f"cov_ex={rows[-1]['cov95_ex_break']:.3f}  ({time.time()-t0:.0f}s)", flush=True)
    return rows

def run_S2():
    rng = np.random.default_rng(2)
    N, T = 100, 60
    A = er_graph(N, 4, rng); W = row_normalise(A)
    out = {}
    for regime, (b1, b2) in {"stable": (0.30, 0.35), "near-unit": (0.55, 0.50)}.items():
        Y = np.zeros((T + 1, N)); Y[0] = rng.normal(0, 1, N)
        for t in range(1, T + 1):
            Y[t] = 0.1 + b1 * (W @ Y[t - 1]) + b2 * Y[t - 1] + rng.normal(0, 1, N)
        B = b1 * W + b2 * np.eye(N); Bnorm = op_norm(B); rho = float(np.max(np.abs(np.linalg.eigvals(B))))
        Wc = np.ones((N, N)) / (N - 1); np.fill_diagonal(Wc, 0)
        res = []
        for alpha in [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:
            Wt = (1 - alpha) * W + alpha * Wc
            dW = op_norm(Wt - W)
            Bt = b1 * Wt + b2 * np.eye(N)
            errs = {}
            for h in [1, 2, 4, 8]:
                diff = []
                for s in range(20, T - h):
                    mu = Y[s].copy(); mut = Y[s].copy()
                    for k in range(h):
                        mu = 0.1 + B @ mu; mut = 0.1 + Bt @ mut
                    diff.append(np.sqrt(np.mean((mu - mut) ** 2)))
                errs[h] = float(np.mean(diff))
            res.append(dict(alpha=alpha, dW=dW, rmse=errs))
        out[regime] = dict(opnorm_B=float(Bnorm), rho_B=rho, b1=b1, b2=b2, rows=res)
        print(f"[S2] {regime}: ||B||={Bnorm:.3f} rho={rho:.3f} " + "; ".join(f"a={r['alpha']}:dW={r['dW']:.2f}:" + ",".join(f"h{h}={v:.3f}" for h, v in r['rmse'].items()) for r in res), flush=True)
    return out

def run_S3():
    """Poisson NSSM, N = 10 nodes, fixed coefficients theta = (-0.30, 0.25, 0.20) (calibrated to the
    Chicago fits on the raw scale), all nodes start at y_0 = 3.  For each specification we simulate S
    paths h = 1,...,5 steps ahead and record the Monte Carlo mean of lambda_{t+h} (average over nodes),
    its median across paths, and the fraction of paths in which some node has lambda > 1e6.
    '+state' adds random-walk coefficient innovations with variance q = 0.01 per step."""
    rng = np.random.default_rng(3)
    N = 10; A = er_graph(N, 3, rng); W = row_normalise(A)
    b0, b1, b2 = -0.30, 0.25, 0.20
    y0 = np.full(N, 3.0)
    out = {}
    for spec in ["raw", "raw+state", "log1p", "log1p+state"]:
        f = np.log1p if spec.startswith("log1p") else (lambda y: y)
        q = 0.01 if spec.endswith("+state") else 0.0
        res = {}
        for S in [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6]:
            r = np.random.default_rng(11)
            y = np.repeat(y0[None, :], S, 0)
            th = np.tile([b0, b1, b2], (S, 1)).astype(float)
            means, medians, tail = {}, {}, {}
            for h in range(1, 6):
                if h > 1 and q > 0:
                    th = th + r.normal(0, np.sqrt(q), (S, 3))
                fy = f(y)
                eta = th[:, [0]] + th[:, [1]] * (fy @ W.T) + th[:, [2]] * fy
                lam = np.exp(np.clip(eta, -700, 700))
                with np.errstate(over="ignore"):
                    means[h] = float(lam.mean())
                medians[h] = float(np.median(lam.mean(1))); tail[h] = float(np.mean(lam.max(1) > 1e6))
                if h < 5:
                    y = r.poisson(np.minimum(lam, 1e15)).astype(float)
            res[S] = dict(mean=means, median=medians, tail=tail)
            print(f"[S3] {spec:12s} S={S:7d} " + " ".join(f"h{h}:{means[h]:.4g}/{medians[h]:.3g}({tail[h]:.1e})" for h in means), flush=True)
        out[spec] = {str(k): v for k, v in res.items()}
    return out

which = sys.argv[1:] or ["S1", "S2", "S3"]
path = f"{OUT}/sim_results.json"
res = json.load(open(path)) if os.path.exists(path) else {}
if "S1" in which: res["S1"] = run_S1()
if "S2" in which: res["S2"] = run_S2()
if "S3" in which: res["S3"] = run_S3()
res["runtime"] = res.get("runtime", 0.0) + (time.time() - t0)
json.dump(res, open(path, "w"))
print("done", round(time.time() - t0), "s")
