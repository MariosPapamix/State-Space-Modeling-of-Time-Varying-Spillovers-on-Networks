import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 9, "legend.fontsize": 7.5, "figure.dpi": 150})
OUT = _os.path.join(_HERE, "out"); FIG = _os.path.join(_ROOT, "paper", "figs"); os.makedirs(FIG, exist_ok=True)
S = json.load(open(f"{OUT}/sim_results.json")); C = json.load(open(f"{OUT}/chicago_results.json"))

# ---------------- Figure 1: identification scaling and filter accuracy (S1) -----------------
fig, ax = plt.subplots(1, 3, figsize=(9.2, 2.8))
fams = ["sparse (d=4)", "dense (d=N/4)", "complete"]
mk = {"sparse (d=4)": "o-", "dense (d=N/4)": "s--", "complete": "^:"}
for f in fams:
    rows = [r for r in S["S1"] if r["family"] == f]
    N = [r["N"] for r in rows]
    ax[0].plot(N, [max(r["info"], 1e-2) for r in rows], mk[f], label=f, ms=4)
    ax[0].plot(N, [max(r["tau_W"], 1e-2) for r in rows], color="grey", lw=0.9, ls="-", alpha=0.7)
    ax[1].plot(N, [r["mse_beta1"] for r in rows], mk[f], label=f, ms=4)
    ax[2].plot(N, [r["cov95_ex_break"] for r in rows], mk[f], label=f, ms=4)
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("$N$"); ax[0].set_ylabel(r"$\mathcal{I}_t(\beta_1)$ (time average)")
ax[0].set_title(r"(a) Profile information (grey: $\tau_W$)")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("$N$"); ax[1].set_ylabel(r"MSE of $\hat\beta_{1,t}$"); ax[1].set_title("(b) Filter MSE for $\\beta_{1,t}$")
ax[2].set_xscale("log"); ax[2].axhline(0.95, color="k", lw=0.7, ls="--"); ax[2].set_ylim(0.8, 1.0)
ax[2].set_xlabel("$N$"); ax[2].set_ylabel("coverage of 95% intervals"); ax[2].set_title("(c) Coverage (excl. break window)")
ax[0].legend(frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_sim_identification.pdf"); plt.close(fig)

# ---------------- Figure 2: multi-step plug-in error (S2) -----------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))
for i, (reg, d) in enumerate(S["S2"].items()):
    for h, m in zip([1, 2, 4, 8], ["o-", "s-", "^-", "d-"]):
        ax[i].plot([r["dW"] for r in d["rows"]], [r["rmse"][str(h)] for r in d["rows"]], m, label=f"$h={h}$", ms=4)
    ax[i].set_xlabel(r"$\|W-\widetilde W\|_{\mathrm{op}}$"); ax[i].set_ylabel("RMS discrepancy of $h$-step mean")
    ax[i].set_title(f"({'ab'[i]}) {reg}: $\\|B\\|_{{\\mathrm{{op}}}}$={d['opnorm_B']:.2f}, $\\rho(B)$={d['rho_B']:.2f}")
ax[0].legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_sim_plugin.pdf"); plt.close(fig)

# ---------------- Figure 3: Chicago coefficient paths and identification diagnostic -----------------
b = np.array(C["betas"]["net"]["m"]); sd = np.array(C["betas"]["net"]["sd"])
months = np.arange(1, C["data"]["T"])       # month index 1..71 (Feb 2010 .. Dec 2015)
xt = [i for i, m in enumerate(months) if m % 12 == 0]
lab = [str(2010 + m // 12) for m in months]
fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.7))
ax[0].plot(months, b[:, 1], color="C0"); ax[0].fill_between(months, b[:, 1] - 1.96 * sd[:, 1], b[:, 1] + 1.96 * sd[:, 1], color="C0", alpha=0.25)
ax[0].axvline(60, color="k", ls=":", lw=0.8); ax[0].set_title(r"(a) $\hat\beta_{1,t}$ (network spillover), 95% band"); ax[0].set_ylabel(r"$\hat\beta_{1,t}$")
ax[1].plot(months, b[:, 1] + b[:, 2], color="C3"); ax[1].axhline(1, color="k", ls="--", lw=0.8)
ax[1].axvline(60, color="k", ls=":", lw=0.8); ax[1].set_title(r"(b) $\hat\beta_{1,t}+\hat\beta_{2,t}$ (log-scale persistence)")
I = C["ident"]
ax[2].plot(I["month"], I["info_net"], color="C0", label="observed network $W$")
ax[2].plot(I["month"], np.maximum(I["info_complete"], 1e-3), color="C1", ls="--", label="complete graph $W_c$ (floor $10^{-3}$)")
ax[2].set_yscale("log"); ax[2].set_title(r"(c) Profile information $\mathcal{I}_t(\beta_1)$"); ax[2].legend(frameon=False, loc="center right")
for a in ax:
    a.set_xticks([months[i] for i in xt]); a.set_xticklabels([lab[i] for i in xt]); a.set_xlabel("month")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_chicago_paths.pdf"); plt.close(fig)

# ---------------- Figure 4: Chicago forecast comparison and PIT -----------------
fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.7))
H = C["horizons"]
for key, labl, m, off in [("net_vs_nonet", "vs. no-network DGLM", "o", -0.15), ("net_vs_complete", "vs. complete-graph DGLM", "s", 0.15)]:
    d = C["pairs"][key]
    y = [d[str(h)]["dls"] for h in H]; lo = [d[str(h)]["dls_ci"][0] for h in H]; hi = [d[str(h)]["dls_ci"][1] for h in H]
    ax[0].errorbar(np.array(H) + off, y, yerr=[np.array(y) - lo, np.array(hi) - y], fmt=m + "-", capsize=3, label=labl, ms=4)
    y = [d[str(h)]["dmae"] for h in H]; lo = [d[str(h)]["dmae_ci"][0] for h in H]; hi = [d[str(h)]["dmae_ci"][1] for h in H]
    ax[1].errorbar(np.array(H) + off, y, yerr=[np.array(y) - lo, np.array(hi) - y], fmt=m + "-", capsize=3, label=labl, ms=4)
ax[0].axhline(0, color="k", lw=0.7); ax[0].set_xlabel("horizon $h$"); ax[0].set_ylabel(r"$\Delta$ log score (network $-$ baseline)"); ax[0].set_title("(a) Log-score gain, 95% block-bootstrap CI"); ax[0].legend(frameon=False)
ax[1].axhline(0, color="k", lw=0.7); ax[1].set_xlabel("horizon $h$"); ax[1].set_ylabel(r"$\Delta$MAE (network $-$ baseline)"); ax[1].set_title("(b) MAE difference, 95% CI")
for key, labl, c in [("net", "network DGLM", "C0"), ("nonet", "no-network DGLM", "C1")]:
    cnt = np.array(C["pit_tests"][key]["counts"]); ax[2].step(np.linspace(0, 1, 11), np.r_[cnt, cnt[-1]] / cnt.sum() * 10, where="post", color=c, label=labl)
ax[2].axhline(1, color="k", ls="--", lw=0.7); ax[2].set_xlabel("randomised PIT"); ax[2].set_ylabel("relative frequency"); ax[2].set_title("(c) PIT histogram, $h=1$"); ax[2].set_ylim(0.7, 1.4); ax[2].legend(frameon=False, loc="upper center", ncol=2)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_chicago_forecast.pdf"); plt.close(fig)
print("figures written")
