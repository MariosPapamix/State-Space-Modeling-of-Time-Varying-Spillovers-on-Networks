"""Write LaTeX table fragments (paper/tables/*.tex) from out/chicago_results.json and out/sim_results.json."""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
import json, os, numpy as np
OUT = _os.path.join(_HERE, "out"); TAB = _os.path.join(_ROOT, "paper", "tables"); os.makedirs(TAB, exist_ok=True)
C = json.load(open(f"{OUT}/chicago_results.json")); S = json.load(open(f"{OUT}/sim_results.json"))

def w(name, body):
    body = body.rstrip()
    if body.endswith(r"\addlinespace"):
        body = body[: -len(r"\addlinespace")].rstrip()
    open(f"{TAB}/{name}.tex", "w").write(body + "\n")

def sci(x, d=1):
    if not np.isfinite(x): return r"$\infty$"
    if x == 0: return "0"
    e = int(np.floor(np.log10(abs(x)))); m = x / 10 ** e
    if abs(e) < 3: return f"{x:.{d}f}" if abs(x) < 100 else f"{x:.0f}"
    return f"${m:.{d}f}\\times 10^{{{e}}}$"

# ---------------- Table: S1 ----------------
rows = []
for fam, lab in [("sparse (d=4)", r"sparse ($d=4$)"), ("dense (d=N/4)", r"dense ($d=N/4$)"), ("complete", "complete")]:
    first = True
    for r in [x for x in S["S1"] if x["family"] == fam]:
        fz = lambda x: ("$<10^{-10}$" if abs(x) < 1e-10 else f"{x:.2f}")
        tau_s = f"{r['tau_W']:.3f}" if r['tau_W'] < 1 else f"{r['tau_W']:.2f}"
        rows.append(f"{lab if first else ''} & {r['N']} & {tau_s} & {fz(r['info'])} & {fz(r['lam_min'])} & "
                    f"{1e3*r['mse_beta1']:.2f} & {r['post_sd_beta1']:.3f} & {r['cov95_ex_break']:.3f} & {r['cov95']:.3f} \\\\")
        first = False
    rows.append(r"\addlinespace")
w("tab_s1", "\n".join(rows))

# ---------------- Table: S2 (appendix) ----------------
rows = []
for reg, d in S["S2"].items():
    rows.append(f"\\multicolumn{{6}}{{l}}{{\\emph{{{reg}}}: $\\beta_1={d['b1']:.2f}$, $\\beta_2={d['b2']:.2f}$, "
                f"$\\|B\\|_{{\\mathrm{{op}}}}={d['opnorm_B']:.3f}$, $\\rho(B)={d['rho_B']:.3f}$}}\\\\")
    for r in d["rows"]:
        rows.append(f"{r['alpha']:.2f} & {r['dW']:.3f} & " + " & ".join(f"{r['rmse'][str(h)]:.4f}" for h in [1, 2, 4, 8]) + r" \\")
    rows.append(r"\addlinespace")
w("tab_s2", "\n".join(rows))

# ---------------- Table: S3 ----------------
rows = []
labs = {"raw": r"raw counts, fixed $\theta$", "raw+state": r"raw counts, random-walk $\theta$ ($q=0.01$)",
        "log1p": r"$\log(1+y)$, fixed $\theta$", "log1p+state": r"$\log(1+y)$, random-walk $\theta$ ($q=0.01$)"}
for spec in ["raw", "raw+state", "log1p", "log1p+state"]:
    rows.append(f"\\multicolumn{{7}}{{l}}{{\\emph{{{labs[spec]}}}}}\\\\")
    for Sn, v in S["S3"][spec].items():
        cells = []
        for h in ["1", "2", "3", "4", "5"]:
            m, md, t = v["mean"][h], v["median"][h], v["tail"][h]
            if spec.startswith("raw"):
                cells.append(f"{sci(m, 2)} ({md:.2f})")
            else:
                cells.append(f"{m:.3f}")
        tail4 = v["tail"]["4"]
        rows.append(f"$10^{{{int(np.log10(int(Sn)))}}}$ & " + " & ".join(cells) + f" & {tail4:.3f} \\\\")
    rows.append(r"\addlinespace")
w("tab_s3", "\n".join(rows))

# ---------------- Table: Chicago accuracy ----------------
names = [("net", r"NSSM, observed network $W$"), ("nonet", "DGLM, no network"), ("complete", r"NSSM, complete graph $W_c$"),
         ("static_net", "static GLM with network"), ("static_nonet", "static GLM, no network"), ("static_season", "static region $\\times$ month")]
rows = []
for key, lab in names:
    s = C["summary"][key]
    cells = []
    for h in ["1", "2", "4", "8"]:
        cells.append(f"{s[h]['mae']:.3f} & {s[h]['ls']:.1f} & {s[h]['cov']:.3f}")
    rows.append(f"{lab} & " + " & ".join(cells) + r" \\")
w("tab_chicago_acc", "\n".join(rows))

# ---------------- Table: paired comparisons ----------------
comp = [("net_vs_nonet", "no-network DGLM"), ("net_vs_complete", "complete-graph NSSM"), ("net_vs_static_net", "static GLM with network"),
        ("net_vs_static_nonet", "static GLM, no network"), ("net_vs_static_season", "static region $\\times$ month")]
rows = []
for key, lab in comp:
    d = C["pairs"][key]
    first = True
    for h in ["1", "2", "4", "8"]:
        v = d[h]
        star = lambda lo, hi: r"$^{*}$" if (lo > 0 or hi < 0) else ""
        rows.append(f"{lab if first else ''} & {h} & ${v['dmae']:+.4f}$ & [{v['dmae_ci'][0]:+.4f}, {v['dmae_ci'][1]:+.4f}]{star(*v['dmae_ci'])} & {v['dmae_dm']:+.2f} & "
                    f"${v['dls']:+.1f}$ & [{v['dls_ci'][0]:+.1f}, {v['dls_ci'][1]:+.1f}]{star(*v['dls_ci'])} & {v['dls_dm']:+.2f} & {v['prop_targets_better_ls']:.2f} \\\\")
        first = False
    rows.append(r"\addlinespace")
w("tab_chicago_pairs", "\n".join(rows))

# ---------------- Table: stress / placebo ----------------
labs = {"original": r"observed $W$", "mix_uniform_0.25": r"$0.75\,W+0.25\,W_c$", "mix_uniform_0.5": r"$0.50\,W+0.50\,W_c$",
        "mix_uniform_0.75": r"$0.25\,W+0.75\,W_c$", "mix_uniform_1.0": r"$W_c$ (complete graph)", "edge_delete_0.2": "20\\% of edges deleted",
        "rewire_degseq": "degree-preserving rewiring", "permute_labels_0": "label permutation 1", "permute_labels_1": "label permutation 2",
        "permute_labels_2": "label permutation 3"}
rows = []
for key, lab in labs.items():
    v = C["stress"][key]
    star = lambda ci: r"$^{*}$" if (ci[0] > 0 or ci[1] < 0) else ""
    info = v["mean_profile_info"]; info_s = f"{info:.1f}" if info > 0.05 else sci(info, 1)
    rows.append(f"{lab} & {v['opnorm_diff']:.2f} & {v['tau']:.1f} & {info_s} & ${v['mean_beta1']:+.2f}$ ({v['mean_sd_beta1']:.3f}) & "
                f"${v['dls_h1']:+.2f}$ [{v['dls_h1_ci'][0]:+.1f}, {v['dls_h1_ci'][1]:+.1f}]{star(v['dls_h1_ci'])} & "
                f"${v['dls_h2']:+.2f}$ [{v['dls_h2_ci'][0]:+.1f}, {v['dls_h2_ci'][1]:+.1f}]{star(v['dls_h2_ci'])} \\\\")
w("tab_chicago_stress", "\n".join(rows))

# ---------------- Table: S-scan raw vs log1p ----------------
rows = []
for spec, lab in [("raw", r"raw counts"), ("log1p", r"$\log(1+y)$")]:
    rows.append(f"\\multicolumn{{7}}{{l}}{{\\emph{{{lab}}}}}\\\\")
    for h in ["1", "2", "3", "4"]:
        d = C["sscan"][spec][h]
        cells = [f"{sci(d[Sn]['mean_total'], 2)}" for Sn in ["100", "1000", "5000"]]
        med = d["5000"]["median_total"]; mx = d["5000"]["max_total"]; fr = d["5000"]["frac_gt_1e4"]
        rows.append(f"{h} & " + " & ".join(cells) + f" & {med:.0f} & {sci(mx, 1)} & {fr:.3f} \\\\")
    rows.append(r"\addlinespace")
w("tab_chicago_sscan", "\n".join(rows))

# raw per-horizon summary line (MAE etc.) for the text
r = C["summary"]["raw_net"]
w("raw_summary", f"% raw_net: h1 MAE {r['1']['mae']:.3f} LS {r['1']['ls']:.1f}; h2 MAE {r['2']['mae']:.3f} LS {r['2']['ls']:.1f} maxlam {r['2']['max_lam']:.0f}; h4 MAE {r['4']['mae']:.3g} MSE {r['4']['mse']:.3g} explosive {r['4']['explosive']:.3f} maxlam {r['4']['max_lam']:.2g}\n")

# ---------------- Table: tuning (appendix) ----------------
rows = []
qs = ["0.0001", "0.0003", "0.001", "0.003", "0.01", "0.03"]
for key, lab in [("net", r"NSSM, observed $W$ ($\log(1+y)$)"), ("nonet", "DGLM, no network"), ("complete", r"NSSM, complete graph"),
                 ("raw", r"NSSM, observed $W$ (raw counts)")]:
    tab = C["tune_tables"][key] if key != "raw" else C["tune_table_raw"]
    best = max(tab, key=lambda q: tab[q])
    cells = [(f"\\textbf{{{tab[q]:.0f}}}" if q == best else f"{tab[q]:.0f}") for q in qs]
    rows.append(f"{lab} & " + " & ".join(cells) + r" \\")
w("tab_chicago_tuning", "\n".join(rows))

# ---------------- Table: PIT counts (appendix) ----------------
rows = []
for key, lab in [("net", r"NSSM, observed $W$"), ("nonet", "DGLM, no network"), ("complete", r"NSSM, complete graph")]:
    v = C["pit_tests"][key]
    rows.append(f"{lab} & " + " & ".join(str(c) for c in v["counts"]) + f" & {v['chi2']:.1f} \\\\")
w("tab_chicago_pit", "\n".join(rows))


# ---------------- Table: raw-count model accuracy (appendix) ----------------
rows = []
for key, lab in [("net", r"$\log(1+y)$ feedback"), ("raw_net", "raw-count feedback")]:
    s = C["summary"][key]
    for h in ["1", "2", "4"]:
        v = s[h]
        rows.append(f"{lab if h == '1' else ''} & {h} & {sci(v['mae'], 2) if v['mae'] > 100 else f'{v[chr(109)+chr(97)+chr(101)]:.3f}'} & {v['ls']:.1f} & {v['explosive']:.3f} & {sci(v['max_lam'], 1)} \\\\")
    rows.append(r"\addlinespace")
w("tab_chicago_raw", "\n".join(rows))

# ---------------- numbers used inline (macro file) ----------------
d = C["data"]; I = C["ident"]; b = np.array(C["betas"]["net"]["m"]); sd = np.array(C["betas"]["net"]["sd"]); rho = b[:, 1] + b[:, 2]
rw = C["regionwise"]
macros = {
    "NumRegions": d["N"], "NumMonths": d["T"], "NumEdges": d["n_edges"], "DegMean": f"{d['deg_mean']:.2f}", "DegMax": int(d["deg_max"]), "DegMin": int(d["deg_min"]),
    "MeanCount": f"{d['mean_count']:.2f}", "ZeroFrac": f"{100*d['zero_frac']:.1f}", "MaxCount": d["max_count"],
    "TotalFirst": d["monthly_total_first"], "TotalLast": d["monthly_total_last"],
    "TauW": f"{d['tau_W']:.1f}", "TauWc": f"{d['tau_Wc']:.4f}", "HarmonicSum": f"{d['tau_regular_bound']:.1f}",
    "InfoNetMean": f"{np.mean(I['info_net']):.1f}", "InfoNetMin": f"{np.min(I['info_net']):.1f}", "InfoNetMax": f"{np.max(I['info_net']):.1f}",
    "InfoCompMax": sci(np.max(I["info_complete"]), 1), "RoughMean": f"{np.mean(I['roughness']):.2f}",
    "LamMinMean": f"{np.mean(I['lam_min_net']):.1f}", "LamMinMin": f"{np.min(I['lam_min_net']):.1f}",
    "BetaOneFirst": f"{b[0,1]:.2f}", "BetaOneLast": f"{b[-1,1]:.2f}", "BetaOneMin": f"{b[:,1].min():.2f}", "BetaOneMax": f"{b[:,1].max():.2f}",
    "BetaTwoMean": f"{b[12:,2].mean():.2f}", "BetaZeroMean": f"{b[12:,0].mean():.2f}",
    "SdBetaOneMean": f"{sd[12:,1].mean():.3f}", "SdBetaOneMax": f"{sd[12:,1].max():.3f}",
    "RhoMin": f"{rho.min():.2f}", "RhoMax": f"{rho.max():.2f}", "RhoFirstAboveOne": int(np.argmax(rho > 1)) + 1,
    "RhoLastMean": f"{rho[-12:].mean():.2f}",
    "RegQten": f"{rw['q10']:+.3f}", "RegQtf": f"{rw['q25']:+.3f}", "RegMed": f"{rw['median']:+.4f}", "RegQsf": f"{rw['q75']:+.3f}", "RegQninety": f"{rw['q90']:+.3f}",
    "RegPropBetter": f"{100*rw['prop_better']:.1f}",
    "PitN": C["pit_tests"]["net"]["n"], "PitChiNet": f"{C['pit_tests']['net']['chi2']:.0f}", "PitChiNonet": f"{C['pit_tests']['nonet']['chi2']:.0f}", "PitChiComp": f"{C['pit_tests']['complete']['chi2']:.0f}",
    "QNet": C["tuned_q"]["net"], "QNonet": C["tuned_q"]["nonet"], "QComp": C["tuned_q"]["complete"], "QRaw": C["q_raw"],
    "RuntimeChicago": f"{C['runtime_sec']/60:.1f}", "RuntimeSim": f"{S['runtime']/60:.1f}",
    "SdBetaCompMean": f"{np.array(C['betas']['complete']['sd'])[12:,1].mean():.3f}",
    "BetaRawMean": ", ".join(f"{x:.2f}" for x in np.array(C["betas"]["raw_net"]["m"])[12:].mean(0)),
}
# Theorem 4(d) margin for the Chicago fit: kappa_t and 2 kappa lambda_max(P)
Lt = np.log1p(d["max_count"]); ct = np.sqrt(3) * (1 + Lt); qbar = C["tuned_q"]["net"]
kappa = max(3.0, np.sqrt(3) * ct + 2 * qbar * ct ** 2)
trP = (sd[12:] ** 2).sum(1).max() + qbar          # upper bound on lambda_max(P_{t+1|t}) via trace
macros.update({"KappaChicago": f"{kappa:.1f}", "LamMaxPChicago": f"{trP:.4f}", "MarginChicago": f"{2*kappa*trP:.2f}", "LtChicago": f"{Lt:.2f}"})
lines = [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in macros.items()]
w("numbers", "\n".join(lines) + "\n")
print("tables written:", sorted(os.listdir(TAB)))
print(macros)
