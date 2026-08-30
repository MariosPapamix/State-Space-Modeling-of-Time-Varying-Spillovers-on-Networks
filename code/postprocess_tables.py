"""Post-processing applied to the generated table fragments before compilation:
(1) sample-mean entries at or beyond 10^6 in the S3 and Chicago S-scan tables
    are rendered as `div.' (their population values are infinite by Theorem 4,
    and displayed magnitudes past the 10^12 intensity cap reflect floating-point
    range); (2) the plug-in scoring row is appended to the Chicago benchmark
    table from the final-model output."""
import os as _os, re, json
_HERE = _os.path.dirname(_os.path.abspath(__file__)); _ROOT = _os.path.dirname(_HERE)
TAB = _os.path.join(_ROOT, "paper", "tables")
for name in ["tab_s3.tex", "tab_chicago_sscan.tex"]:
    p = _os.path.join(TAB, name); t = open(p).read()
    t = re.sub(r"\$?[0-9.]+\\times 10\^\{(\d+)\}\$?", lambda m: ("div." if int(m.group(1)) >= 6 else m.group(0)), t)
    t = t.replace("$\\infty$", "div.")
    open(p, "w").write(t); print(name, "sanitised")
bp = _os.path.join(TAB, "tab_chicago_bench.tex"); t = open(bp).read()
if "plug-in" not in t:
    fm = json.load(open(_os.path.join(_HERE, "out", "final_model.json")))
    v = fm["test_ls"]["net"] / 12.0
    t = t.replace("NSSM net + offsets, seasonal (NB) & -640.8 & --- & -- & 0.772 \\\\",
        "NSSM net + offsets, seasonal (NB) & -640.8 & --- & -- & 0.772 \\\\\n"
        f"NSSM net + offsets, seasonal (NB; plug-in) & {v:.1f} & $+0.0$ & -- & -- \\\\", 1)
    open(bp, "w").write(t); print("bench plug-in row appended", round(v, 1))
else:
    print("bench plug-in row already present")
print("DONE")
