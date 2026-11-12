"""Figure 2: MR forest plot of the selected (highlighted) method for each main and negative-control analysis.
Input : results/tables/MR_estimates_ST9_ST11_ST12.tsv (written by code/R/02_collect_results.R)
Output: figures/figure_2_mr_forest_plot.png (300 dpi) and .pdf
Usage : python code/python/figure2_mr_forest.py
"""
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / ".mplconfig"))
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PKG = Path(__file__).resolve().parents[2]
est = pd.read_csv(PKG / "results/tables/MR_estimates_ST9_ST11_ST12.tsv", sep="\t")
est = est[est["analysis_type"].isin(["main", "negative_control"])]
sel = est[est["selected"] == True].set_index("analysis")
assert sel.index.is_unique and len(sel) == 12

COL = {"Major Depressive Disorder": "#7570B3", "Short sleep duration": "#1B9E77",
       "Long sleep duration": "#D95F02", "Left handedness": "#66A61E",
       "Blond hair color": "#A6761D"}
left = {"Short sleep duration": [("Major Depressive Disorder", "short_to_MDD"),
                                 ("Left handedness", "short_to_lefthand"), ("Blond hair color", "short_to_blond")],
        "Long sleep duration": [("Major Depressive Disorder", "long_to_MDD"),
                                ("Left handedness", "long_to_lefthand"), ("Blond hair color", "long_to_blond")]}
right = {"Major Depressive Disorder": [("Short sleep duration", "MDD_to_short"), ("Long sleep duration", "MDD_to_long"),
                                       ("Left handedness", "MDD_to_lefthand_set_short"),   # MDD->short IV set
                                       ("Blond hair color", "MDD_to_blond_set_short")]}    # MDD->short IV set
FS = 12
rows = []


def panel(ax, groups, outcomes_order, width=0.2):
    for gi, (expo, items) in enumerate(groups.items()):
        n = len(items)
        for k, (out, a) in enumerate(items):
            r = sel.loc[a]
            x = gi + (k - (n - 1) / 2) * width
            ax.vlines(x, r.OR_lower95, r.OR_upper95, color=COL[out], lw=1.6)
            ax.plot(x, r.OR, "o", color=COL[out], ms=8)
            rows.append(dict(exposure=expo, outcome=out, analysis=a, method=r.method, nsnp=int(r.nsnp),
                             OR=r.OR, OR_lower95=r.OR_lower95, OR_upper95=r.OR_upper95, p_value=r.p_value))
    ax.axhline(1, ls="--", color="black", lw=1.2)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(list(groups.keys()), fontsize=FS)
    ax.set_xlim(-0.6, len(groups) - 0.4) if len(groups) > 1 else ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("Exposure", fontsize=FS)
    ax.set_ylabel("OR per 1-unit increase in exposure log-odds", fontsize=FS)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=FS)
    handles = [Line2D([0], [0], color=COL[o], marker="o", lw=1.6, ms=8) for o in outcomes_order]
    ax.legend(handles, outcomes_order, title="Outcome", frameon=False, fontsize=FS,
              title_fontsize=FS, loc="center left", bbox_to_anchor=(1.02, 0.6), alignment="left")


fig, axes = plt.subplots(1, 2, figsize=(20, 8), gridspec_kw={"width_ratios": [1, 1]})
panel(axes[0], left, ["Major Depressive Disorder", "Left handedness", "Blond hair color"])
panel(axes[1], right, ["Short sleep duration", "Long sleep duration", "Left handedness", "Blond hair color"], width=0.125)
for ax, lab in zip(axes, ["a", "b"]):
    ax.text(-0.12, 1.0, lab, transform=ax.transAxes, fontsize=24, fontweight="bold", va="center", ha="left")
plt.tight_layout(w_pad=6)
out = PKG / "figures"
out.mkdir(exist_ok=True)
plt.savefig(out / "figure_2_mr_forest_plot.png", dpi=300, bbox_inches="tight")
plt.savefig(out / "figure_2_mr_forest_plot.pdf", bbox_inches="tight")
pd.DataFrame(rows).to_csv(out / "figure_2_mr_forest_plot_values.tsv", sep="\t", index=False)
print("Figure 2 written:", out / "figure_2_mr_forest_plot.png")
