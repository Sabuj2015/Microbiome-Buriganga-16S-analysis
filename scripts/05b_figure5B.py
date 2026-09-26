#!/usr/bin/env python3
"""Figure 5B: extended error-bar plot, sediment vs water (Welch's two-sided t-test).

The 65 genera and their order are those of the published panel (data/figure5B_genera.csv).
For each genus, the mean relative abundance in each habitat, the difference between means with
its 95% confidence interval (Welch) and the p-value are calculated from per-sample relative
abundances (normalized OTU table if present, otherwise the OTU table). The p-values are not
adjusted for multiple testing; adjusted results are in Additional file 1: Table S2.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIG = HERE.parent / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SED, WAT = "#f4a261", "#5b7bb5"


def welch(s, w):
    t = stats.ttest_ind(s, w, equal_var=False)
    vs, vw = s.var(ddof=1) / len(s), w.var(ddof=1) / len(w)
    df = (vs + vw) ** 2 / (vs ** 2 / (len(s) - 1) + vw ** 2 / (len(w) - 1))
    diff = s.mean() - w.mean()
    half = stats.t.ppf(0.975, df) * np.sqrt(vs + vw)
    return diff, diff - half, diff + half, t.pvalue


def main():
    norm_file = DATA / "otu_table_normalized.csv"
    otu = pd.read_csv(norm_file if norm_file.exists() else DATA / "otu_table.csv")
    meta = pd.read_csv(DATA / "sample_metadata.csv").set_index("Sample_ID")
    samples = list(meta.index)
    genus = otu.Genus.fillna("").replace("", "Unassigned")
    prop = otu[samples].groupby(genus).sum() / otu[samples].sum()           # proportions, as plotted
    pub = pd.read_csv(DATA / "figure5B_genera.csv")
    sed, wat = meta.index[meta.Habitat == "Sediment"], meta.index[meta.Habitat == "Water"]
    rows = []
    for g in pub.Genus:
        s, w = prop.loc[g, sed].to_numpy(), prop.loc[g, wat].to_numpy()
        d, lo, hi, p = welch(s, w)
        rows.append(dict(Genus=g, Mean_sediment=s.mean(), Mean_water=w.mean(), Difference=d, CI_low=lo, CI_high=hi, p=p))
    res = pd.DataFrame(rows)
    res["p_display"] = [f"{p:.3f}" if p >= 0.001 else "<0.001" for p in res.p]
    res.to_csv(FIG / "Figure_5B_welch_recomputed.csv", index=False)

    n = len(res)
    fig = plt.figure(figsize=(9, 0.22 * n + 1.5))
    ax1 = fig.add_axes([0.28, 0.06, 0.2, 0.9])
    ax2 = fig.add_axes([0.52, 0.06, 0.3, 0.9])
    y = np.arange(n)[::-1]
    for i, yy in enumerate(y):
        if i % 2 == 0:
            for ax in (ax1, ax2):
                ax.axhspan(yy - 0.5, yy + 0.5, color="#e6e6e6", zorder=0)
    ax1.barh(y + 0.18, res.Mean_sediment, 0.36, color=SED, edgecolor="#333333", linewidth=0.5, label="Sediment")
    ax1.barh(y - 0.18, res.Mean_water, 0.36, color=WAT, edgecolor="#333333", linewidth=0.5, label="Water")
    ax1.set_yticks(y)
    ax1.set_yticklabels(res.Genus, fontsize=7)
    ax1.set_xlim(0, 0.06)
    ax1.set_xlabel("Means in groups")
    ax1.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, frameon=False, fontsize=7)
    colours = [SED if d > 0 else WAT for d in res.Difference]
    ax2.errorbar(res.Difference, y, xerr=[res.Difference - res.CI_low, res.CI_high - res.Difference], fmt="none",
                 ecolor="#333333", elinewidth=0.8, capsize=2)
    ax2.scatter(res.Difference, y, c=colours, edgecolors="#333333", s=22, zorder=3)
    ax2.axvline(0, color="#333333", linestyle="--", linewidth=0.7)
    ax2.set_xlim(-0.1, 0.1)
    ax2.set_yticks([])
    ax2.set_xlabel("Difference between groups")
    ax2.set_title("95% confidence intervals", fontsize=7)
    for yy, lab in zip(y, res.p_display):
        ax2.text(0.105, yy, lab, va="center", fontsize=7)
    for ax in (ax1, ax2):
        ax.set_ylim(-0.5, n - 0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.text(0.9, 0.5, "p_value", rotation=270, va="center", fontsize=7)
    fig.savefig(FIG / "Figure_5B_reproduced.png", dpi=300)
    plt.close(fig)
    print(f"Figure 5B: {n} genera written")


if __name__ == "__main__":
    main()
