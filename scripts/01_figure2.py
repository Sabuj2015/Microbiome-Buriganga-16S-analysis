#!/usr/bin/env python3
"""Figure 2: sequencing output (A), shared OTUs (B-D) and rarefaction curves (E).

A  Tags per sample and OTUs per sample (sequencing statistics).
B-D Venn diagrams of OTUs among zones (sediment, water) and between habitats.
E  Rarefaction curves (expected number of OTUs) of the six zone x habitat groups (mean +/- SD).

B-E use the normalized (rarefied) OTU table of the original analysis,
data/otu_table_normalized.csv (same layout as otu_table.csv). If that file is absent, the
OTU table is rarefied here to the smallest library (seed 2023), which gives values close to,
but not identical with, the published panels.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle
from scipy.special import gammaln

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIG = HERE.parent / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})

ORDER = ["DeulS", "NbbS", "AbbS", "BsbS", "MdcS", "KrcS", "BsgS", "JjfS", "BbbS", "CmbS", "PgbS", "HnbS", "PaglaS",
         "FatulS", "DeulW", "NbbW", "AbbW", "BsbW", "MdcW", "KrcW", "BsgW", "JjfW", "BbbW", "CmbW", "PgbW", "HnbW",
         "PaglaW", "FatulW"]
GROUP_LABEL = {("Upstream", "Sediment"): "UpstrmS", ("Midstream", "Sediment"): "MidstrmS",
               ("Downstream", "Sediment"): "DownstmS", ("Upstream", "Water"): "UpstrmW",
               ("Midstream", "Water"): "MidstrmW", ("Downstream", "Water"): "DownstmW"}


def panel_a(ax, seq):
    s = seq.set_index("Sample_ID").loc[ORDER]
    bars = [("Total_tags", "Total Tags", "#d7263d"), ("Taxon_tags", "Taxon Tags", "#1f4e9c"),
            ("Unclassified_tags", "Unclassified Tags", "#20b2aa"), ("Unique_tags", "Unique Tags", "#f39c12")]
    x = np.arange(len(ORDER))
    w = 0.17
    for k, (col, lab, colour) in enumerate(bars):
        v = s[col].to_numpy()
        ax.bar(x + (k - 2) * w, v, w, color=colour, label=f"{lab}(avg:{v.mean():.0f})")
        for xi, vi in zip(x + (k - 2) * w, v):
            ax.text(xi, vi + 500, f"{vi:.0f}", rotation=90, ha="center", va="bottom", fontsize=4.5)
    ax2 = ax.twinx()
    v = s["OTUs"].to_numpy()
    ax2.bar(x + 2 * w, v, w, color="#8e7cc3", label=f"OTUs(avg:{v.mean():.0f})")
    for xi, vi in zip(x + 2 * w, v):
        ax2.text(xi, vi + 30, f"{vi:.0f}", rotation=90, ha="center", va="bottom", fontsize=4.5)
    ax.set_ylim(0, 80000)
    ax2.set_ylim(0, 5000)
    ax.set_ylabel("Tags Number")
    ax2.set_ylabel("OTUs Number")
    ax.set_xticks(x)
    ax.set_xticklabels(ORDER, rotation=45, ha="right")
    ax.set_xlabel("Sample Name")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, ncol=3, frameon=False, loc="upper left", bbox_to_anchor=(0, 1.18), fontsize=7)


def rarefy(counts, depth, seed=2023):
    rng = np.random.default_rng(seed)
    out = np.zeros_like(counts)
    for j in range(counts.shape[1]):
        out[:, j] = rng.multivariate_hypergeometric(counts[:, j], depth, method="marginals")
    return out


def expected_richness(c, d):
    c = c[c > 0].astype(float)
    n = c.sum()
    if d >= n:
        return float(len(c))
    with np.errstate(invalid="ignore"):
        la = np.where(n - c >= d, gammaln(n - c + 1) - gammaln(d + 1) - gammaln(n - c - d + 1)
                      - (gammaln(n + 1) - gammaln(d + 1) - gammaln(n - d + 1)), -np.inf)
    return float((1 - np.exp(la)).sum())


def venn3(ax, sets, labels, colours):
    a, b, c = sets
    regions = {"100": len(a - b - c), "010": len(b - a - c), "001": len(c - a - b), "110": len(a & b - c),
               "101": len(a & c - b), "011": len(b & c - a), "111": len(a & b & c)}
    centres = [(-0.45, 0.35), (0.45, 0.35), (0.0, -0.4)]
    for (x, y), col in zip(centres, colours):
        ax.add_patch(Circle((x, y), 0.75, facecolor=col, alpha=0.55, edgecolor="white", linewidth=1))
    pos = {"100": (-0.75, 0.55), "010": (0.75, 0.55), "001": (0, -0.85), "110": (0, 0.75), "101": (-0.5, -0.25),
           "011": (0.5, -0.25), "111": (0, 0.1)}
    for k, (x, y) in pos.items():
        ax.text(x, y, f"{regions[k]}", ha="center", va="center", fontsize=7, family="serif")
    for (x, y), lab, col, dy in zip(centres, labels, colours, (0.95, 0.95, -0.95)):
        ax.text(x * 1.6, y + dy, lab, ha="center", va="center", fontsize=7, color=col, weight="bold", family="serif")
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    return regions


def venn2(ax, a, b, labels, colours):
    regions = {"10": len(a - b), "01": len(b - a), "11": len(a & b)}
    for (x, col) in zip((-0.4, 0.4), colours):
        ax.add_patch(Circle((x, 0), 0.75, facecolor=col, alpha=0.55, edgecolor="white", linewidth=1))
    for k, x in (("10", -0.75), ("11", 0), ("01", 0.75)):
        ax.text(x, 0, f"{regions[k]}", ha="center", va="center", fontsize=7, family="serif")
    for x, lab, col in zip((-1.05, 1.05), labels, colours):
        ax.text(x, 0.7, lab, ha="center", fontsize=7, color=col, weight="bold", family="serif")
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    return regions


def main():
    seq = pd.read_csv(DATA / "sequencing_stats.csv")
    meta = pd.read_csv(DATA / "sample_metadata.csv").set_index("Sample_ID")
    norm_file = DATA / "normalized_otu_table.csv"
    otu = pd.read_csv(norm_file if norm_file.exists() else DATA / "otu_table.csv")
    samples = list(meta.index)
    counts = otu[samples].to_numpy()
    table = counts if norm_file.exists() else rarefy(counts, int(counts.sum(axis=0).min()))
    present = pd.DataFrame(table > 0, columns=samples, index=otu.OTU_ID)

    def otus(zone=None, habitat=None):
        ids = [s for s in samples if (zone is None or meta.loc[s, "Zone"] == zone)
               and (habitat is None or meta.loc[s, "Habitat"] == habitat)]
        return set(present.index[present[ids].any(axis=1)])

    fig = plt.figure(figsize=(12, 7))
    ax_a = fig.add_axes([0.05, 0.56, 0.9, 0.34])
    panel_a(ax_a, seq)
    colours = ["#8faadc", "#70ad47", "#d48bab"]
    ax_b = fig.add_axes([0.02, 0.02, 0.22, 0.38])
    rb = venn3(ax_b, [otus(z, "Sediment") for z in ("Upstream", "Midstream", "Downstream")],
               ["UpstrmS", "MidstrmS", "DownstmS"], colours)
    ax_c = fig.add_axes([0.25, 0.02, 0.22, 0.38])
    rc = venn2(ax_c, otus(habitat="Sediment"), otus(habitat="Water"), ["Soil", "Water"], colours[:2])
    ax_d = fig.add_axes([0.48, 0.02, 0.22, 0.38])
    rd = venn3(ax_d, [otus(z, "Water") for z in ("Upstream", "Midstream", "Downstream")],
               ["UpstrmW", "MidstrmW", "DownstmW"], colours)

    # E: rarefaction curves (expected OTUs), group means +/- SD
    ax_e = fig.add_axes([0.75, 0.07, 0.23, 0.33])
    depth = int(counts.sum(axis=0).min())
    depths = np.linspace(0, depth, 7).round().astype(int)
    styles = {"UpstrmS": ("#c0392b", "o"), "MidstrmS": ("#1f3a93", "^"), "DownstmS": ("#16a085", "+"),
              "UpstrmW": ("#e67e22", "x"), "MidstrmW": ("#8e44ad", "D"), "DownstmW": ("#27ae60", "v")}
    for (zone, hab), lab in GROUP_LABEL.items():
        ids = [s for s in samples if meta.loc[s, "Zone"] == zone and meta.loc[s, "Habitat"] == hab]
        curves = np.array([[expected_richness(counts[:, samples.index(s)], d) for d in depths] for s in ids])
        col, mk = styles[lab]
        ax_e.errorbar(depths, curves.mean(axis=0), yerr=curves.std(axis=0, ddof=1), color=col, marker=mk,
                      markersize=3, linewidth=1, capsize=2, label=lab)
    ax_e.set_xlabel("Sequences Number")
    ax_e.set_ylabel("Observed Species Number")
    ax_e.legend(fontsize=5, frameon=True)
    for ax, lab, x, y in ((ax_a, "A", 0.01, 0.93), (ax_b, "B", 0.02, 0.38), (ax_c, "C", 0.26, 0.38),
                          (ax_d, "D", 0.49, 0.38), (ax_e, "E", 0.72, 0.38)):
        fig.text(x, y, lab, fontsize=16)
    fig.savefig(FIG / "Figure_2_reproduced.png", dpi=300)
    fig.savefig(FIG / "Figure_2_reproduced.tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    pd.DataFrame([dict(Panel="B (sediment zones)", **rb), dict(Panel="D (water zones)", **rd)]).to_csv(
        FIG / "Figure_2BD_venn_counts.csv", index=False)
    pd.DataFrame([dict(Panel="C (sediment vs water)", **rc)]).to_csv(FIG / "Figure_2C_venn_counts.csv", index=False)
    print("Figure 2: shared sediment-water OTUs =", rc["11"],
          "(normalized OTU table)" if norm_file.exists() else f"(rarefied to {depth} reads, seed 2023)")


if __name__ == "__main__":
    main()
