#!/usr/bin/env python3
"""Additional file 2: UpSet plot of OTU membership across the six zone x habitat groups.

An OTU belongs to a group if it has at least one read in any sample of that group (OTU table,
17,664 OTUs, each OTU counted once). The 30 largest intersections are shown.
"""

from itertools import compress
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIG = HERE.parent / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SETS = [("Up-Sed", "Upstream", "Sediment", "#8b4513"), ("Mid-Sed", "Midstream", "Sediment", "#cd853f"),
        ("Down-Sed", "Downstream", "Sediment", "#deb887"), ("Up-Wat", "Upstream", "Water", "#1f4e89"),
        ("Mid-Wat", "Midstream", "Water", "#3a87c8"), ("Down-Wat", "Downstream", "Water", "#8ec9ea")]
KIND = {"sed": ("Sediment-only", "#8b4513"), "wat": ("Water-only", "#3a87c8"), "both": ("Shared (sediment + water)", "#6a5acd")}


def main():
    otu = pd.read_csv(DATA / "otu_table.csv")
    meta = pd.read_csv(DATA / "sample_metadata.csv").set_index("Sample_ID")
    member = pd.DataFrame({name: (otu[list(meta.index[(meta.Zone == z) & (meta.Habitat == h)])] > 0).any(axis=1)
                           for name, z, h, _ in SETS})
    combos = member.value_counts().rename("OTUs").reset_index()
    combos = combos[combos[[s[0] for s in SETS]].any(axis=1)].head(30).reset_index(drop=True)
    combos.to_csv(FIG / "Additional_file_2_upset_counts.csv", index=False)
    names = [s[0] for s in SETS]
    fig = plt.figure(figsize=(16, 8))
    ax_bar = fig.add_axes([0.13, 0.42, 0.85, 0.45])
    ax_dot = fig.add_axes([0.13, 0.05, 0.85, 0.33], sharex=ax_bar)
    x = np.arange(len(combos))
    colours, kinds = [], []
    for _, r in combos.iterrows():
        sed_any = any(r[n] for n in names[:3])
        wat_any = any(r[n] for n in names[3:])
        k = "both" if sed_any and wat_any else ("sed" if sed_any else "wat")
        kinds.append(k)
        colours.append(KIND[k][1])
    ax_bar.bar(x, combos.OTUs, color=colours, width=0.65)
    for xi, v in zip(x, combos.OTUs):
        ax_bar.text(xi, v + combos.OTUs.max() * 0.01, f"{v:,}", ha="center", va="bottom", fontsize=6, weight="bold")
    ax_bar.set_ylabel("OTU count")
    ax_bar.grid(axis="y", linestyle="--", alpha=0.5)
    for side in ("top", "right"):
        ax_bar.spines[side].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in KIND.values()]
    ax_bar.legend(handles, [l for l, _ in KIND.values()], loc="upper right")
    ax_bar.set_title("UpSet Plot: OTU membership across zone × habitat combinations\n"
                     "Buriganga–Turag River System  |  January 2023  |  top 30 intersections", weight="bold")
    plt.setp(ax_bar.get_xticklabels(), visible=False)
    for i, (name, _, _, col) in enumerate(SETS):
        yy = len(SETS) - 1 - i
        if i % 2 == 0:
            ax_dot.axhspan(yy - 0.5, yy + 0.5, color="#f2f2f2", zorder=0)
        for xi, r in combos.iterrows():
            ax_dot.scatter(xi, yy, s=60, color=col if r[name] else "#dddddd", zorder=3)
        tr = ax_dot.get_yaxis_transform()
        ax_dot.text(-0.075, yy, name, transform=tr, ha="right", va="center", color=col, weight="bold")
        ax_dot.text(-0.012, yy, f"{int(member[name].sum()):,}", transform=tr, ha="right", va="center", color=col, weight="bold")
    for xi, r in combos.iterrows():
        ys = [len(SETS) - 1 - i for i, n in enumerate(names) if r[n]]
        if len(ys) > 1:
            ax_dot.plot([xi, xi], [min(ys), max(ys)], color="#555555", linewidth=1.5, zorder=2)
    ax_dot.set_ylim(-0.6, len(SETS) - 0.4)
    ax_dot.set_xlim(-0.6, len(combos) - 0.4)
    ax_dot.axis("off")
    fig.text(0.12, 0.385, "Set\nsize", ha="center", fontsize=8, style="italic")
    fig.savefig(FIG / "Additional_file_2_UpSet.png", dpi=300)
    # TIFF without an alpha channel (RGB, LZW compression, 300 dpi) for journal submission systems
    Image.open(FIG / "Additional_file_2_UpSet.png").convert("RGB").save(
        FIG / "Additional_file_2_UpSet.tiff", dpi=(300, 300), compression="tiff_lzw")
    core = int(member.all(axis=1).sum())
    print(f"UpSet: {len(member[member.any(axis=1)])} OTUs; shared by all six groups: {core}; set sizes:",
          {n: int(member[n].sum()) for n in names})


if __name__ == "__main__":
    main()
