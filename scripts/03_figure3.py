#!/usr/bin/env python3
"""Figure 3: alpha diversity (A) and LEfSe cladograms (B-D).

A  Box plots of the per-sample alpha diversity indices (ACE, Chao1, Shannon, Simpson)
   in sediment and water (data/alpha_diversity.csv).
B-D Cladograms from 02_lefse.sh (LEfSe; LDA score > 4).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
RES = HERE.parent / "results"
FIG = RES / "figures"
FIG.mkdir(parents=True, exist_ok=True)
Image.MAX_IMAGE_PIXELS = None


def panel_a(path):
    alpha = pd.read_csv(DATA / "alpha_diversity.csv").set_index("Sample_ID")
    meta = pd.read_csv(DATA / "sample_metadata.csv").set_index("Sample_ID")
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    for ax, (col, lab) in zip(axes, [("ACE", "ACE"), ("chao1", "chao1"), ("shannon", "shannon"), ("simpson", "simpson")]):
        data = [alpha.loc[meta.index[meta.Habitat == h], col] for h in ("Sediment", "Water")]
        bp = ax.boxplot(data, widths=0.75, patch_artist=True, showfliers=True,
                        flierprops=dict(marker="o", markersize=3, markerfacecolor="black", markeredgecolor="black"),
                        medianprops=dict(color="#333333"))
        for patch, colour in zip(bp["boxes"], ("#dc143c", "#0000ff")):
            patch.set_facecolor(colour)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Sediment", "Water"], rotation=30)
        ax.set_xlabel("Group")
        ax.set_ylabel(lab)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main():
    panel_a(FIG / "Figure_3A_reproduced.png")
    top = Image.open(FIG / "Figure_3A_reproduced.png").convert("RGB")
    clado = [Image.open(RES / "lefse" / f"cladogram_{k}.png").convert("RGB") for k in "BCD"]
    h = max(c.height for c in clado)
    row = Image.new("RGB", (sum(c.width for c in clado), h), "white")
    x = 0
    for c in clado:
        row.paste(c, (x, 0))
        x += c.width
    row = row.resize((top.width, round(row.height * top.width / row.width)), Image.LANCZOS)
    canvas = Image.new("RGB", (top.width, top.height + row.height), "white")
    canvas.paste(top, (0, 0))
    canvas.paste(row, (0, top.height))
    canvas.save(FIG / "Figure_3_reproduced.tiff", dpi=(300, 300), compression="tiff_lzw")
    canvas.thumbnail((2400, 2400))
    canvas.save(FIG / "Figure_3_reproduced.png")
    print("Figure 3 written")


if __name__ == "__main__":
    main()
