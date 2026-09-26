#!/usr/bin/env python3
"""Write LEfSe input tables (Fig. 3B-D) from the normalized OTU table (or the OTU table if absent).

Each feature is a taxon at one rank, written as its full lineage (k__...|p__...|...), with its
relative abundance in each sample (reads of the taxon / reads in the OTU table). Unclassified
ranks end the lineage. Three comparisons are written, as in the manuscript:
  B  sediment vs water (28 samples)
  C  sediment samples by zone (14 samples)
  D  water samples, upstream vs downstream (8 samples)
"""

import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = HERE.parent / "results" / "lefse"
OUT.mkdir(parents=True, exist_ok=True)
RANKS = [("k", "Kingdom"), ("p", "Phylum"), ("c", "Class"), ("o", "Order"), ("f", "Family"), ("g", "Genus"), ("s", "Species")]


def lineage_table(otu, samples):
    rel = otu[samples] / otu[samples].sum()
    rows = {}
    for depth in range(1, len(RANKS) + 1):
        names = []
        for _, r in otu[[c for _, c in RANKS[:depth]]].iterrows():
            parts = []
            for (prefix, col) in RANKS[:depth]:
                v = r[col]
                if not isinstance(v, str) or not v or v == "Unassigned":
                    parts = None
                    break
                parts.append(f"{prefix}__" + re.sub(r"[^A-Za-z0-9_.]", "_", v))   # LEfSe/R-safe names
            names.append("|".join(parts) if parts else None)
        key = pd.Series(names, index=otu.index)
        grouped = rel[key.notna()].groupby(key[key.notna()]).sum()
        rows.update({k: v for k, v in grouped.iterrows()})
    return pd.DataFrame(rows).T[samples]


def write(table, meta, samples, label_col, name):
    header = pd.DataFrame([meta.loc[samples, label_col].tolist(), samples], index=["class", "subject"], columns=samples)
    out = pd.concat([header, table[samples].map(lambda x: f"{x:.8f}")])
    out.to_csv(OUT / f"{name}.txt", sep="\t", header=False)


def main():
    norm_file = DATA / "otu_table_normalized.csv"
    otu = pd.read_csv(norm_file if norm_file.exists() else DATA / "otu_table.csv")
    meta = pd.read_csv(DATA / "sample_metadata.csv").set_index("Sample_ID")
    meta["Group"] = meta.Zone.str.slice(0, 4) + "_" + meta.Habitat
    table = lineage_table(otu, list(meta.index))
    sed = list(meta.index[meta.Habitat == "Sediment"])
    wat_ud = list(meta.index[(meta.Habitat == "Water") & meta.Zone.isin(["Upstream", "Downstream"])])
    write(table, meta, list(meta.index), "Habitat", "lefse_B_sediment_vs_water")
    write(table, meta, sed, "Zone", "lefse_C_sediment_zones")
    write(table, meta, wat_ud, "Zone", "lefse_D_water_upstream_vs_downstream")
    print("LEfSe input written:", len(table), "features")


if __name__ == "__main__":
    main()
