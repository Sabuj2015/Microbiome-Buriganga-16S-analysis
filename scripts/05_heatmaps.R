#!/usr/bin/env Rscript
# Figure 5A and Figure 6A-B: heatmaps drawn with pheatmap.
#  5A: group means of the 35 genera with the highest mean relative abundance in any of the
#      six zone x habitat groups; rows scaled (z-scores), rows clustered (Euclidean, complete).
#  6A/6B: habitat means of the KEGG level-2 / level-3 categories (PICRUSt2 v2.3.0);
#      rows scaled across the two habitat means, so every z-score is +/-0.707 (1/sqrt(2)):
#      the colours show which habitat is higher, not the size of the difference.
suppressPackageStartupMessages(library(pheatmap))
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
data <- file.path(here, "..", "data"); out <- file.path(here, "..", "results", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Figure 5A
g <- read.csv(file.path(data, "genus_group_means.csv"), check.names = FALSE, row.names = 1)
g <- g[rownames(g) != "Others", c("UpstrmS", "MidstrmS", "DownstmS", "UpstrmW", "MidstrmW", "DownstmW")]
top <- head(rownames(g)[order(apply(g, 1, max), decreasing = TRUE)], 35)
otu <- read.csv(file.path(data, "otu_table.csv"), check.names = FALSE)
phylum <- tapply(otu$Phylum, otu$Genus, function(x) names(sort(table(x), decreasing = TRUE))[1])
ann <- data.frame(Phylum = phylum[top], row.names = top)
pheatmap(as.matrix(g[top, ]), scale = "row", cluster_cols = FALSE, clustering_distance_rows = "euclidean",
         clustering_method = "complete", annotation_row = ann, cellwidth = 40, cellheight = 14,
         filename = file.path(out, "Figure_5A_reproduced.png"), width = 9, height = 9)

# ---- Figure 6
meta <- read.csv(file.path(data, "sample_metadata.csv"))
for (lev in c(2, 3)) {
  k <- read.csv(file.path(data, sprintf("picrust2_kegg_level_%d.csv", lev)), check.names = FALSE, row.names = 1)
  m <- cbind(Sediment = rowMeans(k[, meta$Sample_ID[meta$Habitat == "Sediment"]]),
             Water = rowMeans(k[, meta$Sample_ID[meta$Habitat == "Water"]]))
  z <- t(scale(t(m)))
  cat(sprintf("Figure 6%s: z-scores range %.4f to %.4f\n", c("", "", "A", "B")[lev + 1], min(z), max(z)))
  pheatmap(m, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, cellwidth = 60, cellheight = 14,
           filename = file.path(out, sprintf("Figure_6%s_reproduced.png", c("", "", "A", "B")[lev + 1])),
           width = 6, height = 9)
}
