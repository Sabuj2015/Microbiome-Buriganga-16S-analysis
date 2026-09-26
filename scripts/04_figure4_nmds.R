#!/usr/bin/env Rscript
# Figure 4B-C: NMDS of Bray-Curtis dissimilarities (vegan::metaMDS with default settings,
# i.e. square-root and Wisconsin double standardization) of the normalized OTU table
# (data/otu_table_normalized.csv). If that file is absent, the OTU table is rarefied here to
# the smallest library, which gives a configuration close to, but not identical with, the
# published panels. Panel A (phylogenetic tree of dominant genera, from the OTU representative
# sequences) is included as data/panels/Figure_4A_phylogenetic_tree.png.
suppressPackageStartupMessages({library(vegan); library(ggplot2)})
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
data <- file.path(here, "..", "data"); out <- file.path(here, "..", "results", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

otu  <- read.csv(file.path(data, "otu_table.csv"), check.names = FALSE)
meta <- read.csv(file.path(data, "sample_metadata.csv"))
X <- t(as.matrix(otu[, meta$Sample_ID])); rownames(X) <- meta$Sample_ID
norm_file <- file.path(data, "otu_table_normalized.csv")
if (file.exists(norm_file)) {
  nt <- read.csv(norm_file, check.names = FALSE)
  R <- t(as.matrix(nt[, meta$Sample_ID])); rownames(R) <- meta$Sample_ID
} else {
  set.seed(2023)
  R <- rrarefy(X, min(rowSums(X)))
}
R <- R[, colSums(R) > 0]                      # OTUs lost by rarefaction
set.seed(1)
nm <- metaMDS(R, distance = "bray", k = 2, trace = 0)
sc <- as.data.frame(scores(nm, display = "sites")); sc$Sample_ID <- rownames(sc)
sc <- merge(sc, meta, by = "Sample_ID")
lab <- c(Upstream = "Upstrm", Midstream = "Midstrm", Downstream = "Downstm")
sc$Group <- factor(paste0(lab[sc$Zone], ifelse(sc$Habitat == "Sediment", "S", "W")),
                   levels = c("UpstrmS", "MidstrmS", "DownstmS", "UpstrmW", "MidstrmW", "DownstmW"))
write.csv(sc, file.path(out, "Figure_4BC_nmds_coordinates.csv"), row.names = FALSE)
stress <- sprintf("Stress = %.3f", nm$stress)
theme_nmds <- theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"),
                                 legend.title = element_blank())
pB <- ggplot(sc, aes(NMDS1, NMDS2, colour = Habitat, shape = Habitat)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  geom_point(size = 3) + stat_ellipse(level = 0.95, show.legend = FALSE) +
  scale_colour_manual(values = c(Sediment = "#d62d3a", Water = "#232a7d")) +
  scale_shape_manual(values = c(Sediment = 15, Water = 16)) +
  annotate("text", x = Inf, y = Inf, label = stress, hjust = 1.1, vjust = 1.5) +
  labs(title = "NMDS Plot", x = "MDS1", y = "MDS2") + theme_nmds
pC <- ggplot(sc, aes(NMDS1, NMDS2, colour = Group, shape = Group)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(UpstrmS = "#d62d3a", MidstrmS = "#232a7d", DownstmS = "#1b86a8",
                                 UpstrmW = "#f39c12", MidstrmW = "#8e6fa6", DownstmW = "#7cc59b")) +
  scale_shape_manual(values = c(UpstrmS = 15, MidstrmS = 16, DownstmS = 17, UpstrmW = 18, MidstrmW = 20, DownstmW = 7)) +
  annotate("text", x = Inf, y = Inf, label = stress, hjust = 1.1, vjust = 1.5) +
  labs(title = "NMDS Plot", x = "MDS1", y = "MDS2") + theme_nmds
ggsave(file.path(out, "Figure_4B_reproduced.png"), pB, width = 6, height = 5, dpi = 300)
ggsave(file.path(out, "Figure_4C_reproduced.png"), pC, width = 6, height = 5, dpi = 300)
cat("NMDS", stress, "\n")
