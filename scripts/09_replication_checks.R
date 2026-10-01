#!/usr/bin/env Rscript
# Checks reported in the manuscript and the response to reviewers (second revision):
#  1. PERMANOVA of sediment vs water (Bray-Curtis, relative abundances), with permutations
#     either unrestricted or restricted within sites. Each site provided one sediment and one
#     water sample, so the restricted test compares the two habitats within each site.
#  2. Kruskal-Wallis tests of the relative abundance of the potentially pathogenic genera in
#     water among the three zones (the five genera of Table 3, summed and individually).
# Relative abundance = reads of an OTU or genus / reads of the sample in the OTU table (as in Table 3).
suppressPackageStartupMessages(library(vegan))
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
data <- file.path(here, "..", "data"); out <- file.path(here, "..", "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

otu  <- read.csv(file.path(data, "otu_table.csv"), check.names = FALSE)
meta <- read.csv(file.path(data, "sample_metadata.csv"))
X <- t(as.matrix(otu[, meta$Sample_ID])); rownames(X) <- meta$Sample_ID
d <- vegdist(X / rowSums(X), method = "bray")

set.seed(2023)
free <- adonis2(d ~ Habitat, data = meta, permutations = 9999)
set.seed(2023)
within_site <- adonis2(d ~ Habitat, data = meta,
                       permutations = how(nperm = 9999, blocks = factor(meta$Site_number)))

genus <- ifelse(is.na(otu$Genus) | otu$Genus == "", "Unassigned", otu$Genus)
rel <- rowsum(otu[, meta$Sample_ID], genus)
rel <- sweep(rel, 2, colSums(otu[, meta$Sample_ID]), "/") * 100
pathogens <- c("Ralstonia", "Arcobacter", "Acinetobacter", "Cloacibacterium", "Aeromonas")
w <- meta$Habitat == "Water"
zone <- factor(meta$Zone[w], levels = c("Upstream", "Midstream", "Downstream"))
kw <- function(v) kruskal.test(v ~ zone)
summed <- colSums(rel[pathogens, meta$Sample_ID[w]])

rows <- list(
  data.frame(test = "PERMANOVA sediment vs water, permutations unrestricted (9999)", n = 28,
             statistic = sprintf("R2 = %.3f, F = %.2f", free$R2[1], free$F[1]), p = free$`Pr(>F)`[1]),
  data.frame(test = "PERMANOVA sediment vs water, permutations restricted within sites (9999)", n = 28,
             statistic = sprintf("R2 = %.3f, F = %.2f", within_site$R2[1], within_site$F[1]),
             p = within_site$`Pr(>F)`[1]),
  data.frame(test = sprintf("Kruskal-Wallis, water, zones: sum of %s (zone means %s %%)",
                            paste(pathogens, collapse = "+"),
                            paste(sprintf("%s %.2f", levels(zone), tapply(summed, zone, mean)), collapse = "; ")),
             n = sum(w), statistic = sprintf("chi2 = %.2f, df = 2", kw(summed)$statistic), p = kw(summed)$p.value))
for (g in pathogens) {
  k <- kw(unlist(rel[g, meta$Sample_ID[w]]))
  rows[[length(rows) + 1]] <- data.frame(test = paste("Kruskal-Wallis, water, zones:", g), n = sum(w),
                                         statistic = sprintf("chi2 = %.2f, df = 2", k$statistic), p = k$p.value)
}
res <- do.call(rbind, rows)
res$p <- signif(res$p, 3)
write.csv(res, file.path(out, "replication_checks.csv"), row.names = FALSE)
print(res, right = FALSE, row.names = FALSE)
