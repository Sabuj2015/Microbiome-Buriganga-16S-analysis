# Microbiome-Buriganga-16S-analysis

Data and code for the manuscript

> **Sediment and water microbiomes of the Buriganga and Turag Rivers along an anthropogenic pollution gradient in Dhaka, Bangladesh** (BMC Microbiology, under review)

Raw reads: NCBI Sequence Read Archive, BioProject [PRJNA1399629](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1399629). The run accession of each sample is listed in `data/sra_runs.csv`.

## Study design

- 14 sites along the Turag (upstream, 3 sites) and Buriganga (midstream, 6 sites; downstream, 5 sites) Rivers, sampled in January 2023.
- At each site, three water and three sediment sub-samples were collected within a 5–8 m radius and pooled into one composite water sample and one composite sediment sample (28 samples).
- The sites are the replicates in all statistical comparisons: 14 per habitat, and 3, 6 and 5 per zone.

## Figures and tables

| Figure or table | Script | Input (`data/`) |
| --- | --- | --- |
| Figure 2A (tags and OTUs per sample) | `01_figure2.py` | `sequencing_stats.csv` |
| Figure 2B–E (Venn diagrams, rarefaction curves) | `01_figure2.py` | normalized OTU table * |
| Figure 3A (alpha diversity) | `03_figure3.py` | `alpha_diversity.csv` |
| Figure 3B–D (LEfSe cladograms; LDA score > 4) | `02_lefse.sh`, `02_lefse_input.py` | normalized OTU table * |
| Figure 4A (phylogenetic tree of dominant genera) | – | `panels/Figure_4A_phylogenetic_tree.png` |
| Figure 4B–C (NMDS, Bray–Curtis) | `04_figure4_nmds.R` | normalized OTU table * |
| Figure 5A (heatmap of 35 dominant genera) | `05_heatmaps.R` | `genus_group_means.csv` |
| Figure 5B (Welch's t-test, sediment vs water) | `05b_figure5B.py` | normalized OTU table *, `figure5B_genera.csv` |
| Figure 6A–B (PICRUSt2 v2.3.0, KEGG categories) | `05_heatmaps.R` | `picrust2_kegg_level_2.csv`, `picrust2_kegg_level_3.csv` |
| Tables 1–3; Additional file 1 (Tables S1–S3) | `06_additional_file_1.py` | as listed in the script |
| Additional file 2 (UpSet plot) | `07_additional_file_2_upset.py` | `otu_table.csv` |
| Figures 4–6 assembled | `08_compose_figures.py` | – |
| Sediment–water PERMANOVA with permutations restricted within sites; Kruskal–Wallis tests of potentially pathogenic genera among water zones | `09_replication_checks.R` | `otu_table.csv`, `sample_metadata.csv` |

\* These panels use the normalized (rarefied) OTU table of the original analysis, `data/otu_table_normalized.csv`, which has the same layout as `otu_table.csv`. If that file is absent, the scripts rarefy `otu_table.csv` to the smallest library (seed 2023), and the values are then close to, but not identical with, the published panels.

## Repository structure

```text
data/
  otu_table.csv                     17,664 OTUs x 28 samples with SILVA 138 taxonomy
  sample_metadata.csv               site number and name, river, zone, habitat, coordinates
  sra_runs.csv                      SRA run, experiment and BioSample accessions of each sample
  sequencing_stats.csv              per-sample read counts and quality statistics
  alpha_diversity.csv               per-sample alpha diversity indices (QIIME v1.9.1)
  alpha_group_means.csv             group means (Table 2)
  alpha_wilcoxon_sediment_vs_water.csv, alpha_kruskal_pairwise.csv
  adonis_bray_curtis.csv            PERMANOVA (Adonis, Bray-Curtis) results
  genus_*_means.csv, phylum_*_means.csv
  picrust2_kegg_level_1/2/3.csv     predicted KEGG categories per sample (PICRUSt2 v2.3.0)
  figure5B_genera.csv               genera shown in Figure 5B
  panels/                           Figure 4A
scripts/                            00-09 as listed above
results/                            regenerated figures, LEfSe output, Additional files 1 and 2,
                                    replication_checks.csv
run_all.sh                          runs scripts 01-09
requirements.txt
```

## Reproducing the figures

Requirements:

- Python 3 with the packages in `requirements.txt`;
- R with the `vegan`, `ggplot2` and `pheatmap` packages;
- optionally, LEfSe 1.1.2 (bioconda) for Figure 3B–D. Without LEfSe, the results already in `results/lefse/` are used.

```bash
pip install -r requirements.txt
Rscript -e 'install.packages(c("vegan", "ggplot2", "pheatmap"))'
./run_all.sh
```

With rpy2 ≥ 3.5, LEfSe 1.1.2 needs a two-line fix in `lefse/lefse.py`. Replace `scal = robjects.r('wfinal <- w.unit * effect.size')` with `robjects.r('wfinal <- w.unit * effect.size'); scal = robjects.r('wfinal')`, and make the same change for `rres = robjects.r('mm <- z$means')`.

`data/` can be recreated from the results workbook with `scripts/00_extract_results_workbook.py`.

## Repository history

The repository was reorganized in September 2026. Earlier files have been superseded and removed; use only the files listed above.

## Contact

Muhammad Shahdat Hossain (corresponding author), National Institute of Biotechnology, Savar, Dhaka, Bangladesh — <shahdatfbd@nib.gov.bd>

Repository maintained by Sabuj Biswas (GitHub: Sabuj2015).
