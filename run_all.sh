#!/usr/bin/env bash
# Regenerate the figures and supplementary tables of the manuscript from the files in data/.
# Requirements: Python 3 (requirements.txt), R with vegan, ggplot2 and pheatmap, and (for Fig. 3B-D) LEfSe.
set -euo pipefail
cd "$(dirname "$0")/scripts"

python3 01_figure2.py                 # Figure 2
if command -v lefse_run.py >/dev/null 2>&1; then
  ./02_lefse.sh                       # Figure 3B-D (LEfSe)
else
  echo "LEfSe not found: using the LEfSe results already in results/lefse (see README)"
fi
python3 03_figure3.py                 # Figure 3
Rscript 04_figure4_nmds.R             # Figure 4B-C
Rscript 05_heatmaps.R                 # Figures 5A and 6
python3 05b_figure5B.py               # Figure 5B
python3 06_additional_file_1.py       # Additional file 1 (Tables S1-S3, Tables 1-3)
python3 07_additional_file_2_upset.py # Additional file 2 (UpSet plot)
python3 08_compose_figures.py         # Figures 4-6 assembled
echo "Done. Outputs are in results/."
