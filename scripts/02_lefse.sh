#!/usr/bin/env bash
# LEfSe for Fig. 3B-D (LDA score > 4; alpha = 0.05 for the Kruskal-Wallis and Wilcoxon tests).
# Requires LEfSe 1.1.2 (bioconda: micromamba create -n lefse -c conda-forge -c bioconda lefse).
# With rpy2 >= 3.5, LEfSe needs a two-line fix in lefse/lefse.py (see README).
set -euo pipefail
cd "$(dirname "$0")"
python3 02_lefse_input.py
cd ../results/lefse
for x in B:lefse_B_sediment_vs_water C:lefse_C_sediment_zones D:lefse_D_water_upstream_vs_downstream; do
  k=${x%%:*}; f=${x#*:}
  lefse_format_input.py "$f.txt" "$k.in" -c 1 -s -1 -u 2 -o 1000000
  lefse_run.py "$k.in" "$k.res" -l 4
  lefse_plot_cladogram.py "$k.res" "cladogram_$k.png" --format png --dpi 300 --title "Cladogram" \
      --class_legend_font_size 9 --labeled_stop_lev 7 --abrv_stop_lev 7
  awk -F'\t' 'BEGIN{OFS="\t"; print "Taxon","Group","LDA_log10","p_KW"} $3!=""{print $1,$3,$4,$5}' "$k.res" \
      | sort -t$'\t' -k2,2 -k3,3gr > "lefse_significant_$k.tsv"
done
