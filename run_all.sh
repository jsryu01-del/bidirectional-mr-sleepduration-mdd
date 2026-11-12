#!/usr/bin/env bash
# Reproduce all main and negative-control MR results, Supplementary Tables 8-12 values and the figures.
# Usage: bash run_all.sh [n_cores]
#   Environment variables RSCRIPT and PYTHON may point to specific interpreters (default: Rscript, python3).
#   Runtime: MR-PRESSO (10,000 simulations) takes ~30 min for the 165-170-SNP analyses; with >=12 cores all
#   analyses run in parallel (~35 min total). Each analysis uses one core.
set -euo pipefail
cd "$(dirname "$0")"
NC=${1:-4}; RSCRIPT=${RSCRIPT:-Rscript}; PYTHON=${PYTHON:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results/logs
"$RSCRIPT" --vanilla code/R/00_check_environment.R
tail -n +2 data/manifest.tsv | cut -f1 | \
  xargs -P "$NC" -I{} sh -c "\"$RSCRIPT\" --vanilla code/R/01_run_mr_analysis.R {} > results/logs/{}.log 2>&1 || echo 'FAILED: {} (see results/logs/{}.log)'"
for a in $(tail -n +2 data/manifest.tsv | cut -f1); do
  [ -f "results/runs/published_seed/$a/status.txt" ] || { echo "Analysis $a did not finish"; exit 1; }
done
"$RSCRIPT" --vanilla code/R/02_collect_results.R
"$RSCRIPT" --vanilla code/R/03_compare_with_published.R
"$PYTHON" code/python/figure2_mr_forest.py
"$PYTHON" code/python/figure_main_diagnostics.py
echo "Done. Tables: results/tables/  Comparison: results/comparison_summary.tsv  Figures: figures/"
