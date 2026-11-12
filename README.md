# Reproducibility package: sleep duration and major depressive disorder — main, negative-control and leave-many-out MR

This package contains the final harmonized datasets used for every main, negative-control and leave-many-out
Mendelian randomization (MR) analysis and code that reproduces, from those datasets alone:

- the instrument-strength values of Supplementary Table 8 (main and negative-control rows),
- all estimates and diagnostics of Supplementary Tables 9 (main), 11 (leave-many-out) and 12 (negative controls)
  (IVW, MR-Egger, weighted median, MR-PRESSO raw/outlier-corrected, MR-RAPS for main and leave-many-out analyses;
  Cochran's Q; MR-Egger intercept; BH-adjusted P values; highlighted method),
- the leave-many-out exclusions of Supplementary Table 10, from the complete predictor look-up table,
- Figure 2 (MR forest plot) and the four main-analysis diagnostic figures (scatter, funnel, single-SNP forest,
  leave-one-out) for short→MDD, long→MDD, MDD→short and MDD→long,

and verifies them automatically against the printed supplementary tables (ST8, ST9 and ST11 at full precision;
ST12, which is printed rounded, at its printed precision).
Table numbers refer to the manuscript's Supplementary Tables 6–12 (ST6–ST12). The comparison values extracted from
those tables are included in `data/published_values/published_values.tsv`; the workbook itself is
not required for this workflow. Reproduction starts from the supplied harmonized SNP-level data.
Raw-GWAS extraction, LD clumping, proxy selection and harmonization are not rerun here.

## Supplementary table scope

This package uses Supplementary Tables 6–12 (ST6–ST12): ST6–ST7 document the
SNP-level instruments supplied as harmonized inputs; ST8 reports instrument strength;
ST9 reports primary MR results; ST10 lists leave-many-out exclusions; ST11 reports
leave-many-out MR results; and ST12 reports negative-control MR results.
Automated numerical comparisons cover ST8–ST12; ST6–ST7 provide the input-data context.

## Included files

| Path | Purpose |
|---|---|
| `run_all.sh` | Run the complete workflow locally using Bash |
| `code/R/00_check_environment.R` | Check packages and record versions |
| `code/R/01_run_mr_analysis.R` | Run one MR analysis |
| `code/R/02_collect_results.R` | Collect estimates, diagnostics, instrument strength and figure inputs |
| `code/R/03_compare_with_published.R` | Compare results with the included reference values |
| `code/python/figure2_mr_forest.py` | Generate the main MR forest plot |
| `code/python/figure_main_diagnostics.py` | Generate primary-analysis diagnostic figures |
| `data/manifest.tsv` | Analysis definitions, input paths, instrument counts and seeds |
| `data/harmonized/` | Twelve SNP-level input datasets |
| `data/leave_many_out/` | Predictor panels and complete lookup results |
| `data/published_values/` | Reference values for automated comparisons |
| `environment/install_R_packages.R` | Install R dependencies |
| `environment/requirements_python.txt` | Python dependency versions |
| `environment/package_versions_original_analysis.tsv` | R/package versions used in the original analysis |

## Analyses

| analysis | type | exposure → outcome | instruments |
|---|---|---|---|
| short_to_MDD | main | short sleep → MDD | 26 / 26 |
| long_to_MDD | main | long sleep → MDD | 7 / 7 |
| MDD_to_short | main | MDD → short sleep | 170 / 184 |
| MDD_to_long | main | MDD → long sleep | 168 / 184 |
| short_to_lefthand | negative control | short sleep → left-handedness | 26 / 26 |
| short_to_blond | negative control | short sleep → blond hair | 25 / 26 |
| long_to_lefthand | negative control | long sleep → left-handedness | 7 / 7 |
| long_to_blond | negative control | long sleep → blond hair | 6 / 7 |
| MDD_to_lefthand_set_short | negative control | MDD → left-handedness (MDD→short IV set) | 170 / 170 |
| MDD_to_blond_set_short | negative control | MDD → blond hair (MDD→short IV set) | 167 / 170 |
| MDD_to_lefthand_set_long | negative control | MDD → left-handedness (MDD→long IV set) | 168 / 168 |
| MDD_to_blond_set_long | negative control | MDD → blond hair (MDD→long IV set) | 165 / 168 |
| LMO_short_to_MDD | leave-many-out | short sleep → MDD | 14 / 26 |
| LMO_MDD_to_long | leave-many-out | MDD → long sleep | 143 / 168 |

Candidates are the harmonised, LD-clumped instruments before Steiger filtering (main) or the main-analysis
instrument set (negative controls). Rows with `included_in_MR = TRUE` are exactly the instruments analysed and
printed in Supplementary Tables 6–7.
Leave-many-out analyses use the main harmonized file (`short_to_MDD.tsv`, `MDD_to_long.tsv`) excluding rows with
`LMO_excluded = TRUE`, i.e. instruments associated at P ≤ 5 × 10⁻⁸ with any candidate predictor in
`data/leave_many_out/predictor_lookup_all_pairs.tsv` (OpenGWAS look-ups; 18 predictors for short sleep → MDD, 9 for
MDD → long sleep).

## Requirements and installation

- R 4.4.2 with data.table 1.15.4, TwoSampleMR 0.6.15, MRPRESSO 1.0, mr.raps 0.4.2
  (`Rscript environment/install_R_packages.R` installs the exact GitHub commits used).
- Python 3.13 with pandas 2.2.3, numpy 2.2.2, matplotlib 3.10.0 (`pip install -r environment/requirements_python.txt`).

Version differences can affect numerical results, particularly Monte Carlo estimates. The workflow records
R/package versions in its generated outputs. The R installer pins the three GitHub packages to specific
commits, but installs data.table and remotes from the configured repository; it does not pin data.table to 1.15.4.

Install R and Python first. On Windows, use a Bash-compatible environment with both interpreters
available (for example, WSL). From the extracted package root, install dependencies:

```bash
Rscript environment/install_R_packages.R
python3 -m pip install -r environment/requirements_python.txt
```

## Run

```bash
bash run_all.sh 12                       # 12 = number of parallel processes (one per analysis)
# or with specific interpreters:
RSCRIPT=/path/to/Rscript PYTHON=/path/to/python bash run_all.sh 12
```


MR-PRESSO with 10,000 simulations takes about 30 minutes for the 165–170-instrument analyses; with 12 cores the
full run was approximately 35 minutes in the original environment; runtime depends on hardware. A single analysis: `Rscript code/R/01_run_mr_analysis.R MDD_to_long`.

## Outputs

- `results/tables/IV_strength_ST8.tsv`, `MR_estimates_ST9_ST11_ST12.tsv`, `Heterogeneity_Q_ST9_ST11_ST12.tsv`,
  `Egger_intercept_ST9_ST11_ST12.tsv`, `MR_PRESSO_tests.tsv`, `Selected_method.tsv`, `Steiger_recheck_main.tsv`,
  `Leave_many_out_screening_summary.tsv`
- `results/comparison_summary.tsv` (and per-value `comparison_with_published.tsv`, `comparison_selected_method.tsv`,
  `comparison_ST10_excluded_associations.tsv`): the run stops with an error if any reproduced value differs from the
  printed tables. ST12 values match when the reproduced value rounds to the printed value (half a unit of the last
  printed digit).
- `results/runs/published_seed/<analysis>/`: per-analysis inputs, raw fits, MR-PRESSO outlier tests, warnings.
- `figures/figure_2_mr_forest_plot.{png,pdf}` and `figures/<analysis>_diagnostics.{png,pdf}` for the four main
  directions.

## Methods implemented

IVW: multiplicative random effects (`TwoSampleMR::mr_ivw`). MR-Egger: `mr_egger_regression` (t(n−2) CI; the
intercept is on the outcome log-odds scale and is not exponentiated). Weighted median: `mr_weighted_median`
(10,000 bootstraps). MR-RAPS (main and leave-many-out analyses): `mr.raps::mr.raps.mle(over.dispersion = TRUE, loss.function = 'huber', pruning = FALSE, se.method = 'sandwich')`,
i.e. additive overdispersion (τ²) estimated jointly with the causal effect, robust Huber loss, no instrument pruning and
sandwich (robust) standard errors; the package's shrinkage wrapper `mr.raps()` and the no-overdispersion
`mr.raps.simple()` are not used.
MR-PRESSO: `MRPRESSO::mr_presso` (10,000 simulations, α = 0.05, outlier and distortion tests; t(n_used−1) CI).
Effects are reported as OR per unit increase in exposure log-odds; SE of OR by the delta method.
BH adjustment within method across the 4 main, 8 negative-control or 2 leave-many-out analyses. Highlighted method: Cochran's Q
P < 0.05 and Egger intercept P < 0.05 → MR-Egger; Q P < 0.05 only → weighted median; otherwise MR-PRESSO
outlier-corrected if estimated, else IVW. Instrument strength: F = (β/SE)²; standard-normal liability R²
(formulas in `code/R/02_collect_results.R`).


## Input definitions and interpretation

- Each harmonized row represents a SNP. `beta.exposure`, `se.exposure`, `beta.outcome` and `se.outcome`
  contain the effect estimates and standard errors; allele and EAF fields describe their orientation.
- `units.exposure` and `units.outcome` identify effect scales. BOLT-LMM sleep effects were converted to
  approximate log-odds effects using beta/[P(1-P)]. Left-handedness beta and SE were approximated from
  Z statistics, EAF and effective N (194,198 cases and 1,534,836 controls). Their OR scale depends on
  this approximation.
- `samplesize.*` and `N_type.*` distinguish study totals from variant-specific sample sizes;
  `effective_n.*` records effective N. Effective N = 4*n_case*n_control/(n_case+n_control).
- `rsq.exposure` and `rsq.outcome` record liability-scale R-squared; `F_exposure` is per-variant F.
  `steiger_dir` and `steiger_pval` describe the Steiger assessment. Primary filtering uses the
  standard-normal liability scale; negative controls do not undergo additional Steiger filtering.
- `included_in_MR` identifies eligible rows and `mr_row_order` preserves analysis order.
  `exclusion_reason` records exclusions; `LMO_excluded` identifies leave-many-out exclusions.
  `presso_outlier_original_seed` records reference-run outlier status; the workflow reruns MR-PRESSO.
- `SNP_alias` records alternative rsIDs, including rs386477048/rs2532281 for GRCh37 17:44,245,171 (G/T).
  `selection`, `index_snp`, `proxy_r2` and `proxy_distance_bp` describe instrument/proxy selection.
- Predictor pairs marked `NOT_FOUND` were not screened and do not establish absence of association.
  For BMI (ieu-b-40), 12 of 26 short-sleep instruments and 101 of 168 MDD instruments were unavailable;
  leave-many-out exclusions may therefore be incomplete.
- MR-PRESSO outlier-corrected BH adjustments retain the full analysis-family size even when some
  analyses have no corrected estimate.

## Reproducibility

Random seeds in `data/manifest.tsv` and SNP row order affect weighted-median bootstrap SEs and
MR-PRESSO results. Different seeds or numerical environments can change Monte Carlo estimates and
outlier sets. An alternative seed can be explored with
`Rscript code/R/01_run_mr_analysis.R <analysis> <seed_offset>`.

The workflow generates `results/` and `figures/` when run; precomputed outputs are not required.
Automated comparisons use the reference TSV included in this package, with full-precision comparisons
for ST8, ST9 and ST11, printed-precision comparisons for ST12, and separate checks for selected methods
and ST10 excluded associations. Consult the generated comparison reports to assess reproducibility.

## Data sources and reuse

GWAS sources and sample definitions are described in the manuscript and its supplementary materials.
The inputs here are SNP-level summary data, not participant-level records. Consult the original GWAS
providers for data-use and redistribution terms; a repository code license does not override those terms.
