# ASD and imprinting control region (ICR) methylation: CpG-level regression

Code for the CpG-level association analyses of the manuscript
*"Associations Between ASD Symptom Burden and Epigenetic Variation at Imprinting Control Regions"*
(Charifa Abdrahamane), using the Illumina Human Imprintome array and the
[TDHIA](https://github.com/hoyo-lab-ncsu/hoyo_published_analyses/tree/main/2026_JimaCorliss_tdhia_package) R package (Jima and Corliss).

## What is analyzed

| Set | Samples | Outcome | Model (TDHIA `analyze_association`) |
|---|---|---|---|
| Continuous ASD score | 387 | ASD symptom score (BITSEA ASD total score) | CpG-level Gaussian (linear) regression |
| Continuous ASD score, top 20% vs rest | 387 (all kept) | `asd_total_score >= 80th percentile` (=1) vs the rest (=0) | CpG-level binomial (logistic) regression |
| Confirmed ASD case-control | 36 (12 cases, 24 matched controls) | case status | CpG-level binomial (logistic) regression |

Covariates: continuous set = maternal age at delivery, paternal age, infant sex, maternal pre-pregnancy BMI,
maternal smoking during pregnancy; case-control set = maternal age at delivery and infant sex.

Statistical details (implemented inside TDHIA 0.1.0):
* One model per CpG. Missing values (betas set to missing by the detection p-value filter, and missing covariates) are
  imputed with multiple imputation by chained equations (`mice` 3.19.0, default methods, 20 iterations,
  number of imputations = % of missing cells in the model, rounded up, minimum 5; no imputation if nothing is missing).
  The random seed is fixed inside TDHIA (seed 0), so results are reproducible. Estimates are pooled with Rubin's rules.
* Multiple testing: Benjamini-Hochberg FDR over the number of ICRs tested (`n_p_adj` = number of rows of the ICR matrix:
  1,058 in the continuous set, 1,054 in the case-control set).

## Data

* **Raw IDAT files, processed data and sample information:** deposited in NCBI GEO as two series (continuous set and case-control set). The accession numbers will be added once the data are released.
* Each GEO series contains the IDAT files, `processed_beta_values.csv` (probe-level betas of all 22,819 array probes, replicate probes **not** merged) and a metadata spreadsheet with the sample information.
* Participant-level data of the Newborn Epigenetics STudy (NEST) beyond what is deposited in GEO are not public (participant privacy).

## How to run

1. Install R 4.4.3, the TDHIA package (0.1.0, see the TDHIA repository) and the packages listed below.
2. Open `ASD_ICR_cpg_regression.R` in RStudio and set the working directory to its folder.
3. Edit **SECTION 0** (`BUILD_MATRICES`, `tdhia_dir`) and run the script section by section.

The script has separate sections:

| Section | What it does | Output (in `output/`) |
|---|---|---|
| 2 | IDAT -> probe betas -> quality filters -> CpG matrix -> ICR matrix (TDHIA); slow, reads ~850 IDAT files | `beta_matrices/` |
| 3 | **Analysis A**: continuous ASD score, linear regression, FDR < 0.05 | `results/A_*` |
| 4 | **Analysis B**: binarized score (top 20% vs rest), logistic regression, FDR < 0.05 | `results/B_*` |
| 5 | **Analysis C**: confirmed ASD case-control, logistic regression, FDR < 0.05 and 0.15 | `results/C_*` |
| 6 | R and package versions | `results/sessionInfo.txt` |

Put the IDAT files in `data/idat/` and create `data/metadata_continuous.csv` and
`data/metadata_case_control.csv` (required columns below).

* continuous set: `patient_id` (IDAT base name, e.g. `207448290029_R09C02`), `asd_total_score`, `mom_age_delv`,
  `paternal_age`, `baby_gender`, `bmi_lmp_kgm2`, `maternal_smoking2`
* case-control set: `patient_id`, `case` (1 = case, 0 = control), `baby_gender`, `mom_age_delv`

### Required R packages

`tdhia` (0.1.0), `sesame` (1.24.0) and `mice` (3.19.0) are used by TDHIA; plus `dplyr`, `tibble` (`devtools` only if `tdhia_dir` is used). The versions used by the authors are in
`sessionInfo_authors_run.txt`.

## Notes for reproducibility

* **Probe merging.** In the analyses, TDHIA merged replicate probes that target the same CpG before the beta calculation
  (`merge_replicates = "pre_beta"`): fluorescence is averaged per color channel, and for each sample and CpG the merged or the
  single-probe value with the lowest detection p-value is kept. GEO requires one row per physical probe, so
  `processed_beta_values.csv` was made without merging. Its values are therefore **similar but not identical** to the matrices
  used for the paper. Section 2 starts from the IDAT files and applies the merging, as in the paper. For CpGs with a single
  probe the two files are identical.
* Row names of merged probes (suffix `_MERG`, or the probe chosen for a sample) can differ slightly between runs because they
  are taken from the first sample; this does not change the CpG-level values.
* **Case-control set.** The CpG and ICR matrices were built from the 37 samples with available methylation data and the
  analysis then used 36 children (one case without matched controls was removed). The GEO series contains only these 36 children. Because the quality filters depend on the sample set, building
  the matrices from only 36 samples can give very slightly different filtered CpG/ICR lists.
* Results were produced with R 4.4.3.

## Files in this repository

`ASD_ICR_cpg_regression.R` (the analysis script), `sessionInfo_authors_run.txt` (R and package versions used by the authors), `README.md`.

## Citation

Manuscript citation and Zenodo DOI will be added.
