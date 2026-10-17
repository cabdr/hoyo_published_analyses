# =============================================================================
#  ASD and imprinting control region (ICR) methylation
#  CpG-level regression with the TDHIA package
#
#  Paper : "Associations Between ASD Symptom Burden and Epigenetic Variation at
#          Imprinting Control Regions"
#  Data  : NCBI GEO (accession numbers will be added once the data are released)
#
#  This ONE script has separate, clearly marked sections:
#    SECTION 0  Settings
#    SECTION 2  IDAT files -> probe, CpG and ICR beta matrices (TDHIA)
#    SECTION 3  Analysis A: continuous ASD score           (linear regression)
#    SECTION 4  Analysis B: binarized ASD score, top 20% vs rest (logistic regression)
#    SECTION 5  Analysis C: confirmed ASD case-control     (logistic regression)
#    SECTION 6  Save session information
#
#  How to run: open this file in RStudio, set the working directory to the
#  folder that contains it (Session > Set Working Directory > To Source File
#  Location), edit SECTION 0, then run the sections in order.
#
#  Statistics (done inside TDHIA 0.1.0, function analyze_association):
#   * one regression per CpG: outcome ~ CpG beta + covariates
#   * missing values (betas set to missing by the detection p-value filter and
#     missing covariates) are imputed with MICE (package 'mice' 3.19.0, default
#     methods, 20 iterations; number of imputations = % of missing cells in the
#     model rounded up, minimum 5; none if nothing is missing). The random seed
#     is fixed inside TDHIA (seed = 0), so results are reproducible.
#     Estimates are pooled with Rubin's rules.
#   * multiple testing: Benjamini-Hochberg FDR over the number of ICRs tested
#     (n_p_adj = number of rows of the ICR matrix)
# =============================================================================


# =============================================================================
# SECTION 0  SETTINGS
# =============================================================================

# ---- Which steps to run ------------------------------------------------------
BUILD_MATRICES    <- TRUE    # TRUE: run SECTION 2 (slow, reads ~850 IDAT files).
                             # FALSE: reuse the matrices saved by an earlier run

# ---- Folders -------------------------------------------------------------------
data_dir    <- "data"                                 # input data
idat_dir    <- file.path(data_dir, "idat")            # IDAT files of BOTH sets
beta_dir    <- file.path("output", "beta_matrices")   # probe / CpG / ICR matrices
results_dir <- file.path("output", "results")         # regression results

# One sample-information file per set. Create these files with the columns
# described in the README.
meta_file <- c(
  continuous   = file.path(data_dir, "metadata_continuous.csv"),
  case_control = file.path(data_dir, "metadata_case_control.csv")
)

# ---- TDHIA package ---------------------------------------------------------------
# NULL    : TDHIA is installed as an R package -> library(tdhia)
# "path"  : folder of a local copy of TDHIA   -> devtools::load_all(path)
tdhia_dir <- NULL

# ---- Model settings used in the paper (do not change to reproduce the results) ----
covars_continuous   <- c("mom_age_delv", "paternal_age", "baby_gender",
                         "bmi_lmp_kgm2", "maternal_smoking2")
covars_case_control <- c("mom_age_delv", "baby_gender")
top_fraction_cut    <- 0.80    # "top 20%" = ASD score >= 80th percentile

# ---- Packages and folders ----------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tibble) })
for (d in c(data_dir, beta_dir, results_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

if (is.null(tdhia_dir)) {
  suppressPackageStartupMessages(library(tdhia))
} else {
  suppressPackageStartupMessages(devtools::load_all(tdhia_dir, quiet = TRUE))
}
message("TDHIA version: ", as.character(utils::packageVersion("tdhia")))


# =============================================================================
# Helper functions (used by the sections below)
# =============================================================================

# Run TDHIA on one outcome, annotate the CpGs, and save the result tables.
#   R    : data frame, one column (the outcome), row names = sample IDs
#   P    : data frame, samples x CpGs (beta values)
#   Pe   : data frame, samples x covariates
#   fdr_cutoffs : FDR thresholds for which "significant" tables are written
run_tdhia_cpg_regression <- function(R, P, Pe, family, n_p_adj, max_p_val, prefix, fdr_cutoffs) {
  stopifnot(identical(rownames(R), rownames(P)), identical(rownames(P), rownames(Pe)))
  message("\n--- ", prefix, ": ", nrow(P), " samples, ", ncol(P), " CpGs, family = ", family, " ---")

  fit <- tdhia::analyze_association(
    R         = R,
    P         = P,
    Pe        = Pe,
    family    = family,
    n_p_adj   = n_p_adj,                 # number of ICRs tested
    db_flag   = FALSE,
    rm.na.all = FALSE,
    verbose   = TRUE,
    impute_na = TRUE,                    # MICE (seed fixed inside TDHIA)
    max_p_val = max_p_val,               # only used for the printed summary
    n.cores   = max(parallel::detectCores() - 2, 1)
  )
  saveRDS(fit, file.path(results_dir, paste0(prefix, "_tdhia_fit.rds")))

  # One row per CpG; add ICR / gene annotation from TDHIA
  res <- fit$imp_site
  res$cpg_id <- res$Variable
  ann <- tdhia::add_metadata_to_imp_sites(res$cpg_id, imp_type = "cpg")
  res_annot <- dplyr::bind_cols(res, dplyr::select(ann, -dplyr::any_of("cpg_id"))) |>
    dplyr::select(-dplyr::any_of("Formula"))
  write.csv(res_annot, file.path(results_dir, paste0(prefix, "_all_CpGs_annotated.csv")), row.names = FALSE)

  # Significant CpGs, and the ICRs they belong to, at each FDR cut-off
  for (thr in fdr_cutoffs) {
    sig <- res_annot |> dplyr::filter(ADJ_P_VAL < thr)
    sig_icr <- sig |>
      dplyr::filter(!is.na(icr_id)) |>
      dplyr::group_by(icr_id, icr_conf, is_icr_zinc, Genomic.Coordinates,
                      Nearest.Transcript, Distance.to.Nearest.Transcript) |>
      dplyr::summarise(n_cpg = dplyr::n_distinct(cpg_id),
                       mean_estimate = mean(Estimate, na.rm = TRUE),
                       frac_positive = mean(Estimate > 0, na.rm = TRUE),
                       best_cpg_fdr  = min(ADJ_P_VAL, na.rm = TRUE),
                       cpg_list = paste(unique(cpg_id), collapse = "; "),
                       .groups = "drop") |>
      dplyr::arrange(best_cpg_fdr)
    write.csv(sig,     file.path(results_dir, paste0(prefix, "_sig_CpGs_FDR", thr, ".csv")), row.names = FALSE)
    write.csv(sig_icr, file.path(results_dir, paste0(prefix, "_sig_ICRs_FDR", thr, ".csv")), row.names = FALSE)
    message("FDR < ", thr, ": ", nrow(sig), " CpGs in ", nrow(sig_icr), " ICRs")
  }
  invisible(res_annot)
}


# =============================================================================
# SECTION 2  IDAT FILES -> PROBE, CpG AND ICR BETA MATRICES  (TDHIA)
#
# Steps (settings of the paper):
#   1. IDAT -> probe-level betas and detection p-values. Replicate probes that
#      target the same CpG are merged before the beta calculation
#      (merge_replicates = "pre_beta"). No quantile normalization.
#   2. Quality filters: betas with detection p > 0.20 are set to missing; probes
#      failing in more than 20% of samples and unmapped probes are removed.
#   3. Probes -> CpGs (no smoothing across adjacent CpGs).
#   4. CpGs -> ICRs (ICRs failing in more than 20% of samples are removed).
# Output: output/beta_matrices/<set>_tdhia.rda (object 'data_beta') + CSV copies.
# =============================================================================
if (BUILD_MATRICES) {
  for (set in names(meta_file)) {
    message("\n=== SECTION 2: building matrices for ", set, " ===")
    meta <- read.csv(meta_file[[set]], check.names = FALSE)
    stopifnot(!anyDuplicated(meta$patient_id))

    data_beta <- list()

    data_beta$probe_beta <- tdhia::load_idata_to_probes(
      idat_dir_paths     = idat_dir,
      idat_basenames     = meta$patient_id,
      quantile_norm      = FALSE,
      db_flag            = FALSE,
      multicore          = TRUE,
      merge_replicates   = "pre_beta",
      enforce_req_idats  = TRUE,
      enforce_idat_names = FALSE
    )

    data_beta$filt_probe_beta <- tdhia::filter_probes(
      probe_beta              = data_beta$probe_beta,
      discard_unmapped_probes = TRUE,
      max_sig_pval            = 0.20,
      set_failed_betas_na     = TRUE,
      max_probe_fail_rate     = 0.20,
      discard_failed_probes   = TRUE,
      min_design_score        = NA,
      db_flag                 = FALSE
    )

    data_beta$cpg_beta <- tdhia::convert_probes_to_cpgs(
      data_beta$filt_probe_beta,
      quantile_norm   = FALSE,
      db_flag         = FALSE,
      smooth_adj_cpgs = FALSE
    )

    data_beta$icr_beta <- tdhia::convert_cpgs_to_icrs(
      data_beta$cpg_beta,
      max_icr_fail_rate = 0.20
    )

    cpg <- data_beta$cpg_beta$cpg_beta_df
    icr <- data_beta$icr_beta$icr_beta_df
    message("Samples ", nrow(meta), " | CpG matrix ", nrow(cpg), " x ", ncol(cpg),
            " | ICR matrix ", nrow(icr), " x ", ncol(icr))

    save(data_beta, file = file.path(beta_dir, paste0(set, "_tdhia.rda")))
    write.csv(cbind(ID_REF = rownames(cpg), cpg), file.path(beta_dir, paste0(set, "_cpg_beta.csv")), row.names = FALSE)
    write.csv(cbind(ID_REF = rownames(icr), icr), file.path(beta_dir, paste0(set, "_icr_beta.csv")), row.names = FALSE)
  }
}

# Load the saved matrices and sample information of both sets
load(file.path(beta_dir, "continuous_tdhia.rda"));   cont_beta <- data_beta
load(file.path(beta_dir, "case_control_tdhia.rda")); cc_beta   <- data_beta
rm(data_beta)
meta_cont <- read.csv(meta_file[["continuous"]],   check.names = FALSE)
meta_cc   <- read.csv(meta_file[["case_control"]], check.names = FALSE)


# =============================================================================
# SECTION 3  ANALYSIS A: CONTINUOUS ASD SCORE  (n = 387)
#   Outcome   : ASD symptom score (BITSEA ASD total score), continuous
#   Model     : CpG-level linear (Gaussian) regression
#   Covariates: maternal age at delivery, paternal age, infant sex,
#               maternal pre-pregnancy BMI, maternal smoking in pregnancy
#   Reported  : FDR < 0.05
# =============================================================================
cpg_cont <- cont_beta$cpg_beta$cpg_beta_df      # rows = CpGs, columns = samples
icr_cont <- cont_beta$icr_beta$icr_beta_df      # rows = ICRs

# Put samples in the same order in the metadata and in the CpG matrix
ids_cont  <- intersect(colnames(cpg_cont), meta_cont$patient_id)
meta_c    <- meta_cont[match(ids_cont, meta_cont$patient_id), ]
rownames(meta_c) <- NULL
P_cont    <- as.data.frame(t(cpg_cont[, ids_cont, drop = FALSE]))     # samples x CpGs
covars_c  <- covars_continuous[covars_continuous %in% names(meta_c)]
Pe_cont   <- meta_c |> select(patient_id, all_of(covars_c)) |> column_to_rownames("patient_id")
message("Continuous set: ", nrow(P_cont), " samples (expected 387)")

R_cont <- meta_c |> select(patient_id, asd_total_score) |> column_to_rownames("patient_id")

res_A <- run_tdhia_cpg_regression(
  R = R_cont, P = P_cont, Pe = Pe_cont,
  family      = "gaussian",
  n_p_adj     = nrow(icr_cont),
  max_p_val   = 0.05,
  prefix      = "A_continuous_linear",
  fdr_cutoffs = 0.05
)


# =============================================================================
# SECTION 4  ANALYSIS B: BINARIZED ASD SCORE, TOP 20% vs REST  (n = 387)
#   Outcome   : 1 = ASD score >= 80th percentile (top 20%), 0 = everyone else
#               (all 387 children are kept)
#   Model     : CpG-level logistic (binomial) regression
#   Covariates: same as Analysis A
#   Reported  : FDR < 0.05
# =============================================================================
cut_value <- as.numeric(quantile(meta_c$asd_total_score, probs = top_fraction_cut, na.rm = TRUE))
meta_c$asd_binary_all <- ifelse(is.na(meta_c$asd_total_score), NA,
                                ifelse(meta_c$asd_total_score >= cut_value, 1, 0))
message("Cut-off (", top_fraction_cut * 100, "th percentile) = ", cut_value,
        " | top 20%: ", sum(meta_c$asd_binary_all == 1, na.rm = TRUE),
        " | rest: ", sum(meta_c$asd_binary_all == 0, na.rm = TRUE))

R_bin <- meta_c |> select(patient_id, asd_binary_all) |> column_to_rownames("patient_id")

res_B <- run_tdhia_cpg_regression(
  R = R_bin, P = P_cont, Pe = Pe_cont,       # same CpGs and covariates as Analysis A
  family      = "binomial",
  n_p_adj     = nrow(icr_cont),
  max_p_val   = 0.05,
  prefix      = "B_continuous_top20_logistic",
  fdr_cutoffs = 0.05
)


# =============================================================================
# SECTION 5  ANALYSIS C: CONFIRMED ASD CASE-CONTROL  (12 cases + 24 controls = 36)
#   Outcome   : ASD case (1) vs matched control (0)
#   Model     : CpG-level logistic (binomial) regression
#   Covariates: maternal age at delivery, infant sex (small sample)
#   Reported  : FDR < 0.05, and FDR < 0.15 as exploratory threshold
# =============================================================================
cpg_cc <- cc_beta$cpg_beta$cpg_beta_df
icr_cc <- cc_beta$icr_beta$icr_beta_df

ids_cc <- intersect(colnames(cpg_cc), meta_cc$patient_id)
meta_k <- meta_cc[match(ids_cc, meta_cc$patient_id), ]

# Each case was matched to 2 controls (12 cases, 24 controls). If the metadata has a
# matched-set column, keep only complete matched sets (1 case + 2 controls).
set_col <- intersect(c("matched_set", "case_mom_nestid"), names(meta_k))[1]
if (!is.na(set_col)) {
  grp  <- as.character(meta_k[[set_col]])
  keep <- vapply(seq_len(nrow(meta_k)), function(i) {
    g <- meta_k$case[grp == grp[i]]
    length(g) == 3 && sum(g == 1, na.rm = TRUE) == 1 && sum(g == 0, na.rm = TRUE) == 2
  }, logical(1))
  meta_k <- meta_k[keep, , drop = FALSE]
}
rownames(meta_k) <- NULL

meta_k$case         <- as.integer(meta_k$case)
meta_k$baby_gender  <- as.factor(meta_k$baby_gender)
meta_k$mom_age_delv <- as.numeric(meta_k$mom_age_delv)
stopifnot(sum(meta_k$case == 1) == 12, sum(meta_k$case == 0) == 24)
message("Case-control set: ", nrow(meta_k), " samples (12 cases, 24 controls)")

P_cc   <- as.data.frame(t(cpg_cc[, meta_k$patient_id, drop = FALSE]))
R_cc   <- meta_k |> select(patient_id, case) |> column_to_rownames("patient_id")
covars_k <- covars_case_control[covars_case_control %in% names(meta_k)]
Pe_cc  <- meta_k |> select(patient_id, all_of(covars_k)) |> column_to_rownames("patient_id")

res_C <- run_tdhia_cpg_regression(
  R = R_cc, P = P_cc, Pe = Pe_cc,
  family      = "binomial",
  n_p_adj     = nrow(icr_cc),
  max_p_val   = 0.15,
  prefix      = "C_case_control_logistic",
  fdr_cutoffs = c(0.05, 0.15)
)


# =============================================================================
# SECTION 6  SESSION INFORMATION (R and package versions)
# =============================================================================
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
message("\nDone. Results are in ", results_dir)
