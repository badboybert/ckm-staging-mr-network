# --- portable roots -----------------------------------------------------------------------------
# Injected by scripts/build_repo.py; see the note in the Python header. Set CKM_P4_BASE to override.
P4_BASE <- Sys.getenv("CKM_P4_BASE")
if (!nzchar(P4_BASE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .self <- sub("^--file=", "", .args[grep("^--file=", .args)])
  .here <- if (length(.self)) dirname(.self[1]) else "."
  P4_BASE <- normalizePath(file.path(.here, ".."), winslash = "/", mustWork = FALSE)
}
P4_ROOT <- Sys.getenv("CKM_P4_ROOT"); if (!nzchar(P4_ROOT)) P4_ROOT <- dirname(P4_BASE)
CKM_ROOT <- Sys.getenv("CKM_ROOT");   if (!nzchar(CKM_ROOT)) CKM_ROOT <- dirname(P4_ROOT)
SHARED_LIB <- Sys.getenv("CKM_SHARED_LIB")
if (!nzchar(SHARED_LIB)) SHARED_LIB <- file.path(CKM_ROOT, "_shared")
# ------------------------------------------------------------------------------------------------
# Step 53 (H5): complementary pleiotropy-robust methods for the reverse feedback edges CAD->T2D and
# HF->T2D, to reduce dependence on CAUSE (which failed its HDL->CAD negative control; reviewer M2).
# Methods that do NOT share CAUSE's genome-wide-sharing assumptions: weighted median, weighted mode,
# penalised weighted median, MR-Egger (+ intercept), and radial-IVW outlier scan (RadialMR).
# (mr.raps / MendelianRandomization::mr_conmix were not installable in this environment; the mode-based
#  and radial estimators provide the complementary, assumption-diverse check.)
# Run: Rscript 53_h5_reverse_pleiotropy.R  ->  results/h5_reverse_pleiotropy.txt
suppressMessages({library(TwoSampleMR); library(RadialMR)})
set.seed(20260718)
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised")
sink(file.path(BASE, "results/h5_reverse_pleiotropy.txt"), split = TRUE)
cat("=== H5: pleiotropy-robust methods for the reverse feedback edges (independent of CAUSE) ===\n\n")

edges <- list(c("CAD", "T2D"), c("HF", "T2D"))
for (ed in edges) {
  exp <- ed[1]; out <- ed[2]
  ef <- file.path(INSTR, paste0(exp, ".clumped.tsv"))
  of <- file.path(HARM, paste0(exp, "__", out, ".outcome.tsv"))
  e <- format_data(read.delim(ef, stringsAsFactors = FALSE), type = "exposure",
        snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele",
        other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval",
        samplesize_col = "N", phenotype_col = "Phenotype")
  o <- format_data(read.delim(of, stringsAsFactors = FALSE), type = "outcome",
        snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele",
        other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval", phenotype_col = "Phenotype")
  h <- harmonise_data(e, o, action = 2); h <- h[h$mr_keep, ]
  cat(sprintf("--- %s -> %s : %d instruments ---\n", exp, out, nrow(h)))
  m <- mr(h, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                             "mr_weighted_mode", "mr_penalised_weighted_median"))
  for (i in seq_len(nrow(m)))
    cat(sprintf("   %-28s b=%+.4f se=%.4f p=%.3g\n", m$method[i], m$b[i], m$se[i], m$pval[i]))
  pl <- tryCatch(mr_pleiotropy_test(h), error = function(z) NULL)
  if (!is.null(pl))
    cat(sprintf("   MR-Egger intercept = %+.4f (p=%.3g)  [directional pleiotropy]\n",
                pl$egger_intercept, pl$pval))
  rad <- tryCatch({
    rin <- format_radial(h$beta.exposure, h$beta.outcome, h$se.exposure, h$se.outcome, h$SNP)
    ivw_radial(rin, alpha = 0.05)
  }, error = function(z) NULL)
  if (!is.null(rad)) {
    no <- tryCatch(nrow(rad$outliers), error = function(z) 0)
    cat(sprintf("   radial-IVW b=%+.4f; Q-outliers flagged=%s\n",
                rad$coef[1, 1], ifelse(is.null(no) || is.na(no), "0", no)))
  }
  cat("\n")
}
cat("READ: an edge is robust to CAUSE's failure mode if the mode/median-based estimators (which assume\n")
cat("the plurality/majority of instruments are valid) AGREE with IVW in sign and remain non-null, with\n")
cat("no dominant Egger-intercept. Report these as the CAUSE-independent check on CAD->T2D and HF->T2D.\n")
sink()
cat("wrote results/h5_reverse_pleiotropy.txt\n")
