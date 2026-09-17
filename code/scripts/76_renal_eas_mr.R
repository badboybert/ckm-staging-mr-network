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
# Step 76b: the two East-Asian renal edges (SBP->CKD, BMI->CKD) through the PRIMARY estimator.
# Same battery, same harmonise action and same seed as 06_mr.R / 51_ckd_mr.R, so the renal contrast
# is estimated the way every other edge in the paper is (round-7 C018; PI decision Q3).
# Input : data/harmonised_renal/<exp>_eas__CKD.outcome.tsv  (written by 76_renal_eas_extract.py)
# Output: results/renal_eas_primary.csv
suppressMessages({library(TwoSampleMR)})
set.seed(20260718)
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised_renal")
RES   <- file.path(BASE, "results")

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

rows <- list()
for (exp in c("SBP", "BMI")) {
  of <- file.path(HARM, paste0(exp, "_eas__CKD.outcome.tsv"))
  ef <- file.path(INSTR, paste0(exp, ".eas.clumped.tsv"))
  stopifnot(file.exists(of), file.exists(ef))
  e <- fmt_exp(ef); o <- fmt_out(of)
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep, ]
  stopifnot(nrow(h) >= 3)
  m  <- mr(h, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  pl <- tryCatch(mr_pleiotropy_test(h), error=function(z) data.frame(egger_intercept=NA, pval=NA))
  het<- tryCatch(mr_heterogeneity(h, method_list="mr_ivw"), error=function(z) data.frame(Q=NA, Q_pval=NA))
  st <- tryCatch(directionality_test(h), error=function(z) data.frame(correct_causal_direction=NA, steiger_pval=NA))
  gv <- function(mm, meth, col) { v <- mm[mm$method == meth, col]; if (length(v)) v[1] else NA }
  rows[[exp]] <- data.frame(
    ancestry="EAS", exposure=exp, outcome="CKD", nsnp=gv(m, "Inverse variance weighted", "nsnp"),
    ivw_b=gv(m, "Inverse variance weighted", "b"), ivw_se=gv(m, "Inverse variance weighted", "se"),
    ivw_p=gv(m, "Inverse variance weighted", "pval"),
    egger_b=gv(m, "MR Egger", "b"), egger_p=gv(m, "MR Egger", "pval"),
    egger_intercept=pl$egger_intercept[1], egger_intercept_p=pl$pval[1],
    wm_b=gv(m, "Weighted median", "b"), wm_p=gv(m, "Weighted median", "pval"),
    Q=het$Q[1], Q_p=het$Q_pval[1],
    steiger_correct=st$correct_causal_direction[1], steiger_p=st$steiger_pval[1],
    stringsAsFactors=FALSE)
  cat(sprintf("EAS %-4s -> CKD : nSNP=%3d  IVW b=%+.4f se=%.4f p=%.3e  Q_p=%.2e\n",
      exp, rows[[exp]]$nsnp, rows[[exp]]$ivw_b, rows[[exp]]$ivw_se, rows[[exp]]$ivw_p, rows[[exp]]$Q_p))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES, "renal_eas_primary.csv"), row.names=FALSE)
cat("\nWrote results/renal_eas_primary.csv with", nrow(res), "edges\n")
