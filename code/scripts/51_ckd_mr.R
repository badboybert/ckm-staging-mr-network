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
# Step 51: CKD->Y), same estimators as 06_mr.R.
# CKD instruments = CKD.clumped.tsv; other exposures = <exp>.clumped.tsv. Outcome files in harmonised_ckd/.
# Run: Rscript 41_cac_mr.R  ->  results/network_ckd.csv
suppressMessages({library(TwoSampleMR)})
set.seed(20260718)
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised_ckd")
RES   <- file.path(BASE, "results")

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

ofiles <- list.files(HARM, pattern="__.*\\.outcome\\.tsv$", full.names=TRUE)
rows <- list()
for (of in ofiles) {
  base <- sub("\\.outcome\\.tsv$", "", basename(of))
  parts <- strsplit(base, "__")[[1]]; exp <- parts[1]; out <- parts[2]
  ef <- file.path(INSTR, paste0(exp, ".clumped.tsv"))
  if (!file.exists(ef)) next
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL)
  o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e) || is.null(o) || nrow(o) == 0) next
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep, ]
  if (nrow(h) < 3) { cat(sprintf("%-6s -> %-6s : <3 SNPs\n", exp, out)); next }
  m  <- mr(h, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  pl <- tryCatch(mr_pleiotropy_test(h), error=function(z) data.frame(egger_intercept=NA, pval=NA))
  het <- tryCatch(mr_heterogeneity(h, method_list="mr_ivw"), error=function(z) data.frame(Q=NA, Q_pval=NA))
  st <- tryCatch(directionality_test(h), error=function(z) data.frame(correct_causal_direction=NA, steiger_pval=NA))
  gv <- function(mm, meth, col) { v <- mm[mm$method == meth, col]; if (length(v)) v[1] else NA }
  rows[[base]] <- data.frame(
    exposure=exp, outcome=out, nsnp=gv(m, "Inverse variance weighted", "nsnp"),
    ivw_b=gv(m, "Inverse variance weighted", "b"), ivw_se=gv(m, "Inverse variance weighted", "se"),
    ivw_p=gv(m, "Inverse variance weighted", "pval"),
    egger_b=gv(m, "MR Egger", "b"), egger_p=gv(m, "MR Egger", "pval"),
    egger_intercept=pl$egger_intercept[1], egger_intercept_p=pl$pval[1],
    wm_b=gv(m, "Weighted median", "b"), wm_p=gv(m, "Weighted median", "pval"),
    Q=het$Q[1], Q_p=het$Q_pval[1],
    steiger_correct=st$correct_causal_direction[1], steiger_p=st$steiger_pval[1],
    stringsAsFactors=FALSE)
  cat(sprintf("%-6s -> %-6s : nSNP=%3d  IVW b=%+.4f p=%.2e  Steiger=%s\n",
      exp, out, rows[[base]]$nsnp, rows[[base]]$ivw_b, rows[[base]]$ivw_p, rows[[base]]$steiger_correct))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES, "network_ckd.csv"), row.names=FALSE)
cat("\nWrote results/network_ckd.csv with", nrow(res), "CKD edges\n")
