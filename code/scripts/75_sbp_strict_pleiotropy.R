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
# Step 75 (round-4 Tier-1 item 1): MR-Egger, weighted median and Cochran's Q for every X->SBP edge
# on the STRICT allele-matched instrument set.
#
# WHY THIS EXISTS. When step 73 promoted the strict X->SBP extraction to PRIMARY, it BLANKED the
# Egger/WM/Q columns for those rows, correctly, because the shipped values described the old
# position-only instrument set. Step 63 had only ever computed IVW on the strict set, so nothing
# refilled them. Two things then went wrong silently:
#
#   1. The native-scale ledger (step 57) classifies a reverse edge as "pleiotropy-caveated" when its
#      reverse-arm Egger intercept P < 0.05. For CAD->SBP that field was NaN, and `nan < 0.05` is
#      False in every language that follows IEEE-754 - so a MISSING statistic silently became
#      "no pleiotropy detected" and the transition was recorded as plain DISCORDANT. The ledger
#      output therefore read 3 reverse-supported + 1 caveated, while the Results and the Figure 1
#      legend read 2 + 2. That contradiction is what the round-4 review caught.
#   2. The manuscript still quoted "Egger-intercept P = 6 x 10^-3" for CAD->SBP, which is the
#      POSITION-ONLY value, beside a primary point estimate (+1.48) from the STRICT set - a
#      pleiotropy statistic and an effect estimate computed on different instrument sets.
#
# Both are resolved by computing the missing statistics on the primary set, which is what this does.
# Run: Rscript 75_sbp_strict_pleiotropy.R
# Out: results/sbp_strict_pleiotropy.csv / .txt

suppressMessages({library(TwoSampleMR)})
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised")
RES   <- file.path(BASE, "results")

EXPS <- c("BMI","TG","HDL","TC","LDL","HbA1c","CAD","HF","T2D","eGFR","Stroke")

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N",
  phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N",
  phenotype_col="Phenotype")

rows <- list()
for (e in EXPS) {
  ef <- file.path(INSTR, paste0(e, ".clumped.tsv"))
  of <- file.path(HARM, paste0(e, "__SBP.strict.outcome.tsv"))
  if (!file.exists(ef) || !file.exists(of)) { cat(sprintf("%-12s : no strict file\n", paste0(e,"->SBP"))); next }
  ex <- tryCatch(fmt_exp(ef), error=function(z) NULL)
  ou <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(ex) || is.null(ou)) next
  h <- harmonise_data(ex, ou, action=2); h <- h[h$mr_keep, ]
  if (nrow(h) < 4) { cat(sprintf("%-12s : <4 harmonised SNPs\n", paste0(e,"->SBP"))); next }

  m  <- mr(h, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  pl <- mr_pleiotropy_test(h)          # Egger intercept + its P
  q  <- mr_heterogeneity(h, method_list=c("mr_ivw"))
  g  <- function(nm) { v <- m[m$method == nm, ]; if (nrow(v)) v[1, ] else NULL }
  ivw <- g("Inverse variance weighted"); eg <- g("MR Egger"); wm <- g("Weighted median")

  rows[[e]] <- data.frame(
    exposure=e, outcome="SBP", nsnp=ivw$nsnp,
    ivw_b=ivw$b, ivw_se=ivw$se, ivw_p=ivw$pval,
    egger_b=if (is.null(eg)) NA else eg$b, egger_p=if (is.null(eg)) NA else eg$pval,
    egger_intercept=pl$egger_intercept, egger_intercept_se=pl$se, egger_intercept_p=pl$pval,
    wm_b=if (is.null(wm)) NA else wm$b, wm_p=if (is.null(wm)) NA else wm$pval,
    Q=q$Q[1], Q_p=q$Q_pval[1], stringsAsFactors=FALSE)
  cat(sprintf("%-12s : n=%-4d IVW b=%+.4f  Egger-intercept=%+.4f (P=%.3g)\n",
              paste0(e,"->SBP"), ivw$nsnp, ivw$b, pl$egger_intercept, pl$pval))
}

out <- do.call(rbind, rows)
write.csv(out, file.path(RES, "sbp_strict_pleiotropy.csv"), row.names=FALSE)

sink(file.path(RES, "sbp_strict_pleiotropy.txt"))
cat("=== MR-Egger / weighted median / Cochran Q on the STRICT X->SBP instrument set ===\n\n")
cat("These are the PRIMARY-set values. Step 73 blanked the position-only Egger/WM/Q columns and\n")
cat("nothing refilled them, so the ledger read a NaN intercept as 'no pleiotropy' and downgraded\n")
cat("CAD->SBP from pleiotropy-caveated to plain discordant.\n\n")
print(out, row.names=FALSE)
cat("\nCAD->SBP is the row that decides the ledger classification:\n")
cad <- out[out$exposure == "CAD", ]
if (nrow(cad)) {
  cat(sprintf("  Egger intercept %+.4f (SE %.4f), P = %.4g  ->  %s\n",
              cad$egger_intercept, cad$egger_intercept_se, cad$egger_intercept_p,
              ifelse(cad$egger_intercept_p < 0.05,
                     "directional pleiotropy: reverse-PLEIOTROPY-CAVEATED",
                     "no directional pleiotropy detected: reverse-SUPPORTED")))
  cat(sprintf("  IVW %+.4f vs weighted median %+.4f (pleiotropy-robust)\n", cad$ivw_b, cad$wm_b))
}
sink()
cat("\nwrote results/sbp_strict_pleiotropy.csv / .txt\n")
