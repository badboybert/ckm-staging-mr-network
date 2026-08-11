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
suppressMessages({library(TwoSampleMR)})
BASE <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised_eas")
fmt <- function(f,ty) format_data(read.delim(f,stringsAsFactors=FALSE), type=ty, snp_col="SNP",
  beta_col="beta", se_col="se", effect_allele_col="effect_allele", other_allele_col="other_allele",
  eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
E <- fmt(file.path(INSTR,"KoGES_DM.eas.clumped.tsv"),"exposure")
O <- fmt(file.path(HARM,"KoGES_DM__KoGES_BMI.outcome.tsv"),"outcome")
h <- harmonise_data(E,O,action=2); h <- h[h$mr_keep,]

cat("=== LOO drop MC4R rs6567160 alone (random-effects IVW, TwoSampleMR default) ===\n")
h1 <- h[h$SNP!="rs6567160",]
print(mr(h1, method_list=c("mr_ivw","mr_ivw_fe"))[,c("method","nsnp","b","se","pval")], digits=4, row.names=FALSE)
cat("Q for that subset:\n"); print(mr_heterogeneity(h1, method_list="mr_ivw")[,c("Q","Q_df","Q_pval")], digits=4, row.names=FALSE)

cat("\n=== Positive control orientation: KoGES BMI->DM ===\n")
E2 <- fmt(file.path(INSTR,"KoGES_BMI.eas.clumped.tsv"),"exposure")
O2 <- fmt(file.path(HARM,"KoGES_BMI__KoGES_DM.outcome.tsv"),"outcome")
h2 <- harmonise_data(E2,O2,action=2); h2 <- h2[h2$mr_keep,]
print(mr(h2, method_list=c("mr_ivw","mr_weighted_median"))[,c("method","nsnp","b","se","pval")], digits=4, row.names=FALSE)

cat("\n=== F-stat: on 9 mr_keep vs all 10 instruments ===\n")
cat(sprintf("mean F all-10 instruments = %.1f\n", mean((E$beta.exposure/E$se.exposure)^2)))
cat(sprintf("mean F 9 mr_keep          = %.1f\n", mean((h$beta.exposure/h$se.exposure)^2)))
