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
# Step 30b (ITEM 3 robustness): MR for the two-sample KoGES-DM -> TWB-Waist edge (built by 30_koges_dm_to_twb_waist.py).
# Run: Rscript 30b_twb_edge_mr.R
suppressMessages({library(TwoSampleMR)})
BASE <- P4_BASE
fmt <- function(f,ty) format_data(read.delim(f,stringsAsFactors=FALSE), type=ty, snp_col="SNP",
  beta_col="beta", se_col="se", effect_allele_col="effect_allele", other_allele_col="other_allele",
  eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
E <- fmt(file.path(BASE,"data/instruments/KoGES_DM.eas.clumped.tsv"),"exposure")
O <- fmt(file.path(BASE,"data/harmonised_eas/KoGES_DM__TWB_Waist.outcome.tsv"),"outcome")
h <- harmonise_data(E,O,action=2); h <- h[h$mr_keep,]
cat("harmonised SNPs:", nrow(h), "\n")
w<-data.frame(SNP=h$SNP,bx=h$beta.exposure,by=h$beta.outcome,wald=h$beta.outcome/h$beta.exposure)
print(w[order(w$wald),],row.names=FALSE,digits=3)
m <- mr(h, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))
cat("\n"); print(m[,c("method","nsnp","b","se","pval")],row.names=FALSE,digits=4)
cat("\nEgger intercept p:", tryCatch(mr_pleiotropy_test(h)$pval[1],error=function(z)NA),
    "| Q p:", tryCatch(mr_heterogeneity(h)$Q_pval[1],error=function(z)NA),"\n")
