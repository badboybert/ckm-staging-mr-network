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
# Step 29 (ITEM 3): KoGES within-population MR for the ->BMI triangulation.
# Run: Rscript 29_koges_mr.R
suppressMessages({library(TwoSampleMR)})
BASE <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised_eas")
RES <- file.path(BASE,"results")
fmt <- function(f,ty) format_data(read.delim(f,stringsAsFactors=FALSE), type=ty, snp_col="SNP",
  beta_col="beta", se_col="se", effect_allele_col="effect_allele", other_allele_col="other_allele",
  eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

edges <- list(c("KoGES_DM","KoGES_BMI"), c("KoGES_DM","KoGES_WAIST"),
              c("KoGES_BMI","KoGES_DM"), c("KoGES_WAIST","KoGES_DM"))
rows <- list()
for (e in edges) {
  exp <- e[1]; out <- e[2]
  ef <- file.path(INSTR, paste0(exp,".eas.clumped.tsv"))
  of <- file.path(HARM, paste0(exp,"__",out,".outcome.tsv"))
  E <- fmt(ef,"exposure"); O <- fmt(of,"outcome")
  h <- harmonise_data(E,O,action=2); h <- h[h$mr_keep,]
  m <- mr(h, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  ei <- tryCatch(mr_pleiotropy_test(h), error=function(z) data.frame(egger_intercept=NA,pval=NA))
  het <- tryCatch(mr_heterogeneity(h), error=function(z) NULL)
  st <- tryCatch(directionality_test(h), error=function(z) data.frame(correct_causal_direction=NA))
  gv <- function(meth,col){v <- m[m$method==meth,col]; if(length(v)) v[1] else NA}
  # mean-F of instruments (exposure): (beta/se)^2 averaged
  Fstat <- mean((E$beta.exposure/E$se.exposure)^2, na.rm=TRUE)
  Qp <- if(!is.null(het)) het$Q_pval[het$method=="Inverse variance weighted"][1] else NA
  rows[[paste0(exp,"__",out)]] <- data.frame(
    exposure=exp, outcome=out, nsnp=gv("Inverse variance weighted","nsnp"),
    ivw_b=gv("Inverse variance weighted","b"), ivw_se=gv("Inverse variance weighted","se"),
    ivw_p=gv("Inverse variance weighted","pval"), wm_b=gv("Weighted median","b"),
    wm_p=gv("Weighted median","pval"), egger_b=gv("MR Egger","b"),
    egger_int=ei$egger_intercept[1], egger_int_p=ei$pval[1], Q_p=Qp,
    meanF=Fstat, steiger=st$correct_causal_direction[1], stringsAsFactors=FALSE)
  cat(sprintf("%-12s -> %-12s : nSNP=%2d IVW b=%+.4f se=%.4f p=%.2e | WM b=%+.4f p=%.2e | eggInt_p=%.3f meanF=%.0f Steiger=%s\n",
    exp,out,rows[[paste0(exp,"__",out)]]$nsnp, rows[[paste0(exp,"__",out)]]$ivw_b,
    rows[[paste0(exp,"__",out)]]$ivw_se, rows[[paste0(exp,"__",out)]]$ivw_p,
    rows[[paste0(exp,"__",out)]]$wm_b, rows[[paste0(exp,"__",out)]]$wm_p,
    ifelse(is.na(ei$pval[1]),NA,ei$pval[1]), Fstat, rows[[paste0(exp,"__",out)]]$steiger))
}
res <- do.call(rbind, rows); write.csv(res, file.path(RES,"koges_triangulation.csv"), row.names=FALSE)
cat("\nWrote results/koges_triangulation.csv with", nrow(res), "KoGES edges\n")
