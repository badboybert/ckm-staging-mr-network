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
# Step 26: MR over the FULL bidirectional TPMI edge set (forward RF->disease + reverse disease->RF).
# Auto-discovers harmonised_tpmi/<exp>__<out>.outcome.tsv paired with <exp>.tpmi.clumped.tsv.
# Run: Rscript 26_tpmi_full_mr.R  ->  results/network_tpmi_full.csv
suppressMessages({library(TwoSampleMR)})
BASE  <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised_tpmi")
RES   <- file.path(BASE,"results")
fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

ofiles <- list.files(HARM, pattern="__.*\\.outcome\\.tsv$", full.names=TRUE)
rows <- list()
for (of in ofiles) {
  base <- sub("\\.outcome\\.tsv$","",basename(of)); parts <- strsplit(base,"__")[[1]]
  exp<-parts[1]; out<-parts[2]
  ef <- file.path(INSTR, paste0(exp,".tpmi.clumped.tsv"))
  if (!file.exists(ef)) next
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL); o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e)||is.null(o)||nrow(o)==0) next
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep,]
  if (nrow(h)<3) next
  m  <- mr(h, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  pl <- tryCatch(mr_pleiotropy_test(h), error=function(z) data.frame(egger_intercept=NA,pval=NA))
  het<- tryCatch(mr_heterogeneity(h, method_list="mr_ivw"), error=function(z) data.frame(Q=NA,Q_pval=NA))
  st <- tryCatch(directionality_test(h), error=function(z) data.frame(correct_causal_direction=NA,steiger_pval=NA))
  gv <- function(mm,meth,col){ v<-mm[mm$method==meth,col]; if(length(v)) v[1] else NA }
  rows[[base]] <- data.frame(exposure=exp, outcome=out, nsnp=gv(m,"Inverse variance weighted","nsnp"),
    ivw_b=gv(m,"Inverse variance weighted","b"), ivw_se=gv(m,"Inverse variance weighted","se"),
    ivw_p=gv(m,"Inverse variance weighted","pval"),
    egger_b=gv(m,"MR Egger","b"), egger_intercept=pl$egger_intercept[1], egger_intercept_p=pl$pval[1],
    wm_b=gv(m,"Weighted median","b"), Q=het$Q[1], Q_p=het$Q_pval[1],
    steiger_correct=st$correct_causal_direction[1], steiger_p=st$steiger_pval[1], stringsAsFactors=FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES,"network_tpmi_full.csv"), row.names=FALSE)
cat("Wrote network_tpmi_full.csv with", nrow(res), "edges (forward + reverse)\n")
# quick tally
cat("forward RF->disease sig:", sum(res$ivw_p<0.05 & res$exposure %in% c("BMI","SBP","LDL","HDL","TC","TG","HbA1c"), na.rm=TRUE),
    "| reverse disease->RF sig:", sum(res$ivw_p<0.05 & res$exposure %in% c("T2D","CAD","HF","Stroke","CKD"), na.rm=TRUE), "\n")
