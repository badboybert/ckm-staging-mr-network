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
# EAS MR battery over harmonised_eas edges. Run (PowerShell): Rscript 06_mr_eas.R
suppressMessages({library(TwoSampleMR)})
BASE <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised_eas"); RES <- file.path(BASE,"results")
fmt <- function(f,ty) format_data(read.delim(f,stringsAsFactors=FALSE), type=ty, snp_col="SNP", beta_col="beta",
  se_col="se", effect_allele_col="effect_allele", other_allele_col="other_allele", eaf_col="eaf",
  pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
rows<-list()
for (of in list.files(HARM, pattern="__.*\\.outcome\\.tsv$", full.names=TRUE)) {
  base<-sub("\\.outcome\\.tsv$","",basename(of)); pr<-strsplit(base,"__")[[1]]; exp<-pr[1]; out<-pr[2]
  ef<-file.path(INSTR,paste0(exp,".eas.clumped.tsv")); if(!file.exists(ef)) next
  e<-tryCatch(fmt(ef,"exposure"),error=function(z)NULL); o<-tryCatch(fmt(of,"outcome"),error=function(z)NULL)
  if(is.null(e)||is.null(o)||nrow(o)==0){cat(sprintf("%-5s->%-5s: no data\n",exp,out));next}
  h<-harmonise_data(e,o,action=2); h<-h[h$mr_keep,]
  if(nrow(h)<3){cat(sprintf("%-5s->%-5s: <3 SNPs (%d)\n",exp,out,nrow(h)));next}
  m<-mr(h,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  st<-tryCatch(directionality_test(h),error=function(z)data.frame(correct_causal_direction=NA))
  gv<-function(meth,col){v<-m[m$method==meth,col];if(length(v))v[1] else NA}
  rows[[base]]<-data.frame(exposure=exp,outcome=out,nsnp=gv("Inverse variance weighted","nsnp"),
    ivw_b=gv("Inverse variance weighted","b"),ivw_se=gv("Inverse variance weighted","se"),
    ivw_p=gv("Inverse variance weighted","pval"),wm_b=gv("Weighted median","b"),
    egger_b=gv("MR Egger","b"),steiger=st$correct_causal_direction[1],stringsAsFactors=FALSE)
  cat(sprintf("%-5s -> %-5s : nSNP=%3d IVW b=%+.4f p=%.2e Steiger=%s\n",exp,out,rows[[base]]$nsnp,rows[[base]]$ivw_b,rows[[base]]$ivw_p,rows[[base]]$steiger))
}
res<-do.call(rbind,rows); write.csv(res,file.path(RES,"network_eas_edges.csv"),row.names=FALSE)
cat("\nWrote network_eas_edges.csv with",nrow(res),"EAS edges\n")
