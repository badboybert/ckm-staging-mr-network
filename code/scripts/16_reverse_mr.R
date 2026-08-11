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
# Step 16b: MR for the reverse-arm edges that close the 5 ledger gaps.
# Same battery/output schema as 06_mr.R (forward_local_edges.csv) but for the new edges only,
# written to results/reverse_arms_edges.csv. Run (PowerShell): Rscript 16_reverse_mr.R
suppressMessages({library(TwoSampleMR)})
BASE  <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised")
RES   <- file.path(BASE,"results")

EDGES <- list(c("Stroke","BMI"), c("Stroke","T2D"), c("Stroke","CAD"), c("Stroke","HF"),
              c("CAD","SBP"), c("HF","SBP"), c("Stroke","SBP"))

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

rows <- list()
for (ed in EDGES) {
  exp <- ed[1]; out <- ed[2]; base <- paste0(exp,"__",out)
  ef <- file.path(INSTR, paste0(exp,".clumped.tsv")); of <- file.path(HARM, paste0(base,".outcome.tsv"))
  if (!file.exists(ef) || !file.exists(of)) { cat(sprintf("%-14s : missing files\n",base)); next }
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL); o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e)||is.null(o)||nrow(o)==0) { cat(sprintf("%-14s : no data\n",base)); next }
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep,]
  if (nrow(h)<3) { cat(sprintf("%-14s : <3 SNPs (n=%d)\n",base,nrow(h))); next }
  m  <- mr(h, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  pl <- tryCatch(mr_pleiotropy_test(h), error=function(z) data.frame(egger_intercept=NA,pval=NA))
  het<- tryCatch(mr_heterogeneity(h, method_list="mr_ivw"), error=function(z) data.frame(Q=NA,Q_pval=NA))
  st <- tryCatch(directionality_test(h), error=function(z) data.frame(correct_causal_direction=NA,steiger_pval=NA))
  gv <- function(mm,meth,col){ v<-mm[mm$method==meth,col]; if(length(v)) v[1] else NA }
  rows[[base]] <- data.frame(
    exposure=exp, outcome=out, nsnp=gv(m,"Inverse variance weighted","nsnp"),
    ivw_b=gv(m,"Inverse variance weighted","b"), ivw_se=gv(m,"Inverse variance weighted","se"),
    ivw_p=gv(m,"Inverse variance weighted","pval"),
    egger_b=gv(m,"MR Egger","b"), egger_p=gv(m,"MR Egger","pval"),
    egger_intercept=pl$egger_intercept[1], egger_intercept_p=pl$pval[1],
    wm_b=gv(m,"Weighted median","b"), wm_p=gv(m,"Weighted median","pval"),
    Q=het$Q[1], Q_p=het$Q_pval[1],
    steiger_correct=st$correct_causal_direction[1], steiger_p=st$steiger_pval[1], stringsAsFactors=FALSE)
  cat(sprintf("%-14s nSNP=%3d IVW b=%+.4f p=%.2e Steiger=%s\n", base,
      rows[[base]]$nsnp, rows[[base]]$ivw_b, rows[[base]]$ivw_p, rows[[base]]$steiger_correct))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES,"reverse_arms_edges.csv"), row.names=FALSE)
cat("\nWrote reverse_arms_edges.csv with", nrow(res), "edges\n")
