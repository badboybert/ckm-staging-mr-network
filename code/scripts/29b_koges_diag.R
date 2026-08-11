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
# Step 29b: diagnostics for the decisive KoGES DM->BMI edge (per-SNP Wald, LOO, WM-SE sanity, median methods).
suppressMessages({library(TwoSampleMR)})
BASE <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised_eas")
fmt <- function(f,ty) format_data(read.delim(f,stringsAsFactors=FALSE), type=ty, snp_col="SNP",
  beta_col="beta", se_col="se", effect_allele_col="effect_allele", other_allele_col="other_allele",
  eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

E <- fmt(file.path(INSTR,"KoGES_DM.eas.clumped.tsv"),"exposure")
O <- fmt(file.path(HARM,"KoGES_DM__KoGES_BMI.outcome.tsv"),"outcome")
h <- harmonise_data(E,O,action=2); h <- h[h$mr_keep,]
cat("=== KoGES DM->BMI per-SNP Wald ratios (nearest gene guess) ===\n")
w <- data.frame(SNP=h$SNP, bx=h$beta.exposure, by=h$beta.outcome,
                wald=h$beta.outcome/h$beta.exposure,
                wse=abs(h$se.outcome/h$beta.exposure))
w <- w[order(w$wald),]
print(w, row.names=FALSE, digits=3)
cat(sprintf("\nInstruments: %d | DM-raising alleles oriented; wald = per-log-OR-DM effect on BMI(SD)\n", nrow(h)))
cat(sprintf("median wald = %+.4f | mean wald = %+.4f\n", median(w$wald), mean(w$wald)))

cat("\n=== Full method battery ===\n")
m <- mr(h, method_list=c("mr_ivw","mr_ivw_fe","mr_weighted_median","mr_simple_median",
                         "mr_egger_regression","mr_weighted_mode"))
print(m[,c("method","nsnp","b","se","pval")], row.names=FALSE, digits=4)

cat("\n=== Heterogeneity ===\n")
print(mr_heterogeneity(h)[,c("method","Q","Q_df","Q_pval")], row.names=FALSE, digits=4)

cat("\n=== Leave-one-out (IVW) ===\n")
loo <- mr_leaveoneout(h)
loo <- loo[,c("SNP","b","se","p")]; loo <- loo[order(loo$b),]
print(loo, row.names=FALSE, digits=4)

cat("\n=== WM-SE sanity: manual weighted-median bootstrap SE, seed-free ===\n")
# manual weighted median (Bowden 2016) with 1000 parametric bootstraps
bx<-h$beta.exposure; by<-h$beta.outcome; sx<-h$se.exposure; sy<-h$se.outcome
wm_point <- function(bx,by){ r<-by/bx; wt<-(bx^2)/(sy^2); o<-order(r); r<-r[o]; wt<-wt[o]
  cw<-cumsum(wt)/sum(wt); k<-which(cw>=0.5)[1]
  if(k==1) return(r[1]); pl<-(cw[k]-0.5)/(cw[k]-cw[k-1]); r[k]-(r[k]-r[k-1])*pl }
pt <- wm_point(bx,by)
set.seed(42); B<-2000; ests<-numeric(B)
for(i in 1:B){ bxi<-rnorm(length(bx),bx,sx); byi<-rnorm(length(by),by,sy); ests[i]<-wm_point(bxi,byi) }
cat(sprintf("manual WM point=%+.4f bootSE=%.4f z=%.2f p=%.2e\n", pt, sd(ests), pt/sd(ests), 2*pnorm(-abs(pt/sd(ests)))))
