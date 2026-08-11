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
## Independent adversarial reproduction of KoGES DM->BMI (ITEM 3)
## Recomputes WM point, WM bootstrap SE, IVW fixed/random, Q, I2, sign test,
## heterogeneity-inflated WM SE. NO reliance on the committed CSV.
suppressMessages({library(TwoSampleMR)})
BASE <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised_eas")
fmt <- function(f,ty) format_data(read.delim(f,stringsAsFactors=FALSE), type=ty, snp_col="SNP",
  beta_col="beta", se_col="se", effect_allele_col="effect_allele", other_allele_col="other_allele",
  eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

E <- fmt(file.path(INSTR,"KoGES_DM.eas.clumped.tsv"),"exposure")
O <- fmt(file.path(HARM,"KoGES_DM__KoGES_BMI.outcome.tsv"),"outcome")
h <- harmonise_data(E,O,action=2)
cat("=== harmonise: rows total", nrow(h), " mr_keep", sum(h$mr_keep), "===\n")
cat("dropped by mr_keep:", paste(h$SNP[!h$mr_keep], collapse=","), "\n")
cat("palindromic flags:\n"); print(h[,c("SNP","palindromic","ambiguous","mr_keep")])
h <- h[h$mr_keep,]

bx<-h$beta.exposure; by<-h$beta.outcome; sx<-h$se.exposure; sy<-h$se.outcome
r <- by/bx; wr_se <- abs(sy/bx)                # first-order Wald ratio SE (ignores sx)
cat("\n=== per-SNP Wald ratios (sorted) ===\n")
ord<-order(r); tab<-data.frame(SNP=h$SNP[ord], bx=bx[ord], by=by[ord], wald=r[ord], wr_se=wr_se[ord])
print(tab, digits=4, row.names=FALSE)
cat(sprintf("n negative wald = %d / %d ; median wald = %+.4f\n", sum(r<0), length(r), median(r)))

## ---- IVW fixed & random by hand ----
w <- 1/wr_se^2
ivw_fe <- sum(w*r)/sum(w); ivw_fe_se <- sqrt(1/sum(w))
resid <- r - ivw_fe; Q <- sum(w*resid^2); df<-length(r)-1
I2 <- max(0,(Q-df)/Q)
# multiplicative random-effects (TwoSampleMR default): inflate SE by sqrt(Q/df) if >1
phi <- max(1, Q/df); ivw_re_se <- ivw_fe_se*sqrt(phi)
cat(sprintf("\nIVW fixed  b=%+.4f se=%.4f p=%.2e\n", ivw_fe, ivw_fe_se, 2*pnorm(-abs(ivw_fe/ivw_fe_se))))
cat(sprintf("IVW random b=%+.4f se=%.4f p=%.3f  (Q=%.1f df=%d Q_p=%.2e I2=%.2f)\n",
    ivw_fe, ivw_re_se, 2*pnorm(-abs(ivw_fe/ivw_re_se)), Q, df, pchisq(Q,df,lower.tail=FALSE), I2))

## ---- weighted median point (Bowden 2016) ----
wm_point <- function(bx,by,sy){ rr<-by/bx; wt<-(bx^2)/(sy^2); o<-order(rr); rr<-rr[o]; wt<-wt[o]
  cw<-cumsum(wt)/sum(wt); k<-which(cw>=0.5)[1]; if(k==1) return(rr[1])
  pl<-(cw[k]-0.5)/(cw[k]-cw[k-1]); rr[k]-(rr[k]-rr[k-1])*pl }
pt <- wm_point(bx,by,sy)
set.seed(1); B<-10000; ests<-numeric(B)
for(i in 1:B){ bxi<-rnorm(length(bx),bx,sx); byi<-rnorm(length(by),by,sy); ests[i]<-wm_point(bxi,byi,sy) }
wm_se <- sd(ests)
cat(sprintf("\nManual WM point=%+.4f  bootSE=%.4f  z=%.2f  p=%.2e (B=%d, seed1)\n",
    pt, wm_se, pt/wm_se, 2*pnorm(-abs(pt/wm_se)), B))

## ---- TwoSampleMR's own weighted median (its bootstrap) ----
mm <- mr(h, method_list=c("mr_ivw","mr_ivw_fe","mr_weighted_median","mr_simple_median",
                          "mr_weighted_mode","mr_egger_regression"))
cat("\n=== TwoSampleMR method battery ===\n")
print(mm[,c("method","nsnp","b","se","pval")], digits=4, row.names=FALSE)

## ---- sign / binomial test: is 'majority negative' beyond chance? ----
cat(sprintf("\nBinomial sign test (%d neg of %d, H0 p=0.5): p=%.4f\n",
    sum(r<0), length(r), binom.test(sum(r<0), length(r), 0.5)$p.value))

## ---- heterogeneity-aware WM SE: scale bootstrap SE by sqrt(Q/df) ----
cat(sprintf("\nHeterogeneity-inflated WM: se=%.4f*sqrt(%.1f/%d)=%.4f -> z=%.2f p=%.3f\n",
    wm_se, Q, df, wm_se*sqrt(phi), pt/(wm_se*sqrt(phi)), 2*pnorm(-abs(pt/(wm_se*sqrt(phi))))))

## ---- drop the 2 positive adiposity outliers (MC4R rs6567160, ALDH2 rs2074356): what remains? ----
keep <- !(h$SNP %in% c("rs6567160","rs2074356"))
bx2<-bx[keep]; by2<-by[keep]; sx2<-sx[keep]; sy2<-sy[keep]; r2<-by2/bx2; w2<-1/(abs(sy2/bx2))^2
cat(sprintf("\nDrop MC4R+ALDH2 (n=%d, all beta-cell): IVW-fixed b=%+.4f se=%.4f p=%.2e ; Q=%.1f\n",
    length(bx2), sum(w2*r2)/sum(w2), sqrt(1/sum(w2)),
    2*pnorm(-abs((sum(w2*r2)/sum(w2))/sqrt(1/sum(w2)))),
    sum(w2*(r2-sum(w2*r2)/sum(w2))^2)))
