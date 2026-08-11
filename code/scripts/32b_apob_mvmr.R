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
# Step 32b (ITEM 1): fit the ApoB atherogenic-axis MVMR. CAD ~ ApoB + LDL + HDL + TG (+BMI+SBP).
# Direct effects (MVMR-IVW) + conditional F + horizontal-pleiotropy Q, contrasted with the
# univariable (marginal) total effect. KEY: does ApoB retain a direct CAD effect while HDL dissolves?
# Run: Rscript 32b_apob_mvmr.R
suppressMessages({library(MVMR); library(TwoSampleMR)})
BASE  <- P4_BASE
MVDIR <- file.path(BASE,"data/mvmr"); RES <- file.path(BASE,"results"); INSTR <- file.path(BASE,"data/instruments")

# --- univariable totals: LDL/HDL/TG/BMI/SBP->CAD from forward_local_edges.csv ---
uni <- read.csv(file.path(RES,"forward_local_edges.csv"), stringsAsFactors=FALSE)
uni_b <- function(e){ v<-uni[uni$exposure==e & uni$outcome=="CAD",]; if(nrow(v)) c(v$ivw_b[1],v$ivw_p[1]) else c(NA,NA) }

# --- univariable ApoB->CAD (own 214 instruments) ---
fmt <- function(f,ty) format_data(read.delim(f,stringsAsFactors=FALSE), type=ty, snp_col="SNP",
  beta_col="beta", se_col="se", effect_allele_col="effect_allele", other_allele_col="other_allele",
  eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
Eab <- fmt(file.path(INSTR,"ApoB.clumped.tsv"),"exposure")
Oab <- fmt(file.path(RES,"../data/harmonised_eas/ApoB__CAD_uni.outcome.tsv"),"outcome")
hab <- harmonise_data(Eab,Oab,action=2); hab <- hab[hab$mr_keep,]
mab <- mr(hab, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))
apob_uni_b <- mab$b[mab$method=="Inverse variance weighted"][1]
apob_uni_p <- mab$pval[mab$method=="Inverse variance weighted"][1]
cat(sprintf("[univariable ApoB->CAD] nSNP=%d IVW b=%+.4f p=%.2e | WM b=%+.4f | Egger b=%+.4f\n",
    nrow(hab), apob_uni_b, apob_uni_p, mab$b[mab$method=="Weighted median"][1], mab$b[mab$method=="MR Egger"][1]))
UNI <- function(e) if(e=="ApoB") c(apob_uni_b,apob_uni_p) else uni_b(e)

# sub-models derived by COLUMN-SUBSETTING the shared 124-SNP model_ApoB_CAD matrix (identical harmonised
# values, perfectly comparable, no re-extraction). Tests ApoB AS the atherogenic axis (no collinear LDL competitor).
D4 <- read.delim(file.path(MVDIR,"model_ApoB_CAD.tsv"), stringsAsFactors=FALSE)
submodel <- function(expos) {
  cols <- c("SNP", as.vector(rbind(paste0("bx_",expos), paste0("se_",expos))), "by","sey")
  D4[, cols]
}
SUB <- list(
  ApoB_HDL_TG = c("ApoB","HDL","TG"),   # ApoB represents atherogenic axis; does HDL dissolve?
  ApoB_LDL    = c("ApoB","LDL"),        # head-to-head (both collinear/weak, expected inconclusive)
  LDL_HDL_TG  = c("LDL","HDL","TG")     # reference without ApoB (the existing cascade story)
)

fit_matrix <- function(name, d) {
  bx_cols <- grep("^bx_", names(d), value=TRUE); se_cols <- grep("^se_", names(d), value=TRUE)
  expos <- sub("^bx_","",bx_cols)
  BXGs <- as.matrix(d[,bx_cols,drop=FALSE]); seBXGs <- as.matrix(d[,se_cols,drop=FALSE])
  r_in <- format_mvmr(BXGs=BXGs, BYG=d$by, seBXGs=seBXGs, seBYG=d$sey, RSID=d$SNP)
  est <- ivw_mvmr(r_in)
  Fst <- tryCatch(strength_mvmr(r_in, gencov=0), error=function(z) NULL)
  Qpl <- tryCatch(pleiotropy_mvmr(r_in, gencov=0), error=function(z) NULL)
  cat(sprintf("\n=== MODEL %s : CAD ~ %s   (nSNP=%d) ===\n", name, paste(expos,collapse=" + "), nrow(d)))
  rr <- list()
  for (i in seq_along(expos)) {
    e <- expos[i]; dir_b<-est[i,1]; dir_se<-est[i,2]; dir_p<-est[i,4]
    tot <- UNI(e); cf <- if(!is.null(Fst)) as.numeric(Fst[1,i]) else NA
    cat(sprintf("  %-5s direct b=%+.4f se=%.4f p=%.2e | univariable total b=%+.4f p=%.2e | condF=%s\n",
        e, dir_b, dir_se, dir_p, tot[1], tot[2], ifelse(is.na(cf),"NA",sprintf("%.1f",cf))))
    rr[[e]] <- data.frame(model=name, outcome="CAD", exposure=e, nsnp=nrow(d),
      direct_b=dir_b, direct_se=dir_se, direct_p=dir_p, total_b=tot[1], total_p=tot[2], cond_F=cf,
      stringsAsFactors=FALSE)
  }
  if(!is.null(Qpl)) cat(sprintf("  [horizontal pleiotropy] Q=%.1f  Qp=%.2e\n",
      as.numeric(Qpl$Qstat), as.numeric(Qpl$Qpval)))
  do.call(rbind, rr)
}

allrows <- list()
for (name in c("ApoB_CAD","ApoB_CAD_full")) {
  f <- file.path(MVDIR, paste0("model_",name,".tsv"))
  if (!file.exists(f)) { cat("MISSING", f, "\n"); next }
  allrows[[name]] <- fit_matrix(name, read.delim(f, stringsAsFactors=FALSE))
}
for (nm in names(SUB)) allrows[[nm]] <- fit_matrix(nm, submodel(SUB[[nm]]))

res <- do.call(rbind, allrows); write.csv(res, file.path(RES,"mvmr_apob.csv"), row.names=FALSE)
cat("\nWrote results/mvmr_apob.csv with", nrow(res), "rows\n")
