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
# Component A (MVMR): cascade-vs-common-driver via multivariable IVW + conditional F.
# Consumes data/mvmr/model_*.tsv (harmonised by 13a). Run: Rscript 13b_mvmr.R
suppressMessages({library(MVMR)})
BASE  <- P4_BASE
MVDIR <- file.path(BASE,"data/mvmr"); RES <- file.path(BASE,"results")

# univariable IVW (for the cascade contrast: total vs direct effect)
uni <- read.csv(file.path(RES,"forward_local_edges.csv"), stringsAsFactors=FALSE)
uni_b <- function(e,o){ v<-uni[uni$exposure==e & uni$outcome==o,]; if(nrow(v)) c(v$ivw_b[1],v$ivw_se[1],v$ivw_p[1]) else c(NA,NA,NA) }

files <- list.files(MVDIR, pattern="^model_.*\\.tsv$", full.names=TRUE)
allrows <- list()
for (f in files) {
  name <- sub("^model_","",sub("\\.tsv$","",basename(f)))
  d <- read.delim(f, stringsAsFactors=FALSE)
  bx_cols <- grep("^bx_", names(d), value=TRUE)
  se_cols <- grep("^se_", names(d), value=TRUE)
  expos   <- sub("^bx_","",bx_cols)
  outc    <- name  # outcome encoded in model name prefix; recover from mapping
  outmap  <- c(CAD_full="CAD", CAD_BMI_T2D="CAD", HF_via_CAD="HF", Stroke_full="Stroke")
  outc    <- outmap[[name]]
  BXGs  <- as.matrix(d[,bx_cols,drop=FALSE]); seBXGs <- as.matrix(d[,se_cols,drop=FALSE])
  BYG   <- d$by; seBYG <- d$sey
  r_in  <- format_mvmr(BXGs=BXGs, BYG=BYG, seBXGs=seBXGs, seBYG=seBYG, RSID=d$SNP)
  est   <- ivw_mvmr(r_in)                       # direct effects
  Fst   <- tryCatch(strength_mvmr(r_in, gencov=0), error=function(z) NULL)   # conditional F
  Qpl   <- tryCatch(pleiotropy_mvmr(r_in, gencov=0), error=function(z) NULL) # horizontal pleiotropy Q
  cat(sprintf("\n=== MODEL %s : %s ~ %s   (nSNP=%d) ===\n", name, outc, paste(expos,collapse=" + "), nrow(d)))
  for (i in seq_along(expos)) {
    e <- expos[i]
    dir_b <- est[i,1]; dir_se <- est[i,2]; dir_p <- est[i,4]
    tot   <- uni_b(e,outc)
    cf    <- if(!is.null(Fst)) as.numeric(Fst[1,i]) else NA
    cat(sprintf("  %-6s direct b=%+.4f se=%.4f p=%.2e | univariable total b=%+.4f p=%.2e | condF=%s\n",
        e, dir_b, dir_se, dir_p, tot[1], tot[3], ifelse(is.na(cf),"NA",sprintf("%.1f",cf))))
    allrows[[paste(name,e)]] <- data.frame(model=name, outcome=outc, exposure=e, nsnp=nrow(d),
      direct_b=dir_b, direct_se=dir_se, direct_p=dir_p,
      total_b=tot[1], total_p=tot[3], cond_F=cf, stringsAsFactors=FALSE)
  }
  if(!is.null(Qpl)) cat(sprintf("  [horizontal pleiotropy] Q=%.1f  Qp=%.2e\n",
      as.numeric(Qpl$Qstat), as.numeric(Qpl$Qpval)))
}
res <- do.call(rbind, allrows)
write.csv(res, file.path(RES,"mvmr_cascade.csv"), row.names=FALSE)
cat("\nWrote mvmr_cascade.csv with", nrow(res), "exposure-rows\n")
