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
# EAS clumping (1000G-EAS merged panel). Usage (PowerShell): Rscript 04_clump_eas.R BMI CAD T2D eGFR
suppressMessages({library(ieugwasr); library(dplyr)})
BASE <- P4_BASE
PLINK <- file.path(CKM_ROOT, "lhcmr_eas_refs", "tools", "plink.exe")
BFILE <- file.path(BASE,"data/g1000_eas_merged")
args <- commandArgs(trailingOnly=TRUE); if(!length(args)) args <- c("BMI","CAD","T2D","eGFR")
for (node in args) {
  f <- file.path(BASE,"data/instruments",paste0(node,".eas.sig.tsv"))
  if (!file.exists(f)) { cat(node,": no sig file\n"); next }
  d <- read.delim(f, stringsAsFactors=FALSE)
  d2 <- data.frame(rsid=d$SNP, pval=as.numeric(d$pval), id=node, stringsAsFactors=FALSE)
  d2 <- d2[is.finite(d2$pval),]
  cl <- tryCatch(ieugwasr::ld_clump(d2, plink_bin=PLINK, bfile=BFILE, clump_r2=0.001, clump_kb=10000, clump_p=1),
                 error=function(e){cat(node,"CLUMP ERR:",conditionMessage(e),"\n"); NULL})
  if (is.null(cl)) next
  out <- d[d$SNP %in% cl$rsid, ]
  write.table(out, file.path(BASE,"data/instruments",paste0(node,".eas.clumped.tsv")), sep="\t", row.names=FALSE, quote=FALSE)
  cat(sprintf("[clump-EAS] %-6s %6d sig -> %4d instruments\n", node, nrow(d), nrow(out)))
}
