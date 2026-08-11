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
# Step 21: LD-clump TPMI instruments on the 1000G-EAS panel (preserves MarkerID column).
# Run (PowerShell): Rscript 21_tpmi_clump.R BMI SBP LDL HDL TC TG HbA1c
suppressMessages({library(ieugwasr); library(dplyr)})
BASE  <- P4_BASE
PLINK <- file.path(CKM_ROOT, "lhcmr_eas_refs", "tools", "plink.exe")
BFILE <- file.path(P4_BASE, "data", "g1000_eas_merged")
args <- commandArgs(trailingOnly=TRUE); if(!length(args)) args <- c("BMI","SBP","LDL","HDL","TC","TG","HbA1c")
for (node in args) {
  f <- file.path(BASE,"data/instruments",paste0(node,".tpmi.sig.tsv"))
  if (!file.exists(f)) { cat(node,": no sig file\n"); next }
  d <- read.delim(f, stringsAsFactors=FALSE); d <- d[is.finite(as.numeric(d$pval)),]
  d <- d[!duplicated(d$SNP),]
  cl <- tryCatch(ieugwasr::ld_clump(dplyr::tibble(rsid=d$SNP, pval=as.numeric(d$pval), id=node),
                   plink_bin=PLINK, bfile=BFILE, clump_r2=0.001, clump_kb=10000, clump_p=1),
                 error=function(e){cat(node,"CLUMP ERR:",conditionMessage(e),"\n"); NULL})
  if (is.null(cl)) next
  out <- d[d$SNP %in% cl$rsid, ]
  outf <- file.path(BASE,"data/instruments",paste0(node,".tpmi.clumped.tsv"))
  write.table(out, outf, sep="\t", row.names=FALSE, quote=FALSE)
  cat(sprintf("[tpmi clump] %-6s %5d sig -> %4d independent -> %s\n", node, nrow(d), nrow(out), basename(outf)))
}
