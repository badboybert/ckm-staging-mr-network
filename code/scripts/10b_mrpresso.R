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
# Step 10b: MR-PRESSO (global pleiotropy + outlier-corrected estimate) on headline edges.
# Background job. Run (PowerShell): Rscript 10b_mrpresso.R
#
# AUTHORITATIVE WRITER of results/mrpresso_headline.csv (settled 2026-07-23).
#   scripts/10_harden.R used to write the same file with different settings and no instrument cap;
#   it is now deprecated and refuses to run. Only this script may write mrpresso_headline.csv.
#   Provenance of the committed file, from its own contents:
#     - nb_dist = 1000 in every row  -> matches the NbDistribution = 1000 hardcoded below.
#       (10_harden.R used NbDistribution = min(20000, max(2000, ceiling(nSNP/0.05))), i.e. >= 2000,
#       so it cannot have produced a file whose nb_dist column reads 1000.)
#     - 5 rows only, exactly the edges of EDGES below with <= 130 harmonised instruments
#       (HF->CAD 45, HF->T2D 33, T2D->HF 115, T2D->CAD 116, CAD->T2D 128); 10_harden.R has no
#       >130 skip and would have written all 16 edges.
# PARAMETERS OF RECORD for results/mrpresso_headline.csv:
#   NbDistribution = 1000; SignifThreshold = 0.05; seed = 42; OUTLIER + DISTORTION tests on.
#   Instrument-count window actually analysed: >= 6 and <= 130 harmonised SNPs (edges outside the
#   window are skipped with a printed reason). Realised range in the committed file: 33-128 SNPs.
suppressMessages({library(TwoSampleMR); library(MRPRESSO)})
set.seed(42)
BASE  <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised")
RES   <- file.path(BASE,"results")

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

EDGES <- list(
  c("LDL","CAD"), c("CAD","HF"), c("HF","CAD"), c("BMI","CAD"), c("BMI","T2D"),
  c("T2D","CAD"), c("T2D","HF"), c("SBP","Stroke"), c("CAD","T2D"), c("HF","T2D"),
  c("HbA1c","T2D"), c("LDL","HF"), c("BMI","HF"), c("BMI","Stroke"), c("SBP","CAD"), c("SBP","HF"))

rows <- list()
for (ed in EDGES) {
  exp <- ed[1]; out <- ed[2]; key <- paste0(exp,"__",out)
  ef <- file.path(INSTR, paste0(exp,".clumped.tsv")); of <- file.path(HARM, paste0(key,".outcome.tsv"))
  if (!file.exists(ef) || !file.exists(of)) next
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL); o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e)||is.null(o)||nrow(o)==0) next
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep,]
  if (nrow(h)<6) { cat(sprintf("%-14s : <6 SNPs, skip PRESSO\n",key)); next }
  # MR-PRESSO (even the global test) is O(NbDist * nSNP^2); it is intractable for the
  # large forward edges (LDL/BMI/SBP with 300-500 SNPs), which already have Egger + WM +
  # Cochran-Q robustness. Run PRESSO only on the modest (<=130-SNP) edges — which are
  # exactly the load-bearing backward/asymmetry edges where outlier-correction matters.
  if (nrow(h) > 130) { cat(sprintf("%-14s : %d SNPs >130, skip PRESSO (Egger/WM/Q cover it)\n",key,nrow(h))); next }
  dfp <- data.frame(by=h$beta.outcome, bx=h$beta.exposure, sy=h$se.outcome, sx=h$se.exposure)
  do_outlier <- TRUE
  nb  <- 1000
  pr  <- tryCatch(mr_presso(BetaOutcome="by", BetaExposure="bx", SdOutcome="sy", SdExposure="sx",
            data=dfp, OUTLIERtest=do_outlier, DISTORTIONtest=do_outlier, NbDistribution=nb, SignifThreshold=0.05, seed=42),
            error=function(z){cat(sprintf("%-14s PRESSO fail: %s\n",key,conditionMessage(z))); NULL})
  if (is.null(pr)) next
  mm <- pr$`Main MR results`
  gv <- function(k,c) { v <- mm[mm$`MR Analysis`==k, c]; if(length(v)) v[1] else NA }
  gt <- pr$`MR-PRESSO results`$`Global Test`
  outl <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  n_out <- if (is.null(outl) || (is.character(outl) && length(outl)==1)) 0 else length(outl)
  dist_p <- pr$`MR-PRESSO results`$`Distortion Test`$Pvalue
  rows[[key]] <- data.frame(exposure=exp, outcome=out, nsnp=nrow(h), nb_dist=nb,
    raw_b=gv("Raw","Causal Estimate"), raw_p=gv("Raw","P-value"),
    corrected_b=gv("Outlier-corrected","Causal Estimate"), corrected_p=gv("Outlier-corrected","P-value"),
    global_RSSobs=if(!is.null(gt$RSSobs)) gt$RSSobs else NA,
    global_p=if(!is.null(gt$Pvalue)) gt$Pvalue else NA,
    n_outliers=n_out, distortion_p=ifelse(is.null(dist_p),NA,dist_p), stringsAsFactors=FALSE)
  cat(sprintf("%-14s raw=%+.4f(%.1e) corr=%s(%s) glob_p=%s outliers=%d dist_p=%s\n",
      key, rows[[key]]$raw_b, rows[[key]]$raw_p,
      ifelse(is.na(rows[[key]]$corrected_b),"NA",sprintf("%+.4f",rows[[key]]$corrected_b)),
      ifelse(is.na(rows[[key]]$corrected_p),"NA",sprintf("%.1e",rows[[key]]$corrected_p)),
      ifelse(is.na(rows[[key]]$global_p),"NA",format(rows[[key]]$global_p,digits=2)),
      n_out, ifelse(is.na(rows[[key]]$distortion_p),"NA",format(rows[[key]]$distortion_p,digits=2))))
  flush.console()
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES,"mrpresso_headline.csv"), row.names=FALSE)
cat("\nWrote mrpresso_headline.csv with", nrow(res), "edges\n")
