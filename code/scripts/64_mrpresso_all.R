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
# Step 64: MR-PRESSO across EVERY Bonferroni-significant multi-instrument edge (reviewer MODERATE-16).
#
# The shipped analysis runs MR-PRESSO on 5 edges (10b_mrpresso.R caps at <=130 instruments because
# the outlier/distortion tests are O(NbDistribution * nSNP^2) and are intractable for the 300-500-SNP
# forward edges). The reviewer's point is that selective deep analysis is defensible but must be
# described precisely, and asks for the full sweep "where computationally feasible".
#
# This runs the sweep honestly, in two tiers, and records which tier each edge got:
#   nSNP <= 130 : full MR-PRESSO (global + outlier + distortion tests)   [as shipped]
#   nSNP >  130 : GLOBAL TEST ONLY (OUTLIERtest=FALSE), which is the tractable part
# so that every significant edge has a global pleiotropy test and the two tiers are labelled rather
# than blurred.
#
# RESUMABLE BY DESIGN: each edge is appended to the CSV as soon as it finishes and a re-run skips
# edges already present, so a kill costs at most one edge.
#
# Run: Rscript 64_mrpresso_all.R
# Out: results/mrpresso_all_edges.csv (incremental) and .txt (summary, written at the end)

suppressMessages({library(TwoSampleMR); library(MRPRESSO)})
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised")
RES   <- file.path(BASE, "results")
OUTC  <- file.path(RES, "mrpresso_all_edges.csv")

NB <- 1000
FULL_MAX <- 130          # above this, global test only

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

net  <- read.csv(file.path(RES, "forward_local_edges.csv"), stringsAsFactors=FALSE)
BONF <- 0.05 / nrow(net)
sig  <- net[net$ivw_p < BONF, c("exposure", "outcome", "nsnp", "ivw_b", "ivw_p")]
sig  <- sig[order(sig$nsnp), ]
cat(sprintf("Bonferroni-significant edges to sweep: %d (threshold %.3e)\n", nrow(sig), BONF))

done <- character(0)
if (file.exists(OUTC)) {
  prev <- read.csv(OUTC, stringsAsFactors=FALSE)
  done <- paste0(prev$exposure, "__", prev$outcome)
  cat(sprintf("resuming: %d edge(s) already complete\n", length(done)))
}

hdr <- c("exposure","outcome","tier","nsnp_harmonised","raw_b","raw_p","corrected_b","corrected_p",
         "global_rssobs","global_p","n_outliers","distortion_p")
if (!file.exists(OUTC)) {
  write.table(t(hdr), OUTC, sep=",", row.names=FALSE, col.names=FALSE, qmethod="double")
}

for (i in seq_len(nrow(sig))) {
  exp <- sig$exposure[i]; out <- sig$outcome[i]; key <- paste0(exp, "__", out)
  if (key %in% done) next
  ef <- file.path(INSTR, paste0(exp, ".clumped.tsv")); of <- file.path(HARM, paste0(key, ".outcome.tsv"))
  if (!file.exists(ef) || !file.exists(of)) { cat(sprintf("%-16s : missing files\n", key)); next }
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL); o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e) || is.null(o) || nrow(o) == 0) { cat(sprintf("%-16s : no data\n", key)); next }
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep, ]
  if (nrow(h) < 6) { cat(sprintf("%-16s : <6 SNPs, skip\n", key)); next }

  full <- nrow(h) <= FULL_MAX
  dfp  <- data.frame(by=h$beta.outcome, bx=h$beta.exposure, sy=h$se.outcome, sx=h$se.exposure)
  t0 <- Sys.time()
  pr <- tryCatch(mr_presso(BetaOutcome="by", BetaExposure="bx", SdOutcome="sy", SdExposure="sx",
          data=dfp, OUTLIERtest=full, DISTORTIONtest=full, NbDistribution=NB,
          SignifThreshold=0.05, seed=42),
          error=function(z){cat(sprintf("%-16s PRESSO fail: %s\n", key, conditionMessage(z))); NULL})
  if (is.null(pr)) next
  mm <- pr$`Main MR results`
  gv <- function(k, c) { v <- mm[mm$`MR Analysis` == k, c]; if (length(v)) v[1] else NA }
  gt <- pr$`MR-PRESSO results`$`Global Test`
  outl <- if (full) pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices` else NULL
  dp   <- if (full) pr$`MR-PRESSO results`$`Distortion Test`$Pvalue else NA

  row <- data.frame(
    exposure=exp, outcome=out, tier=if (full) "full" else "global_only",
    nsnp_harmonised=nrow(h),
    raw_b=gv("Raw","Causal Estimate"), raw_p=gv("Raw","P-value"),
    corrected_b=gv("Outlier-corrected","Causal Estimate"),
    corrected_p=gv("Outlier-corrected","P-value"),
    global_rssobs=if (!is.null(gt)) gt$RSSobs else NA,
    global_p=if (!is.null(gt)) gt$Pvalue else NA,
    n_outliers=if (is.null(outl)) NA_integer_ else length(outl),
    distortion_p=if (is.null(dp)) NA_real_ else dp,
    stringsAsFactors=FALSE)
  write.table(row, OUTC, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE, qmethod="double")
  cat(sprintf("%-16s [%-11s] n=%4d  global P=%-8s outliers=%-4s  (%.0fs)  %d/%d\n",
      key, row$tier, nrow(h),
      if (is.na(row$global_p)) "NA" else format(row$global_p, digits=3),
      if (is.na(row$n_outliers)) "-" else as.character(row$n_outliers),
      as.numeric(difftime(Sys.time(), t0, units="secs")), i, nrow(sig)))
  flush.console()
}

# ---------------------------------------------------------------- summary
res <- read.csv(OUTC, stringsAsFactors=FALSE)
# MR-PRESSO returns the LITERAL STRING "<0.001" for a global P below 1/NbDistribution, so the
# column comes back as character and `res$global_p < 0.05` would silently do a STRING comparison.
# Parse it to a numeric bound before any thresholding.
.numP <- function(v) {
  v <- trimws(as.character(v))
  out <- suppressWarnings(as.numeric(v))
  cens <- !is.na(v) & grepl("^<", v)
  out[cens] <- suppressWarnings(as.numeric(sub("^<", "", v[cens]))) / 2   # midpoint of (0, bound)
  out
}
res$global_p_num     <- .numP(res$global_p)
res$distortion_p_num <- .numP(res$distortion_p)
res$global_p_censored <- grepl("^<", trimws(as.character(res$global_p)))
sink(file.path(RES, "mrpresso_all_edges.txt"))
cat("=== MR-PRESSO across every Bonferroni-significant edge ===\n\n")
cat(sprintf("Network Bonferroni = 0.05/%d = %.3e; %d edges significant.\n", nrow(net), BONF, nrow(sig)))
cat(sprintf("Swept: %d edges (%d full outlier+distortion, %d global test only).\n",
            nrow(res), sum(res$tier == "full"), sum(res$tier == "global_only")))
cat(sprintf("Tiering rule: full MR-PRESSO for <=%d harmonised instruments; above that the outlier and\n", FULL_MAX))
cat("distortion tests are O(NbDistribution * nSNP^2) and only the global test is run.\n\n")
gsig <- res[!is.na(res$global_p_num) & res$global_p_num < 0.05, ]
cat(sprintf("Global pleiotropy test P < 0.05 in %d/%d edges.\n", nrow(gsig), sum(!is.na(res$global_p))))
cat("(A significant global test indicates heterogeneity/outliers, not necessarily a biased estimate.)\n\n")
full <- res[res$tier == "full" & !is.na(res$n_outliers), ]
if (nrow(full)) {
  cat(sprintf("Among the %d fully-tested edges: %d have >=1 outlier; distortion P<0.05 in %d.\n",
              nrow(full), sum(full$n_outliers > 0), sum(!is.na(full$distortion_p_num) & full$distortion_p_num < 0.05)))
  dist <- full[!is.na(full$distortion_p_num) & full$distortion_p_num < 0.05, ]
  if (nrow(dist)) {
    cat("\nEdges whose estimate is DISTORTED by outliers (distortion P<0.05):\n")
    for (i in seq_len(nrow(dist)))
      cat(sprintf("  %s->%s  raw %+.4f -> corrected %+.4f  (%d outliers, distortion P=%s)\n",
          dist$exposure[i], dist$outcome[i], dist$raw_b[i], dist$corrected_b[i],
          dist$n_outliers[i], as.character(dist$distortion_p[i])))
  }
}
cat("\nFULL TABLE (sorted by global-test P):\n")
o <- res[order(res$global_p_num), ]
cat(sprintf("%-16s %-12s %6s %12s %12s %9s %9s\n","edge","tier","nSNP","raw b","corrected b","globalP","distortP"))
for (i in seq_len(nrow(o))) {
  cat(sprintf("%-16s %-12s %6d %12s %12s %9s %9s\n",
      paste0(o$exposure[i],"->",o$outcome[i]), o$tier[i], o$nsnp_harmonised[i],
      formatC(o$raw_b[i], format="f", digits=4),
      if (is.na(o$corrected_b[i])) "-" else formatC(o$corrected_b[i], format="f", digits=4),
      if (is.na(o$global_p_num[i])) "-" else as.character(o$global_p[i]),
      if (is.na(o$distortion_p_num[i])) "-" else as.character(o$distortion_p[i])))
}
sink()
cat(sprintf("\nwrote results/mrpresso_all_edges.{csv,txt} (%d edges)\n", nrow(res)))
