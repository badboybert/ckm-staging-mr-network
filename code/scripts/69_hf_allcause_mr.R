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
# Step 69: MR for the ALL-AETIOLOGY HFpEF/HFrEF edges (Enzan 2025), the round-2 BLOCKER-2 replacement
# for the non-ischemic specificity analysis. Same estimators as 06_mr.R / 44_subtype_mr.R.
#
# Writes a SEPARATE file (network_hfsubtypes_allcause.csv) so the shipped non-ischemic
# network_hfsubtypes.csv is untouched. Reads only *__HFpEF/__HFrEF outcome files (step 68).
#
# The load-bearing question the reviewer raised: in the NON-ischemic analysis, CAD and T2D "not
# reaching" HF subtypes is partly built into the case definition (ischemia excluded). Here CAD is NOT
# excluded, so this re-tests the specificity claim honestly. Run: Rscript 69_hf_allcause_mr.R

suppressMessages({library(TwoSampleMR)})
set.seed(20260718)
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised_hfsub")
RES   <- file.path(BASE, "results")

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

ofiles <- list.files(HARM, pattern="__(HFpEF|HFrEF)\\.outcome\\.tsv$", full.names=TRUE)
rows <- list()
for (of in ofiles) {
  base <- sub("\\.outcome\\.tsv$", "", basename(of)); parts <- strsplit(base, "__")[[1]]
  exp <- parts[1]; out <- parts[2]
  ef <- file.path(INSTR, paste0(exp, ".clumped.tsv")); if (!file.exists(ef)) next
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL); o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e) || is.null(o) || nrow(o) == 0) next
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep, ]
  if (nrow(h) < 3) { cat(sprintf("%-6s -> %-6s : <3 SNPs\n", exp, out)); next }
  m  <- mr(h, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  pl <- tryCatch(mr_pleiotropy_test(h), error=function(z) data.frame(egger_intercept=NA, pval=NA))
  het <- tryCatch(mr_heterogeneity(h, method_list="mr_ivw"), error=function(z) data.frame(Q=NA, Q_pval=NA))
  gv <- function(mm, meth, col) { v <- mm[mm$method == meth, col]; if (length(v)) v[1] else NA }
  rows[[base]] <- data.frame(
    exposure=exp, outcome=out, nsnp=gv(m, "Inverse variance weighted", "nsnp"),
    ivw_b=gv(m, "Inverse variance weighted", "b"), ivw_se=gv(m, "Inverse variance weighted", "se"),
    ivw_p=gv(m, "Inverse variance weighted", "pval"),
    egger_b=gv(m, "MR Egger", "b"), egger_intercept=pl$egger_intercept[1], egger_intercept_p=pl$pval[1],
    wm_b=gv(m, "Weighted median", "b"), Q=het$Q[1], Q_p=het$Q_pval[1], stringsAsFactors=FALSE)
  cat(sprintf("%-6s -> %-6s : nSNP=%3d  IVW b=%+.4f p=%.2e\n",
      exp, out, rows[[base]]$nsnp, rows[[base]]$ivw_b, rows[[base]]$ivw_p))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES, "network_hfsubtypes_allcause.csv"), row.names=FALSE)

# ------------------------------------------------------------------ specificity contrast + summary
sink(file.path(RES, "hf_allcause_subtypes.txt"))
cat("=== All-aetiology HFpEF / HFrEF MR (Enzan 2025, PMID 41184235) ===\n")
cat("HFpEF = GCST90654629 (N=483,263; 26,743 cases: 19,589 EUR + 7,154 EAS);\n")
cat("HFrEF = GCST90654628 (N=480,269; 23,749 cases: 19,495 EUR + 4,254 EAS). CROSS-ANCESTRY, build37.\n")
cat("This REPLACES the non-ischemic specificity analysis (Henry 2025, 3,590/4,975 cases) that\n")
cat("reviewer BLOCKER-2 called circular: CAD is NOT excluded from these all-aetiology outcomes, so a\n")
cat("null CAD->HFpEF edge is evidence, not case-definition.\n")
cat("PROVENANCE: EUR+EAS meta-analyses (NOT EUR-only). EUR instruments on a ~42%-EAS outcome attenuate\n")
cat("rather than inflate, so flips are conservative; EAS controls overlap the paper's BBJ HF outcome\n")
cat("GCST90668009. Report DIRECTIONS, not magnitudes.\n\n")
BONF <- 0.05 / nrow(res)
res$sig <- res$ivw_p < BONF
cat(sprintf("Within-family Bonferroni = 0.05/%d = %.3e\n\n", nrow(res), BONF))
cat(sprintf("%-8s %-7s %5s %10s %10s %10s  %s\n","exposure","subtype","nSNP","IVW b","IVW se","IVW P","sig"))
for (sub in c("HFpEF","HFrEF")) {
  s <- res[res$outcome==sub, ]; s <- s[order(-abs(s$ivw_b)), ]
  for (i in seq_len(nrow(s)))
    cat(sprintf("%-8s %-7s %5d %+10.4f %10.4f %10.2e  %s\n",
        s$exposure[i], sub, s$nsnp[i], s$ivw_b[i], s$ivw_se[i], s$ivw_p[i],
        ifelse(s$sig[i],"*","")))
  cat("\n")
}
w <- function(e,o){ r<-res[res$exposure==e & res$outcome==o,]; if(nrow(r)) r else NULL }
cat("SPECIFICITY CONTRAST (the claim under review):\n")
for (sub in c("HFpEF","HFrEF")) {
  for (e in c("BMI","T2D","CAD","SBP")) {
    r <- w(e,sub); if (is.null(r)) next
    cat(sprintf("  %-4s -> %-6s : b=%+.4f  P=%.2e  %s\n", e, sub, r$ivw_b, r$ivw_p,
        ifelse(r$ivw_p<BONF, "REACHES", "not significant")))
  }
}
cat("\nREAD: with CAD included in the case definition, whether CAD and T2D reach the HF subtypes is\n")
cat("now an empirical result rather than a definitional one. Compare each edge's direction and\n")
cat("significance with the non-ischemic analysis (network_hfsubtypes.csv) and with all-cause HF\n")
cat("(forward_local_edges.csv) to read adiposity's obesity-cardiomyopathy route honestly.\n")
sink()
cat("\nWrote results/network_hfsubtypes_allcause.csv (", nrow(res), "edges) and hf_allcause_subtypes.txt\n")
