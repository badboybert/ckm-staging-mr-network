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
# Step 71: MR for unadjusted WHR -> {CAD, HF, T2D, HFpEF, HFrEF}, and the adiposity-distribution
# triangulation table (reviewer MAJOR-7 / Tier-2 #10).
#
# Puts three adiposity exposures side by side against each outcome:
#   BMI        = overall mass                        (from forward_local_edges.csv + step 69)
#   WHRadjBMI  = central distribution, BMI-residualised (shipped network_whradjbmi.csv + subtype MR)
#   WHR_unadj  = central distribution WITH mass       (this step)
# so the claim "mass, not distribution, drives HF" can be read from evidence rather than asserted
# from a single residualised null.
#
# Run: Rscript 71_whr_triangulation.R   (after 70 sig -> clump -> 70 out)

suppressMessages({library(TwoSampleMR)})
set.seed(20260718)
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised_whr")
RES   <- file.path(BASE, "results")

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

ef <- file.path(INSTR, "WHR_unadj.clumped.tsv")
stopifnot(file.exists(ef))
e <- fmt_exp(ef)
OUTS <- c("CAD", "HF", "T2D", "HFpEF", "HFrEF")
rows <- list()
for (o in OUTS) {
  of <- file.path(HARM, paste0("WHR_unadj__", o, ".outcome.tsv"))
  if (!file.exists(of)) { cat(sprintf("WHR_unadj -> %-6s : missing\n", o)); next }
  oo <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(oo) || nrow(oo) == 0) { cat(sprintf("WHR_unadj -> %-6s : no data\n", o)); next }
  h <- harmonise_data(e, oo, action=2); h <- h[h$mr_keep, ]
  if (nrow(h) < 3) { cat(sprintf("WHR_unadj -> %-6s : <3 SNPs\n", o)); next }
  m  <- mr(h, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  pl <- tryCatch(mr_pleiotropy_test(h), error=function(z) data.frame(egger_intercept=NA, pval=NA))
  gv <- function(meth, col) { v <- m[m$method == meth, col]; if (length(v)) v[1] else NA }
  rows[[o]] <- data.frame(exposure="WHR_unadj", outcome=o,
    nsnp=gv("Inverse variance weighted", "nsnp"),
    ivw_b=gv("Inverse variance weighted", "b"), ivw_se=gv("Inverse variance weighted", "se"),
    ivw_p=gv("Inverse variance weighted", "pval"),
    egger_b=gv("MR Egger", "b"), egger_intercept_p=pl$pval[1], wm_b=gv("Weighted median", "b"),
    stringsAsFactors=FALSE)
  cat(sprintf("WHR_unadj -> %-6s : nSNP=%3d  IVW b=%+.4f p=%.2e\n",
      o, rows[[o]]$nsnp, rows[[o]]$ivw_b, rows[[o]]$ivw_p))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES, "network_whr_unadj.csv"), row.names=FALSE)

# ------------------------------------------------------------------ triangulation table
get1 <- function(csv, exp, out, bcol="ivw_b", pcol="ivw_p") {
  if (!file.exists(csv)) return(c(NA, NA))
  d <- read.csv(csv, stringsAsFactors=FALSE)
  r <- d[d$exposure == exp & d$outcome == out, ]
  if (!nrow(r)) return(c(NA, NA))
  c(as.numeric(r[[bcol]][1]), as.numeric(r[[pcol]][1]))
}
FLE  <- file.path(RES, "forward_local_edges.csv")     # BMI-> and WHRadjBMI is NOT here
WADJ <- file.path(RES, "network_whradjbmi.csv")        # shipped WHRadjBMI edges
WADJS<- file.path(RES, "network_hfsubtypes.csv")       # WHRadjBMI->niHF* not here; subtype = BMI/CAD/T2D
ALLC <- file.path(RES, "network_hfsubtypes_allcause.csv")

sink(file.path(RES, "whr_triangulation.txt"))
cat("=== Adiposity distribution vs mass: WHR triangulation ===\n\n")
cat("Three adiposity exposures against each outcome (IVW b [P]):\n")
cat("  BMI       = overall mass\n")
cat("  WHRadjBMI = central distribution, BMI-residualised (can be collider-biased; reviewer MAJOR-7)\n")
cat("  WHR_unadj = central distribution retaining the mass component (this step)\n\n")
cat(sprintf("%-8s | %-22s | %-22s | %-22s\n", "outcome", "BMI", "WHRadjBMI", "WHR_unadj"))
cat(strrep("-", 82), "\n")
fmt <- function(v) if (any(is.na(v))) "        -            " else sprintf("%+7.4f (P=%.1e)", v[1], v[2])
for (o in c("CAD", "HF", "T2D", "HFpEF", "HFrEF")) {
  bmi <- get1(FLE, "BMI", o)
  if (o %in% c("HFpEF","HFrEF")) bmi <- get1(ALLC, "BMI", o)
  wadj <- get1(WADJ, "WHRadjBMI", o)
  wun  <- get1(file.path(RES,"network_whr_unadj.csv"), "WHR_unadj", o)
  cat(sprintf("%-8s | %-22s | %-22s | %-22s\n", o, fmt(bmi), fmt(wadj), fmt(wun)))
}
cat("\nREAD (interpret against the reviewer's collider caveat, do not over-read):\n")
cat("- If WHR_unadj reaches HF/HFpEF/HFrEF while WHRadjBMI does not, the HF signal is carried by the\n")
cat("  overall-mass component that BMI-adjustment removes - consistent with 'mass, not distribution',\n")
cat("  but WHR_unadj and BMI share that mass component, so this triangulates rather than proves.\n")
cat("- If WHR_unadj is ALSO null on HF, then central adiposity (with or without mass) does not reach\n")
cat("  HF on this evidence, and the 'spares' language must be dropped entirely.\n")
cat("- WHR_unadj includes UK Biobank, overlapping UKB-containing outcome GWAS; report the direction,\n")
cat("  not the magnitude, and note the overlap.\n")
sink()
cat("\nWrote results/network_whr_unadj.csv and whr_triangulation.txt\n")
