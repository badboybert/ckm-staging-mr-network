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
# Step 61c: Steiger direction margins for the CONTINUOUS->CONTINUOUS staging-graph edges.
#
# WHY THIS EXISTS. Step 61 sweeps assumed PREVALENCE, so by construction it only covers edges with
# at least one binary endpoint -- 34 of the staging graph's 54. Step 61b then wrote "Smallest margin
# in the whole staging graph", and the S20 caption, Methods and STROBE 12a inherited that scope
# claim. The 20 uncovered edges are all continuous->continuous, and they include the reciprocal
# lipid pairs (TC<->LDL, TG<->LDL, ...) where the two traits are near-collinear and the margin is
# therefore expected to sit close to 1 -- i.e. exactly the edges most likely to be fragile were the
# ones excluded from the fragility analysis.
#
# Prevalence cannot enter a continuous->continuous call, so these margins are prevalence-invariant:
# margin_lo == margin == margin_hi. That is a property worth recording, not a gap.
#
# Method is deliberately IDENTICAL to 10a_steiger.R: same format_data columns, same harmonise_data
# action=2, same mr_keep filter, same get_r_from_pn(p, n) * sign(b) per side, margin = r2x / r2y.
#
# Run: Rscript 61c_steiger_margin_continuous.R
# Out: results/steiger_margin_continuous.csv
suppressMessages({library(TwoSampleMR)})
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised")
RES   <- file.path(BASE, "results"); DAT <- file.path(BASE, "figures/data")

BIN <- c("CAD", "HF", "T2D", "Stroke", "CKD")
is_bin <- function(x) x %in% BIN

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors = FALSE), type = "exposure",
  snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele",
  other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval",
  samplesize_col = "N", phenotype_col = "Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors = FALSE), type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele",
  other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval",
  samplesize_col = "N", phenotype_col = "Phenotype")

r_cont <- function(h, side) {
  get_r_from_pn(p = h[[paste0("pval.", side)]], n = h[[paste0("samplesize.", side)]]) *
    sign(h[[paste0("beta.", side)]])
}

# The edge list is DERIVED from the shipped staging graph, never typed: every edge in
# fig1_edges.csv with two continuous endpoints, i.e. exactly those step 61 cannot reach.
g <- read.csv(file.path(DAT, "fig1_edges.csv"), stringsAsFactors = FALSE)
sel <- !is_bin(g$exp) & !is_bin(g$out)
cat(sprintf("staging graph: %d edges; continuous->continuous: %d\n", nrow(g), sum(sel)))

rows <- list()
for (i in which(sel)) {
  exp <- g$exp[i]; out <- g$out[i]; key <- paste0(exp, "__", out)
  ef <- file.path(INSTR, paste0(exp, ".clumped.tsv"))
  of <- file.path(HARM, paste0(key, ".outcome.tsv"))
  if (!file.exists(ef) || !file.exists(of)) {
    cat(sprintf("%-16s : MISSING FILES\n", key)); next
  }
  e <- tryCatch(fmt_exp(ef), error = function(z) NULL)
  o <- tryCatch(fmt_out(of), error = function(z) NULL)
  if (is.null(e) || is.null(o) || nrow(o) == 0) { cat(sprintf("%-16s : no data\n", key)); next }
  h <- harmonise_data(e, o, action = 2); h <- h[h$mr_keep, ]
  if (nrow(h) < 4) { cat(sprintf("%-16s : <4 SNPs\n", key)); next }
  rx <- r_cont(h, "exposure"); ry <- r_cont(h, "outcome")
  r2x <- sum(rx^2, na.rm = TRUE); r2y <- sum(ry^2, na.rm = TRUE)
  margin <- if (r2y > 0) r2x / r2y else Inf
  rows[[key]] <- data.frame(exposure = exp, outcome = out, cross = g$cross[i],
    direction = g$direction[i], nsnp = nrow(h), r2_exp = r2x, r2_out = r2y,
    margin = round(margin, 4), margin_lo = round(margin, 4), margin_hi = round(margin, 4),
    prevalence_invariant = TRUE, fragile = margin < 1.5, stringsAsFactors = FALSE)
  cat(sprintf("%-16s nsnp=%3d r2exp=%.5f r2out=%.5f margin=%6.2f%s\n",
      key, nrow(h), r2x, r2y, margin, ifelse(margin < 1.5, "   <- FRAGILE", "")))
}

res <- do.call(rbind, rows)
stopifnot(nrow(res) > 0)
write.csv(res, file.path(RES, "steiger_margin_continuous.csv"), row.names = FALSE)

cat(sprintf("\n%d of %d continuous->continuous edges computed\n", nrow(res), sum(sel)))
cat(sprintf("margins: min %.2f (%s->%s), max %.2f\n", min(res$margin),
    res$exposure[which.min(res$margin)], res$outcome[which.min(res$margin)], max(res$margin)))
cat(sprintf("fragile (margin < 1.5): %d\n", sum(res$fragile)))
if (any(res$fragile)) {
  f <- res[res$fragile, ]
  cat("  ", paste(sprintf("%s->%s %.2f (%s)", f$exposure, f$outcome, f$margin,
      ifelse(f$cross == "TRUE" | f$cross == TRUE, "CROSS-STAGE", "intra-stage")),
      collapse = "; "), "\n")
}
cat("\nwrote results/steiger_margin_continuous.csv\n")
