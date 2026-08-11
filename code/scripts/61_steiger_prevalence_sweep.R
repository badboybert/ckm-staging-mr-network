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
# Step 61: liability-scale Steiger sensitivity across plausible prevalences, for EVERY binary-
# involving edge that enters the staging graph.
#
# Why: Steiger direction is a HARD graph-entry gate, so the headline staging count (25 of 27
# cross-stage edges) is conditional on it. The shipped analysis gates the network on the
# QUANTITATIVE approximation and recomputes the liability scale on only 18 load-bearing edges at a
# SINGLE fixed prevalence each. get_r_from_lor is prevalence-dependent, so a direction call can in
# principle flip across clinically plausible prevalences. This script tests whether any actually does.
#
# Design: harmonise each edge once, compute the per-SNP r for both sides, then sweep the prevalence
# of whichever side(s) are binary over a clinically plausible grid, recording the Steiger direction
# (sum r^2_exposure > sum r^2_outcome) at every grid point.
#
# Run: Rscript 61_steiger_prevalence_sweep.R
# Out: results/steiger_prevalence_sweep.csv  (one row per edge x prevalence combination)
#      results/steiger_prevalence_sweep.txt  (summary: which edges are direction-stable)

suppressMessages({library(TwoSampleMR)})
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised")
RES   <- file.path(BASE, "results"); FIG <- file.path(BASE, "figures/data")

# Sample sizes as in 10a_steiger.R; prevalence becomes a swept parameter rather than a constant.
BIN <- list(
  CAD    = list(ncase=181522, ncontrol=1047847, base=0.06, grid=c(0.03, 0.045, 0.06, 0.08, 0.10)),
  HF     = list(ncase=139533, ncontrol=1054983, base=0.02, grid=c(0.01, 0.015, 0.02, 0.035, 0.05)),
  T2D    = list(ncase=74124,  ncontrol=824006,  base=0.10, grid=c(0.05, 0.075, 0.10, 0.125, 0.15)),
  Stroke = list(ncase=73652,  ncontrol=1234808, base=0.05, grid=c(0.02, 0.035, 0.05, 0.065, 0.08)),
  CKD    = list(ncase=41395,  ncontrol=439303,  base=0.09, grid=c(0.05, 0.07, 0.09, 0.12, 0.15)))
is_bin <- function(x) x %in% names(BIN)

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

# r for one side at a GIVEN prevalence (continuous traits ignore prev)
r_side_at <- function(h, side, trait, prev) {
  b <- h[[paste0("beta.", side)]]; p <- h[[paste0("pval.", side)]]
  n <- h[[paste0("samplesize.", side)]]; eaf <- h[[paste0("eaf.", side)]]
  if (is_bin(trait)) {
    m <- BIN[[trait]]; af <- eaf; af[is.na(af)] <- 0.5
    get_r_from_lor(lor=b, af=af, ncase=m$ncase, ncontrol=m$ncontrol, prevalence=prev) * sign(b)
  } else {
    get_r_from_pn(p=p, n=n) * sign(b)
  }
}

# every edge in the staging graph that involves at least one binary trait
edges <- read.csv(file.path(FIG, "fig1_edges.csv"), stringsAsFactors=FALSE)
edges <- edges[is_bin(edges$exp) | is_bin(edges$out), c("exp", "out", "cross", "direction")]
cat(sprintf("staging-graph edges involving a binary trait: %d\n", nrow(edges)))

rows <- list(); k <- 0
for (i in seq_len(nrow(edges))) {
  exp <- edges$exp[i]; out <- edges$out[i]; key <- paste0(exp, "__", out)
  ef <- file.path(INSTR, paste0(exp, ".clumped.tsv"))
  of <- file.path(HARM, paste0(key, ".outcome.tsv"))
  if (!file.exists(ef) || !file.exists(of)) { cat(sprintf("%-16s : missing files\n", key)); next }
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL)
  o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e) || is.null(o) || nrow(o) == 0) { cat(sprintf("%-16s : no data\n", key)); next }
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep, ]
  if (nrow(h) < 4) { cat(sprintf("%-16s : <4 SNPs\n", key)); next }

  ge <- if (is_bin(exp)) BIN[[exp]]$grid else NA
  go <- if (is_bin(out)) BIN[[out]]$grid else NA
  for (pe in ge) for (po in go) {
    rx <- r_side_at(h, "exposure", exp, pe)
    ry <- r_side_at(h, "outcome",  out, po)
    r2x <- sum(rx^2, na.rm=TRUE); r2y <- sum(ry^2, na.rm=TRUE)
    k <- k + 1
    rows[[k]] <- data.frame(
      exposure=exp, outcome=out, cross=edges$cross[i], direction=edges$direction[i],
      nsnp=nrow(h),
      prev_exp=ifelse(is.na(pe), NA_real_, pe), prev_out=ifelse(is.na(po), NA_real_, po),
      at_base=( (is.na(pe) || isTRUE(all.equal(pe, BIN[[exp]]$base))) &&
                (is.na(po) || isTRUE(all.equal(po, BIN[[out]]$base))) ),
      r2_exp=r2x, r2_out=r2y, correct_dir=(r2x > r2y), stringsAsFactors=FALSE)
  }
  # combinations just written for THIS edge = length(ge) * length(go), where a continuous side
  # contributes the single NA level. Counting only the non-NA levels gave 0 for one-binary-side
  # edges and printed a misleading "grid= 1"; the written rows were always correct.
  n_comb <- length(ge) * length(go)
  fl <- do.call(rbind, rows[(k - n_comb + 1):k])
  cat(sprintf("%-16s nSNP=%4d  grid=%2d  direction TRUE in %d/%d combinations%s\n",
      key, nrow(h), nrow(fl), sum(fl$correct_dir), nrow(fl),
      if (length(unique(fl$correct_dir)) > 1) "   <-- FLIPS" else ""))
}

res <- do.call(rbind, rows)
write.csv(res, file.path(RES, "steiger_prevalence_sweep.csv"), row.names=FALSE)

# ---- summary -------------------------------------------------------------------------------
sink(file.path(RES, "steiger_prevalence_sweep.txt"))
cat("=== Liability-scale Steiger sensitivity to assumed prevalence ===\n")
cat("Every staging-graph edge involving a binary trait, swept over a clinically plausible\n")
cat("prevalence grid for each binary side. Steiger direction = sum(r^2 exposure) > sum(r^2 outcome).\n\n")
cat("Prevalence grids swept:\n")
# Only print grids that were actually EXERCISED. CKD is defined in BIN but no staging-graph edge
# involves it (fig1_edges.csv has 12 nodes and no CKD), so advertising a CKD grid here invited a
# Methods sentence claiming the renal edges were covered by this sensitivity. They were not.
.swept <- unique(c(res$exposure, res$outcome))
for (nm in names(BIN)) cat(sprintf("  %-7s base %.3f   grid %s%s
", nm, BIN[[nm]]$base,
                                   paste(sprintf("%.3f", BIN[[nm]]$grid), collapse=", "),
                                   if (nm %in% .swept) "" else
                                     "   [NOT EXERCISED - no staging-graph edge involves this trait]"))
cat("\n")
agg <- aggregate(correct_dir ~ exposure + outcome + cross + direction, data=res,
                 FUN=function(v) c(n=length(v), t=sum(v)))
agg <- do.call(data.frame, agg)
names(agg)[5:6] <- c("n_comb", "n_true")
agg$stable <- agg$n_true == agg$n_comb | agg$n_true == 0
agg$flips  <- !agg$stable
cat(sprintf("edges swept: %d   |   direction-STABLE across the whole grid: %d   |   FLIPPING: %d\n\n",
            nrow(agg), sum(agg$stable), sum(agg$flips)))
cat(sprintf("%-16s %-6s %-9s %7s %8s  %s\n", "edge", "cross", "direction", "n_comb", "n_TRUE", "verdict"))
for (i in order(agg$flips, agg$exposure, decreasing=c(TRUE, FALSE), method="radix")) {
  cat(sprintf("%-16s %-6s %-9s %7d %8d  %s\n",
      paste0(agg$exposure[i], "->", agg$outcome[i]), agg$cross[i], agg$direction[i],
      agg$n_comb[i], agg$n_true[i], if (agg$flips[i]) "FLIPS ACROSS GRID" else "stable"))
}
cat("\nBASE-PREVALENCE CALLS (the values the shipped analysis uses):\n")
b <- res[res$at_base, ]
cat(sprintf("  %d edges at base prevalence; Steiger-correct in %d\n", nrow(b), sum(b$correct_dir)))
bad <- b[!b$correct_dir, ]
if (nrow(bad)) {
  cat("  Edges the LIABILITY scale calls wrong-direction at base prevalence (the shipped graph\n")
  cat("  admitted them on the QUANTITATIVE approximation):\n")
  for (i in seq_len(nrow(bad)))
    cat(sprintf("    %s->%s  (cross=%s, %s)  r2exp=%.5f r2out=%.5f\n",
        bad$exposure[i], bad$outcome[i], bad$cross[i], bad$direction[i], bad$r2_exp[i], bad$r2_out[i]))
} else {
  cat("  none\n")
}
sink()
cat("\nwrote results/steiger_prevalence_sweep.{csv,txt}\n")
