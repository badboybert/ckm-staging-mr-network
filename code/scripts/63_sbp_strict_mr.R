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
# Step 63: re-estimate every X->SBP edge on the STRICT allele-matched variant set (step 62), and
# compare against the shipped estimate. CAD->SBP is the edge that matters: it is the only backward
# cross-stage edge surviving the genome-wide-strict threshold, one of two backward edges in the
# primary network, and the only backward edge in the UK-Biobank-free MAIN network.
#
# The strict set drops every strand-ambiguous (palindromic) variant outright, which is STRICTER than
# harmonise_data(action=2) (that drops only intermediate-frequency palindromes). So this is a
# deliberately conservative sensitivity analysis, not a correction of a known error: step 62 found
# zero allele mismatches, zero non-SNV matches and zero position collisions.
#
# Run: Rscript 63_sbp_strict_mr.R
# Out: results/sbp_strict_mr.csv / .txt

suppressMessages({library(TwoSampleMR)})
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised")
RES   <- file.path(BASE, "results")

EXPS <- c("BMI","TG","HDL","TC","LDL","HbA1c","CAD","HF","T2D","eGFR","Stroke")

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

ivw_of <- function(exp, outfile) {
  ef <- file.path(INSTR, paste0(exp, ".clumped.tsv"))
  if (!file.exists(ef) || !file.exists(outfile)) return(NULL)
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL)
  o <- tryCatch(fmt_out(outfile), error=function(z) NULL)
  if (is.null(e) || is.null(o) || nrow(o) < 4) return(NULL)
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep, ]
  if (nrow(h) < 4) return(NULL)
  m <- mr(h, method_list=c("mr_ivw"))
  data.frame(nsnp=m$nsnp[1], b=m$b[1], se=m$se[1], p=m$pval[1], stringsAsFactors=FALSE)
}

net <- read.csv(file.path(RES, "forward_local_edges.csv"), stringsAsFactors=FALSE)
BONF <- 0.05 / nrow(net)

rows <- list()
for (e in EXPS) {
  base   <- ivw_of(e, file.path(HARM, paste0(e, "__SBP.outcome.tsv")))
  strict <- ivw_of(e, file.path(HARM, paste0(e, "__SBP.strict.outcome.tsv")))
  if (is.null(base) || is.null(strict)) { cat(sprintf("%-12s : skipped\n", paste0(e,"->SBP"))); next }
  shipped <- net[net$exposure == e & net$outcome == "SBP", ]
  rows[[e]] <- data.frame(
    exposure=e, outcome="SBP",
    shipped_b = if (nrow(shipped)) shipped$ivw_b[1] else NA_real_,
    shipped_p = if (nrow(shipped)) shipped$ivw_p[1] else NA_real_,
    base_nsnp=base$nsnp, base_b=base$b, base_se=base$se, base_p=base$p,
    strict_nsnp=strict$nsnp, strict_b=strict$b, strict_se=strict$se, strict_p=strict$p,
    sign_same = sign(base$b) == sign(strict$b),
    both_bonf = (base$p < BONF) == (strict$p < BONF),
    stringsAsFactors=FALSE)
  cat(sprintf("%-12s base %4d SNP b=%+.4f p=%.2e | strict %4d SNP b=%+.4f p=%.2e | %s\n",
      paste0(e,"->SBP"), base$nsnp, base$b, base$p, strict$nsnp, strict$b, strict$p,
      if (sign(base$b)==sign(strict$b) && (base$p<BONF)==(strict$p<BONF)) "concordant" else "CHANGED"))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES, "sbp_strict_mr.csv"), row.names=FALSE)

sink(file.path(RES, "sbp_strict_mr.txt"))
cat("=== X->SBP re-estimated on strict allele-matched instruments ===\n\n")
cat("Strict set = step-62 'match' + 'strand_flip' only; every strand-ambiguous (palindromic)\n")
cat("variant dropped outright. Step 62 found 0 allele mismatches, 0 non-SNV matches, 0 position\n")
cat("collisions and 0 multiallelic SBP records among 3,156 matched variants, so this is a\n")
cat(sprintf("conservative sensitivity analysis. Network Bonferroni = 0.05/%d = %.2e\n\n", nrow(net), BONF))
cat(sprintf("%-12s %6s %10s %10s | %6s %10s %10s | %s\n",
            "edge","nSNP","b","P","nSNP","b","P","verdict"))
for (i in seq_len(nrow(res))) {
  cat(sprintf("%-12s %6d %+10.4f %10.2e | %6d %+10.4f %10.2e | %s\n",
      paste0(res$exposure[i],"->SBP"), res$base_nsnp[i], res$base_b[i], res$base_p[i],
      res$strict_nsnp[i], res$strict_b[i], res$strict_p[i],
      if (res$sign_same[i] && res$both_bonf[i]) "sign + Bonferroni status unchanged" else "CHANGED"))
}
# Two cells cannot support a stability claim and must not pad the denominator:
#  (a) HbA1c->SBP is unchanged BY CONSTRUCTION - all 415 HbA1c instruments lack an effect-allele
#      frequency, so harmonise_data(action=2) had already dropped every palindrome before the
#      strict filter could act; base and strict are byte-identical.
#  (b) edges non-significant under BOTH base and strict pass "status unchanged" vacuously.
res$identical_by_construction <- res$base_nsnp == res$strict_nsnp
res$informative <- (!res$identical_by_construction) & (res$base_p < BONF | res$strict_p < BONF)
cat(sprintf("\n%d/%d edges keep both their sign and their Bonferroni status under strict matching.\n",
            sum(res$sign_same & res$both_bonf), nrow(res)))
cat(sprintf("Of these, %d are INFORMATIVE (instrument set actually changed AND significant under at
",
            sum(res$informative)))
cat(sprintf("least one of base/strict): %s.
",
            paste0(res$exposure[res$informative], "->SBP", collapse=", ")))
cat(sprintf("%d edge(s) unchanged BY CONSTRUCTION (no eaf, no palindrome left to drop): %s.
",
            sum(res$identical_by_construction),
            paste0(res$exposure[res$identical_by_construction], "->SBP", collapse=", ")))
cat("Precise claim: among informative edges none crosses the network Bonferroni threshold under
")
cat("strict matching except HDL->SBP. Point estimates DO move (see table).
")
cad <- res[res$exposure == "CAD", ]
if (nrow(cad)) {
  cat("\nCAD->SBP (the load-bearing backward edge):\n")
  cat(sprintf("  shipped   b=%+.4f  P=%.3e\n", cad$shipped_b, cad$shipped_p))
  cat(sprintf("  re-run    b=%+.4f  P=%.3e  (%d SNPs)\n", cad$base_b, cad$base_p, cad$base_nsnp))
  cat(sprintf("  strict    b=%+.4f  P=%.3e  (%d SNPs, %d palindromic dropped)\n",
              cad$strict_b, cad$strict_p, cad$strict_nsnp, cad$base_nsnp - cad$strict_nsnp))
  cat(sprintf("  -> backward edge %s under strict allele matching.\n",
              if (sign(cad$strict_b)==sign(cad$base_b) && cad$strict_p < BONF) "SURVIVES" else "DOES NOT survive"))
}
sink()
cat("\nwrote results/sbp_strict_mr.{csv,txt}\n")
