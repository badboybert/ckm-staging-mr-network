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
# Step 56: leave-one-out IVW for the four headline edges (Supplementary Figure S4).
# Regenerates results/leaveoneout_headline.csv.
# Estimation path mirrors 06_mr.R exactly: same format_data() column mapping,
# harmonise_data(action=2), subset to mr_keep, then TwoSampleMR::mr_leaveoneout().
# Run (PowerShell): Rscript 56_leaveoneout_headline.R
suppressMessages({library(TwoSampleMR)})
# mr_leaveoneout is deterministic (IVW only), but seed anyway so the script is
# byte-identical to 06_mr.R's environment.
set.seed(20260718)
BASE <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised")
RES <- file.path(BASE,"results"); dir.create(RES, showWarnings=FALSE, recursive=TRUE)

fmt_exp <- function(f, name) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

# Edge order matches the committed CSV.
EDGES <- list(c("CAD","HF"), c("BMI","CAD"), c("LDL","CAD"), c("BMI","T2D"))

rows <- list()
for (ed in EDGES) {
  exp <- ed[1]; out <- ed[2]
  ef <- file.path(INSTR, paste0(exp,".clumped.tsv"))
  of <- file.path(HARM, sprintf("%s__%s.outcome.tsv", exp, out))
  stopifnot(file.exists(ef), file.exists(of))
  e <- fmt_exp(ef, exp)
  o <- fmt_out(of)
  h <- harmonise_data(e, o, action=2)
  h <- h[h$mr_keep,]
  loo <- mr_leaveoneout(h)
  loo$edge <- sprintf("%s -> %s", exp, out)
  loo$exposure <- exp; loo$outcome <- out
  rows[[paste0(exp,"__",out)]] <- loo[, c("edge","exposure","outcome","SNP","b","se","p","samplesize")]
  cat(sprintf("%-4s -> %-4s : %d LOO rows (%d SNPs + All)\n", exp, out, nrow(loo), nrow(loo)-1L))
}
res <- do.call(rbind, rows)
outfile <- file.path(RES,"leaveoneout_headline.REGEN.csv")
write.csv(res, outfile, row.names=FALSE)
cat("\nWrote", outfile, "with", nrow(res), "rows\n")
