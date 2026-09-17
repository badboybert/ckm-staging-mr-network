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
# Step 78: what changes if the European HbA1c exposure carries its TRUE sample size?
#
# Round-7 provenance finding (2026-09-16). scripts/02_download.sh pulls the European HbA1c file from
# GWAS Catalog accession GCST90014006, which the Catalog attributes to Mbatchou et al. 2021
# (PMID 34017140), UK Biobank, 389,889 European participants. The package described it as MAGIC /
# Chen et al. 2021 with n = 146,806, and lib_formats.py assigns that n as a constant because the
# harmonised file carries no sample-size column. The filename had become the fact.
#
# The betas and standard errors are unaffected - they are whatever the file holds - so no MR point
# estimate moves. What DOES depend on n is the Steiger directionality test, which is the gate that
# decides which edges enter the 53-edge staging graph. This script re-runs Steiger for every edge
# with HbA1c on either side under both sample sizes and prints any call that changes.
#
# Output: results/hba1c_provenance_steiger.csv
suppressMessages({library(TwoSampleMR)})
set.seed(20260718)
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised")
RES   <- file.path(BASE, "results")
N_OLD <- 146806      # what the package assumed (MAGIC / Chen 2021)
N_NEW <- 389889      # what GCST90014006 actually is (Mbatchou 2021, UK Biobank)

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

ofiles <- list.files(HARM, pattern="__.*\\.outcome\\.tsv$", full.names=TRUE)
rows <- list()
for (of in ofiles) {
  base  <- sub("\\.outcome\\.tsv$", "", basename(of))
  parts <- strsplit(base, "__")[[1]]; exp <- parts[1]; out <- parts[2]
  if (exp != "HbA1c" && out != "HbA1c") next
  ef <- file.path(INSTR, paste0(exp, ".clumped.tsv"))
  if (!file.exists(ef)) next
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL)
  o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e) || is.null(o) || nrow(o) == 0) next
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep, ]
  if (nrow(h) < 3) next
  one <- function(n_hba1c) {
    hh <- h
    if (exp == "HbA1c") hh$samplesize.exposure <- n_hba1c
    if (out == "HbA1c") hh$samplesize.outcome  <- n_hba1c
    st <- tryCatch(directionality_test(hh),
                   error=function(z) data.frame(correct_causal_direction=NA, steiger_pval=NA))
    c(as.character(st$correct_causal_direction[1]), st$steiger_pval[1])
  }
  a <- one(N_OLD); b <- one(N_NEW)
  rows[[base]] <- data.frame(exposure=exp, outcome=out, nsnp=nrow(h),
                             steiger_old=a[1], steiger_p_old=as.numeric(a[2]),
                             steiger_new=b[1], steiger_p_new=as.numeric(b[2]),
                             changed=(a[1] != b[1]), stringsAsFactors=FALSE)
  cat(sprintf("%-6s -> %-6s n=%3d | n=%s: %-5s (p=%.2e) | n=%s: %-5s (p=%.2e) | %s\n",
      exp, out, nrow(h), format(N_OLD, big.mark=","), a[1], as.numeric(a[2]),
      format(N_NEW, big.mark=","), b[1], as.numeric(b[2]),
      if (a[1] != b[1]) "*** CALL CHANGES ***" else "unchanged"))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES, "hba1c_provenance_steiger.csv"), row.names=FALSE)
cat(sprintf("\n%d HbA1c edges re-tested; %d Steiger call(s) change with the corrected sample size\n",
            nrow(res), sum(res$changed)))
