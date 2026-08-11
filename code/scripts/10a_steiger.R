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
# Step 10a: PROPER binary Steiger via get_r_from_lor (pre-computed r columns).
# directionality_test() silently falls back to the quantitative approximation unless
# r.exposure/r.outcome are already present, so we compute them explicitly.
# Run (PowerShell): Rscript 10a_steiger.R
suppressMessages({library(TwoSampleMR)})
BASE  <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised")
RES   <- file.path(BASE,"results")

BIN <- list(
  CAD   =list(ncase=181522, ncontrol=1047847, prev=0.06),
  HF    =list(ncase=139533, ncontrol=1054983, prev=0.02),
  T2D   =list(ncase=74124,  ncontrol=824006,  prev=0.10),
  Stroke=list(ncase=73652,  ncontrol=1234808, prev=0.05),
  CKD   =list(ncase=41395,  ncontrol=439303,  prev=0.09))
is_bin <- function(x) x %in% names(BIN)

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

# per-SNP correlation for one side of a harmonised table
r_side <- function(h, side, trait) {
  b   <- h[[paste0("beta.",side)]]; p <- h[[paste0("pval.",side)]]
  n   <- h[[paste0("samplesize.",side)]]; eaf <- h[[paste0("eaf.",side)]]
  if (is_bin(trait)) {
    m <- BIN[[trait]]
    af <- eaf; af[is.na(af)] <- 0.5           # neutral fallback if EAF missing
    r <- get_r_from_lor(lor=b, af=af, ncase=m$ncase, ncontrol=m$ncontrol, prevalence=m$prev)
    r * sign(b)
  } else {
    get_r_from_pn(p=p, n=n) * sign(b)
  }
}

EDGES <- list(
  c("LDL","CAD"), c("CAD","HF"), c("HF","CAD"), c("BMI","CAD"), c("BMI","T2D"),
  c("T2D","CAD"), c("T2D","HF"), c("SBP","Stroke"), c("CAD","T2D"), c("HF","T2D"),
  c("HbA1c","T2D"), c("LDL","HF"), c("BMI","HF"), c("BMI","Stroke"), c("SBP","CAD"), c("SBP","HF"),
  c("T2D","BMI"), c("CAD","BMI"))

rows <- list()
for (ed in EDGES) {
  exp <- ed[1]; out <- ed[2]; key <- paste0(exp,"__",out)
  ef <- file.path(INSTR, paste0(exp,".clumped.tsv"))
  of <- file.path(HARM, paste0(key,".outcome.tsv"))
  if (!file.exists(ef) || !file.exists(of)) { cat(sprintf("%-14s : missing files\n",key)); next }
  e <- tryCatch(fmt_exp(ef), error=function(z) NULL)
  o <- tryCatch(fmt_out(of), error=function(z) NULL)
  if (is.null(e)||is.null(o)||nrow(o)==0) { cat(sprintf("%-14s : no data\n",key)); next }
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep,]
  if (nrow(h)<4) { cat(sprintf("%-14s : <4 SNPs\n",key)); next }
  h$r.exposure <- r_side(h,"exposure",exp)
  h$r.outcome  <- r_side(h,"outcome", out)
  dt <- tryCatch(directionality_test(h), error=function(z){cat("  dt fail:",conditionMessage(z),"\n"); NULL})
  # manual snp_r2 as cross-check (independent of directionality_test internals)
  r2x <- sum(h$r.exposure^2, na.rm=TRUE); r2y <- sum(h$r.outcome^2, na.rm=TRUE)
  if (!is.null(dt)) {
    rows[[key]] <- data.frame(exposure=exp, outcome=out, nsnp=nrow(h),
      snp_r2_exp=dt$snp_r2.exposure, snp_r2_out=dt$snp_r2.outcome,
      correct_dir=dt$correct_causal_direction, steiger_p=dt$steiger_pval,
      manual_r2x=r2x, manual_r2y=r2y, manual_correct=(r2x>r2y),
      exp_binary=is_bin(exp), out_binary=is_bin(out), stringsAsFactors=FALSE)
    cat(sprintf("%-14s r2exp=%.4f r2out=%.4f correct=%-5s p=%.2e  [binary exp=%s out=%s]\n",
        key, dt$snp_r2.exposure, dt$snp_r2.outcome, as.character(dt$correct_causal_direction),
        dt$steiger_pval, is_bin(exp), is_bin(out)))
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES,"steiger_lor.csv"), row.names=FALSE)
cat("\nWrote steiger_lor.csv with", nrow(res), "edges\n")
