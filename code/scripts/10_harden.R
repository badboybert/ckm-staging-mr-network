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
# Step 10: HARDENING — proper binary Steiger (get_r_from_lor) + MR-PRESSO on headline edges.
# Reuses the same harmonisation as 06_mr.R.
#
# DEPRECATED 2026-07-23. This script wrote BOTH results/steiger_lor.csv AND
# results/mrpresso_headline.csv, and it produced NEITHER of the committed files:
#   - steiger_lor.csv (committed) carries the columns snp_r2_exp/manual_r2x/manual_r2y/
#     manual_correct written by 10a_steiger.R; this script writes snp_r2.exposure/snp_r2.outcome
#     and no manual_* columns.
#   - mrpresso_headline.csv (committed) has nb_dist = 1000 in all 5 rows and contains only the
#     edges with <= 130 instruments; this script uses NbDistribution = min(20000, max(2000,
#     ceiling(nSNP/0.05))) (always >= 2000) and has no instrument cap, so it would write 16 rows.
#     It was written by 10b_mrpresso.R.
# Keeping two writers for the same file meant the shipped parameters depended on run order.
# Authoritative writers now: 10a_steiger.R -> steiger_lor.csv; 10b_mrpresso.R -> mrpresso_headline.csv.
# Kept only as build history. Running it would silently overwrite both committed files.
if (!("--i-know-this-is-deprecated" %in% commandArgs(trailingOnly = TRUE))) {
  stop("10_harden.R is DEPRECATED and must not write results/. Run scripts/10a_steiger.R for ",
       "steiger_lor.csv and scripts/10b_mrpresso.R for mrpresso_headline.csv (the authoritative ",
       "writers). It duplicated both outputs with different MR-PRESSO settings (NbDistribution ",
       ">=2000, no <=130-instrument cap), so whichever script ran last decided the shipped file.",
       call. = FALSE)
}
suppressMessages({library(TwoSampleMR); library(MRPRESSO)})
set.seed(42)
BASE  <- P4_BASE
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised")
RES   <- file.path(BASE,"results")

# Binary-trait metadata for liability-scale Steiger (published case/control; prevalence approx).
# CAD=Aragam2022 GCST90132314; HF=Henry2025 GCST90728695; T2D=Mahajan2018 GCST006867;
# Stroke=GIGASTROKE GCST90104539 (EUR any); CKD=Wuttke2019 GCST008065.
BIN <- list(
  CAD   =list(ncase=181522, ncontrol=1047847, prev=0.06),
  HF    =list(ncase=139533, ncontrol=1054983, prev=0.02),
  T2D   =list(ncase=74124,  ncontrol=824006,  prev=0.10),
  Stroke=list(ncase=73652,  ncontrol=1234808, prev=0.05),
  CKD   =list(ncase=41395,  ncontrol=439303,  prev=0.09)
)
is_bin <- function(x) x %in% names(BIN)

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

annotate_bin <- function(h, trait, which=c("exposure","outcome")) {
  which <- match.arg(which); m <- BIN[[trait]]
  h[[paste0("units.",which)]]      <- "log odds"
  h[[paste0("ncase.",which)]]      <- m$ncase
  h[[paste0("ncontrol.",which)]]   <- m$ncontrol
  h[[paste0("prevalence.",which)]] <- m$prev
  h[[paste0("samplesize.",which)]] <- m$ncase + m$ncontrol
  h
}

EDGES <- list(
  c("LDL","CAD"), c("CAD","HF"), c("HF","CAD"), c("BMI","CAD"), c("BMI","T2D"),
  c("T2D","CAD"), c("T2D","HF"), c("SBP","Stroke"), c("CAD","T2D"), c("HF","T2D"),
  c("HbA1c","T2D"), c("LDL","HF"), c("BMI","HF"), c("BMI","Stroke"), c("SBP","CAD"), c("SBP","HF")
)

st_rows <- list(); pr_rows <- list()
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
  # annotate units for proper Steiger
  if (is_bin(exp)) h <- annotate_bin(h, exp, "exposure") else h$units.exposure <- "SD"
  if (is_bin(out)) h <- annotate_bin(h, out, "outcome")  else h$units.outcome  <- "SD"
  dt <- tryCatch(directionality_test(h), error=function(z) NULL)
  if (!is.null(dt)) {
    st_rows[[key]] <- data.frame(exposure=exp, outcome=out, nsnp=nrow(h),
      snp_r2.exposure=dt$snp_r2.exposure, snp_r2.outcome=dt$snp_r2.outcome,
      correct_dir=dt$correct_causal_direction, steiger_p=dt$steiger_pval,
      exp_binary=is_bin(exp), out_binary=is_bin(out), stringsAsFactors=FALSE)
    cat(sprintf("%-14s Steiger-lor: r2exp=%.4f r2out=%.4f correct=%s p=%.2e\n",
        key, dt$snp_r2.exposure, dt$snp_r2.outcome, dt$correct_causal_direction, dt$steiger_pval))
  }
  # MR-PRESSO
  dfp <- data.frame(by=h$beta.outcome, bx=h$beta.exposure, sy=h$se.outcome, sx=h$se.exposure)
  nb  <- min(20000, max(2000, ceiling(nrow(h)/0.05)))
  pr  <- tryCatch(mr_presso(BetaOutcome="by", BetaExposure="bx", SdOutcome="sy", SdExposure="sx",
            data=dfp, OUTLIERtest=TRUE, DISTORTIONtest=TRUE, NbDistribution=nb, SignifThreshold=0.05, seed=42),
            error=function(z){cat("  PRESSO fail:",conditionMessage(z),"\n"); NULL})
  if (!is.null(pr)) {
    mm <- pr$`Main MR results`
    raw_b  <- mm$`Causal Estimate`[mm$`MR Analysis`=="Raw"][1]
    raw_p  <- mm$`P-value`[mm$`MR Analysis`=="Raw"][1]
    cor_b  <- mm$`Causal Estimate`[mm$`MR Analysis`=="Outlier-corrected"][1]
    cor_p  <- mm$`P-value`[mm$`MR Analysis`=="Outlier-corrected"][1]
    gt     <- pr$`MR-PRESSO results`$`Global Test`
    glob_p <- if(!is.null(gt$Pvalue)) gt$Pvalue else NA
    glob_r <- if(!is.null(gt$RSSobs)) gt$RSSobs else NA
    outl   <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    n_out  <- if(is.null(outl) || (length(outl)==1 && is.character(outl))) 0 else length(outl)
    dist_p <- pr$`MR-PRESSO results`$`Distortion Test`$Pvalue
    pr_rows[[key]] <- data.frame(exposure=exp, outcome=out, nsnp=nrow(h), nb_dist=nb,
      raw_b=raw_b, raw_p=raw_p, corrected_b=cor_b, corrected_p=cor_p,
      global_RSSobs=glob_r, global_p=glob_p, n_outliers=n_out,
      distortion_p=ifelse(is.null(dist_p),NA,dist_p), stringsAsFactors=FALSE)
    cat(sprintf("%-14s PRESSO: raw=%+.4f(p=%.2e) corr=%s(p=%s) global_p=%s outliers=%d\n",
        key, raw_b, raw_p, ifelse(is.na(cor_b),"NA",sprintf("%+.4f",cor_b)),
        ifelse(is.na(cor_p),"NA",sprintf("%.2e",cor_p)),
        ifelse(is.na(glob_p),"NA",format(glob_p,digits=3)), n_out))
  }
}
st <- do.call(rbind, st_rows); pr <- do.call(rbind, pr_rows)
write.csv(st, file.path(RES,"steiger_lor.csv"), row.names=FALSE)
write.csv(pr, file.path(RES,"mrpresso_headline.csv"), row.names=FALSE)
cat("\nWrote steiger_lor.csv (", nrow(st), ") + mrpresso_headline.csv (", nrow(pr), ")\n")
