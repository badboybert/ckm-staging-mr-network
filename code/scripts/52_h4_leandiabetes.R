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
# Step 52 (H4): robust lean-diabetes suite for the EAS T2D->BMI edge (reviewer M6).
# (1) Hartung-Knapp random-effects meta of BBJ + TPMI (metafor, knha=TRUE) vs DL/FE.
# (2) Robust MR on the population KoGES DM->BMI: IVW/Egger/weighted-median/weighted-mode/penalised-WM,
#     radial-IVW outlier scan (RadialMR), leave-one-out, and leave-one-locus-out for 12q24 + MC4R.
# Run: Rscript 52_h4_leandiabetes.R  ->  results/h4_leandiabetes.txt
suppressMessages({library(TwoSampleMR); library(RadialMR); library(metafor)})
set.seed(20260718)
BASE  <- P4_BASE
INSTR <- file.path(BASE, "data/instruments"); HARM <- file.path(BASE, "data/harmonised_eas")
sink(file.path(BASE, "results/h4_leandiabetes.txt"), split = TRUE)

cat("=== H4: robust lean-diabetes suite (EAS T2D -> BMI) ===\n\n")

## (1) Hartung-Knapp meta of the two hospital cohorts
## The two inputs are READ from the committed network CSVs (not hardcoded): BBJ T2D->BMI from
## results/network_eas_edges.csv, TPMI T2D->BMI from results/network_tpmi_full.csv.
## NOTE: the committed run of this script used the estimates ROUNDED TO 4 DECIMALS
## (BBJ -0.0977/0.0188; TPMI -0.0456/0.0169), so the round() below is deliberate and required to
## reproduce results/h4_leandiabetes.txt exactly. Using full CSV precision instead shifts the
## reported CIs by ~1e-4 (e.g. HK [-0.4023,+0.2603] p=0.2241 vs the committed [-0.4019,+0.2599]
## p=0.2238). Do not remove the round() without re-deriving every downstream quoted number.
read_edge <- function(path, exp, out) {
  d <- read.csv(file.path(BASE, path), stringsAsFactors = FALSE)
  r <- d[d$exposure == exp & d$outcome == out, ]
  stopifnot(nrow(r) == 1)
  c(b = round(as.numeric(r$ivw_b), 4), se = round(as.numeric(r$ivw_se), 4))
}
bbj  <- read_edge("results/network_eas_edges.csv",  "T2D", "BMI")
tpmi <- read_edge("results/network_tpmi_full.csv", "T2D", "BMI")
b  <- c(BBJ = unname(bbj["b"]),  TPMI = unname(tpmi["b"]))
se <- c(unname(bbj["se"]), unname(tpmi["se"]))
# assert the parsed+rounded values still equal the literals the committed run used (drift guard)
stopifnot(all(abs(b  - c(-0.0977, -0.0456)) < 1e-4),
          all(abs(se - c( 0.0188,  0.0169)) < 1e-4))
cat(sprintf("(1) EAS hospital-cohort meta of T2D->BMI (BBJ b=%.3f se=%.3f; TPMI b=%.3f se=%.3f):\n",
            b[1], se[1], b[2], se[2]))
fe  <- rma(yi = b, sei = se, method = "FE")
dl  <- rma(yi = b, sei = se, method = "DL")
hk  <- rma(yi = b, sei = se, method = "REML", test = "knha")     # Hartung-Knapp
for (nm in c("FE", "DL", "HK")) {
  m <- get(tolower(nm))
  cat(sprintf("   %-3s: b=%+.4f  95%% CI [%+.4f, %+.4f]  p=%.4f  (I2=%.0f%%)\n",
              nm, m$beta, m$ci.lb, m$ci.ub, m$pval, ifelse(is.na(m$I2), 0, m$I2)))
}
cat("   -> Hartung-Knapp (the conservative RE inference for k=2) is the honest headline.\n\n")

## (2) robust MR on population KoGES DM->BMI
ef <- file.path(INSTR, "KoGES_DM.eas.clumped.tsv")
of <- file.path(HARM, "KoGES_DM__KoGES_BMI.outcome.tsv")
e <- read.delim(ef, stringsAsFactors = FALSE)
e <- format_data(e, type = "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se",
                 effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                 eaf_col = "eaf", pval_col = "pval", samplesize_col = "N", phenotype_col = "Phenotype")
o <- read.delim(of, stringsAsFactors = FALSE)
o <- format_data(o, type = "outcome", snp_col = "SNP", beta_col = "beta", se_col = "se",
                 effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                 eaf_col = "eaf", pval_col = "pval", phenotype_col = "Phenotype")
h <- harmonise_data(e, o, action = 2); h <- h[h$mr_keep, ]
cat(sprintf("(2) Population KoGES DM->BMI: %d harmonised instruments\n", nrow(h)))
m <- mr(h, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                           "mr_weighted_mode", "mr_penalised_weighted_median"))
for (i in seq_len(nrow(m)))
  cat(sprintf("   %-28s b=%+.4f se=%.4f p=%.4f\n", m$method[i], m$b[i], m$se[i], m$pval[i]))

## radial-IVW outlier scan
rad <- tryCatch({
  rin <- format_radial(h$beta.exposure, h$beta.outcome, h$se.exposure, h$se.outcome, h$SNP)
  ivw_radial(rin, alpha = 0.05)
}, error = function(z) NULL)
if (!is.null(rad)) {
  no <- tryCatch(nrow(rad$outliers), error = function(z) 0)
  cat(sprintf("   radial-IVW: b=%+.4f  Q-outliers flagged=%s\n",
              rad$coef[1, 1], ifelse(is.null(no) || is.na(no), "0", no)))
}

## leave-one-out + leave-one-locus-out (12q24=rs2074356, MC4R=rs6567160)
cat("\n   Leave-one-out (IVW, drop each SNP):\n")
loo <- mr_leaveoneout(h)
for (i in which(loo$SNP != "All"))
  cat(sprintf("     drop %-12s b=%+.4f p=%.4f\n", loo$SNP[i], loo$b[i], loo$p[i]))
allb <- loo$b[loo$SNP == "All"][1]
cat(sprintf("   full IVW = %+.4f; LOO range [%+.4f, %+.4f]\n",
            allb, min(loo$b[loo$SNP != "All"]), max(loo$b[loo$SNP != "All"])))
for (loc in c("rs2074356", "rs6567160")) {
  if (loc %in% h$SNP) {
    h2 <- h[h$SNP != loc, ]
    mm <- mr(h2, method_list = "mr_ivw")
    cat(sprintf("   drop %s (%s): IVW b=%+.4f p=%.4f (n=%d)\n",
                loc, ifelse(loc == "rs2074356", "12q24/ALDH2", "MC4R"), mm$b, mm$pval, nrow(h2)))
  }
}
cat("\nREAD: report the Hartung-Knapp meta as the conservative EAS hospital estimate; the population\n")
cat("KoGES estimate under weighted-mode/median + radial (pleiotropy-robust) and its sensitivity to\n")
cat("12q24/MC4R govern whether the negative edge is a robust signature or an outlier-driven artifact.\n")
sink()
cat("wrote results/h4_leandiabetes.txt\n")
