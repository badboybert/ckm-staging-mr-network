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
# Step 36: refit the FULL-DENSITY lipid MVMR (LDL + HDL + TG -> CAD, no ApoB-coverage restriction)
# and COMMIT it as data.
#
# Why this script exists: the 630-instrument model is the source of a main-text number — HDL retains
# a direct effect of -0.165 after conditioning, which is what makes "HDL attenuates ~50% but a
# residual survives" true rather than "HDL dissolves". That model was only ever printed into
# results/apob_mvmr_summary.txt as prose; it was never written to a results CSV, so Supplementary Figure S7b/S7d
# carried it as hardcoded literals and the submission package shipped no source data for it.
# Same method as 32b_apob_mvmr.R (MVMR::format_mvmr -> ivw_mvmr -> strength_mvmr).
# Run: Rscript 36_lipid_fulldensity_mvmr.R
suppressMessages({library(MVMR)})

BASE  <- P4_BASE
MVDIR <- file.path(BASE, "data/mvmr")
RES   <- file.path(BASE, "results")

d <- read.delim(file.path(MVDIR, "model_LDL_HDL_TG_fulldens.tsv"), stringsAsFactors = FALSE)
cat("SNPs in full-density model:", nrow(d), "\n")

EXP   <- c("LDL", "HDL", "TG")
BXGs   <- as.matrix(d[, paste0("bx_", EXP)])
seBXGs <- as.matrix(d[, paste0("se_", EXP)])

r_in <- format_mvmr(BXGs = BXGs, BYG = d$by, seBXGs = seBXGs, seBYG = d$sey, RSID = d$SNP)
est  <- ivw_mvmr(r_in)
Fst  <- tryCatch(strength_mvmr(r_in, gencov = 0), error = function(z) NULL)

out <- data.frame(
  model      = "LDL_HDL_TG_630",
  outcome    = "CAD",
  exposure   = EXP,
  nsnp       = nrow(d),
  direct_b   = as.numeric(est[, 1]),
  direct_se  = as.numeric(est[, 2]),
  direct_p   = as.numeric(est[, 4]),
  cond_F     = if (!is.null(Fst)) as.numeric(unlist(Fst)[seq_along(EXP)]) else NA_real_,
  stringsAsFactors = FALSE
)
print(out)

# The committed analysis output (results/apob_mvmr_summary.txt) records this model as
# LDL +0.499, HDL -0.165, TG +0.140. A refit that does not reproduce those is a red flag, not a
# new result — fail loudly rather than silently ship a different number under the same label.
expect <- c(LDL = 0.499, HDL = -0.165, TG = 0.140)
for (e in EXP) {
  got <- out$direct_b[out$exposure == e]
  if (abs(got - expect[[e]]) > 0.005)
    stop(sprintf("full-density MVMR did not reproduce %s: refit %+.4f vs committed %+.3f",
                 e, got, expect[[e]]))
}
cat("reproduces the committed LDL/HDL/TG direct effects within 0.005\n")

write.csv(out, file.path(RES, "mvmr_lipid_fulldensity.csv"), row.names = FALSE)
cat("wrote results/mvmr_lipid_fulldensity.csv\n")
