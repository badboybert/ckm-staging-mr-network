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
# Step 17: MVMR robustness — MVMR-Egger (InSIDE/directional-pleiotropy) + MVMR qhet.
# MendelianRandomization pkg has a broken xfun::attr import on this box, so MVMR-Egger is
# implemented directly: orient every SNP by the primary exposure's sign, then a 1/sey^2-weighted
# multivariable regression WITH an intercept -> the intercept is the directional-pleiotropy test
# and the slopes are the pleiotropy-adjusted direct effects (lm SEs carry the over-dispersion
# factor, matching ivw_mvmr). Consumes data/mvmr/model_*.tsv. Run: Rscript 17_mvmr_robust.R
suppressMessages({library(MVMR)})
BASE  <- P4_BASE
MVDIR <- file.path(BASE,"data/mvmr"); RES <- file.path(BASE,"results")
outmap <- c(CAD_full="CAD", CAD_BMI_T2D="CAD", HF_via_CAD="HF", Stroke_full="Stroke")

rows <- list()
for (f in list.files(MVDIR, pattern="^model_.*\\.tsv$", full.names=TRUE)) {
  name <- sub("^model_","",sub("\\.tsv$","",basename(f)))
  d <- read.delim(f, stringsAsFactors=FALSE)
  bx_cols <- grep("^bx_", names(d), value=TRUE); se_cols <- grep("^se_", names(d), value=TRUE)
  expos <- sub("^bx_","",bx_cols); outc <- outmap[[name]]
  BX <- as.matrix(d[,bx_cols,drop=FALSE]); BXse <- as.matrix(d[,se_cols,drop=FALSE])
  BY <- d$by; BYse <- d$sey; w <- 1/BYse^2

  ## --- MVMR-IVW (no intercept), over-dispersion SEs via lm ---
  ivwfit <- lm(BY ~ 0 + BX, weights=w); ivws <- summary(ivwfit)$coefficients

  ## --- MVMR-Egger: orient by primary exposure (col 1), fit WITH intercept ---
  orient <- sign(BX[,1]); orient[orient==0] <- 1
  BXo <- BX * orient; BYo <- BY * orient
  eggfit <- lm(BYo ~ BXo, weights=w); eggs <- summary(eggfit)$coefficients
  egg_int  <- eggs[1,1]; egg_intp <- eggs[1,4]     # intercept = directional pleiotropy

  ## --- MVMR qhet (heterogeneity-robust); pcor approximated by instrument-beta correlation ---
  pcor <- tryCatch(cor(BX), error=function(e) diag(length(expos)))
  rin  <- format_mvmr(BXGs=BX, BYG=BY, seBXGs=BXse, seBYG=BYse, RSID=d$SNP)
  qh   <- tryCatch(qhet_mvmr(rin, pcor, CI=FALSE, iterations=50), error=function(e) NULL)

  cat(sprintf("\n=== %s : %s ~ %s  (nSNP=%d) ===\n", name, outc, paste(expos,collapse=" + "), nrow(d)))
  cat(sprintf("  MVMR-Egger intercept=%+.4f p=%.3g  (InSIDE/directional-pleiotropy test)\n", egg_int, egg_intp))
  for (i in seq_along(expos)) {
    e <- expos[i]
    ivw_b <- ivws[i,1]; ivw_p <- ivws[i,4]
    egg_b <- eggs[i+1,1]; egg_p <- eggs[i+1,4]      # +1 for intercept row
    qh_b  <- if(!is.null(qh)) as.numeric(qh[i,1]) else NA
    cat(sprintf("  %-6s ivw=%+.3f(p=%.2e)  egger=%+.3f(p=%.2e)  qhet=%s\n",
        e, ivw_b, ivw_p, egg_b, egg_p, ifelse(is.na(qh_b),"NA",sprintf("%+.3f",qh_b))))
    rows[[paste(name,e)]] <- data.frame(model=name, outcome=outc, exposure=e, nsnp=nrow(d),
      ivw_b=ivw_b, ivw_p=ivw_p, egger_b=egg_b, egger_p=egg_p,
      egger_intercept=egg_int, egger_intercept_p=egg_intp, qhet_b=qh_b, stringsAsFactors=FALSE)
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(RES,"mvmr_robust.csv"), row.names=FALSE)
cat("\nWrote mvmr_robust.csv with", nrow(res), "rows\n")
