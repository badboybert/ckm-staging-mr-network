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
# Step 18: CAUSE negative-control (specificity) — the verifier's requested binary->binary check.
# Take the cached CAD->HF merged data and PERMUTE the outcome SNP-effects across variants
# (paired beta_hat_2/seb2/p2), destroying the true SNP->outcome relationship while keeping the
# exposure instruments and marginal distributions intact. A correct pipeline must NOT return
# CAUSAL on the sham. Run (PowerShell, rtools on PATH): Rscript 18_cause_negcontrol.R
suppressMessages({library(data.table); library(cause)})
MV  <- file.path(P4_BASE, "data", "cause")
OUT <- file.path(P4_BASE, "results", "cause_negcontrol.csv")
L <- readRDS(file.path(MV,"dbg_X.rds")); X<-L$X; pruned<-L$pruned; X$p1<-NULL

extract <- function(r){
  s <- summary(r); qc <- s$quants[[2]]
  data.table(z=as.numeric(s$z), p=as.numeric(s$p),
             gamma_med=qc[1,"gamma"], gamma_lo=qc[2,"gamma"], gamma_hi=qc[3,"gamma"],
             q_med=qc[1,"q"],
             verdict=ifelse(as.numeric(s$p)<0.05 && as.numeric(s$z)<0 && ((qc[2,"gamma"]>0)==(qc[3,"gamma"]>0)),
                            "CAUSAL","SHARING/NULL"))
}

rows <- list()
for(seed in c(101, 202, 303)){
  set.seed(seed)
  Xs <- copy(X)
  idx <- sample(nrow(Xs))                       # permute outcome effects across SNPs
  Xs$beta_hat_2 <- X$beta_hat_2[idx]; Xs$seb2 <- X$seb2[idx]; Xs$p2 <- X$p2[idx]
  params <- est_cause_params(Xs, sample(Xs$snp, min(2e5, nrow(Xs))))
  r <- tryCatch(cause(X=Xs, variants=pruned, param_ests=params), error=function(e){cat("fail",seed,conditionMessage(e),"\n");NULL})
  if(is.null(r)) next
  row <- extract(r); row$sham_seed <- seed
  rows[[as.character(seed)]] <- row
  cat(sprintf("SHAM seed=%d: z=%+.2f p=%.3g gamma=%.3f[%.3f,%.3f] q=%.2f -> %s\n",
      seed, row$z, row$p, row$gamma_med, row$gamma_lo, row$gamma_hi, row$q_med, row$verdict)); flush.console()
}
res <- rbindlist(rows)
fwrite(res, OUT)
cat("\n[compare] REAL CAD->HF was CAUSAL z=-9.10 p=4.7e-20. Above are the outcome-permuted shams.\n")
cat("Wrote", OUT, "\n")
