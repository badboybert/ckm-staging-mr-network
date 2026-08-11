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
# Step 15: CAUSE (correlated-pleiotropy vs causal) on headline edges.
# CAUSE needs GENOME-WIDE stats for both traits; we substitute PLINK LD-clumping (1000G EUR)
# for cause::ld_prune (no authors' LD reference on disk). Checkpointed per edge -> resumable.
# Run (PowerShell, with rtools on PATH for gzip): Rscript 15_cause.R
suppressMessages({library(data.table); library(cause)})
set.seed(42)
ROOT <- CKM_ROOT
BASE <- file.path(ROOT,"paper 4/independent_build")
RES  <- file.path(BASE,"results"); MV <- file.path(BASE,"data/cause"); dir.create(MV,showWarnings=FALSE,recursive=TRUE)
PLINK <- file.path(ROOT,"lhcmr_eas_refs/tools/plink.exe")
PANEL <- file.path(ROOT,"paper 6/analysis/g1000_eur/g1000_eur")
OUT   <- file.path(RES,"cause_headline.csv")

# node -> (file, snp, ea, oa, beta, se, p)
NODES <- list(
 CAD  =list(f="paper 6/analysis/data/outcomes/CAD_GCST90132314.h.tsv.gz", snp="rsid", ea="effect_allele", oa="other_allele", b="beta", se="standard_error", p="p_value"),
 HF   =list(f="paper 6/analysis/data/outcomes/HF_GCST90728695.tsv.gz",    snp="rsid", ea="effect_allele", oa="other_allele", b="beta", se="standard_error", p="p_value"),
 T2D  =list(f="paper 6/analysis/data/outcomes/T2D_GCST006867.h.tsv.gz",   snp="hm_rsid", ea="effect_allele", oa="other_allele", b="beta", se="standard_error", p="p_value"),
 LDL  =list(f="paper 6/analysis/data/outcomes/LDL_GLGC2021_EUR.gz",       snp="rsID", ea="ALT", oa="REF", b="EFFECT_SIZE", se="SE", p="pvalue"),
 BMI  =list(f="paper 4/independent_build/data/raw/BMI_Yengo2018_EUR.txt.gz", snp="SNP", ea="Tested_Allele", oa="Other_Allele", b="BETA", se="SE", p="P"),
 HbA1c=list(f="paper 4/independent_build/data/raw/HbA1c_MAGIC_EUR.h.tsv.gz", snp="hm_rsid", ea="effect_allele", oa="other_allele", b="beta", se="standard_error", p="p_value")
)
# ordered discriminating-first (backward/asymmetry edges use medium CAD/HF/T2D files; huge LDL last)
EDGES <- list(c("CAD","HF"), c("HF","CAD"), c("CAD","T2D"), c("HF","T2D"), c("T2D","CAD"),
              c("T2D","HF"), c("BMI","T2D"), c("HbA1c","T2D"), c("BMI","CAD"), c("LDL","CAD"))

CACHE <- new.env()
read_node <- function(n){
  if(!is.null(CACHE[[n]])) return(CACHE[[n]])
  cfg <- NODES[[n]]; path <- file.path(ROOT,cfg$f)
  dt <- fread(cmd=sprintf('gzip -dc "%s"', path),
              select=c(cfg$snp,cfg$ea,cfg$oa,cfg$b,cfg$se,cfg$p), showProgress=FALSE)
  setnames(dt, c(cfg$snp,cfg$ea,cfg$oa,cfg$b,cfg$se,cfg$p), c("snp","A1","A2","beta","se","p"))
  dt <- dt[!is.na(beta) & !is.na(se) & se>0 & grepl("^rs",snp)]
  dt[, A1:=toupper(A1)][, A2:=toupper(A2)]
  dt <- unique(dt, by="snp")
  CACHE[[n]] <- dt; dt
}

done <- if(file.exists(OUT)) fread(OUT) else NULL
if(!is.null(done)) done_edges <- paste(done$exposure, done$outcome) else done_edges <- character(0)

rows <- list()
for(ed in EDGES){
  E <- ed[1]; O <- ed[2]; key <- paste(E,O)
  if(key %in% done_edges){ cat(key,"already done, skip\n"); next }
  cat("\n=== CAUSE", E, "->", O, "===\n"); flush.console()
  res_row <- tryCatch({
    x1 <- read_node(E); x2 <- read_node(O)
    X <- gwas_merge(x1, x2, snp_name_cols=c("snp","snp"), beta_hat_cols=c("beta","beta"),
                    se_cols=c("se","se"), A1_cols=c("A1","A1"), A2_cols=c("A2","A2"))
    cat("  merged SNPs:", nrow(X), "\n"); flush.console()
    # nuisance params from a random subset (200k is ample; params are stable well below 1M)
    vars <- X$snp; nsub <- min(2e5, length(vars))
    params <- est_cause_params(X, sample(vars, nsub))
    # LD-prune via PLINK on exposure p (p1 from beta_hat_1/seb1).
    # Keep p1 as a VECTOR — do NOT add it as an X column (pass pristine cause_data to cause()).
    p1vec <- 2*pnorm(-abs(X$beta_hat_1/X$seb1))
    assoc <- file.path(MV, sprintf("assoc_%s_%s.txt", E,O))
    fwrite(data.table(SNP=X$snp, P=p1vec), assoc, sep="\t")
    pref <- file.path(MV, sprintf("clump_%s_%s", E,O))
    system2(PLINK, c("--bfile",shQuote(PANEL),"--clump",shQuote(assoc),"--clump-p1","1e-4",
      "--clump-r2","0.01","--clump-kb","10000","--clump-snp-field","SNP","--clump-field","P",
      "--out",shQuote(pref)), stdout=FALSE, stderr=FALSE)
    cl <- fread(paste0(pref,".clumped"), header=TRUE)
    pruned <- cl$SNP[grepl("^rs",cl$SNP)]
    cat("  LD-pruned instruments (p1<1e-4, r2<0.01):", length(pruned), "\n"); flush.console()
    r <- cause(X=X, variants=pruned, param_ests=params)
    saveRDS(r, file.path(MV, sprintf("cause_obj_%s_%s.rds", E,O)))   # insurance: re-extract later w/o re-run
    s <- summary(r)
    # summary$z,$p = the sharing-vs-causal test (z<0 => causal model preferred; p one-sided)
    z <- as.numeric(s$z); pcz <- as.numeric(s$p)
    qc <- s$quants[[2]]   # causal-model params; rows = [median, 2.5%, 97.5%] (unnamed)
    gamma_med <- qc[1,"gamma"]; gamma_lo <- qc[2,"gamma"]; gamma_hi <- qc[3,"gamma"]
    eta_med   <- qc[1,"eta"];   q_med    <- qc[1,"q"]
    gamma_sig <- !is.na(gamma_lo) && ((gamma_lo>0)==(gamma_hi>0))
    verdict <- if(!is.na(pcz) && pcz<0.05 && z<0 && gamma_sig) "CAUSAL" else "SHARING/PLEIOTROPY"
    cat(sprintf("  sharing->causal z=%.2f p=%.3g | gamma=%.3f [%.3f,%.3f] eta=%.3f q=%.2f -> %s\n",
        z, pcz, gamma_med, gamma_lo, gamma_hi, eta_med, q_med, verdict)); flush.console()
    data.table(exposure=E, outcome=O, n_merged=nrow(X), n_pruned=length(pruned),
      z_sharing_vs_causal=z, p=pcz, gamma_med=gamma_med, gamma_lo=gamma_lo, gamma_hi=gamma_hi,
      eta_med=eta_med, q_med=q_med, verdict=verdict)
  }, error=function(e){ cat("  CAUSE fail:", conditionMessage(e), "\n"); NULL })
  if(!is.null(res_row)){
    fwrite(res_row, OUT, append=file.exists(OUT))   # checkpoint per edge
    rows[[key]] <- res_row
  }
  gc()
}
cat("\nDone. Results ->", OUT, "\n")
