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
# Step 19: PROPER CAUSE negative control — real confounded-but-not-causal pairs.
# HDL->CAD and HDL->T2D: HDL-C is genetically correlated with CAD/T2D through the LDL/TG lipid
# axis (univariable MR gives a spurious signal) but is NOT causal after conditioning on LDL/ApoB
# (Voight 2012, PMID 22607825). CAUSE must attribute these to correlated pleiotropy => SHARING.
# This tests CAUSE's actual specificity property (unlike the permutation sham) — real shared
# architecture present, no causal link. Contrast: LDL->CAD = CAUSAL positive control (already run).
# Per-EDGE checkpoint (resumable). Run (PowerShell, rtools on PATH for gzip): Rscript 19_cause_negcontrol_pair.R
suppressMessages({library(data.table); library(cause)})
set.seed(42)
ROOT <- CKM_ROOT
BASE <- file.path(ROOT,"paper 4/independent_build")
RES  <- file.path(BASE,"results"); MV <- file.path(BASE,"data/cause")
PLINK <- file.path(ROOT,"lhcmr_eas_refs/tools/plink.exe")
PANEL <- file.path(ROOT,"paper 6/analysis/g1000_eur/g1000_eur")
OUT   <- file.path(RES,"cause_negcontrol_pairs.csv")

NODES <- list(
 CAD=list(f="paper 6/analysis/data/outcomes/CAD_GCST90132314.h.tsv.gz", snp="rsid", ea="effect_allele", oa="other_allele", b="beta", se="standard_error", p="p_value"),
 T2D=list(f="paper 6/analysis/data/outcomes/T2D_GCST006867.h.tsv.gz",   snp="hm_rsid", ea="effect_allele", oa="other_allele", b="beta", se="standard_error", p="p_value"),
 HDL=list(f="paper 4/independent_build/data/raw/HDL_GLGC2021_EUR.gz",    snp="rsID", ea="ALT", oa="REF", b="EFFECT_SIZE", se="SE", p="pvalue")
)
EDGES <- list(c("HDL","CAD"), c("HDL","T2D"))   # both expected SHARING (real confounded-not-causal)

read_node <- function(n){
  cfg <- NODES[[n]]
  dt <- fread(cmd=sprintf('gzip -dc "%s"', file.path(ROOT,cfg$f)),
              select=c(cfg$snp,cfg$ea,cfg$oa,cfg$b,cfg$se,cfg$p), showProgress=FALSE)
  setnames(dt, c(cfg$snp,cfg$ea,cfg$oa,cfg$b,cfg$se,cfg$p), c("snp","A1","A2","beta","se","p"))
  dt <- dt[!is.na(beta)&!is.na(se)&se>0&grepl("^rs",snp)]; dt[,A1:=toupper(A1)][,A2:=toupper(A2)]
  unique(dt, by="snp")
}
done <- if(file.exists(OUT)) fread(OUT) else NULL
done_edges <- if(!is.null(done)) paste(done$exposure,done$outcome) else character(0)

for(ed in EDGES){
  E<-ed[1]; O<-ed[2]; key<-paste(E,O)
  if(key %in% done_edges){ cat(key,"done, skip\n"); next }
  cat("\n=== CAUSE negctrl", E, "->", O, "===\n"); flush.console()
  row <- tryCatch({
    X <- gwas_merge(read_node(E), read_node(O), snp_name_cols=c("snp","snp"),
                    beta_hat_cols=c("beta","beta"), se_cols=c("se","se"), A1_cols=c("A1","A1"), A2_cols=c("A2","A2"))
    cat("  merged:", nrow(X), "\n"); flush.console()
    params <- est_cause_params(X, sample(X$snp, min(2e5, nrow(X))))
    p1vec <- 2*pnorm(-abs(X$beta_hat_1/X$seb1))
    assoc <- file.path(MV, sprintf("nc_assoc_%s_%s.txt", E,O)); fwrite(data.table(SNP=X$snp, P=p1vec), assoc, sep="\t")
    pref <- file.path(MV, sprintf("nc_clump_%s_%s", E,O))
    system2(PLINK, c("--bfile",shQuote(PANEL),"--clump",shQuote(assoc),"--clump-p1","1e-4","--clump-r2","0.01",
      "--clump-kb","10000","--clump-snp-field","SNP","--clump-field","P","--out",shQuote(pref)), stdout=FALSE, stderr=FALSE)
    pruned <- fread(paste0(pref,".clumped"))$SNP; pruned <- pruned[grepl("^rs",pruned)]
    cat("  pruned:", length(pruned), "\n"); flush.console()
    r <- cause(X=X, variants=pruned, param_ests=params)
    saveRDS(r, file.path(MV, sprintf("cause_obj_%s_%s.rds", E,O)))
    s <- summary(r); z<-as.numeric(s$z); pcz<-as.numeric(s$p); qc<-s$quants[[2]]
    verdict <- if(!is.na(pcz)&&pcz<0.05&&z<0&&((qc[2,"gamma"]>0)==(qc[3,"gamma"]>0))) "CAUSAL(fail)" else "SHARING(pass)"
    cat(sprintf("  z=%+.2f p=%.3g gamma=%.3f[%.3f,%.3f] q=%.2f -> %s\n",
        z,pcz,qc[1,"gamma"],qc[2,"gamma"],qc[3,"gamma"],qc[1,"q"],verdict)); flush.console()
    dt <- data.table(exposure=E,outcome=O,n_merged=nrow(X),n_pruned=length(pruned),
      z_sharing_vs_causal=z,p=pcz,gamma_med=qc[1,"gamma"],gamma_lo=qc[2,"gamma"],gamma_hi=qc[3,"gamma"],
      q_med=qc[1,"q"],expected="SHARING",verdict=verdict)
    dt
  }, error=function(e){cat("  fail:",conditionMessage(e),"\n"); NULL})
  if(!is.null(row)){ fwrite(row, OUT, append=file.exists(OUT)); }   # PER-EDGE checkpoint
}
cat("\nDone ->", OUT, "\n")
