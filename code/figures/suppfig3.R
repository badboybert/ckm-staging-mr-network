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
# SupplFig3 — MR scatter + funnel for the headline edges (CAD->HF, BMI->CAD, LDL->CAD).
# Re-harmonises the EXACT pipeline inputs (data/instruments/<EXP>.clumped.tsv +
# data/harmonised/<EXP>__<OUT>.outcome.tsv) with TwoSampleMR::harmonise_data(action=2), as in scripts/06_mr.R.
# Top row: SNP-exposure vs SNP-outcome scatter with IVW / MR-Egger / weighted-median slopes.
# Bottom row: funnel plot (single-SNP Wald ratio vs instrument strength 1/SE) for directional pleiotropy.
source(file.path(P4_BASE, "figures", "fig_setup.R"))
suppressMessages(library(TwoSampleMR))
INSTR <- file.path(BASE,"data/instruments"); HARM <- file.path(BASE,"data/harmonised")

fmt_exp <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")
fmt_out <- function(f) format_data(read.delim(f, stringsAsFactors=FALSE), type="outcome",
  snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="effect_allele",
  other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="N", phenotype_col="Phenotype")

EDGES <- list(c("CAD","HF"), c("BMI","CAD"), c("LDL","CAD"))
raw <- rd("forward_local_edges.csv")

scat <- list(); slp <- list(); fun <- list(); funref <- list(); chk <- list()
for (ed in EDGES) {
  exp <- ed[1]; out <- ed[2]; tag <- sprintf("%s -> %s", exp, out)
  e <- fmt_exp(file.path(INSTR, paste0(exp,".clumped.tsv")))
  o <- fmt_out(file.path(HARM, paste0(exp,"__",out,".outcome.tsv")))
  h <- harmonise_data(e, o, action=2); h <- h[h$mr_keep, ]
  m  <- mr(h, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  pl <- mr_pleiotropy_test(h)
  ss <- mr_singlesnp(h)
  b_ivw <- m$b[m$method=="Inverse variance weighted"]
  b_egg <- m$b[m$method=="MR Egger"]; i_egg <- pl$egger_intercept[1]
  b_wm  <- m$b[m$method=="Weighted median"]
  # cross-check vs the committed edge table
  ref <- raw[raw$exposure==exp & raw$outcome==out, ]
  chk[[tag]] <- sprintf("%-12s IVW re-run b=%+.4f (nSNP=%d) | table b=%+.4f nSNP=%d | dIVW=%.2e",
                        tag, b_ivw, nrow(h), ref$ivw_b[1], ref$nsnp[1], abs(b_ivw-ref$ivw_b[1]))
  # scatter (orient exposure betas positive, as in mr_scatter_plot)
  fl <- sign(h$beta.exposure)
  scat[[tag]] <- data.frame(edge=tag, bx=h$beta.exposure*fl, by=h$beta.outcome*fl,
                            sx=h$se.exposure, sy=h$se.outcome)
  slp[[tag]] <- data.frame(edge=tag,
                           method=c("IVW","MR-Egger","Weighted median"),
                           slope=c(b_ivw, b_egg, b_wm), intercept=c(0, i_egg, 0))
  # funnel: per-SNP Wald ratio (b) vs 1/se
  sp <- ss[!grepl("^All", ss$SNP), ]
  fun[[tag]] <- data.frame(edge=tag, biv=sp$b, inv_se=1/sp$se)
  funref[[tag]] <- data.frame(edge=tag, method=c("IVW","MR-Egger"), x=c(b_ivw, b_egg))
}
cat("== IVW re-harmonise cross-check vs forward_local_edges.csv ==\n"); for (s in chk) cat(s, "\n")

scatD <- do.call(rbind, scat); slpD <- do.call(rbind, slp)
funD  <- do.call(rbind, fun);  funR <- do.call(rbind, funref)
lev <- sapply(EDGES, function(x) sprintf("%s -> %s", x[1], x[2]))
# Source data FIRST, with the ASCII "X -> Y" edge KEY intact: fig3.R matches on it. Only the facet
# strip is typeset (disp_edge), so the printed panel and the machine-readable key can differ safely.
write.csv(scatD, file.path(FIG,"suppfig_data","suppfig3_scatter.csv"), row.names=FALSE)
write.csv(slpD,  file.path(FIG,"suppfig_data","suppfig3_slopes.csv"),  row.names=FALSE)
scatD$edge <- factor(scatD$edge, levels=lev, labels=disp_edge(lev))
slpD$edge  <- factor(slpD$edge,  levels=lev, labels=disp_edge(lev))
funD$edge  <- factor(funD$edge,  levels=lev, labels=disp_edge(lev))
funR$edge  <- factor(funR$edge,  levels=lev, labels=disp_edge(lev))

MCOL <- c("IVW"="#B2182B", "MR-Egger"="#4393C3", "Weighted median"="#E69F00")

p_scat <- ggplot(scatD, aes(bx, by)) +
  geom_hline(yintercept=0, linewidth=0.2, colour="grey80") + geom_vline(xintercept=0, linewidth=0.2, colour="grey80") +
  geom_errorbar(aes(ymin=by-sy, ymax=by+sy), width=0, linewidth=0.2, colour="grey70") +
  geom_errorbarh(aes(xmin=bx-sx, xmax=bx+sx), height=0, linewidth=0.2, colour="grey70") +
  geom_point(size=0.7, colour="grey30", alpha=0.7) +
  geom_abline(data=slpD, aes(slope=slope, intercept=intercept, colour=method), linewidth=0.6) +
  facet_wrap(~edge, scales="free", nrow=1) +
  scale_colour_manual(values=MCOL, name=NULL) +
  labs(x="SNP effect on exposure", y="SNP effect on outcome",
       title="Per-SNP scatter with IVW / MR-Egger / weighted-median slopes") +
  theme_ckm(legend="top")

p_fun <- ggplot(funD, aes(biv, inv_se)) +
  geom_vline(data=funR, aes(xintercept=x, colour=method), linewidth=0.6) +
  geom_point(size=0.8, colour="grey30", alpha=0.7) +
  facet_wrap(~edge, scales="free", nrow=1) +
  scale_colour_manual(values=c("IVW"="#B2182B","MR-Egger"="#4393C3"), name=NULL) +
  # plotmath, not the R variable name "beta_iv", which shipped on this axis through round 4.
  labs(x=expression("Single-SNP Wald ratio (" * beta[IV] * ")"), y="Instrument strength  1 / SE",
       title="Funnel plot: symmetry indicates no directional pleiotropy") +
  theme_ckm(legend="top")

sf3 <- (p_scat / p_fun) +
  plot_annotation(tag_levels="a", theme=theme(plot.tag=element_text(size=FS_TAG, face="bold")))
save_fig(sf3, "SupplFig3", 170, 143)
cat("SF3 done\n")
