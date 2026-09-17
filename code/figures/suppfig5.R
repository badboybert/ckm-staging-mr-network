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
# SupplFig5 — UK-Biobank-free outcome-side sensitivity: headline edges are not sample-overlap artifacts.
# Sources: results/forward_local_edges.csv (PRIMARY, UKB-containing outcome) +
#          results/sensitivity_noukb_edges.csv (UKB-free outcome GWAS + fully-independent FinnGen R12).
# Per edge: IVW beta (95% CI) for PRIMARY vs UKB-free-outcome vs FinnGen. Filled point = passes
# network Bonferroni (IVW p < 0.05/132 = 3.8e-4). Free x per facet (edges span per-mmHg to log-OR).
source(file.path(P4_BASE, "figures", "fig_setup.R"))
suppressMessages(library(dplyr))

prim <- rd("forward_local_edges.csv")
sens <- rd("sensitivity_noukb_edges.csv")
# DERIVE from the committed edge table (was hardcoded 0.05/111, the superseded denominator).
BONF <- 0.05/nrow(read.csv(file.path(RES, "forward_local_edges.csv")))

edges <- data.frame(
  exp = c("BMI","LDL","TC","TG","SBP", "BMI","HbA1c","CAD","HF", "CAD","SBP","BMI", "SBP","BMI"),
  out = c("CAD","CAD","CAD","CAD","CAD", "T2D","T2D","T2D","T2D", "HF","HF","HF", "Stroke","Stroke"),
  stringsAsFactors = FALSE)
edges$edge <- sprintf("%s -> %s", edges$exp, edges$out)

mkrow <- function(exp, out, b, se, p, source){
  data.frame(exp=exp, out=out, edge=sprintf("%s -> %s", exp, out),
             b=b, lo=b-1.96*se, hi=b+1.96*se, p=p, source=source, stringsAsFactors=FALSE)
}
rows <- list()
for (i in seq_len(nrow(edges))) {
  ex <- edges$exp[i]; ou <- edges$out[i]
  rp <- prim[prim$exposure==ex & prim$outcome==ou, ]
  rows[[length(rows)+1]] <- mkrow(ex, ou, rp$ivw_b[1], rp$ivw_se[1], rp$ivw_p[1], "Primary (incl. UKB)")
  ss <- sens[sens$exposure==ex & sens$outcome==ou, ]
  fg <- ss[ss$provider=="FinnGen", ]
  uk <- ss[ss$provider!="FinnGen", ]
  if (nrow(uk)) rows[[length(rows)+1]] <- mkrow(ex, ou, uk$ivw_b[1], uk$ivw_se[1], uk$ivw_p[1], "UKB-free outcome")
  if (nrow(fg)) rows[[length(rows)+1]] <- mkrow(ex, ou, fg$ivw_b[1], fg$ivw_se[1], fg$ivw_p[1], "FinnGen (independent)")
}
d <- do.call(rbind, rows)
d$source <- factor(d$source, levels = rev(c("Primary (incl. UKB)","UKB-free outcome","FinnGen (independent)")))
d$sig    <- d$p < BONF
SRC_COL <- c("Primary (incl. UKB)"="#0072B2", "UKB-free outcome"="#009E73", "FinnGen (independent)"="#CC79A7")

# write source data FIRST, with the ASCII "X->Y" edge KEY intact; the strip label is typeset below
# (disp_edge) and must not leak back into a machine-readable column.
write.csv(d[,c("edge","source","b","lo","hi","p","sig")],
          file.path(FIG,"suppfig_data","suppfig5_source.csv"), row.names=FALSE)
d$edge <- factor(d$edge, levels = edges$edge, labels = disp_edge(edges$edge))

sf5 <- ggplot(d, aes(b, source, colour=source)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(aes(shape=sig), size=1.9, fill="white", stroke=0.7) +
  facet_wrap(~edge, scales="free_x", ncol=4) +
  scale_colour_manual(values=SRC_COL, name=NULL, breaks=c("Primary (incl. UKB)","UKB-free outcome","FinnGen (independent)")) +
  scale_shape_manual(values=c(`TRUE`=16, `FALSE`=1), guide="none") +
  scale_x_continuous(n.breaks=4) +
  # The threshold in the axis title is DERIVED from BONF above. It was hardcoded at the superseded
  # 0.05/111 = 4.5e-4 and disagreed with both the plotted shading and this figure's own legend.
  labs(x=sprintf("IVW causal estimate (95%% CI)  |  filled = passes Bonferroni (p<%.1e)", BONF), y=NULL,
       title="Headline edges hold under outcome-side UKB removal and in independent FinnGen") +
  theme_ckm(legend="top") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.spacing=unit(6,"pt"), strip.text=element_text(size=8))

save_fig(sf5, "SupplFig5", 170, 132)
cat("SF5 done:", nrow(d), "rows across", length(unique(d$edge)), "edges;",
    sum(d$sig), "of", nrow(d), "pass Bonferroni\n")
