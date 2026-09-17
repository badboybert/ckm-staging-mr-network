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
# Figure 4 — Suggestive, not definitive, portability to East Asians; staging clears chance only in EUR.
# a 3-cohort staging concordance · b each cohort's own permutation null · c EUR-vs-EAS effect scatter
# d forward-cascade replication incl. ZERO-OVERLAP BBJ->TPMI · e BBJ anomalies resolved by TPMI
# f H3 dependence-robust portability (edge- vs node-block bootstrap)
#
# 2026-07-22 REVISION (132-edge rewrite):
#   - panel a NO LONGER hardcodes the staging table. EUR is read from data/fig1_staging_stat.csv (the
#     132-edge network: 0.926, exact P=1.5e-3); the previous literal P=0.0018 was the SUPERSEDED 111-edge value.
#     BBJ/TPMI are read from suppfig_data/suppfig6_observed.csv (EAS networks, unaffected by EUR completion).
#   - panel b EUR null now uses the 132-edge null vector (data/fig1_permnull.csv); BBJ/TPMI keep their own.
#   - panel f REPLACES the retired node-netflow/Kendall-tau panel. Kendall tau=0.49 and the 44/65 binomial
#     assumed edge independence and are NO LONGER claimed in the Results; H3's node-block bootstrap is the
#     dependence-robust statistic and it INCLUDES chance. Data: data/fig4_portability*.csv (H3).
#   - panel d adds the zero-overlap BBJ-exposure -> TPMI-outcome series (answers the C5 overlap concern).
source(file.path(P4_BASE, "figures", "fig_setup.R"))
suppressMessages(library(dplyr))
# Round-7 C089: the panels printed raw P values (P=0.001515) while the legend and the Results print
# 1.5 × 10⁻³. One formatter, matching the manuscript's convention exactly (three decimals at or
# above 0.01, otherwise a mantissa to two significant figures and a real superscript exponent), is
# used for every P drawn in this figure, so panel and prose print the SAME string.
.sup <- function(e) {
  d <- c("⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹")
  paste0("⁻", paste(d[as.integer(strsplit(as.character(abs(e)), "")[[1]]) + 1], collapse=""))
}
.pf1 <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p >= 0.01) return(sprintf("%.3f", p))
  e <- floor(log10(p)); m <- p / 10^e
  paste0(sub("\\.0$", "", sprintf("%.1f", m)), " × 10", .sup(e))
}
.pf <- function(p) vapply(p, .pf1, character(1))
D    <- file.path(FIG, "data")
DD6  <- file.path(FIG, "suppfig_data")
ee   <- rd("eur_vs_eas_comparison.csv"); tb <- rd("tpmi_bbj_eur_comparison.csv")
xc   <- rd("network_crosscohort_bbj_tpmi.csv")

# ---------------------------------------------------------------- 4a. 3-cohort staging concordance
stat132 <- read.csv(file.path(D, "fig1_staging_stat.csv"), stringsAsFactors=FALSE)
getv <- function(k) as.numeric(stat132$value[stat132$metric==k])
obs6 <- read.csv(file.path(DD6, "suppfig6_observed.csv"), stringsAsFactors=FALSE)
eas  <- obs6[obs6$cohort %in% c("BBJ","TPMI"), c("cohort","concordance","perm_p")]
st <- rbind(data.frame(cohort="EUR", concordance=getv("observed_concordance"), perm_p=getv("perm_p")), eas)
# EUR staging P must equal the canonical value from fig1_staging_stat.csv (now the EXACT enumerated
# 3/1980 = 1.5e-3), and it must equal the value in suppfig6_observed.csv (same statistic, two files).
stopifnot(abs(st$perm_p[st$cohort=="EUR"] - as.numeric(stat132$value[stat132$metric=="perm_p"])) < 1e-9)
stopifnot(abs(st$perm_p[st$cohort=="EUR"] - obs6$perm_p[obs6$cohort=="EUR"]) < 1e-9)
st$cohort <- factor(st$cohort, levels=c("EUR","BBJ","TPMI"))
st$sig <- ifelse(st$perm_p < 0.05, "beyond chance", "not significant")
p_stage <- ggplot(st, aes(cohort, concordance, fill=sig)) +
  geom_col(width=0.62) +
  geom_hline(yintercept=0.5, linetype="dashed", colour="grey60", linewidth=0.3) +
  annotate("label", x=3.4, y=0.52, label="chance", size=FS/.pt-1.2, colour="grey55", hjust=1,
           fill="white", label.size=0, label.padding=unit(0.5,"pt")) +
  geom_text(aes(label=sprintf("%.3f\nP=%s", concordance, .pf(perm_p))), vjust=-0.3,
            size=FS/.pt-1.3, lineheight=0.9) +
  scale_fill_manual(values=c("beyond chance"="#2E7D32","not significant"=NULL_COL), name=NULL) +
  labs(x=NULL, y="Cross-stage concordance") +
  coord_cartesian(ylim=c(0,1.15)) + theme_ckm(legend="bottom")

# ---------------------------------------------------------------- 4b. each cohort's own permutation null
nulls6 <- rbind(
  data.frame(cohort="EUR",  concordance=read.csv(file.path(D,"fig1_permnull.csv"))$concordance),
  data.frame(cohort="BBJ",  concordance=read.csv(file.path(DD6,"suppfig6_null_BBJ.csv"))$concordance),
  data.frame(cohort="TPMI", concordance=read.csv(file.path(DD6,"suppfig6_null_TPMI.csv"))$concordance))
nulls6$cohort <- factor(nulls6$cohort, levels=c("EUR","BBJ","TPMI"))
ob <- st; ob$lab <- sprintf("obs=%.3f\nP=%s", ob$concordance, .pf(ob$perm_p))
# the plotted EUR null must reproduce the printed EUR P (self-check, mirrors fig1.R)
p_emp <- mean(nulls6$concordance[nulls6$cohort=="EUR"] >= ob$concordance[ob$cohort=="EUR"])
stopifnot(abs(p_emp - ob$perm_p[ob$cohort=="EUR"]) < 1e-9)
COH3 <- c(EUR=COHORT_COL[["EUR"]], BBJ=COHORT_COL[["BBJ"]], TPMI=COHORT_COL[["TPMI"]])
p_null <- ggplot(nulls6, aes(concordance)) +
  geom_histogram(aes(y=after_stat(density)), bins=34, fill="grey82", colour="white", linewidth=0.12) +
  geom_vline(data=ob, aes(xintercept=concordance, colour=cohort), linewidth=0.7) +
  geom_text(data=ob, aes(x=ifelse(concordance>0.6, concordance-0.02, concordance+0.02), y=Inf,
                         label=lab, colour=cohort),
            hjust=ifelse(ob$concordance>0.6,1,0), vjust=1.3, size=FS/.pt-1.6, lineheight=0.9, show.legend=FALSE) +
  facet_wrap(~cohort, ncol=1, strip.position="right") +
  scale_colour_manual(values=COH3, guide="none") +
  labs(x="Cross-stage concordance", y="Null density") +
  coord_cartesian(xlim=c(0.1, 1.04)) + theme_ckm() +
  theme(strip.text.y=element_text(angle=0))

# ---------------------------------------------------------------- 4c. EUR vs EAS effect scatter
# NON-COMPARABLE TRAITS EXCLUDED, on EITHER side of the edge. The E5 scale dictionary (Supp Table S15)
# marks SBP (EUR per-mmHg vs EAS per-SD), HbA1c (EUR %-units vs EAS rank-inverse-normal) and eGFR
# (EUR log-eGFR per-SD vs EAS RINT) as magnitude-non-comparable -> "sign only". An identity line is
# meaningless for any edge touching one of them. This is the SAME exclusion rule the Methods state for
# the comparable-scale effect correlation; it previously dropped only SBP-as-exposure.
# The fitted EUR->EAS global scale offset, read from the file that estimates it (script 67), so this
# panel and Supplementary Figure S6e cannot print different slopes.
.gs <- read.csv(file.path(BASE, "results", "interaction_global_scale.csv"), stringsAsFactors=FALSE)
GLOBAL_SLOPE <- unique(.gs$global_slope)
stopifnot(length(GLOBAL_SLOPE) == 1, is.finite(GLOBAL_SLOPE))

NONCOMP <- c("SBP", "HbA1c", "eGFR")
.touches <- function(edges, traits) {
  parts <- strsplit(edges, "->", fixed=TRUE)
  vapply(parts, function(p) any(p %in% traits), logical(1))
}
sc <- ee[!.touches(ee$edge, NONCOMP), ]
stopifnot(!any(.touches(sc$edge, NONCOMP)))
sc$both <- ifelse(as.character(sc$both_sig) %in% c("TRUE","True","true"), "both significant", "EUR-anchored only")
labpts <- sc[sc$edge %in% c("LDL->CAD","TC->CAD","BMI->HF","BMI->Stroke","T2D->CAD","BMI->T2D","TG->CAD"), ]
p_scatter <- ggplot(sc, aes(EUR_b, EAS_b)) +
  geom_hline(yintercept=0, linewidth=0.25, colour="grey80") + geom_vline(xintercept=0, linewidth=0.25, colour="grey80") +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey75", linewidth=0.3) +
  # The global EUR->EAS scale offset fitted in Supplementary Figure S6e. Drawing it here stops the
  # identity line reading as the expectation: it is not, and the paper says so.
  geom_abline(slope=GLOBAL_SLOPE, intercept=0, linetype="solid", colour=POS_COL, linewidth=0.45) +
  geom_point(aes(colour=both), size=1.9) +
  ggrepel::geom_text_repel(data=labpts, aes(label=disp_edge(edge)), size=FS/.pt-1.5, colour="grey25",
     segment.size=0.2, min.segment.length=0, max.overlaps=20, box.padding=0.2) +
  scale_colour_manual(values=c("both significant"=unname(COHORT_COL["BBJ"]), "EUR-anchored only"="grey62"), name=NULL) +
  labs(x="EUR causal effect", y="EAS (BBJ) causal effect") +
  annotate("text", x=-0.45, y=0.9, label="y = x", colour="grey55", size=FS/.pt-1, angle=32) +
  # Round 4 (T2-7): the solid line was visible but not self-identifying. The label is BUILT from
  # GLOBAL_SLOPE, never typed, so it cannot drift from the line it names or from Supplementary
  # Figure S6e, which reports the same fit. Gate F13 asserts the number reaches Figure4.pdf.
  # It is placed in the EMPTY upper-left quadrant, NOT riding the line: an angled label on the line
  # ran straight through the TG->CAD / BMI->HF labels and the central point cloud (checked in the PNG).
  # plotmath, not "beta_EAS", so the panel prints a real subscript rather than a variable name.
  # Round-7 C029: at 170 mm the two-line note reached the TC->CAD point label. It moves to the
  # empty lower-right quadrant, where no comparable-scale edge falls.
  annotate("text", x=-0.53, y=1.02, hjust=0, colour=POS_COL, size=FS/.pt-1.5,
           label="solid line: fitted global offset") +
  annotate("text", x=-0.53, y=0.92, hjust=0, colour=POS_COL, size=FS/.pt-1.5, parse=TRUE,
           label=sprintf("beta[EAS] == %.3f %%*%% beta[EUR]", GLOBAL_SLOPE)) +
  annotate("text", x=1.03, y=-0.5, label="comparable-scale edges only\n(SBP, HbA1c, eGFR excluded)",
           hjust=1, size=FS/.pt-1.7, colour="grey45", lineheight=0.9) +
  coord_cartesian(xlim=c(-0.55,1.05), ylim=c(-0.55,1.05)) + theme_ckm(legend="bottom")

# ---------------------------------------------------------------- 4d. cascade replication + ZERO-OVERLAP
# HbA1c->T2D is NOT plotted here: the x-axis is a per-SD magnitude axis and HbA1c is rank-inverse-normal
# in East Asians but %-units in Europeans (Supp Table S15 = "sign only"), so its cross-cohort magnitudes
# are not comparable. Its replication is reported on sign in the text, not on this axis.
cas <- c("LDL->CAD","TC->CAD","TG->CAD","BMI->CAD","BMI->HF")
stopifnot(!any(.touches(cas, NONCOMP)))
d <- tb[tb$edge %in% cas, c("edge","EUR_b","BBJ_b","TPMI_b","EUR_p","BBJ_p","TPMI_p")]
L <- do.call(rbind, lapply(c("EUR","BBJ","TPMI"), function(co)
  data.frame(edge=d$edge, cohort=co, b=as.numeric(d[[paste0(co,"_b")]]), p=as.numeric(d[[paste0(co,"_p")]]))))
xc$edge <- paste0(xc$exposure, "->", xc$outcome)
xcs <- xc[xc$edge %in% cas, ]
L <- rbind(L, data.frame(edge=xcs$edge, cohort="BBJ→TPMI", b=xcs$ivw_b, p=xcs$ivw_p))
L$cohort <- factor(L$cohort, levels=c("EUR","BBJ","TPMI","BBJ→TPMI"))
L$sig <- ifelse(!is.na(L$p) & L$p<0.05, "sig", "ns")
L$edge <- factor(L$edge, levels=rev(cas), labels=rev(disp_edge(cas)))
p_casc <- ggplot(L, aes(b, edge, colour=cohort, shape=sig)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_point(size=2, position=position_dodge(width=0.68)) +
  scale_colour_manual(values=c(COHORT_COL[c("EUR","BBJ","TPMI")], "BBJ→TPMI"="#333333"), name=NULL) +
  scale_shape_manual(values=c(sig=16, ns=1), guide="none") +
  labs(x="Causal effect (per SD)", y=NULL) +
  guides(colour=guide_legend(nrow=2)) +
  coord_cartesian(xlim=c(-0.15,1.65)) + theme_ckm(legend="bottom")

# ---------------------------------------------------------------- 4e. BBJ anomalies resolved by TPMI
ano <- c("BMI->T2D","BMI->CAD","HbA1c->HF")
da <- tb[tb$edge %in% ano, c("edge","EUR_b","BBJ_b","TPMI_b","EUR_p","BBJ_p","TPMI_p")]
La <- do.call(rbind, lapply(c("EUR","BBJ","TPMI"), function(co)
  data.frame(edge=da$edge, cohort=co, b=as.numeric(da[[paste0(co,"_b")]]), p=as.numeric(da[[paste0(co,"_p")]]))))
La$cohort <- factor(La$cohort, levels=c("EUR","BBJ","TPMI"))
La$sig <- ifelse(!is.na(La$p) & La$p<0.05, "sig", "ns")
La$edge <- factor(La$edge, levels=rev(ano), labels=rev(disp_edge(ano)))
p_anom <- ggplot(La, aes(b, edge, colour=cohort, shape=sig)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_point(size=2.4, position=position_dodge(width=0.6)) +
  scale_colour_manual(values=COHORT_COL[c("EUR","BBJ","TPMI")], name=NULL) +
  scale_shape_manual(values=c(sig=16, ns=1), name=NULL, labels=c(sig="P<0.05", ns="n.s.")) +
  # Round-7 C037: this annotation used to sit above BMI→T2D — the one row the legend and the
  # Results say is NOT resolved by the second cohort, because its TPMI estimate comes from the
  # overlap-prone within-cohort design. It belongs on BMI→CAD, the anomaly that does resolve.
  annotate("text", x=0.72, y=2.42, label="TPMI recovers the\nEuropean direction", size=FS/.pt-1.5,
           colour="grey30", lineheight=0.9) +
  annotate("text", x=0.72, y=3.42, label="not resolved:\nTPMI estimate is\nwithin-cohort only", size=FS/.pt-1.6,
           colour="grey45", lineheight=0.9) +
  labs(x="Causal effect", y=NULL) +
  coord_cartesian(xlim=c(-0.35,1.15)) + theme_ckm(legend="bottom")

# ---------------------------------------------------------------- 4f. H3 dependence-robust portability
po <- read.csv(file.path(D, "fig4_portability.csv"), stringsAsFactors=FALSE)
mt <- read.csv(file.path(D, "fig4_portability_meta.csv"), stringsAsFactors=FALSE)
gm <- function(k) mt$value[mt$stat==k]
sg <- po[po$stat=="Sign concordance", ]
# 2026-07-25: the NODE-level (dependence-robust) interval is the PRIMARY portability result and the
# edge-level one is secondary, because resampling edges independently is anticonservative when edges
# share nodes. The hierarchy is made visible in the panel rather than left to the legend: node-level
# is drawn on top (last factor level) and both rows carry their rank. verdict is computed from the raw
# block name BEFORE relabelling, so the label text cannot silently change which row is which.
sg$verdict <- ifelse(sg$block=="Node-level", "includes chance", "excludes chance")
sg$block <- factor(ifelse(sg$block=="Node-level", "Node-level\n(primary)", "Edge-level\n(secondary)"),
                   levels=c("Edge-level\n(secondary)", "Node-level\n(primary)"))
stopifnot(nlevels(sg$block)==2, !any(is.na(sg$block)))
p_port <- ggplot(sg, aes(estimate, block, colour=verdict)) +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, fill="grey94") +
  geom_vline(xintercept=0.5, linetype="dashed", linewidth=0.4, colour="grey45") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.7, orientation="y") +
  geom_point(size=2.8) +
  geom_text(aes(label=sprintf("[%.2f, %.2f]", lo, hi)), vjust=-1.25, size=FS/.pt-1.3, colour="black") +
  geom_text(aes(label=verdict), vjust=2.1, size=FS/.pt-1.6, colour="grey35") +
  annotate("text", x=0.5, y=2.62, label="chance", size=FS/.pt-1.5, colour="grey45", hjust=1.08, vjust=0) +
  annotate("text", x=1.02, y=1.42, hjust=1, size=FS/.pt-1.6, colour="grey35", lineheight=0.9,
           label=sprintf("%d/%d EUR-significant edges concordant\nsign-flip P = %.3f",
                         gm("k_concordant"), gm("n_tested"), gm("perm_p"))) +
  scale_colour_manual(values=c("excludes chance"="#2E7D32","includes chance"=POS_COL), guide="none") +
  scale_x_continuous(breaks=seq(0.4,1.0,0.2)) +
  coord_cartesian(xlim=c(0.38,1.02), ylim=c(0.5,2.85), clip="off") +
  labs(x="Edge sign concordance (95% CI)", y=NULL) +   # NB: keep SHORT — this is a half-width panel and a
                                                        # longer label is clipped at the right edge.
  theme_ckm() + theme(axis.text.y=element_text(face=2))

# ---------------------------------------------------------------- assemble
# tag order = a staging, b nulls, c scatter, d cascade(+zero-overlap), e anomalies, f portability
fig4 <- (p_stage | p_null) / (p_scatter | p_casc) / (p_anom | p_port) +
  plot_layout(heights=c(1, 1, 1)) +
  plot_annotation(tag_levels="a", theme=theme(plot.tag=element_text(size=FS_TAG, face="bold")))
save_fig(fig4, "Figure4", 170, 214)
