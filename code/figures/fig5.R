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
# 2026-07-25: this figure MOVED to the supplement (round-2 item B-2 #15). Panels and
# numbers are unchanged; only the output name is, so it renders as SupplFig6.
# Supplementary Figure S6 — Coherent, mechanistically-interpretable ancestry divergences.
# a SBP->CKD EUR-causal / EAS-null (Zheng replication) · b the ->BMI triangulation (hospital vs population)
# c per-SNP mechanism (beta-cell vs adiposity) · d T2D->BMI fixed-vs-random knife-edge · e genuine sign divergences
source(file.path(P4_BASE, "figures", "fig_setup.R"))

# Every panel below reads from figures/suppfig_data/fig5_source.csv, which prep_fig5.py derives
# from results/ and ships with the submission package. This figure previously held all five panels
# as inline literals: the numbers were right, but nothing backed them and nothing stopped them
# drifting from the analysis they came from.
S5 <- read.csv(file.path(FIG, "suppfig_data", "fig5_source.csv"), stringsAsFactors = FALSE)
pan <- function(k) S5[S5$panel == k, ]

# ---------- 5a. Zheng SBP->CKD / BMI->CKD, EUR vs EAS ----------
.a <- pan("a")
z <- data.frame(edge=.a$label, anc=.a$group, b=.a$b, se=.a$se, p=.a$p)
z$lo<-z$b-1.96*z$se; z$hi<-z$b+1.96*z$se
z$row <- factor(paste(z$edge, z$anc), levels=rev(c("SBP → CKD EUR","SBP → CKD EAS","BMI → CKD EUR","BMI → CKD EAS")))
# Category labels are RENDERED into the figure legend, so they must hold calibration invariant 1
# (MR language is "consistent with a causal effect", never a bare "causal"). Label by significance.
z$sig <- ifelse(z$p<0.05, "P < 0.05", "n.s.")
pa <- ggplot(z, aes(b, row, colour=sig)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(aes(shape=anc), size=2.3) +
  geom_text(aes(label=sprintf("P=%s", ifelse(p<0.01, sprintf("%.0e",p), sprintf("%.2f",p)))),
            vjust=-1, size=FS/.pt-1.3, colour="grey30") +
  scale_colour_manual(values=c("P < 0.05"=POS_COL, "n.s."=NULL_COL), name=NULL) +
  scale_shape_manual(values=c(EUR=16, EAS=17), name=NULL) +
  scale_y_discrete(labels=function(x) sub(" (EUR|EAS)$","",x)) +
  labs(x="Causal effect on CKD (log-OR)", y=NULL,
       title="SBP→CKD: EUR-supported, EAS underpowered (interaction n.s.)") +
  coord_cartesian(xlim=c(-0.35,0.72)) + theme_ckm(legend="bottom")

# ---------- 5b. ->BMI triangulation: hospital vs population ----------
.b <- pan("b")
tr <- data.frame(lab=.b$label, grp=.b$group, b=.b$b, se=.b$se)
tr$lo<-tr$b-1.96*tr$se; tr$hi<-tr$b+1.96*tr$se
tr$lab <- factor(tr$lab, levels=rev(tr$lab))
tr$grp <- factor(tr$grp, levels=c("Hospital biobank","Population cohort"))
pb <- ggplot(tr, aes(b, lab, colour=grp)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2.1) +
  scale_colour_manual(values=c("Hospital biobank"=unname(COHORT_COL["Hospital"]),
                               "Population cohort"=unname(COHORT_COL["Population"])),
                      labels=c("Hospital biobank"="Hospital","Population cohort"="Population"), name=NULL) +
  labs(x="T2D / DM → adiposity effect", y=NULL, title="→BMI does not disappear in population cohorts") +
  coord_cartesian(xlim=c(-0.16,0.03)) + theme_ckm(legend="bottom") +
  theme(axis.text.y=element_text(size=8))

# ---------- 5c. per-SNP mechanism (DM instrument effect on BMI) ----------
.c <- pan("c")
ps <- data.frame(snp=.c$row, lab=.c$label, wald=.c$b, se=.c$se, type=.c$group)
ps$lo<-ps$wald-1.96*ps$se; ps$hi<-ps$wald+1.96*ps$se
ps <- ps[order(ps$wald),]; ps$lab <- factor(ps$lab, levels=ps$lab)
pc <- ggplot(ps, aes(wald, lab, colour=type)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_vline(xintercept=-0.090, linetype="dashed", colour="grey55", linewidth=0.3) +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2) +
  scale_colour_manual(values=c("β-cell / insulin-secretion"=NEG_COL, "adiposity"=POS_COL), name=NULL) +
  labs(x="Per-SNP effect of T2D-raising allele on BMI", y=NULL,
       title="Lean-diabetes architecture (WM = −0.09)") +
  annotate("text", x=-0.09, y=9.45, vjust=1, label="weighted\nmedian", size=FS/.pt-1.6,
           colour="grey50", lineheight=0.85) +
  coord_cartesian(xlim=c(-0.25,0.6)) + theme_ckm(legend="bottom") +
  theme(axis.text.y=element_text(size=8))

# ---------- 5d. T2D->BMI fixed vs random knife-edge ----------
.d <- pan("d")
ke <- data.frame(model=.d$label, b=.d$b, se=.d$se, p=.d$p, surv=.d$group)
ke$lo<-ke$b-1.96*ke$se; ke$hi<-ke$b+1.96*ke$se
ke$model <- factor(ke$model, levels=rev(ke$model))
pd <- ggplot(ke, aes(b, model, colour=surv)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2.6) +
  geom_text(aes(label=ifelse(p<0.01, sprintf("P=%.1e", p), sprintf("P=%.2f", p))),
            vjust=-1.1, size=FS/.pt-1.2, colour="black") +
  # 2026-07-26 (round-4 Tier 1) shortened the axis title and the legend labels because "the
  # composite width is fixed by the journal column, so length is the only lever". It was not enough:
  # round 5 (M3) found the legend STILL clipped, with the PDF text layer literally ending at
  # "Bonferroni su". Shortening is a fragile lever — it depends on a string length nobody re-measures
  # — so the two keys are STACKED instead. The legend box then needs the width of one label rather
  # than two, which removes the dependence on label length altogether.
  scale_colour_manual(values=c("survives"=POS_COL,"fails"=NULL_COL),
                      labels=c("Bonferroni fails","Bonferroni survives"), name=NULL) +
  guides(colour=guide_legend(nrow=2, byrow=TRUE)) +
  labs(x="Meta-analysed T2D→BMI", y=NULL,
       title="T2D→BMI meta: Hartung–Knapp CI spans zero (I²=0.76)") +
  coord_cartesian(xlim=c(-0.42,0.28)) + theme_ckm(legend="bottom")

# ---------- 5e. ancestry differences vs a single global scaling offset ----------
# Replaces the retired "genuine sign reversals" slopegraph. Over the 67-edge eligible family, EAS
# effects are ~0.436x EUR (a single global attenuation, dashed); the identity line is solid grey. If
# the differences were edge-specific interactions the points would scatter off the offset line, not
# lie along it. The two edges significant in BOTH ancestries (T2D->HDL, TG->LDL) are the only ones
# carried as hypotheses. All values READ from results/interaction_global_scale.csv (no literals).
gs <- read.csv(file.path(FIG, "..", "results", "interaction_global_scale.csv"), stringsAsFactors=FALSE)
.slope <- gs$global_slope[1]
gs$edge <- paste0(gs$exposure, "→", gs$outcome)
gs$class <- ifelse(gs$sig_raw=="True" & gs$p_raw < gs$bonferroni & (gs$b_eur_used>0)==(gs$b_eas>0),
                   "same-sign difference",
                   ifelse(gs$sig_raw=="True", "opposite-sign difference", "not significant"))
gs$both <- gs$exposure %in% c("T2D","TG") & gs$outcome %in% c("HDL","LDL") &
           ((gs$exposure=="T2D"&gs$outcome=="HDL")|(gs$exposure=="TG"&gs$outcome=="LDL"))
.lim <- max(abs(c(gs$b_eur_used, gs$b_eas))) * 1.05
pe <- ggplot(gs, aes(b_eur_used, b_eas)) +
  geom_hline(yintercept=0, linewidth=0.3, colour="grey85") +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey85") +
  geom_abline(slope=1, intercept=0, linewidth=0.35, colour="grey55") +
  geom_abline(slope=.slope, intercept=0, linetype="dashed", linewidth=0.5, colour=POS_COL) +
  geom_point(aes(colour=class), size=1.7, alpha=0.85) +
  geom_point(data=gs[gs$both,], shape=21, size=3, stroke=0.7, colour="black", fill=NA) +
  ggrepel::geom_text_repel(data=gs[gs$both,], aes(label=edge), size=FS/.pt-1.5,
     min.segment.length=0, segment.size=0.2, box.padding=0.5) +
  annotate("text", x=.lim*0.62, y=.lim*0.62, label="identity", size=FS/.pt-1.7, colour="grey55",
           angle=45, vjust=-0.4) +
  annotate("text", x=.lim*0.78, y=.lim*0.78*.slope, label=sprintf("EAS = %.3f × EUR", .slope),
           size=FS/.pt-1.7, colour=POS_COL, angle=atan(.slope)*180/pi, vjust=1.5) +
  scale_colour_manual(values=c("same-sign difference"="#6699CC",
                               "opposite-sign difference"=NEG_COL, "not significant"="grey75"), name=NULL) +
  labs(x="European effect", y="East Asian effect",
       title="Ancestry differences track one global offset, not edge-specific interaction") +
  coord_cartesian(xlim=c(-.lim,.lim), ylim=c(-.lim,.lim)) + theme_ckm(legend="bottom")

# ---------- 5f. adiposity DISTRIBUTION vs MASS: WHRadjBMI spares heart failure ----------
# Central fat distribution (WHRadjBMI, Pulit 2019) reaches CAD and T2D but NOT heart failure, whereas
# overall mass (BMI) reaches all three. Round 4 (Tier 1): report this as NO DETECTED WHRadjBMI-HF
# association only. WHRadjBMI is conditioned on BMI and the estimates are heterogeneous, so it cannot
# establish that HF risk is carried by mass rather than distribution — do not restore that contrast.
# Both series are read from results/ (no literals). Collider + UKB-overlap caveats live in the legend.
wh  <- rd("network_whradjbmi.csv"); fw <- rd("forward_local_edges.csv")
OUT3 <- c("CAD","T2D","HF")
w3 <- wh[wh$exposure=="WHRadjBMI" & wh$outcome %in% OUT3, c("outcome","ivw_b","ivw_se","ivw_p")]
b3 <- fw[fw$exposure=="BMI"       & fw$outcome %in% OUT3, c("outcome","ivw_b","ivw_se","ivw_p")]
w3$series <- "Central distribution (WHRadjBMI)"; b3$series <- "Overall mass (BMI)"
ad <- rbind(w3, b3)
ad$lo <- ad$ivw_b - 1.96*ad$ivw_se; ad$hi <- ad$ivw_b + 1.96*ad$ivw_se
ad$sig <- ifelse(ad$ivw_p < 0.05, "P < 0.05", "n.s.")
ad$outcome <- factor(ad$outcome, levels=rev(OUT3))
ad$series  <- factor(ad$series, levels=c("Overall mass (BMI)","Central distribution (WHRadjBMI)"))
pf <- ggplot(ad, aes(ivw_b, outcome, colour=series, shape=sig)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.5,
                position=position_dodge(width=0.55), orientation="y") +
  geom_point(size=2.3, position=position_dodge(width=0.55)) +
  # NB: half-width panel — long legend labels + a second (shape) legend overflow the canvas and are
  # clipped. Labels are kept short and the shape guide is dropped; open = n.s. is stated in the legend.
  scale_colour_manual(values=c("Overall mass (BMI)"=POS_COL,
                               "Central distribution (WHRadjBMI)"=NEG_COL),
                      labels=c("Overall mass (BMI)"="BMI (overall mass)",
                               "Central distribution (WHRadjBMI)"="WHRadjBMI (central)"), name=NULL) +
  scale_shape_manual(values=c("P < 0.05"=16, "n.s."=1), guide="none") +
  annotate("text", x=0.62, y=1.34, hjust=0.5, size=FS/.pt-1.6, colour="grey35", lineheight=0.9,
           label="no WHRadjBMI–HF\nassociation") +
  guides(colour=guide_legend(nrow=2)) +
  labs(x="Causal effect (per SD)", y=NULL) +
  coord_cartesian(xlim=c(-0.1, 1.15)) + theme_ckm(legend="bottom") +
  theme(axis.text.y=element_text(face=2))

# Panel e is a 2-position slopegraph: at FULL width it stretched to ~3:1 and flattened the EUR→EAS lines.
# It is now paired with panel f, which gives it ~91 mm — close to the ~100 mm the previous spacer row gave it
# — so the spacer hack is retired. DO NOT restore a patchwork `design` grid ("#" cells): it forces ONE shared
# column grid across all rows and clipped panel d's overhanging legend at the canvas edge. Equally, do not put
# a large plot.margin on pe — patchwork propagates panel alignment and squeezes every other panel to a sliver.
fig5 <- (pa | pb) / (pc | pd) / (pe | pf) +
  plot_layout(heights=c(1, 1.05, 0.92)) +
  plot_annotation(tag_levels="a", theme=theme(plot.tag=element_text(size=FS_TAG, face="bold")))
save_fig(fig5, "SupplFig6", 183, 232)
