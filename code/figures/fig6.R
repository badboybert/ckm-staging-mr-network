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
# numbers are unchanged; only the output name is, so it renders as SupplFig7.
# Supplementary Figure S7 — The atherogenic lipoprotein axis, and a tiered synthesis.
# a atherogenic-axis MVMR (ApoB survives) · b HDL is largely axis-confounding (marginal->conditional)
# c ApoB/LDL not separable (conditional-F collapse) · d full-density lipid MVMR · e tiered-evidence synthesis
source(file.path(P4_BASE, "figures", "fig_setup.R"))
ma <- rd("mvmr_apob.csv")
getd <- function(model, exp, col){ v <- ma[ma$model==model & ma$exposure==exp, col]; if(length(v)) v[1] else NA }

# ---------- 6a. atherogenic-axis MVMR: CAD ~ ApoB + HDL + TG ----------
ax <- data.frame(exp=c("ApoB","HDL","TG"),
                 b=c(getd("ApoB_HDL_TG","ApoB","direct_b"), getd("ApoB_HDL_TG","HDL","direct_b"), getd("ApoB_HDL_TG","TG","direct_b")),
                 se=c(getd("ApoB_HDL_TG","ApoB","direct_se"), getd("ApoB_HDL_TG","HDL","direct_se"), getd("ApoB_HDL_TG","TG","direct_se")),
                 cf=c(getd("ApoB_HDL_TG","ApoB","cond_F"), getd("ApoB_HDL_TG","HDL","cond_F"), getd("ApoB_HDL_TG","TG","cond_F")),
                 p=c(getd("ApoB_HDL_TG","ApoB","direct_p"), getd("ApoB_HDL_TG","HDL","direct_p"), getd("ApoB_HDL_TG","TG","direct_p")))
ax$lo<-ax$b-1.96*ax$se; ax$hi<-ax$b+1.96*ax$se
ax$exp <- factor(ax$exp, levels=rev(c("ApoB","HDL","TG")))
ax$surv <- ifelse(ax$p<0.05, "direct effect (P<0.05)", "not significant")
pa <- ggplot(ax, aes(b, exp, colour=surv)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2.4) +
  geom_text(aes(label=sprintf("%.2f (F=%.0f)", b, cf)), vjust=-1, size=FS/.pt-1.2, colour="black") +
  scale_colour_manual(values=c("direct effect (P<0.05)"=POS_COL, "not significant"=NULL_COL), name=NULL) +
  labs(x="Direct effect on CAD", y=NULL, title="ApoB represents the atherogenic axis (CAD ~ ApoB+HDL+TG)") +
  coord_cartesian(xlim=c(-0.2,0.45)) + theme_ckm(legend="bottom")

# ---------- 6b. HDL is largely (not fully) axis-confounding ----------
# Read every value; none is typed here. The 630-SNP row is the source of the main-text claim that a
# residual HDL effect survives conditioning, and it previously existed only as a literal.
fdz  <- rd("mvmr_lipid_fulldensity.csv")
getf <- function(exp, col){ v <- fdz[fdz$exposure==exp, col]; if(length(v)) as.numeric(v[1]) else NA }
.mb  <- as.numeric(getd("ApoB_HDL_TG","HDL","total_b"))
.mp  <- as.numeric(getd("ApoB_HDL_TG","HDL","total_p"))
.mse <- abs(.mb) / qnorm(.mp/2, lower.tail=FALSE)     # marginal SE back from the reported b and P
hd <- data.frame(
  model=c("Univariable\n(marginal)","MVMR 124-SNP\n(ApoB-restricted)","MVMR 630-SNP\n(full density)"),
  b =c(.mb, as.numeric(getd("ApoB_HDL_TG","HDL","direct_b")),  getf("HDL","direct_b")),
  se=c(.mse, as.numeric(getd("ApoB_HDL_TG","HDL","direct_se")), getf("HDL","direct_se")),
  note=c("","lower precision","residual survives"))
hd$lo<-hd$b-1.96*hd$se; hd$hi<-hd$b+1.96*hd$se
hd$model <- factor(hd$model, levels=rev(hd$model))
hd$sig <- ifelse(hd$hi<0, "P<0.05","n.s.")
pb <- ggplot(hd, aes(b, model, colour=sig)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2.4) +
  geom_text(aes(label=note), vjust=-1.1, size=FS/.pt-1.4, colour="grey35") +
  scale_colour_manual(values=c("P<0.05"=NEG_COL, "n.s."=NULL_COL), name=NULL) +
  labs(x="HDL → CAD effect", y=NULL, title="HDL→CAD attenuates ~50%, but does not vanish") +
  coord_cartesian(xlim=c(-0.4,0.05)) + theme_ckm(legend="bottom") +
  theme(axis.text.y=element_text(size=8))

# ---------- 6c. ApoB/LDL not separable (conditional-F collapse) ----------
# Both bar-pairs are now READ from mvmr_apob.csv (already loaded as `ma`). They were hand-typed as
# 57.7 / 9.2 / 13.6 — correct against source, but unprotected. The second pair is the ApoB_CAD model,
# in which ApoB and LDL are co-modelled ALONGSIDE HDL and TG; the old label "ApoB + LDL (joint model)"
# collided with the separate ApoB_LDL model in the same file, so it now names the model it plots.
cf <- data.frame(model=c("ApoB as axis\n(ApoB+HDL+TG)","ApoB + LDL co-modelled\n(+HDL+TG)"),
                 ApoB=c(getd("ApoB_HDL_TG","ApoB","cond_F"), getd("ApoB_CAD","ApoB","cond_F")),
                 LDL =c(NA,                                  getd("ApoB_CAD","LDL","cond_F")))
stopifnot(sum(is.finite(unlist(cf[, c("ApoB","LDL")]))) == 3L)   # 3 bars plotted, one cell is NA by design
cfl <- reshape(cf, direction="long", varying=c("ApoB","LDL"), v.names="condF", times=c("ApoB","LDL"), timevar="exp")
cfl <- cfl[!is.na(cfl$condF),]
cfl$model <- factor(cfl$model, levels=cf$model)   # levels follow the data; a typed copy silently NA'd them
pc <- ggplot(cfl, aes(model, condF, fill=exp)) +
  geom_col(position=position_dodge(width=0.7), width=0.6) +
  geom_hline(yintercept=10, linetype="dashed", colour=POS_COL, linewidth=0.4) +
  annotate("text", x=0.55, y=64, label="F = 10 (weak-instrument threshold)", size=FS/.pt-1.4,
           colour="grey35", hjust=0) +
  geom_text(aes(label=sprintf("%.0f", condF)), position=position_dodge(width=0.7), vjust=-0.4, size=FS/.pt-1.2) +
  scale_fill_manual(values=c(ApoB="#5D4037", LDL="#8D6E63"), name=NULL) +
  labs(x=NULL, y="Conditional F", title="ApoB & LDL collinear → not separable") +
  coord_cartesian(ylim=c(0,66)) + theme_ckm(legend="bottom") +
  theme(axis.text.x=element_text(size=8))

# ---------- 6d. full-density lipid MVMR (630 SNPs) ----------
.e <- c("LDL","HDL","TG")
fd <- data.frame(exp=.e,
                 b =vapply(.e, function(x) getf(x,"direct_b"),  numeric(1)),
                 se=vapply(.e, function(x) getf(x,"direct_se"), numeric(1)),
                 cf=vapply(.e, function(x) getf(x,"cond_F"),    numeric(1)),
                 p =vapply(.e, function(x) getf(x,"direct_p"),  numeric(1)),
                 row.names=NULL)
fd$lo<-fd$b-1.96*fd$se; fd$hi<-fd$b+1.96*fd$se
fd$exp <- factor(fd$exp, levels=rev(c("LDL","HDL","TG")))
fd$sign <- ifelse(fd$b>0,"pos","neg")
pd <- ggplot(fd, aes(b, exp, colour=sign)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2.4) +
  geom_text(aes(label=sprintf("%.2f (F=%.0f)", b, cf)), vjust=-1, size=FS/.pt-1.2, colour="black") +
  scale_colour_manual(values=c(pos=POS_COL, neg=NEG_COL), guide="none") +
  labs(x="Direct effect on CAD", y=NULL, title="Full-density lipid MVMR (630 instruments)") +
  coord_cartesian(xlim=c(-0.25,0.6)) + theme_ckm()

# ---------- 6e. tiered-evidence synthesis ----------
# 2026-07-22: this strip is RENDERED TEXT and had drifted from the analysis. Four rows were retired:
#   "T2D reaches HF only via CAD"          -> E7 shows 68-79% mediated with a wide CI ("predominantly")
#   "P=0.002" for the staging concordance  -> now READ from fig1_staging_stat.csv (exact enumeration;
#                                             1.2e-3 was itself superseded by 1.5e-3 on 2026-07-24)
#   "topology is portable ... (tau=0.49)"  -> H3's node-block CI includes chance; tau/binomial retired,
#                                             and portability moves from ASSERT to DEVELOP
#   "SBP->CKD is a replicated ancestry divergence" -> E6 interaction P=0.43, a power gap, not a divergence
# The tiers below mirror the Results synthesis paragraph exactly.
.st  <- read.csv(file.path(FIG, "data", "fig1_staging_stat.csv"), stringsAsFactors=FALSE)
.pp  <- as.numeric(.st$value[.st$metric=="perm_p"])
# Every statistic in this synthesis strip is now READ, not typed. Only .pp was derived before, so when
# the staging null switched to exact enumeration this panel kept printing P=0.0012 against a canonical
# 0.0015 — and the concordance beside it was a bare "0.93" literal that no check could have caught.
.conc <- as.numeric(.st$value[.st$metric=="observed_concordance"])
.med  <- read.csv(file.path(FIG, "data", "fig2_mediation.csv"), stringsAsFactors=FALSE)
.bmi_direct <- 100 - .med$pm_product[.med$exposure=="BMI"]
# 2026-07-26: this row USED to print a T2D mediated range of 68-79%, taken from the product- and
# difference-method POINT estimates in fig2_mediation.csv. Those are the pre-covariance-sweep values.
# Figure 2c, the Results, the legends and the abstract all now report the covariance ENVELOPE
# (38-142% under independence; 44-116% under strong correlation) and say explicitly that no single
# proportion is reportable — so a second, narrower range on the synthesis panel contradicted every
# other surface. The number is dropped rather than re-derived: the panel's job is the evidence tier,
# and the one honest summary of that quantity is that it is imprecise.
stopifnot(is.finite(.conc), is.finite(.bmi_direct))
# Round-7 Q7 / C014: the graded-synthesis TEXT PANEL that used to be panel (e) has been promoted
# to main-text Table 1 (scripts/build_tables.py), because a six-line summary of the paper's whole
# evidence base does not belong inside a supplementary lipid figure - and the panel and the
# Results paragraph that restated it had drifted into disagreeing about the East Asian staging
# order. This figure is now the four lipid panels.

fig6 <- (pa | pb) / (pc | pd) +
  plot_layout(heights=c(1, 1)) +
  plot_annotation(tag_levels="a", theme=theme(plot.tag=element_text(size=FS_TAG, face="bold")))
save_fig(fig6, "SupplFig7", 170, 150)
