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
# SupplFig1 — Full 132-edge bidirectional causal network heatmap (European ancestry).
# Complete-network companion to main Fig 1. Source: results/forward_local_edges.csv.
# Fill = signed IVW beta (diverging, centred 0, clamped +/-0.6 for legibility; exact value in cell);
# asterisk = Bonferroni-significant (IVW p < 0.05/nrow = 3.8e-4 at 132 edges).
source(file.path(P4_BASE, "figures", "fig_setup.R"))

raw <- rd("forward_local_edges.csv")
# DERIVE the threshold from the edge table rather than typing it: this read 0.05/111 and would have
# kept asterisking against the superseded denominator after the network was completed to 132 edges.
BONF <- 0.05/nrow(raw)

# CKM-role node order (AHA stage 1 -> 2 -> 4, then kidney)
ORD <- c("BMI","SBP","HbA1c","TG","TC","LDL","HDL","T2D","CAD","HF","Stroke","eGFR")

# full 12x12 grid; fill only where an edge was tested
grid <- expand.grid(exposure = ORD, outcome = ORD, stringsAsFactors = FALSE)
grid <- merge(grid, raw[, c("exposure","outcome","ivw_b","ivw_p")],
              by = c("exposure","outcome"), all.x = TRUE)
grid$exposure <- factor(grid$exposure, levels = rev(ORD))   # BMI at top
grid$outcome  <- factor(grid$outcome,  levels = ORD)
grid$sig      <- !is.na(grid$ivw_p) & grid$ivw_p < BONF
grid$lab      <- ifelse(is.na(grid$ivw_b), "", sprintf("%.2f", grid$ivw_b))
grid$star     <- ifelse(grid$sig, "*", "")
grid$diag     <- as.character(grid$exposure) == as.character(grid$outcome)

LIM <- 0.6   # clamp for the diverging fill; exact beta printed in-cell so nothing is lost

sf1 <- ggplot(grid, aes(outcome, exposure)) +
  geom_tile(aes(fill = ivw_b), colour = "grey85", linewidth = 0.3) +
  # grey out the diagonal (self-edges, never tested)
  geom_tile(data = subset(grid, diag), fill = "grey93", colour = "grey85", linewidth = 0.3) +
  geom_text(aes(label = lab), size = FS/.pt - 1.2, colour = "grey15") +
  geom_text(aes(label = star), nudge_x = 0.30, nudge_y = 0.20,
            size = FS/.pt + 1.4, colour = "black", fontface = 2) +
  scale_fill_gradient2(low = NEG_COL, mid = "white", high = POS_COL, midpoint = 0,
                       limits = c(-LIM, LIM), oob = scales::squish, na.value = "grey96",
                       name = "IVW beta\n(signed)", breaks = c(-0.6,-0.3,0,0.3,0.6),
                       labels = c("<=-0.6","-0.3","0","0.3",">=0.6")) +
  scale_x_discrete(position = "top") +
  labs(x = "Outcome", y = "Exposure",
       # DERIVED: this title read "Full 111-edge ..." and survived the completion to 132.
       title = sprintf("Full %d-edge bidirectional causal network (EUR)", nrow(raw))) +
  coord_equal() +
  theme_ckm(legend = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
        panel.grid = element_blank(),
        legend.key.height = unit(16, "pt"))

save_fig(sf1, "SupplFig1", 178, 165)
cat("SF1 done:", nrow(raw), "edges,", sum(grid$sig), "Bonferroni-significant cells\n")
