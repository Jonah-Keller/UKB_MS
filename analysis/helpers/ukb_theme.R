#!/usr/bin/env Rscript
# ukb_theme.R — shared ggplot2 theme and colour palette for UKB_MS figures

# ── Stage identity colours (UMAP, trajectory lines, group labels) ─────────────
COL_PRE  <- "#CC0066"   # hot pink/magenta — pre-onset MS
COL_POST <- "#2B4C7E"   # navy            — post-onset MS (trajectories, UMAP)
COL_COMB <- "#C46200"   # amber           — combined (all) MS (UMAP; amber reads on white)
COL_HC   <- "#888888"   # mid-grey        — healthy controls

# ── Directional colours for DEP visualization (volcanos, barplots, GO, heatmaps)
# All pairs are Okabe-Ito / colorblind-safe. Post-onset avoids green/red.
#
#   Analysis     Up          Down        Notes
#   Pre-onset    pink        sky-blue    established
#   Post-onset   orange      blue        Okabe-Ito canonical colorblind-safe pair
#   Combined     yellow      blue-violet ~160° CCW from yellow
#   Pre vs Post  COL_PRE     COL_POST    stage identity colours

COL_PRE_UP    <- "#CC0066"   # hot pink  (= COL_PRE)
COL_PRE_DOWN  <- "#56B4E9"   # sky blue  (Okabe-Ito)

COL_POST_UP   <- "#E69F00"   # Okabe-Ito orange (H ≈ 42°)
COL_POST_DOWN <- "#0072B2"   # Okabe-Ito blue   (H ≈ 210°) — safe for red-green colorblindness

COL_COMB_UP   <- "#F0C030"   # golden yellow  (H ≈ 47°)
COL_COMB_DOWN <- "#5B50D6"   # blue-violet / indigo (H ≈ 247°)

# Backwards-compatibility alias (used in Figure 2 and older code)
COL_DOWN <- "#56B4E9"

# ── Compartment colours (cell-type heatmap / overrep compartment bands) ──────
COL_VASC   <- "#E67E22"   # vascular (Walchli)
COL_CNS    <- "#27AE60"   # CNS glia + neurons (HPA)
COL_IMMUNE <- "#2980B9"   # peripheral immune (HPA)

# ── Endophenotype cluster colours (Figure 5 / 5S / GMM) ─────────────────────
# Canonical C0/C1/C2 anchor colours; cluster_palette(n) extends to any k.
# Scripts needing a "None" or "All pre-onset" entry extend with c(p, "None"=...).
CLUST_COLS <- c("C0" = "#E6A817", "C1" = "#CC0066", "C2" = "#2B4C7E")

# Generate a CLUST_COLS-style named palette for k clusters.  Anchors the
# first three on the canonical UKB_MS palette (so existing MS figures are
# byte-identical) and extends with a categorical Okabe-Ito sequence.
cluster_palette <- function(n) {
    extra <- c("#117733", "#88CCEE", "#882255", "#999933", "#44AA99",
               "#AA4499", "#332288", "#DDCC77")
    full  <- c(unname(CLUST_COLS), extra)
    if (n > length(full)) {
        full <- c(full, scales::hue_pal()(n - length(full)))
    }
    setNames(full[seq_len(n)], paste0("C", seq_len(n) - 1L))
}

# ── Legacy / auxiliary ────────────────────────────────────────────────────────
UKB_PALETTE    <- c("#1A1A1A", "#E6A817", "#56B4E9", "#CC0066")
UKB_NAVY       <- "#2B4C7E"
VOLCANO_COLORS <- c("Up in MS" = COL_PRE, "Down in MS" = COL_DOWN, "NS" = "grey75")

# ── Base theme ────────────────────────────────────────────────────────────────
theme_ukb <- function(base_size = 11) {
    theme_classic(base_size = base_size) %+replace%
        theme(
            legend.position   = "bottom",
            legend.key.size   = unit(0.9, "lines"),
            strip.background  = element_blank(),
            strip.text        = element_text(face = "bold"),
            plot.title        = element_text(size = base_size, face = "plain",
                                             margin = margin(b = 4)),
            plot.subtitle     = element_text(size = base_size - 2, colour = "grey40",
                                             margin = margin(b = 6)),
            axis.line         = element_line(colour = "black", linewidth = 0.4),
            axis.ticks        = element_line(colour = "black", linewidth = 0.3)
        )
}
