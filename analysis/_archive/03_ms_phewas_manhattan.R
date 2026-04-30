#!/usr/bin/env Rscript
# 03_ms_phewas_manhattan.R — PheWAS Manhattan plot for MS vs HC
#
# Style matches figure5.R panel A exactly (chapter colours, point size by OR,
# text-repel without boxes, inter-chapter gaps, grey-dashed FDR line).
# Extended to cover all 24 ICD-10 chapters now available.
#
# Output: results/comorbidity/ms_phewas_manhattan.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

PHEWAS_FILE <- file.path(PROJ_DIR, "results", "comorbidity", "ms_full_phewas.csv")
OUT_DIR     <- file.path(PROJ_DIR, "results", "comorbidity")

# ── Chapter colours — matched to figure5.R for G–Q, extended for A/B/E/F/R–Z ──
CHAP_COLS <- c(
    A = "#6D8B3A",   B = "#9DB66A",   C = "grey75",    D = "#FF9933",
    E = "#20B2AA",   F = "#7B2D8B",   G = "#2B4C7E",   H = "#5B9BD5",
    I = "#E6A817",   J = "#4CAF50",   K = "#A0522D",   L = "#CC0066",
    M = "#9370DB",   N = "#20B2AA",   O = "#DC143C",   P = "#CDDC39",
    Q = "#808080",   R = "#607D8B",   S = "#795548",   T = "#FF5722",
    U = "#00897B",   V = "#1565C0",   Y = "#BF360C",   Z = "#9E9E9E"
)
CHAP_COLS_ALT <- c(
    A = "#A8C06B",   B = "#C5D98D",   C = "grey85",    D = "#FFCC99",
    E = "#7DCCCC",   F = "#B57EC7",   G = "#5A7FAD",   H = "#9DC4E8",
    I = "#F0C55D",   J = "#88D898",   K = "#C47A5A",   L = "#E87EB9",
    M = "#B89BE5",   N = "#7DCCCC",   O = "#F08080",   P = "#E6EE9C",
    Q = "#AAAAAA",   R = "#90A4AE",   S = "#A1887F",   T = "#FF8A65",
    U = "#4DB6AC",   V = "#5C8FD6",   Y = "#E67635",   Z = "#BDBDBD"
)
CHAP_NAMES <- c(
    A = "Infections",    B = "Parasitic",    C = "Neoplasms",   D = "Blood",
    E = "Endocrine",     F = "Mental",       G = "Neurological",H = "Eye/Ear",
    I = "Circulatory",   J = "Respiratory",  K = "Digestive",   L = "Skin",
    M = "Musculoskeletal",N= "Genitourinary",O = "Obstetric",   P = "Perinatal",
    Q = "Congenital",    R = "Symptoms",     S = "Injury",      T = "Poisoning",
    U = "COVID",         V = "Transport",    Y = "External",    Z = "Health contact"
)

# ── Load ──────────────────────────────────────────────────────────────────────
cat("Loading PheWAS results...\n")
phewas <- fread(PHEWAS_FILE)
cat(sprintf("  %d conditions, %d FDR<0.05\n",
            nrow(phewas), sum(phewas$fdr < 0.05, na.rm = TRUE)))

FDR_THR <- 0.05

# Ordered by standard ICD-10 chapter sequence
chap_order <- names(CHAP_COLS)
chap_order <- chap_order[chap_order %in% unique(phewas$chapter)]
phewas[, chapter := factor(chapter, levels = chap_order)]
phewas <- phewas[!is.na(chapter)]
setorder(phewas, chapter, pval)

# ── x-positions with inter-chapter gaps (4 units) ────────────────────────────
CHAPTER_GAP <- 4L
phewas[, x_pos := seq_len(.N), by = chapter]
offset <- 0L
for (ch in chap_order[chap_order %in% phewas$chapter]) {
    idx <- phewas$chapter == ch
    phewas[idx, x_pos := x_pos + offset]
    offset <- offset + max(phewas[idx, x_pos]) + CHAPTER_GAP
}
phewas[, x_pos := as.numeric(x_pos)]

chap_mids <- phewas[, .(x_mid = mean(x_pos), n = .N), by = chapter]
chap_mids[, chap_name := CHAP_NAMES[as.character(chapter)]]

# ── FDR threshold y-position ──────────────────────────────────────────────────
fdr_line_y <- {
    ord <- phewas[order(fdr)]
    first_ns <- which(ord$fdr >= FDR_THR)[1L]
    if (is.na(first_ns)) NA_real_ else -log10(ord$pval[first_ns])
}
cat(sprintf("  FDR 0.05 threshold: -log10(p) = %.2f\n", fdr_line_y))

# ── Timing category (for bold labels) ────────────────────────────────────────
phewas[, timing_cat := fcase(
    !is.na(pct_pre_onset) & pct_pre_onset >= 50, "Pre-onset dominant",
    !is.na(pct_pre_onset) & pct_pre_onset >= 20, "Mixed",
    default = "Post-onset dominant"
)]

# ── Point styling ─────────────────────────────────────────────────────────────
phewas[, sig      := fdr < FDR_THR]
phewas[, log10p   := -log10(pval)]
phewas[, log2or   := log2(pmax(OR, 0.01))]
phewas[, pt_size  := pmin(abs(log2or), 4) * 0.7 + 0.3]
phewas[, col_idx  := as.integer(chapter) %% 2L]
phewas[, fill_col := ifelse(col_idx == 0L,
                             CHAP_COLS[as.character(chapter)],
                             CHAP_COLS_ALT[as.character(chapter)])]
phewas[sig == FALSE, fill_col := "grey82"]

# ── Labels: top 25 FDR-sig ────────────────────────────────────────────────────
label_dt <- phewas[sig == TRUE][order(pval)][seq_len(min(25L, sum(phewas$sig)))]

# ── Plot ──────────────────────────────────────────────────────────────────────
cat("Building plot...\n")

p <- ggplot(phewas, aes(x = x_pos, y = log10p)) +
    { if (!is.na(fdr_line_y))
        geom_hline(yintercept = fdr_line_y,
                   linewidth = 0.35, linetype = "dashed", colour = "grey45")
    } +
    geom_point(aes(size = pt_size, colour = fill_col),
               alpha = ifelse(phewas$sig, 0.9, 0.40), stroke = 0) +
    geom_point(data = phewas[sig == TRUE],
               aes(size = pt_size, colour = fill_col),
               alpha = 1, stroke = 0) +
    geom_text_repel(
        data        = label_dt,
        aes(label   = icd10),
        size        = 2.1,
        segment.size = 0.25,
        max.overlaps = 30,
        box.padding  = 0.3,
        colour      = "grey10",
        fontface    = ifelse(label_dt$timing_cat == "Pre-onset dominant", "bold", "plain")
    ) +
    scale_colour_identity() +
    scale_size_identity() +
    scale_x_continuous(
        breaks = chap_mids$x_mid,
        labels = chap_mids$chap_name,
        expand = expansion(add = c(2, 2))
    ) +
    annotate("text",
             x = max(phewas$x_pos) * 0.02,
             y = max(phewas$log10p) * 0.97,
             label = sprintf(
                 "n = %d tested | %d FDR<0.05\nBold = pre-onset dominant (>50%% before Dx)",
                 nrow(phewas), sum(phewas$sig)),
             size = 2.4, hjust = 0, vjust = 1, colour = "grey30") +
    labs(
        x        = "ICD-10 chapter",
        y        = expression(-log[10](italic(p))),
        title    = "MS comorbidity PheWAS (ICD-10 first-occurrence, all chapters)",
        subtitle = sprintf("MS n=499, HC n=47,801 | adjusted for age + sex | BH FDR correction")
    ) +
    theme_ukb(base_size = 9) +
    theme(
        axis.text.x        = element_text(size = 7.5, angle = 35, hjust = 1),
        panel.grid.major.x = element_blank()
    )

out_path <- file.path(OUT_DIR, "ms_phewas_manhattan.pdf")
ggsave(out_path, p, width = 13, height = 4.5, device = cairo_pdf)
cat(sprintf("Saved: %s\n", basename(out_path)))
cat("03_ms_phewas_manhattan.R complete.\n")
