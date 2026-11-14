# =============================================================================
# Figure 2 — Diagnosis coverage differs 2.6-fold, treatment coverage 1.25-fold
# Backend: R (ggplot2 + patchwork)   |  Data: outputs/33_std_ci.csv
# -----------------------------------------------------------------------------
# FIGURE CONTRACT
#   Core conclusion : diagnosis coverage among people with pain spans 2.62-fold
#                     (non-overlapping 95% CI); treatment coverage among
#                     diagnosed people with pain spans only 1.25-fold.
#   Evidence chain  : (a) diagnosis, hero panel, 15-80 axis
#                     (b) treatment, contrast panel, 38-66 axis, ELSA absent
#                     (c) fold-span, makes the difference in scale explicit
#   Archetype       : asymmetric mixed-modality figure
#   Export          : 180 x 108 mm; PDF/SVG vector + 600 dpi TIFF
# -----------------------------------------------------------------------------
# DELIBERATE DESIGN NOTE
#   (a) and (b) use different x scales. This is declared in both axis titles and
#   in the caption; panel (c) supplies the un-visualised ratio so the scale
#   difference cannot be mistaken for effect size.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr)
  library(readr); library(scales)
})

invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))

root <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("--file=", "", a[grep("--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f), "..")) else normalizePath("..")
})
src  <- file.path(root, "outputs", "33_std_ci.csv")
outd <- file.path(root, "manuscript", "figures")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

## ---- data -----------------------------------------------------------------
d <- read_csv(src, show_col_types = FALSE)

panel_lab <- c(A = "Frequent pain", B = "Joint pain")
cohort_lab <- c(HRS = "HRS (USA)", SHARE = "SHARE (Europe)", ELSA = "ELSA (UK)",
                MHAS = "MHAS (Mexico)", CHARLS = "CHARLS (China)",
                KLoSA = "KLoSA (Korea)", LASI = "LASI (India)")

# fixed row geometry: 8 data rows + 2 group headers, shared by (a) and (b)
# the gap between MHAS and SHARE(B) hosts the panel-B header
row_y_A <- c(HRS = 8.6, SHARE = 7.6, ELSA = 6.6, MHAS = 5.6)
row_y_B <- c(SHARE = 3.6, CHARLS = 2.6, KLoSA = 1.6, LASI = 0.6)
all_y   <- c(row_y_A, row_y_B)
Y_LIM   <- c(0.15, 9.95)

mk <- function(tier, metric) {
  d |> filter(tier == !!tier, metric == !!metric) |>
    mutate(y = ifelse(tier == "A", row_y_A[cohort], row_y_B[cohort]),
           tier_lab = panel_lab[tier])
}

dx <- bind_rows(mk("A", "dx_given_pain"), mk("B", "dx_given_pain"))
tx <- bind_rows(mk("A", "tx_given_dx_pain"), mk("B", "tx_given_dx_pain")) |>
  mutate(analgesic = cohort == "SHARE")

## ---- theme ----------------------------------------------------------------
pal <- c(A = "#1F4E79", B = "#6E9BC5")
sep <- 4.55  # dotted rule between the two panels, inside the gap

th <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line       = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks      = element_line(linewidth = 0.35, colour = "black"),
    axis.text       = element_text(size = 6.4, colour = "black"),
    axis.title      = element_text(size = 6.8, colour = "black"),
    plot.title      = element_text(size = 7.2, face = "bold", hjust = 0),
    plot.subtitle   = element_text(size = 6.2, colour = "grey25", hjust = 0),
    plot.margin     = margin(4, 6, 4, 4),
    panel.grid      = element_blank(),
    legend.position = "none"
  )

y_labs <- c("HRS (USA)", "SHARE (Europe, w5–8)", "ELSA (UK)", "MHAS (Mexico)",
            "SHARE (Europe, w1/2/4)", "CHARLS (China)", "KLoSA (Korea)", "LASI (India)")

## ---- panel (a) ------------------------------------------------------------
# value column at a fixed x, so labels never touch error-bar caps
XL_A <- c(15, 88); XLAB_A <- 86.5
pa <- ggplot(dx, aes(x = est, y = y, colour = tier)) +
  geom_hline(yintercept = sep, linewidth = 0.25, linetype = "dotted", colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.42) +
  geom_point(size = 1.45) +
  geom_text(aes(x = XLAB_A, label = sprintf("%.1f", est)),
            hjust = 1, size = 2.1, family = "Arial", colour = "grey15") +
  annotate("text", x = XL_A[1] + 0.5, y = 9.55, label = "Frequent pain (panel A)",
           hjust = 0, size = 2.3, fontface = "bold", family = "Arial") +
  annotate("text", x = XL_A[1] + 0.5, y = 4.98, label = "Joint pain (panel B)",
           hjust = 0, size = 2.3, fontface = "bold", family = "Arial") +
  scale_colour_manual(values = pal) +
  scale_y_continuous(breaks = all_y, labels = y_labs,
                     limits = Y_LIM, expand = c(0, 0)) +
  scale_x_continuous(limits = XL_A, breaks = seq(20, 80, 10), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Diagnosed among people with pain (%)   [axis 15–80]", y = NULL,
       title = "a   Diagnosis coverage", subtitle = "weighted, age-standardised") +
  th

## ---- panel (b) ------------------------------------------------------------
XL_B <- c(38, 74); XLAB_B <- 72.5
pb <- ggplot(tx, aes(x = est, y = y, colour = tier)) +
  geom_hline(yintercept = sep, linewidth = 0.25, linetype = "dotted", colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.42) +
  geom_point(aes(shape = analgesic), size = 1.45, stroke = 0.55) +
  geom_text(aes(x = XLAB_B, label = sprintf("%.1f", est)),
            hjust = 1, size = 2.1, family = "Arial", colour = "grey15") +
  annotate("text", x = XL_B[1] + 0.4, y = 6.6, label = "not available",
           hjust = 0, size = 2.05, colour = "grey35", family = "Arial",
           fontface = "italic") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1)) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(breaks = all_y, labels = rep("", 8),
                     limits = Y_LIM, expand = c(0, 0)) +
  scale_x_continuous(limits = XL_B, breaks = seq(40, 65, 5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Treated among diagnosed with pain (%)   [axis 38–66]", y = NULL,
       title = "b   Treatment coverage",
       subtitle = "open symbol: analgesic only") +
  th

## ---- panel (c) ------------------------------------------------------------
span <- function(df, step_name) {
  df |> group_by(tier) |>
    summarise(fold = max(est) / min(est), .groups = "drop") |>
    mutate(step = step_name)
}
fc <- bind_rows(span(dx, "Diagnosis"), span(tx, "Treatment")) |>
  mutate(y = case_when(step == "Diagnosis" & tier == "A" ~ 3.15,
                       step == "Diagnosis" & tier == "B" ~ 2.25,
                       step == "Treatment" & tier == "A" ~ 1.35,
                       step == "Treatment" & tier == "B" ~ 0.45))

pc <- ggplot(fc, aes(y = y)) +
  geom_segment(aes(x = 1, xend = fold, yend = y, colour = step),
               linewidth = 2.0, lineend = "butt") +
  geom_text(aes(x = fold, label = sprintf("%.2f×", fold), colour = step),
            hjust = -0.28, size = 2.25, family = "Arial") +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey35") +
  scale_colour_manual(values = c(Diagnosis = "#1F4E79", Treatment = "#5B8CB8")) +
  scale_y_continuous(
    breaks = fc$y,
    labels = c("Diagnosis · frequent pain", "Diagnosis · joint pain",
               "Treatment · frequent pain", "Treatment · joint pain"),
    limits = c(-0.05, 3.65), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1, 3.05), breaks = c(1, 1.5, 2, 2.5, 3),
                     expand = c(0, 0)) +
  labs(x = "Fold-span (highest ÷ lowest cohort)", y = NULL,
       title = "c   How wide is the gap at each step",
       subtitle = "ratio of the highest to the lowest cohort") +
  th + theme(axis.text.y = element_text(size = 6.0))

## ---- assemble -------------------------------------------------------------
fig <- (pa | pb) / pc + plot_layout(heights = c(3.0, 1.45), widths = c(1.05, 1))

w_mm <- 180; h_mm <- 108
save_pub_r <- function(plot, filename, width_mm = 180, height_mm = 108, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi)
  print(plot); dev.off()
}
save_pub_r(fig, file.path(outd, "Fig2_coverage_forest"), w_mm, h_mm, 600)

ragg::agg_png(file.path(outd, "Fig2_preview.png"),
              width = w_mm / 25.4, height = h_mm / 25.4, units = "in", res = 400)
print(fig); dev.off()

write_csv(bind_rows(
  transmute(dx, panel = tier, step = "diagnosis", cohort, n, est, lo, hi),
  transmute(tx, panel = tier, step = "treatment", cohort, n, est, lo, hi)
), file.path(outd, "Fig2_source_data.csv"))

cat("EXPORTED\n"); print(fc)
