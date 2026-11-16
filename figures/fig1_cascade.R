# =============================================================================
# Figure 1 — Where the rheumatic pain cascade thins
# Backend: R (ggplot2 + patchwork)   |  Data: outputs/42_fig1_funnel.csv
# -----------------------------------------------------------------------------
# FIGURE CONTRACT
#   Core conclusion : of every 100 adults with pain, between 11 (India) and 40
#                     (USA) reach treatment; in every cohort except the USA the
#                     larger absolute loss occurs between pain and diagnosis.
#   Evidence chain  : (a) frequent pain, 3 cohorts
#                     (b) joint pain, 3 cohorts
#                     (c) SHARE twice: same cohort, different item and waves
#   Archetype       : quantitative grid, common-denominator stacked bars
#   Export          : 180 x 121 mm; PDF/SVG vector + 600 dpi TIFF
# -----------------------------------------------------------------------------
# DELIBERATE DESIGN NOTE
#   Every bar starts at 100 people with pain and splits into three mutually
#   exclusive segments, so widths add to 100 and are comparable across rows.
#   This is the only honest way to compare "loss at diagnosis" with "loss at
#   treatment": conditional percentages have different denominators and cannot
#   be chained end to end. SHARE appears only in (c), never twice in one panel.
#   Cohorts are restricted to waves in which a treatment item was fielded
#   (HRS w5-11, MHAS w3-5, CHARLS w1-4, KLoSA w3-8, LASI w1, SHARE w5-8 /
#   w1,2,4); this moves diagnosis coverage by <=1.1 pp (43_diag_txwaves.txt).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr)
  library(tidyr); library(readr); library(svglite); library(ragg)
})

invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))

root <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("--file=", "", a[grep("--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f), "..")) else normalizePath("..")
})
src  <- file.path(root, "outputs", "42_fig1_funnel.csv")
outd <- file.path(root, "outputs", "figures")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

## ---- data -----------------------------------------------------------------
wide <- read_csv(src, show_col_types = FALSE) |>
  select(tier, cohort, metric, est, lo, hi) |>
  pivot_wider(names_from = metric, values_from = c(est, lo, hi)) |>
  # ELSA has no treatment item: the Python step wrote 0 for the yield, which
  # must not be drawn as "nobody reached treatment".
  mutate(across(starts_with(c("est_tx", "lo_tx", "hi_tx")),
                ~ ifelse(cohort == "ELSA", NA_real_, .))) |>
  mutate(dx      = est_dx_given_pain,
         treated = est_tx_per100_pain,
         undx    = 100 - dx,
         dx_untx = ifelse(is.na(treated), NA_real_, dx - treated),
         notmeas = ifelse(is.na(treated), dx, NA_real_))

GEO <- list(
  a = list(hdr = 8.62, y = c(HRS = 7.85, MHAS = 6.85, ELSA = 5.85),
           lab = c(HRS = "HRS (USA)", MHAS = "MHAS (Mexico)", ELSA = "ELSA (UK)"),
           title = "(a) Frequent pain"),
  b = list(hdr = 4.72, y = c(CHARLS = 3.95, KLoSA = 2.95, LASI = 1.95),
           lab = c(CHARLS = "CHARLS (China)", KLoSA = "KLoSA (Korea)", LASI = "LASI (India)"),
           title = "(b) Joint pain"),
  c = list(hdr = 0.82, y = c(SHARE_A = 0.05, SHARE_B = -0.95),
           lab = c(SHARE_A = "Frequent pain (w5–8)", SHARE_B = "Joint pain (w1/2/4)"),
           title = "(c) SHARE: same cohort, two items")
)

rows <- bind_rows(
  wide |> filter(tier == "A", cohort %in% c("HRS", "MHAS", "ELSA")) |>
    mutate(panel = "a", key = cohort),
  wide |> filter(tier == "B", cohort %in% c("CHARLS", "KLoSA", "LASI")) |>
    mutate(panel = "b", key = cohort),
  wide |> filter(cohort == "SHARE") |>
    mutate(panel = "c", key = paste0("SHARE_", tier))
) |>
  rowwise() |>
  mutate(y = GEO[[panel]]$y[[key]],
         row_lab = GEO[[panel]]$lab[[key]]) |>
  ungroup()

## ---- long segments --------------------------------------------------------
seg_order <- c("undx", "notmeas", "dx_untx", "treated")
long <- rows |>
  select(panel, y, row_lab, undx, notmeas, dx_untx, treated) |>
  pivot_longer(all_of(seg_order), names_to = "seg", values_to = "w") |>
  filter(!is.na(w), w > 0) |>
  group_by(y) |>
  arrange(match(seg, seg_order), .by_group = TRUE) |>
  mutate(xmax = cumsum(w), xmin = xmax - w) |>
  ungroup() |>
  mutate(seg = factor(seg, levels = seg_order))

hdrs <- tibble(panel = names(GEO),
               y = vapply(GEO, function(g) g$hdr, numeric(1)),
               title = vapply(GEO, function(g) g$title, character(1)))
divs <- tibble(y = c(5.30, 1.38))

Y_LIM <- c(-2.05, 9.10)

## ---- palette / theme ------------------------------------------------------
pal <- c(undx = "#DCDCDC", notmeas = "#F2F2F2",
         dx_untx = "#6E9BC5", treated = "#1F4E79")

th <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line.x      = element_line(linewidth = 0.35, colour = "black"),
        axis.line.y      = element_blank(),
        axis.ticks.x     = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks.y     = element_blank(),
        axis.text.x      = element_text(size = 6.4, colour = "black"),
        axis.text.y      = element_text(size = 6.3, colour = "black", hjust = 1),
        axis.title.x     = element_text(size = 6.8, colour = "black"),
        panel.grid       = element_blank(),
        plot.margin      = margin(4, 2, 2, 7))

## ---- main panel -----------------------------------------------------------
p_main <- ggplot(long, aes(y = y)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - 0.30, ymax = y + 0.30,
                fill = seg), colour = NA) +
  # outline every bar so the pale blocks still read as one bar
  geom_rect(data = rows, aes(xmin = 0, xmax = 100, ymin = y - 0.30,
                             ymax = y + 0.30),
            inherit.aes = FALSE, fill = NA, colour = "grey45", linewidth = 0.22) +
  # dashed outline marks "treatment not recorded" (ELSA)
  geom_rect(data = long |> filter(seg == "notmeas"),
            aes(xmin = xmin, xmax = xmax, ymin = y - 0.30, ymax = y + 0.30),
            inherit.aes = FALSE, fill = NA, colour = "grey35",
            linetype = "22", linewidth = 0.30) +
  # diagnosis boundary
  geom_segment(data = rows, aes(x = dx, xend = dx, y = y - 0.30, yend = y + 0.30),
               inherit.aes = FALSE, linewidth = 0.32, colour = "grey10") +
  # panel headers and separators
  geom_text(data = hdrs, aes(x = 0, y = y, label = title), hjust = 0,
            size = 2.65, fontface = "bold", family = "Arial", inherit.aes = FALSE) +
  geom_segment(data = divs, aes(x = 0, xend = 100, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.22, colour = "grey70",
               linetype = "22") +
  scale_y_continuous(limits = Y_LIM, breaks = rows$y, labels = rows$row_lab,
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  scale_fill_manual(values = pal, drop = FALSE, name = NULL,
                    breaks = c("treated", "dx_untx", "notmeas", "undx"),
                    labels = c(treated = "Reached treatment",
                               dx_untx = "Diagnosed, untreated",
                               notmeas = "Diagnosed, treatment not recorded",
                               undx = "Undiagnosed")) +
  labs(x = "Per 100 people with pain", y = NULL) +
  coord_cartesian(clip = "off") +
  th

## ---- right-hand numeric columns -------------------------------------------
p_txt <- ggplot(rows) +
  geom_text(aes(x = 0, y = y + 0.13, label = sprintf("%.1f", dx)),
            hjust = 0, size = 2.25, colour = "grey15", family = "Arial") +
  geom_text(aes(x = 1.0, y = y + 0.13,
                label = ifelse(is.na(treated), "—", sprintf("%.1f", treated))),
            hjust = 0, size = 2.4, fontface = "bold", colour = "grey5",
            family = "Arial") +
  geom_text(aes(x = 1.0, y = y - 0.22,
                label = ifelse(is.na(lo_tx_per100_pain), "",
                               sprintf("(%.1f–%.1f)", lo_tx_per100_pain,
                                       hi_tx_per100_pain))),
            hjust = 0, size = 1.8, colour = "grey42", family = "Arial") +
  annotate("text", x = 0, y = 9.08, label = "Diagnosed\n(%)", hjust = 0, vjust = 1,
           size = 2.0, lineheight = 0.92, colour = "grey30", family = "Arial") +
  annotate("text", x = 1.0, y = 9.08, label = "Reached\ntreatment\n(per 100)",
           hjust = 0, vjust = 1, size = 2.0, lineheight = 0.92, colour = "grey30",
           family = "Arial", fontface = "bold") +
  scale_x_continuous(limits = c(0, 2.05), expand = c(0, 0)) +
  scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(4, 2, 2, 0))

fig <- p_main + p_txt +
  plot_layout(widths = c(4.22, 1.22), guides = "collect") &
  theme(legend.position   = "bottom",
        legend.title      = element_blank(),
        legend.text       = element_text(size = 6.2, family = "Arial"),
        legend.key.size   = unit(3.2, "mm"),
        legend.key.height = unit(2.2, "mm"),
        legend.spacing.x  = unit(1.3, "mm"))

## ---- export ---------------------------------------------------------------
## final width 180 mm, final height 121 mm (double column)
cairo_pdf(file.path(outd, "Fig1_cascade.pdf"), width = 180/25.4, height = 121/25.4, bg = "white")
print(fig); dev.off()

svglite::svglite(file.path(outd, "Fig1_cascade.svg"), width = 180/25.4, height = 121/25.4, bg = "white")
print(fig); dev.off()

ragg::agg_tiff(file.path(outd, "Fig1_cascade.tiff"), width = 180, height = 121,
               units = "mm", res = 600, bg = "white", compression = "lzw")
print(fig); dev.off()

agg_png(file.path(outd, "Fig1_preview.png"), width = 180, height = 121,
        units = "mm", res = 400, bg = "white")
print(fig); dev.off()

## ---- source data ----------------------------------------------------------
rows |>
  transmute(panel, cohort, tier, row_lab,
            undiagnosed = undx,
            diagnosed_treatment_not_recorded = notmeas,
            diagnosed_untreated = dx_untx,
            reached_treatment = treated,
            diagnosed_pct = dx,
            reached_treatment_lo = lo_tx_per100_pain,
            reached_treatment_hi = hi_tx_per100_pain) |>
  write_csv(file.path(outd, "Fig1_source_data.csv"))

cat("EXPORTED\n")
