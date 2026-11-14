# =============================================================================
# Figure 3 — Functional limitation across cascade states
# Backend: R (ggplot2 + patchwork)  |  Data: outputs/31_outcome_v2.csv, 31_pooled_v2.csv
# -----------------------------------------------------------------------------
# FIGURE CONTRACT
#   Core conclusion : undiagnosed pain already carries 60-65% of the excess
#                     burden seen in diagnosed pain, and the treated subgroup is
#                     highest in all five cohorts - confounding by indication,
#                     not harm.
#   Evidence chain  : (a) three states vs neither pain nor diagnosis, 8 cohorts
#                         + 2 pooled diamonds; hero = pooled
#                     (b) treated subgroup separated, 5 cohorts with a tx item
#                     (c) diagnosis vs no diagnosis, before/after adjustment for
#                         pain severity, SAME cohorts on both sides
#   Archetype       : quantitative grid, log-scale forest
#   Export          : 180 x 132 mm; PDF/SVG vector + 600 dpi TIFF
# -----------------------------------------------------------------------------
# DELIBERATE DESIGN NOTE
#   Panel (c) recomputes the unadjusted pooled estimate on exactly the cohorts
#   that survive severity adjustment (panel B loses SHARE, which has no severity
#   item in waves 1/2/4). Comparing the published 4-cohort unadjusted estimate
#   with a 3-cohort adjusted one would confound adjustment with sample change.
#   All four state labels come from one shared scale so the legend appears once.
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
outd <- file.path(root, "manuscript", "figures")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

per <- read_csv(file.path(root, "outputs", "31_outcome_v2.csv"), show_col_types = FALSE)
pol <- read_csv(file.path(root, "outputs", "31_pooled_v2.csv"), show_col_types = FALSE)

## ---- DerSimonian-Laird, for the same-cohort recomputation in panel (c) -----
dl <- function(b, se) {
  w  <- 1/se^2; fe <- sum(w*b)/sum(w); k <- length(b)
  Q  <- sum(w*(b-fe)^2); C <- sum(w) - sum(w^2)/sum(w)
  tau2 <- max(0, (Q-(k-1))/C)
  wr <- 1/(se^2+tau2); re <- sum(wr*b)/sum(wr); sere <- sqrt(1/sum(wr))
  c(OR = exp(re), lo = exp(re-1.96*sere), hi = exp(re+1.96*sere),
    I2 = max(0, 100*(Q-(k-1))/Q), k = k)
}

b_adj <- pol |> filter(panel == "B", mode == "M2_adj_sev")
keep  <- strsplit(b_adj$cohorts, ",")[[1]]
b_same <- dl(per |> filter(panel == "B", mode == "M2_no_sev", cohort %in% keep) |> pull(b),
             per |> filter(panel == "B", mode == "M2_no_sev", cohort %in% keep) |> pull(se))
cat(sprintf("panel B, same %d cohorts: unadjusted OR %.2f (%.2f-%.2f), I2=%.1f\n",
            b_same[["k"]], b_same[["OR"]], b_same[["lo"]], b_same[["hi"]], b_same[["I2"]]))

## ---- one shared state vocabulary ------------------------------------------
LEV <- c("undx", "dx_any", "dx_untx", "tx", "legacy")
LAB <- c(undx          = "Pain, undiagnosed",
         dx_any        = "Pain, diagnosed, any treatment",
         dx_untx       = "Pain, diagnosed, untreated",
         tx            = "Pain, diagnosed, treated",
         legacy        = "No pain, previously diagnosed")
SHP <- c(undx = 16, dx_any = 15, dx_untx = 0, tx = 15, legacy = 17)
COL <- c(undx = "#9E9E9E", dx_any = "#6E9BC5", dx_untx = "#6E9BC5",
         tx = "#1F4E79", legacy = "#B08968")

one_scale <- function() {
  list(scale_colour_manual(values = COL, breaks = LEV, labels = LAB, drop = FALSE),
       scale_shape_manual(values = SHP, breaks = LEV, labels = LAB, drop = FALSE))
}

lab_cohort <- c(HRS = "HRS (USA)", SHARE = "SHARE (Europe)", ELSA = "ELSA (UK)",
                MHAS = "MHAS (Mexico)", CHARLS = "CHARLS (China)",
                KLoSA = "KLoSA (Korea)", LASI = "LASI (India)")

th <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text.y      = element_text(size = 6.3, colour = "black", hjust = 1, lineheight = 0.9),
        axis.text.x      = element_text(size = 6.4, colour = "black"),
        axis.title.x     = element_text(size = 6.8, colour = "black"),
        axis.line.x      = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks.x     = element_line(linewidth = 0.35, colour = "black"),
        panel.grid       = element_blank(),
        plot.margin      = margin(3, 4, 3, 6))

XSC <- scale_x_log10(limits = c(0.62, 12.2), breaks = c(1, 2, 4, 8),
                     labels = c("1", "2", "4", "8"), expand = c(0, 0))

## =========================== panel (a) =====================================
yA <- c(HRS = 9.9, SHARE_A = 9.0, ELSA = 8.1, MHAS = 7.2,
        CHARLS = 4.35, SHARE_B = 3.45, KLoSA = 2.55, LASI = 1.65)
yPOOL <- c(A = 6.20, B = 0.85)

a_rows <- per |>
  filter(mode == "three") |>
  mutate(key = case_when(cohort == "SHARE" & panel == "A" ~ "SHARE_A",
                         cohort == "SHARE" & panel == "B" ~ "SHARE_B",
                         TRUE ~ cohort),
         y = yA[key],
         state = factor(c(s1 = "undx", s2 = "dx_any", s3 = "legacy")[term], levels = LEV),
         lab = lab_cohort[cohort])

pool_lab <- pol |>
  filter(mode == "three") |>
  arrange(panel, OR) |>
  group_by(panel) |>
  summarise(txt = paste(sprintf("%.2f", OR), collapse = " · "), .groups = "drop") |>
  transmute(y = yPOOL[panel],
            lab = paste0("Pooled (4 cohorts)\n", txt))

a_brk <- bind_rows(a_rows |> distinct(y, lab), pool_lab) |> arrange(y)

a_pool <- pol |> filter(mode == "three") |>
  mutate(y = yPOOL[panel],
         state = factor(c(s1 = "undx", s2 = "dx_any", s3 = "legacy")[term], levels = LEV)) |>
  arrange(panel, OR) |>
  group_by(panel) |>
  # stagger labels when two pooled estimates sit closer than 12% on the log axis
  mutate(lab_dy = ifelse(row_number() > 1 & OR / lag(OR) < 1.12, -0.46, 0.40)) |>
  ungroup()

p_a <- ggplot() +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey55", linetype = "22") +
  geom_segment(data = a_rows, aes(x = lo, xend = hi, y = y, yend = y, colour = state),
               linewidth = 0.32) +
  geom_point(data = a_rows, aes(x = OR, y = y, colour = state, shape = state),
             size = 1.3, stroke = 0.32) +
  geom_rect(data = a_pool |> distinct(y), aes(xmin = 0.62, xmax = 12.2,
                                              ymin = y - 0.48, ymax = y + 0.48),
            inherit.aes = FALSE, fill = "grey94", colour = NA) +
  geom_segment(data = a_pool, aes(x = lo, xend = hi, y = y, yend = y, colour = state),
               linewidth = 0.45) +
  geom_point(data = a_pool, aes(x = OR, y = y, colour = state, shape = state),
             size = 1.9, stroke = 0.45, fill = "white") +
  annotate("text", x = 0.63, y = 10.85, label = "Frequent pain (panel A)",
           hjust = 0, size = 2.5, fontface = "bold", family = "Arial") +
  annotate("text", x = 0.63, y = 5.28, label = "Joint pain (panel B)",
           hjust = 0, size = 2.5, fontface = "bold", family = "Arial") +
  annotate("segment", x = 0.63, xend = 17, y = 5.62, yend = 5.62,
           linewidth = 0.2, colour = "grey80", linetype = "22") +
  scale_y_continuous(limits = c(0.02, 11.2), breaks = a_brk$y,
                     labels = a_brk$lab, expand = c(0, 0)) +
  XSC + one_scale() +
  labs(x = "Odds ratio for ADL limitation (log scale)", y = NULL) +
  coord_cartesian(clip = "off") + th +
  theme(axis.title.x = element_blank())

## =========================== panel (b) =====================================
yB <- c(HRS = 5.05, MHAS = 4.15, CHARLS = 2.95, KLoSA = 2.05, LASI = 1.15)
b_rows <- per |>
  filter(mode == "four") |>
  mutate(y = yB[cohort],
         state = factor(c(s1 = "undx", s2 = "dx_untx", s3 = "tx")[term], levels = LEV),
         lab = lab_cohort[cohort])
b_brk <- b_rows |> distinct(y, lab) |> arrange(desc(y))

p_b <- ggplot(b_rows) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey55", linetype = "22") +
  geom_segment(aes(x = lo, xend = hi, y = y, yend = y, colour = state),
               linewidth = 0.32, show.legend = FALSE) +
  geom_point(aes(x = OR, y = y, colour = state, shape = state),
             size = 1.3, stroke = 0.35, show.legend = FALSE) +
  annotate("segment", x = 0.63, xend = 17, y = 3.52, yend = 3.52,
           linewidth = 0.2, colour = "grey80", linetype = "22") +
  annotate("text", x = 0.63, y = 5.95,
           label = "(b) Treatment subgroup, five cohorts with a treatment item",
           hjust = 0, size = 2.6, fontface = "bold", family = "Arial") +
  scale_y_continuous(limits = c(0.72, 6.20), breaks = b_brk$y,
                     labels = b_brk$lab, expand = c(0, 0)) +
  XSC + one_scale() +
  labs(x = NULL, y = NULL) + coord_cartesian(clip = "off") + th +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

## =========================== panel (c) =====================================
c_dat <- tibble(
  panel = c("A", "A", "B", "B"),
  y     = c(1.00, 1.00, 0.00, 0.00),
  adj   = c("unadjusted", "adjusted", "unadjusted", "adjusted"),
  OR = c(pol |> filter(panel == "A", mode == "M2_no_sev") |> pull(OR),
         pol |> filter(panel == "A", mode == "M2_adj_sev") |> pull(OR),
         b_same[["OR"]],
         pol |> filter(panel == "B", mode == "M2_adj_sev") |> pull(OR)),
  lo = c(pol |> filter(panel == "A", mode == "M2_no_sev") |> pull(lo),
         pol |> filter(panel == "A", mode == "M2_adj_sev") |> pull(lo),
         b_same[["lo"]],
         pol |> filter(panel == "B", mode == "M2_adj_sev") |> pull(lo)),
  hi = c(pol |> filter(panel == "A", mode == "M2_no_sev") |> pull(hi),
         pol |> filter(panel == "A", mode == "M2_adj_sev") |> pull(hi),
         b_same[["hi"]],
         pol |> filter(panel == "B", mode == "M2_adj_sev") |> pull(hi))
)

c_txt <- c_dat |>
  group_by(panel, y) |>
  summarise(txt = sprintf("%.2f → %.2f (%.2f–%.2f)",
                          OR[adj == "unadjusted"], OR[adj == "adjusted"],
                          lo[adj == "adjusted"], hi[adj == "adjusted"]),
            .groups = "drop")

p_c <- ggplot(c_dat) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey55", linetype = "22") +
  geom_segment(aes(x = lo, xend = hi, y = y + 0.15, yend = y + 0.15),
               linewidth = 0.32, colour = "grey35", show.legend = FALSE) +
  geom_segment(aes(x = lo, xend = hi, y = y - 0.15, yend = y - 0.15),
               linewidth = 0.32, colour = "grey10") +
  geom_point(aes(x = OR, y = y + 0.15), shape = 1, size = 1.45, stroke = 0.38,
             colour = "grey35") +
  geom_point(aes(x = OR, y = y - 0.15), shape = 16, size = 1.35, colour = "grey10") +
  geom_text(data = c_txt, aes(x = 4.35, y = y, label = txt), hjust = 0,
            size = 2.05, colour = "grey5", family = "Arial") +
  annotate("text", x = 0.63, y = 1.92,
           label = "(c) Diagnosis vs no diagnosis, unadjusted → adjusted for pain severity",
           hjust = 0, size = 2.6, fontface = "bold", family = "Arial") +
  scale_y_continuous(limits = c(-0.55, 2.05), breaks = c(1.00, 0.00),
                     labels = c("Frequent pain (4 cohorts)",
                                "Joint pain (3 cohorts)"), expand = c(0, 0)) +
  XSC + labs(x = "Odds ratio for ADL limitation (log scale)", y = NULL) +
  coord_cartesian(clip = "off") + th

## =========================== manual legend =================================
leg <- tibble(
  state = factor(c("undx", "dx_any", "legacy", "dx_untx", "tx"), levels = LEV),
  x = c(0.025, 0.355, 0.685, 0.025, 0.355),
  y = c(1, 1, 1, 0, 0)
) |> mutate(lab = LAB[as.character(state)])

p_leg <- ggplot(leg) +
  geom_point(aes(x = x, y = y, colour = state, shape = state), size = 1.45,
             stroke = 0.4, fill = "white", show.legend = FALSE) +
  geom_text(aes(x = x + 0.017, y = y, label = lab), hjust = 0, size = 2.05,
            family = "Arial", colour = "black") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.55, 1.55), expand = c(0, 0)) +
  one_scale() +
  theme_void(base_family = "Arial") +
  theme(legend.position = "none", plot.margin = margin(2, 4, 2, 6))

## =========================== assemble ======================================
fig <- p_a / p_b / p_c / p_leg +
  plot_layout(heights = c(1.55, 0.82, 0.52, 0.15)) &
  theme(legend.position = "none")

## final width 180 mm, final height 134 mm (double column)
## log axis guard: every odds ratio is exp() of a finite log-OR and therefore
## strictly positive; no zero or negative value can reach scale_x_log10().
cairo_pdf(file.path(outd, "Fig3_outcomes.pdf"), width = 180/25.4, height = 134/25.4, bg = "white")
print(fig); dev.off()

svglite::svglite(file.path(outd, "Fig3_outcomes.svg"), width = 180/25.4, height = 134/25.4, bg = "white")
print(fig); dev.off()

ragg::agg_tiff(file.path(outd, "Fig3_outcomes.tiff"), width = 180, height = 134,
               units = "mm", res = 600, bg = "white", compression = "lzw")
print(fig); dev.off()

agg_png(file.path(outd, "Fig3_preview.png"), width = 180, height = 134,
        units = "mm", res = 400, bg = "white")
print(fig); dev.off()

write_csv(bind_rows(
  a_rows |> transmute(panel_label = paste0("a-", panel), cohort, state, OR, lo, hi, n),
  a_pool |> transmute(panel_label = paste0("a-", panel, "-pooled"), cohort = cohorts,
                      state = "pooled", OR, lo, hi, n = NA),
  b_rows |> transmute(panel_label = "b", cohort, state, OR, lo, hi, n),
  c_dat  |> transmute(panel_label = paste0("c-", panel), cohort = adj,
                      state = "dx_vs_undx", OR, lo, hi, n = NA)
), file.path(outd, "Fig3_source_data.csv"))

cat("EXPORTED\n")
