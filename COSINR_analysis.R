# COSINR AREG analysis
# Run from the project root directory

library(data.table)
library(survival)
library(ggplot2)
library(ggsurvfit)
library(ggpubr)
library(ggsci)
library(forcats)
library(patchwork)

##### constants #####
figdata <- "output/fig_data/"
figpath <- "output/figures/"
dir.create(figpath, showWarnings = FALSE, recursive = TRUE)

tp.levels <- c("Pre-treatment", "Early treatment", "Late treatment",
               "Pre- to early-treatment", "Early- to late-treatment",
               "Pre- to late-treatment change")

tp.colors <- c("Pre-treatment" = "#0073C2FF", "Early treatment" = "#3B3B3BFF", "Late treatment" = "#CD534CFF",
               "Pre- to late-treatment change" = "#00A087B2")

##### functions #####
read_fig <- function(f) {
  dt <- fread(paste0(figdata, f), encoding = "UTF-8")
  if ("id" %in% names(dt)) dt[, id := sprintf("%03d", as.integer(id))]
  chr <- names(dt)[sapply(dt, is.character)]
  dt[, (chr) := lapply(.SD, function(x) ifelse(validUTF8(x), x, iconv(x, "macintosh", "UTF-8"))), .SDcols = chr]
  if ("time.point" %in% names(dt)) dt[, time.point := factor(time.point, levels = tp.levels)]
  dt[]
}

fmt_p <- function(p) {
  fcase(p < .001, "<0.001",
        round(p, 2) == .05, sprintf("%.3f", p),
        p < .01, sprintf("%.3f", p),
        default = sprintf("%.2f", p))
}

km_plot <- function(dt, grp, pval.x = 10, leg.coords = c(0.2, 0.2)) {
  survfit2(as.formula(paste0("Surv(os_time_mo, event_death) ~ ", grp)), data = dt) |>
    ggsurvfit(theme = theme_classic(), linetype_aes = FALSE, linewidth = 0.9) +
    add_censor_mark() +
    add_risktable(risktable_stats = "n.risk") +
    ylim(0, 1) +
    theme(legend.position = "inside",
          legend.position.inside = leg.coords) +
    ylab("Estimated overall survival") +
    scale_color_npg() +
    add_pvalue(location = "annotation", x = pval.x, y = 0.05, size = 3.5)
}

corline <- function(dt, x, y, xttl, yttl) {
  ggplot(dt, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(colour = "black", shape = 16, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey50", alpha = 0.25) +
    stat_cor(method = "spearman", cor.coef.name = "rho") +
    theme_classic() +
    xlab(xttl) +
    ylab(yttl)
}

corline2 <- function(dt, facet.var, y.label, wrap = 20) {
  dt <- copy(dt)[, significant.p := fifelse(cor.test(AREG_NPX, value, method = "spearman", exact = FALSE)$p.value < 0.05,
                                            "p < 0.05", "p >= 0.05"),
                 by = c(facet.var, "time.point")]
  ggplot(dt, aes(x = AREG_NPX, y = value)) +
    geom_point(fill = "grey15", colour = "black", shape = 21, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey50", alpha = 0.25) +
    stat_cor(method = "spearman", cor.coef.name = "rho", size = 2.5, label.x.npc = "right", hjust = 1) +
    geom_rect(aes(color = significant.p), fill = NA, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    scale_color_manual(values = c("p < 0.05" = "red3", "p >= 0.05" = NA), na.value = NA) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.y = element_text(size = 7),
          axis.text.y = element_text(size = 7)) +
    ylab(y.label) +
    xlab("log2 normalized late-treatment AREG expression") +
    facet_grid(as.formula(paste(facet.var, "~ time.point")), scales = "free_y",
               labeller = label_wrap_gen(width = wrap))
}

corbar <- function(dt, lab, xlb) {
  ggplot(copy(dt)[, p.value := fifelse(p.value < 0.05, p.value, NA_real_)],
         aes(x = fct_reorder(.data[[lab]], spearman_rho), y = spearman_rho, fill = p.value)) +
    geom_bar(stat = "identity") +
    xlab(xlb) +
    ylab("Spearman's rho") +
    theme_classic() +
    scale_fill_gradient(limits = c(0.001, 0.05), breaks = c(0.001, 0.05),
                        low = "#B2182B", high = "#FDDBC7", na.value = "grey70",
                        oob = scales::squish) +
    theme(axis.text.x = element_text(size = 5),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          legend.title = element_text(size = 8, vjust = 0.75),
          strip.background = element_blank(),
          legend.position = "bottom") +
    labs(fill = "p-value") +
    coord_flip() +
    facet_grid(. ~ time.point)
}

forestfunc <- function(dt, xlimmin = 0, xlimmax = 3) {
  dt <- copy(dt)[, let(HR.txt  = sprintf("%.2f (%.2f to %.2f)", HR, conf.low, conf.high),
                       p.txt   = fmt_p(p.value),
                       conf.high = pmin(conf.high, xlimmax),
                       metric  = "AREG")]
  hl <- geom_hline(yintercept = c(0.5, 1.5), color = "grey60", linewidth = 0.2)
  dodge <- position_dodge2(width = 1, reverse = TRUE)

  tl <- ggplot(dt, aes(y = metric, group = time.point)) +
    geom_text(aes(x = 0, label = metric), hjust = 0, size = 3) +
    geom_text(aes(x = 1, label = HR.txt), hjust = 0, size = 3, position = dodge) +
    theme_void() +
    coord_cartesian(xlim = c(0, 2)) +
    ggtitle("Assay                 HR") +
    hl

  tf <- ggplot(dt, aes(y = metric, group = time.point)) +
    geom_point(aes(x = HR, color = time.point), shape = 15, size = 1, position = dodge) +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high, color = time.point), position = dodge) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    xlim(c(xlimmin, xlimmax)) +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    hl +
    labs(color = "Time point") +
    xlab("Hazard Ratio") +
    ggtitle("    ") +
    scale_color_manual(values = tp.colors)

  tp <- ggplot(dt, aes(y = metric, group = time.point)) +
    geom_text(aes(x = 0, label = p.txt), hjust = 0, size = 3, position = dodge) +
    theme_void() +
    hl +
    ggtitle("           p-value")

  ly <- c(area(t = 0, l = 0, b = 30, r = 3),
          area(t = 0, l = 4, b = 30, r = 11),
          area(t = 0, l = 11, b = 30, r = 12))

  tl + tf + tp + plot_layout(design = ly)
}

##### Figure 1 #####
fig1a <- forestfunc(read_fig("Fig1a.csv"), 0, 3)
ggsave(paste0(figpath, "Fig1a.png"), fig1a, width = 9, height = 3, dpi = 300)

fig1b <- km_plot(read_fig("Fig1b.csv"), "AREG_high_low")
ggsave(paste0(figpath, "Fig1b.png"), fig1b, width = 7.5, height = 5, dpi = 300)

fig1c <- km_plot(read_fig("Fig1c.csv"), "AREG_change_high_low")
ggsave(paste0(figpath, "Fig1c.png"), fig1c, width = 7.5, height = 5, dpi = 300)

fig1d <- corline(read_fig("Fig1d.csv"), "AREG_NPX", "recist_mm",
                 "Baseline log2 normalized AREG expression", "Baseline tumor burden by RECIST (mm)")
ggsave(paste0(figpath, "Fig1d.png"), fig1d, width = 7.5, height = 6, dpi = 300)

fig1e <- corline(read_fig("Fig1e.csv"), "AREG_NPX_change", "change_recist_mm",
                 "Change in log2 normalized AREG expression", "Change in RECIST measurements (mm)")
ggsave(paste0(figpath, "Fig1e.png"), fig1e, width = 7.5, height = 6, dpi = 300)

##### Figure 2 #####
fig2a <- corbar(read_fig("Fig2a.csv"), "hallmark_label", "ssGSEA Pathway")
ggsave(paste0(figpath, "Fig2a.png"), fig2a, width = 8, height = 5, dpi = 300)

fig2b <- corbar(read_fig("Fig2b.csv"), "metric", "Cell Type by Flow Cytometry")
ggsave(paste0(figpath, "Fig2b.png"), fig2b, width = 8, height = 5, dpi = 300)

fig2c.dat <- read_fig("Fig2c.csv")
fig2c <- corline2(fig2c.dat[full_text == "Total T-cells"], "full_text", "Count") /
  (corline2(fig2c.dat[full_text == "Contracted TCR clones"], "full_text", "Count of TCR clones") + ylim(c(0, 400))) /
  corline2(fig2c.dat[full_text == "Observed Richness"], "full_text", "Count")
ggsave(paste0(figpath, "Fig2c.png"), fig2c, width = 5, height = 7.5, dpi = 300)

##### Supplementary #####
figS1 <- km_plot(read_fig("FigS1.csv"), "lesions_treated", pval.x = 17, leg.coords = c(0.3, 0.2))
ggsave(paste0(figpath, "FigS1.png"), figS1, width = 7.5, height = 5, dpi = 300)

figS2a <- corline2(read_fig("FigS2a.csv"), "hallmark_label", "ssGSEA score")
ggsave(paste0(figpath, "FigS2a.png"), figS2a, width = 5, height = 16, dpi = 300)

figS2b.dat <- read_fig("FigS2b.csv")[, block := cumsum(c(TRUE, diff(as.integer(id)) <= 0)), by = label][
  , time.point := factor(tp.levels[block], levels = tp.levels)]
figS2b <- corline2(figS2b.dat, "label", "Percent of cells*", wrap = 24)
ggsave(paste0(figpath, "FigS2b.png"), figS2b, width = 5, height = 18, dpi = 300)
