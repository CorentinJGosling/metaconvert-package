## Renders the figures used by the redesign mockups.
##
## The view is the one the app opens on for study 01a: confidence interval
## coverage against the biserial population parameter, rho on the horizontal
## axis, split by n, p_exp = 0.3, 1,000 replications per cell. Every number in
## the mockups can therefore be checked against
## data/aggregated/01a_smd_to_cor_r_nrep1000.csv rather than taken on trust.
##
## Run from simulations/:  Rscript app/redesign-2026-08/make-figures.R

library(ggplot2)

DIR <- "app/redesign-2026-08"
MC  <- list(ink = "#1c1a2e", ink2 = "#2c2a47", soft = "#5d5c6e",
            rule = "#e2e2dc", faint = "#83828f")
PAL <- c(lipsey_cooper = "#c0474e", viechtbauer = "#2c2a47")

d <- utils::read.csv("data/aggregated/01a_smd_to_cor_r_nrep1000.csv",
                     stringsAsFactors = FALSE)
d <- d[d$target == "biserial_population" & d$p_exp == 0.3, ]
d$n_lab <- factor(paste0("n = ", d$n), levels = paste0("n = ", sort(unique(d$n))))

base_theme <- theme_minimal(base_size = 12, base_family = "Work Sans") +
  theme(
    text               = element_text(colour = MC$ink),
    axis.title.x       = element_text(size = 10, colour = MC$soft,
                                      margin = margin(t = 10), hjust = 0),
    axis.text          = element_text(size = 9, colour = MC$soft,
                                      family = "Cascadia Code"),
    axis.line.x        = element_line(colour = MC$ink, linewidth = 0.4),
    axis.ticks.x       = element_line(colour = MC$ink, linewidth = 0.4),
    axis.ticks.length  = unit(3, "pt"),
    ## horizontal rules only: rho is a short ordered set of designed values, so
    ## a vertical grid would rule between points that are already marked
    panel.grid.major.y = element_line(colour = MC$rule, linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(1.7, "lines"),
    ## the panel label is a caption, not a button: no boxed strip
    strip.background   = element_blank(),
    strip.text         = element_text(size = 9.5, colour = MC$ink, hjust = 0,
                                      family = "Cascadia Code",
                                      margin = margin(0, 0, 8, 0)),
    legend.position    = "none",
    plot.margin        = margin(2, 6, 2, 2)
  )

## The nominal rate is the thing every line is read against, so it is drawn
## first and underneath, as a reference rule rather than as another series.
p <- ggplot(d, aes(rho, coverage, colour = method, group = method)) +
  geom_hline(yintercept = 0.95, linetype = "22", linewidth = 0.4,
             colour = MC$faint) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.9, stroke = 0) +
  facet_wrap(~ n_lab, nrow = 1) +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75),
                     labels = c("0", ".25", ".50", ".75"),
                     expand = expansion(mult = 0.07)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 0.95),
                     labels = c("0", ".25", ".50", ".75", ".95"),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "true correlation  rho", y = NULL) +
  base_theme

grDevices::svg(file.path(DIR, "fig-coverage.svg"), width = 9.6, height = 2.6,
               bg = "transparent")
print(p)
invisible(grDevices::dev.off())

## Compact single-panel version for the listing layout: n on the horizontal
## axis, one line per rho, so the "gets worse as n grows" direction is the
## axis rather than something the reader has to read across four panels.
d2 <- d[d$method == "lipsey_cooper", ]
d2$rho_lab <- factor(d2$rho)
p2 <- ggplot(d2, aes(n, coverage, group = rho_lab)) +
  geom_hline(yintercept = 0.95, linetype = "22", linewidth = 0.4,
             colour = MC$faint) +
  geom_line(linewidth = 0.6, colour = PAL[["lipsey_cooper"]]) +
  geom_point(size = 1.8, stroke = 0, colour = PAL[["lipsey_cooper"]]) +
  scale_x_continuous(breaks = c(25, 50, 100, 300), trans = "log10") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95),
                     labels = c("0", ".50", ".95")) +
  labs(x = "n (log spacing)", y = NULL) +
  base_theme

grDevices::svg(file.path(DIR, "fig-lipsey-vs-n.svg"), width = 3.2, height = 2.2,
               bg = "transparent")
print(p2)
invisible(grDevices::dev.off())

cat("wrote figures to", DIR, "\n")
