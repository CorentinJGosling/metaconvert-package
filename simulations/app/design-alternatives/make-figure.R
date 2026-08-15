## Renders the figure used by the static design mockups.
##
## It is the app's DEFAULT view of study 05 -- RMSE against rr, split by n,
## target = population, br = 0.01, br_guess = 0.01, p_exp = 0.5 -- so the numbers
## in the mockups can be checked against data/aggregated/05_rr_to_or_nrep1000.csv
## rather than taken on trust. Nothing here is invented for the sake of a layout.
##
## Run from simulations/:  Rscript app/design-alternatives/make-figure.R

library(ggplot2)

OUT <- "app/design-alternatives/figure-05-rmse.svg"

MC <- list(ink = "#211f38", ink2 = "#2c2a47", soft = "#6b6a7d",
           rule = "#e6e6e0", paper = "#f7f7f3")
PAL <- c(dipietrantonj          = "#d66268",
         grant                  = "#2c2a47",
         metaumbrella           = "#8a6712",
         transpose              = "#6b4370",
         vanderweele_rr_squared = "#4a6fa5")

d <- utils::read.csv("data/aggregated/05_rr_to_or_nrep1000.csv",
                     stringsAsFactors = FALSE)
d <- d[d$target == "population" & d$br == 0.01 &
       d$br_guess == 0.01 & d$p_exp == 0.5, ]
d$n_lab <- factor(paste0("n = ", d$n), levels = paste0("n = ", sort(unique(d$n))))

## The risk ratios are 0.25, 0.5, 0.75, 1 and 2. On a linear axis four of the
## five crowd into the left third and their labels collide. A ratio belongs on a
## log axis anyway -- 0.25, 0.5, 1, 2 are equally spaced there -- so the position
## is log2(rr) and the labels are the ratios themselves. The transform is applied
## to the data rather than through scale_x_continuous(trans=/transform=), whose
## argument name changed between ggplot2 versions.
RR <- c(0.25, 0.5, 0.75, 1, 2)
d$x <- log2(d$rr)

p <- ggplot(d, aes(x, rmse, colour = method, group = method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.9, stroke = 0) +
  facet_wrap(~ n_lab, nrow = 1) +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_x_continuous(breaks = log2(RR),
                     labels = c("0.25", "0.5", "0.75", "1", "2"),
                     expand = expansion(mult = 0.06)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
  labs(x = "risk ratio (log spacing)", y = NULL) +
  theme_minimal(base_size = 12, base_family = "Work Sans") +
  theme(
    text              = element_text(colour = MC$ink),
    axis.title.x      = element_text(size = 10.5, colour = MC$soft,
                                     margin = margin(t = 9), hjust = 0),
    axis.text         = element_text(size = 9.5, colour = MC$soft,
                                     family = "Cascadia Code"),
    axis.line.x       = element_line(colour = MC$ink, linewidth = 0.4),
    axis.ticks.x      = element_line(colour = MC$ink, linewidth = 0.4),
    axis.ticks.length = unit(3, "pt"),
    ## horizontal rules only: vertical grid adds nothing when x is a short
    ## ordered set of designed values
    panel.grid.major.y = element_line(colour = MC$rule, linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(1.6, "lines"),
    ## no boxed strip: the panel label is a caption, not a button
    strip.background   = element_blank(),
    strip.text         = element_text(size = 10, colour = MC$ink, hjust = 0,
                                      family = "Cascadia Code",
                                      margin = margin(0, 0, 7, 0)),
    legend.position    = "none",
    plot.margin        = margin(2, 4, 2, 2)
  )

grDevices::svg(OUT, width = 9.4, height = 2.9, bg = "transparent")
print(p)
invisible(grDevices::dev.off())
cat("wrote", OUT, "\n")
