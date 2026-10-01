# Diagnostic plots for the reference fits, from results/reference/<name>/:
#   fit.png          observations against individual predictions, per output
#   convergence.png  each population mean over the iterations
# and results/reference/recovery.png: every estimated population mean
# relative to the typical value in the simulated subjects, all models
# together.
#
# Run from the repository root: Rscript run/03_plots.R [model ...]

suppressMessages(library(ggplot2))
for (f in c("ctl.R", "build.R", "simulate.R")) source(file.path("engine", f))

MODELS <- c("idr1_template", "corticosteroid_hong2007", "mab_mpbpk_tmdd", "adc_tdm1", "cart_stein2019", "protac_kcat")
INK <- "#0b0b0b"
INK2 <- "#52514e"
SURFACE <- "#fcfcfb"
GRID <- "#e6e5e0"
BLUE <- "#2a78d6"

theme_fit <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = SURFACE, colour = NA), panel.background = element_rect(fill = SURFACE, colour = NA),
          panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          text = element_text(colour = INK2), axis.text = element_text(colour = INK2),
          plot.title = element_text(colour = INK, face = "bold", size = 12), plot.title.position = "plot",
          plot.subtitle = element_text(colour = INK2, size = 9.5), strip.text = element_text(colour = INK, hjust = 0, size = 9.5),
          panel.spacing = grid::unit(1.1, "lines"), plot.margin = margin(12, 16, 10, 12))
}
save_png <- function(plot, file, width, height) ggsave(file, plot, width = width, height = height, dpi = 150, bg = SURFACE)
tall <- function(n, ncol) 1.2 + 2.6 * ceiling(n / ncol)

plot_model <- function(name) {
  dir <- file.path("results", "reference", name)
  design <- read_design(name)
  pred <- utils::read.csv(file.path(dir, "predictions.csv"))
  pred$OUTPUT <- factor(design$outputs[pred$CMT], levels = design$outputs)
  # outputs that are positive and span more than two decades read better on log axes
  span <- tapply(pmin(pred$DV, pred$IPRED), pred$OUTPUT, function(z) all(z > 0))
  wide <- tapply(pred$DV, pred$OUTPUT, function(z) max(z) / max(min(z), 1e-300) > 100)
  logged <- names(which(span & wide))
  pred$LOG <- pred$OUTPUT %in% logged
  panels <- lapply(split(pred, pred$OUTPUT, drop = TRUE), function(d) {
    lim <- range(c(d$DV, d$IPRED))
    p <- ggplot(d, aes(IPRED, DV)) +
      geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.4) +
      geom_point(colour = BLUE, alpha = 0.55, size = 1.6, stroke = 0) +
      labs(title = as.character(d$OUTPUT[1]), x = "Individual prediction", y = "Observation") +
      theme_fit() + theme(plot.title = element_text(size = 9.5, face = "plain"))
    if (d$LOG[1]) p + scale_x_log10(limits = lim) + scale_y_log10(limits = lim) else p + coord_cartesian(xlim = lim, ylim = lim)
  })
  ncol <- min(3, length(panels))
  grob <- gridExtra::arrangeGrob(grobs = panels, ncol = ncol, top = grid::textGrob(
    sprintf("%s: observations against individual predictions", design$title), x = 0.012, hjust = 0,
    gp = grid::gpar(fontface = "bold", fontsize = 12, col = INK)))
  png(file.path(dir, "fit.png"), width = 3.3 * ncol + 0.4, height = tall(length(panels), ncol), units = "in", res = 150,
      bg = SURFACE)
  grid::grid.draw(grob)
  dev.off()

  hist <- utils::read.csv(file.path(dir, "iterations.csv"), check.names = FALSE)
  par <- load_model(name)$par
  means <- par$PNAME[par$PTYPE == "P"]
  long <- do.call(rbind, lapply(means, function(k) data.frame(ITER = hist$ITER, PARAMETER = k, VALUE = hist[[k]])))
  long$PARAMETER <- factor(long$PARAMETER, levels = means)
  truth <- data.frame(PARAMETER = factor(means, levels = means), VALUE = design$truth[means])
  burn <- max(par$VARBURN, na.rm = TRUE)
  p <- ggplot(long, aes(ITER, VALUE)) +
    geom_vline(xintercept = burn, colour = GRID, linewidth = 0.6) +
    geom_hline(data = truth, aes(yintercept = VALUE), colour = INK2, linewidth = 0.4, linetype = "22") +
    geom_line(colour = BLUE, linewidth = 0.6) +
    facet_wrap(~PARAMETER, scales = "free_y", ncol = 4) +
    labs(title = sprintf("%s: population means over the iterations", design$title),
         subtitle = sprintf("Dashed: the value the data were simulated from. The variances are held up until iteration %d.", burn),
         x = "Iteration", y = NULL) +
    theme_fit()
  save_png(p, file.path(dir, "convergence.png"), 9.6, 1.5 + 1.9 * ceiling(length(means) / 4))
}

plot_recovery <- function(models) {
  est <- do.call(rbind, lapply(models, function(name) {
    e <- utils::read.csv(file.path("results", "reference", name, "estimates.csv"))
    e <- e[e$KIND == "mean", ]
    data.frame(MODEL = read_design(name)$title, PARAMETER = e$PARAMETER, RATIO = e$ESTIMATE / e$SAMPLE)
  }))
  est$MODEL <- factor(est$MODEL, levels = unique(est$MODEL))
  est$PARAMETER <- factor(est$PARAMETER, levels = rev(unique(est$PARAMETER)))
  p <- ggplot(est, aes(RATIO, PARAMETER)) +
    annotate("rect", xmin = 0.8, xmax = 1.25, ymin = -Inf, ymax = Inf, fill = GRID, alpha = 0.6) +
    geom_vline(xintercept = 1, colour = INK2, linewidth = 0.4) +
    geom_point(colour = BLUE, size = 2.4) +
    facet_wrap(~MODEL, scales = "free_y", ncol = 3) +
    scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 2), limits = c(0.45, 2.2)) +
    labs(title = "Estimated population means relative to the subjects the data were simulated from",
         subtitle = paste("Reference MC-PEM fits. The reference value is the typical (geometric mean) parameter of the",
                          "simulated subjects. 1 is exact recovery; the band is 0.8 to 1.25."),
         x = "Estimate / value in the simulated subjects", y = NULL) +
    theme_fit() + theme(panel.grid.major.y = element_blank())
  save_png(p, file.path("results", "reference", "recovery.png"), 10, 6.4)
}

args <- commandArgs(TRUE)
for (name in if (length(args)) args else MODELS) plot_model(name)
done <- MODELS[file.exists(file.path("results", "reference", MODELS, "estimates.csv"))]
plot_recovery(done)
cat("plots written for", length(if (length(args)) args else MODELS), "model(s); recovery.png covers", length(done), "\n")
