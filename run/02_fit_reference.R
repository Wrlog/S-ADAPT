# Fit every simulated dataset with the reference MC-PEM (engine/mcpem.R),
# starting from the initial estimates in parameter_settings.csv, and compare
# what comes back with what the data were simulated from. Results go to
# results/reference/<name>/:
#   estimates.csv    population means, between-subject variances, residual
#                    parameters: simulated-from, sample, initial, estimate
#   iterations.csv   the history of the estimation, one row per iteration
#   posthoc.csv      each subject's conditional mean parameters
#   predictions.csv  observations with population and individual predictions
#
# The number of iterations is NPOPITER of the settings file, as in S-ADAPT.
# Run from the repository root: Rscript run/02_fit_reference.R [model ...]

for (f in c("ctl.R", "build.R", "simulate.R", "mcpem.R")) source(file.path("engine", f))

MODELS <- c("idr1_template", "corticosteroid_hong2007", "mab_mpbpk_tmdd", "adc_tdm1", "cart_stein2019", "protac_kcat")
NSAMP <- 300

fit_model <- function(name) {
  model <- load_model(name)
  subjects <- read_dataset(name)
  design <- read_design(name)
  ty <- par_types(model)
  par <- model$par[ty$p, ]
  t0 <- Sys.time()
  fit <- mcpem(model, subjects, niter = as.integer(model$settings$NPOPITER), nsamp = NSAMP, verbose = FALSE)
  mins <- as.numeric(Sys.time() - t0, units = "mins")

  # what the subjects in this dataset actually had: the centre and spread of
  # their true parameters on the transformed scale
  indiv <- as.matrix(utils::read.csv(data_paths(name)$truth)[, par$PNAME, drop = FALSE])
  phi <- vapply(seq_len(nrow(par)), function(j) to_phi(par[j, ], indiv[, j]), numeric(nrow(indiv)))
  sample_mean <- from_phi(par, colMeans(phi))[1, ]
  pct <- function(est, ref) round(100 * (est - ref) / ref, 1)
  vn <- model$par$PNAME[ty$v]
  est <- rbind(
    data.frame(PARAMETER = par$PNAME, KIND = "mean", TRANSFORM = par$PTRANSF, SIMULATED = design$truth[par$PNAME],
               SAMPLE = signif(sample_mean, 4), INITIAL = par$PMEAN, ESTIMATE = signif(fit$typical, 4),
               ERROR_PCT = pct(fit$typical, design$truth[par$PNAME])),
    data.frame(PARAMETER = par$PNAME, KIND = "variance", TRANSFORM = par$PTRANSF, SIMULATED = design$bsv[par$PNAME],
               SAMPLE = signif(apply(phi, 2, stats::var), 4), INITIAL = par$PCOV, ESTIMATE = signif(diag(fit$Omega), 4),
               ERROR_PCT = pct(diag(fit$Omega), design$bsv[par$PNAME])),
    data.frame(PARAMETER = vn, KIND = "residual", TRANSFORM = "", SIMULATED = design$truth[vn], SAMPLE = NA,
               INITIAL = model$par$PMEAN[ty$v], ESTIMATE = signif(fit$sigma, 4),
               ERROR_PCT = pct(fit$sigma, design$truth[vn])))
  dir <- file.path("results", "reference", name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  save <- function(x, file) utils::write.csv(x, file.path(dir, file), row.names = FALSE, quote = FALSE)
  save(est, "estimates.csv")
  save(signif(fit$history, 6), "iterations.csv")
  save(data.frame(ID = as.integer(names(subjects)), signif(fit$posthoc, 5)), "posthoc.csv")
  save(data.frame(lapply(mcpem_predictions(model, subjects, fit), function(z) if (is.double(z)) signif(z, 5) else z)),
       "predictions.csv")
  cat(sprintf("\n%s: %d iterations x %d draws, %.1f min, objective %.1f\n", name, fit$settings$niter, NSAMP, mins,
              fit$objective))
  print(est, row.names = FALSE)
  invisible(fit)
}

args <- commandArgs(TRUE)
for (name in if (length(args)) args else MODELS) fit_model(name)
