# One dataset can land above or below the truth by chance. This simulates
# the same study many times, fits each with the reference MC-PEM, and
# reports the bias and spread of the estimates: the simulate-and-fit check
# S-ADAPT offers with its popsets_simulate and popsets_fit commands.
# Results go to results/reference/<name>/replicates.csv.
#
# Run from the repository root:
#   Rscript run/04_replicates.R <model> [replicates] [iterations]

for (f in c("ctl.R", "build.R", "simulate.R", "mcpem.R")) source(file.path("engine", f))

args <- commandArgs(TRUE)
name <- if (length(args) >= 1) args[1] else "protac_kcat"
nrep <- if (length(args) >= 2) as.integer(args[2]) else 16
niter <- if (length(args) >= 3) as.integer(args[3]) else 70

model <- load_model(name)
design <- read_design(name)
ty <- par_types(model)
par <- model$par[ty$p, ]
rows <- lapply(seq_len(nrep), function(r) {
  design$seed <- 1000 + r
  sim <- simulate_dataset(model, design)
  fit <- mcpem(model, split(sim$data, sim$data$ID), niter = niter, nsamp = 300, seed = r, verbose = FALSE)
  indiv <- as.matrix(sim$truth[, par$PNAME, drop = FALSE])
  phi <- vapply(seq_len(nrow(par)), function(j) to_phi(par[j, ], indiv[, j]), numeric(nrow(indiv)))
  cat(sprintf("  replicate %2d  %s\n", r, paste(sprintf("%s=%.4g", par$PNAME, fit$typical), collapse = " ")))
  data.frame(REPLICATE = r, PARAMETER = c(par$PNAME, par$PNAME, names(fit$sigma)),
             KIND = rep(c("mean", "variance", "residual"), c(nrow(par), nrow(par), length(fit$sigma))),
             SIMULATED = c(design$truth[par$PNAME], design$bsv[par$PNAME], design$truth[names(fit$sigma)]),
             SAMPLE = c(from_phi(par, colMeans(phi))[1, ], apply(phi, 2, stats::var), rep(NA, length(fit$sigma))),
             ESTIMATE = c(fit$typical, diag(fit$Omega), fit$sigma))
})
out <- do.call(rbind, rows)
out[c("SAMPLE", "ESTIMATE")] <- lapply(out[c("SAMPLE", "ESTIMATE")], signif, 5)
dir <- file.path("results", "reference", name)
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, file.path(dir, "replicates.csv"), row.names = FALSE, quote = FALSE)

cat(sprintf("\n%s: %d replicates, %d iterations each\n", name, nrep, niter))
cat("mean of log(estimate / reference), with its standard error, in percent\n")
m <- out[out$KIND == "mean", ]
for (k in par$PNAME) {
  z <- m[m$PARAMETER == k, ]
  a <- 100 * log(z$ESTIMATE / z$SIMULATED)
  b <- 100 * log(z$ESTIMATE / z$SAMPLE)
  cat(sprintf("  %-8s against the population value %+5.1f (%.1f)   against the simulated subjects %+5.1f (%.1f)\n", k,
              mean(a), stats::sd(a) / sqrt(nrep), mean(b), stats::sd(b) / sqrt(nrep)))
}
other <- out[out$KIND != "mean", ]
avg <- stats::aggregate(ESTIMATE ~ PARAMETER + KIND + SIMULATED, other, mean)
cat("mean estimate of the variances and residual parameters\n")
print(data.frame(avg[c("PARAMETER", "KIND", "SIMULATED")], MEAN_ESTIMATE = signif(avg$ESTIMATE, 3)), row.names = FALSE)
