# The reference MC-PEM, run briefly on the template dataset from the
# initial estimates of its settings file, must bring the well-determined
# parameters back to what the data were simulated from. A short run is
# enough for these; the full fits are in results/reference/.
#
# Run from the repository root: Rscript tests/test_reference_fit.R

for (f in c("ctl.R", "build.R", "simulate.R", "mcpem.R")) source(file.path("engine", f))

name <- "idr1_template"
model <- load_model(name)
design <- read_design(name)
fit <- mcpem(model, read_dataset(name), niter = 45, nsamp = 200, navg = 5, verbose = FALSE)

fails <- 0
within <- function(what, est, ref, tol) {
  err <- est / ref - 1
  cat(sprintf("  %-28s estimate %8.4g  simulated %8.4g  %+6.1f%%\n", what, est, ref, 100 * err))
  if (!is.finite(err) || abs(err) > tol) fails <<- fails + 1
}
for (k in c("CL", "V1", "BASE")) within(k, fit$typical[[k]], design$truth[[k]], 0.10)
within("Tout12", fit$typical[["Tout12"]], design$truth[["Tout12"]], 0.25)
within("SDsl (residual)", fit$sigma[["SDsl"]], design$truth[["SDsl"]], 0.30)
within("PDsl (residual)", fit$sigma[["PDsl"]], design$truth[["PDsl"]], 0.30)
cat(sprintf("  objective function fell from %.0f to %.0f\n", fit$history$OBJ[1], fit$history$OBJ[45]))
if (fit$history$OBJ[45] >= fit$history$OBJ[1]) fails <- fails + 1

if (fails) stop(fails, " check(s) failed")
cat("the reference MC-PEM recovers the template's parameters\n")
