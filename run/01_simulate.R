# Simulate the dataset of every model from its design (models/<name>/design.R)
# and write it to data/ in the layout S-ADAPT reads:
#   data/<name>.csv        ID, TIME, DV, AMT, RATE, EVID, MDV, CMT
#   data/<name>_cov.csv    ID and GROUP (dose or concentration)
#   data/<name>_truth.csv  the parameters each subject was simulated with
#
# Run from the repository root: Rscript run/01_simulate.R [model ...]

for (f in c("ctl.R", "build.R", "simulate.R")) source(file.path("engine", f))

MODELS <- c("idr1_template", "corticosteroid_hong2007", "mab_mpbpk_tmdd", "adc_tdm1", "cart_stein2019", "protac_kcat")
args <- commandArgs(TRUE)
for (name in if (length(args)) args else MODELS) {
  model <- load_model(name)
  sim <- simulate_dataset(model, read_design(name))
  write_dataset(name, sim)
  obs <- sim$data[sim$data$EVID == 0, ]
  cat(sprintf("%-26s %3d subjects  %5d observations  %s\n", name, nrow(sim$cov), nrow(obs),
              paste(sprintf("Y%d: %.3g to %.3g", sort(unique(obs$CMT)),
                            tapply(obs$DV, obs$CMT, min), tapply(obs$DV, obs$CMT, max)), collapse = "; ")))
}
