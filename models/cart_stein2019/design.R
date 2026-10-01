# Simulated study for the CAR-T model: 60 patients, transgene copies in blood
# over one year, modelled as their natural log. Population means and
# between-patient standard deviations are those of Stein et al. (CPT
# Pharmacometrics Syst Pharmacol 2019;8:285), Table 1. The residual standard
# deviation is chosen for the exercise. There are no dose records: in this
# model the cell level starts from its initial condition.

design <- list(
  title = "CAR-T cellular kinetics",
  outputs = c("Transgene, ln(copies/ug DNA)"),
  seed = 2019,
  truth = c(FOLDX = 3900, TPEAK = 9.3, CPEAK = 24000, KALPHA = 0.16, KBETA = 0.0032, FPERS = 0.0079, SDLN = 0.5),
  bsv = c(FOLDX = 2.4, TPEAK = 0.38, CPEAK = 0.65, KALPHA = 0.91, KBETA = 0.86, FPERS = 0.8)^2,
  subjects = function() {
    days <- c(1, 3, 5, 7, 9, 11, 14, 17, 21, 28, 42, 60, 90, 180, 270, 365)
    rep(list(list(group = 1, rec = subject_rec(obs_rec(days, cmt = 1)))), 60)
  }
)
