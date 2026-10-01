# Simulated experiment for the PROTAC model: twelve independent experiments,
# each exposing cells to six constant degrader concentrations from 1 nM to
# 100 uM and measuring the target protein over 24 hours as percent of
# control. Each experiment is a "subject" with six arms: concentration k is
# given as a bolus into compartment k, which has no elimination, and the
# target in that arm is output k. Population means: kcat 4.6 1/h and a
# 16-hour target half-life, the BTK degrader of Pharmaceutics 2023;15:195.
# Variability between experiments and the assay error are chosen for the
# exercise.

design <- list(
  title = "PROTAC degradation",
  outputs = c("Target at 1 nM (% of control)",
              "Target at 10 nM (% of control)",
              "Target at 100 nM (% of control)",
              "Target at 1 uM (% of control)",
              "Target at 10 uM (% of control)",
              "Target at 100 uM (% of control)"),
  seed = 2023,
  positive = TRUE,
  truth = c(KCAT = 4.6, THALFP = 16, SDIN = 2, SDSL = 0.08),
  bsv = c(KCAT = 0.20, THALFP = 0.15)^2,
  subjects = function() {
    conc <- c(1, 10, 100, 1000, 10000, 100000)                  # nM
    arms <- lapply(seq_along(conc), function(k) rbind(
      dose_rec(0, cmt = k, amt = conc[k]),
      obs_rec(c(0.5, 1, 2, 4, 8, 16, 24), cmt = k)))
    rep(list(list(group = 1, rec = do.call(subject_rec, arms))), 12)
  }
)
