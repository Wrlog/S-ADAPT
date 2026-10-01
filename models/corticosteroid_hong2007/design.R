# Simulated study for the corticosteroid model: 24 subjects given 30 mg of
# methylprednisolone as an intravenous bolus and 8 given placebo, sampled for
# 32 hours. Population means, between-subject variability and residual
# variability are the methylprednisolone estimates of Hong, Mager, Blum &
# Jusko (Pharm Res 2007;24:1088), Tables I-IV. The paper's inter-occasion
# variability is not simulated: each subject is studied once.

design <- list(
  title = "Corticosteroid PK/PD",
  outputs = c("Methylprednisolone (ng/mL)",
              "Cortisol (ng/mL)",
              "Lymphocytes (cells/uL)",
              "Proliferation (% of pre-dose)"),
  seed = 2007,
  positive = TRUE,
  truth = c(CLMPL = 22.8, VMPL = 78.4, KCOUT = 0.300, IC50MC = 1.59, KINL = 1125, KBE = 0.283, IC50CL = 91.1,
            IC50ML = 19.6, KT = 0.217, IC50MW = 0.984,
            SDPK = 0.0537, SDCORT = 0.344, SDLYM = 0.134, SDWBLP = 0.537),
  # the paper reports variability as %CV; the variance of the log is CV^2
  bsv = c(CLMPL = 0.109, VMPL = 0.135, KCOUT = 0.128, IC50MC = 0.773, KINL = 0.350, KBE = 0.282, IC50CL = 0.289,
          IC50ML = 0.775, KT = 0.793, IC50MW = 0.732)^2,
  subjects = function() {
    pd <- c(0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 28, 32)
    active <- list(group = 30, rec = subject_rec(
      dose_rec(0, cmt = 1, amt = 30000),                       # ug, so that concentrations are ng/mL
      obs_rec(c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12), cmt = 1),
      obs_rec(pd, cmt = 2), obs_rec(pd, cmt = 3),
      obs_rec(c(0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 24, 32), cmt = 4)))
    placebo <- list(group = 0, rec = subject_rec(obs_rec(pd, cmt = 2), obs_rec(pd, cmt = 3)))
    c(rep(list(active), 24), rep(list(placebo), 8))
  }
)
