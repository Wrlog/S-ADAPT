# Simulated study for the template model: three dose groups of ten subjects,
# a 1-hour infusion, plasma concentrations (output 1) and the response
# (output 2). The population values are chosen for the exercise; they are
# not from a publication.

design <- list(
  title = "Template: PK with indirect response",
  outputs = c("Drug concentration (mg/L)",
              "Response"),
  seed = 2011,
  positive = TRUE,
  truth = c(CL = 3.5, V1 = 40, Tout12 = 6, Imax = 0.85, BASE = 20, IC50 = 4, HILL = 1.5,
            SDin = 0.05, SDsl = 0.10, PDin = 0.3, PDsl = 0.05),
  # between-subject variance on the transformed scale (log, or logit for Imax)
  bsv = c(CL = 0.09, V1 = 0.04, Tout12 = 0.06, Imax = 0.25, BASE = 0.04, IC50 = 0.16, HILL = 0.01),
  subjects = function() {
    doses <- rep(c(100, 300, 1000), each = 10)                 # mg
    lapply(doses, function(d) list(group = d, rec = subject_rec(
      dose_rec(0, cmt = 1, amt = d, rate = d),                 # over 1 h
      obs_rec(c(0.5, 1, 1.5, 2, 4, 6, 8, 12, 24), cmt = 1),
      obs_rec(c(0, 1, 2, 4, 6, 8, 12, 16, 24, 36, 48, 72), cmt = 2))))
  }
)
