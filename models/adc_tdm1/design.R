# Simulated study for the T-DM1 model: KPL-4 xenograft mice (25 g), a single
# intravenous dose of vehicle or 0.3, 1, 3 or 15 mg/kg, eight mice a group.
# Total trastuzumab (output 1) and T-DM1 (output 2) in plasma over two
# weeks, and tumour volume (output 3) twice a week for four weeks.
# Population means are from Singh & Shah (AAPS J 2017;19:1054), Table I:
# mouse PK, and the KPL-4 doubling time, killing constant and its
# between-animal variability (21%). The other variabilities are chosen for
# the exercise.

design <- list(
  title = "T-DM1 antibody-drug conjugate",
  outputs = c("Total trastuzumab (ug/mL)",
              "T-DM1 conjugate (ug/mL)",
              "Tumour volume (mm3)"),
  seed = 2017,
  positive = TRUE,
  truth = c(CLADC = 0.0934, V1ADC = 0.043, KDEC = 0.241, DTEXP = 8.21, KKILL = 1.96e-3, TV0 = 150,
            SDTT = 0.15, SDADC = 0.15, SDTV = 0.15),
  bsv = c(CLADC = 0.20, V1ADC = 0.15, KDEC = 0.15, DTEXP = 0.15, KKILL = 0.21, TV0 = 0.20)^2,
  subjects = function() {
    mgkg <- rep(c(0, 0.3, 1, 3, 15), each = 8)
    tumour <- obs_rec(c(0, 3, 7, 10, 14, 17, 21, 24, 28), cmt = 3)
    lapply(mgkg, function(d) {
      if (d == 0) return(list(group = d, rec = subject_rec(tumour)))
      nmol <- d * 0.025 * 1e6 / 148500
      pk <- c(0.02, 1, 3, 7, 14)
      # the dose enters total antibody (compartment 1) and conjugate (3)
      list(group = d, rec = subject_rec(dose_rec(0, cmt = 1, amt = nmol), dose_rec(0, cmt = 3, amt = nmol),
                                        obs_rec(pk, cmt = 1), obs_rec(pk, cmt = 2), tumour))
    })
  }
)
