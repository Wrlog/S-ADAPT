# Simulated study for the antibody model: single 90-minute infusions of
# trastuzumab at 0.5, 2, 8 and 20 mg/kg, twelve 70-kg subjects each, serum
# sampled for six weeks. Population means are the trastuzumab estimates of
# Cao & Jusko (J Pharmacokinet Pharmacodyn 2014;41:375), Table 2. The dose
# range is wide on purpose: target-mediated disposition shows at the low
# doses, where the target is not saturated, and linear clearance at the
# high ones. Between-subject and residual variability are chosen for the
# exercise; the reflection coefficients are given almost none, and are
# estimated with the variance burn-then-shrink setting.

design <- list(
  title = "Antibody minimal PBPK with TMDD",
  outputs = c("Trastuzumab in serum (ug/mL)"),
  seed = 2014,
  positive = TRUE,
  truth = c(CLP = 0.00644, SIGMA1 = 0.95, SIGMA2 = 0.512, KSS = 0.99, KSYN = 0.376, SDIN = 0.05, SDSL = 0.12),
  # variance of the log; of the logit for the reflection coefficients
  bsv = c(CLP = 0.09, SIGMA1 = 0.01, SIGMA2 = 0.01, KSS = 0.09, KSYN = 0.09),
  subjects = function() {
    mgkg <- rep(c(0.5, 2, 8, 20), each = 12)
    lapply(mgkg, function(d) {
      nmol <- d * 70 * 1e6 / 148000
      list(group = d, rec = subject_rec(
        dose_rec(0, cmt = 1, amt = nmol, rate = nmol / 1.5),
        obs_rec(c(1.5, 4, 8, 24, 72, 168, 336, 504, 672, 840, 1008), cmt = 1)))
    })
  }
)
