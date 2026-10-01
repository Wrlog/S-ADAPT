# Every control stream, translated to Fortran and solved by the reference
# engine, must agree with an mrgsolve model written independently of it
# (reference/). Four of those are the model files of the PBPK and QSP
# repositories; the other two were written for this check. If they agree,
# the code S-ADAPT will compile computes the published model.
#
# Run from the repository root: Rscript tests/test_models_vs_mrgsolve.R

for (f in c("ctl.R", "build.R", "simulate.R")) source(file.path("engine", f))
if (!requireNamespace("mrgsolve", quietly = TRUE)) stop("mrgsolve is not installed")
suppressMessages(library(mrgsolve))

worst <- c()
report <- function(name, a, b) {
  d <- max(abs(a - b) / pmax(abs(b), 1e-12))
  cat(sprintf("  %-44s %.2e\n", name, d))
  worst[name] <<- d
}
mrg <- function(file, param, ev, times) {
  mod <- suppressMessages(mread(file, project = "reference", quiet = TRUE))
  # mrgsolve starts the clock at its first record, so time zero is always one
  out <- mrgsim_df(param(mod, param), events = if (is.null(ev)) NULL else as.ev(ev), start = 0, end = -1,
                   add = union(0, times), obsonly = TRUE, atol = 1e-12, rtol = 1e-10, maxsteps = 1e5)
  out <- out[!duplicated(out$time, fromLast = TRUE), ]
  out[match(times, out$time), ]
}
setup <- function(name) {
  m <- load_model(name)
  truth <- read_design(name)$truth
  list(m = m, p = truth[m$par$PNAME], truth = as.list(truth))
}
no_dose <- dose_rec(0, 1, 0)[0, ]

# Template: 1000 mg over 1 h
x <- setup("idr1_template")
tt <- c(0.5, 1, 2, 4, 8, 12, 24, 48, 72)
a <- sa_profile(x$m, x$p, dose_rec(0, 1, 1000, 1000), tt)
b <- mrg("idr1_template", with(x$truth, list(CL = CL, V1 = V1, TOUT12 = Tout12, IMAX = Imax, BASE = BASE, IC50 = IC50,
                                             HILL = HILL)), data.frame(time = 0, cmt = 1, amt = 1000, rate = 1000), tt)
report("template concentration", a$Y1, b$CP)
report("template response", a$Y2, b$RESP)

# Corticosteroid: methylprednisolone 30 mg, and placebo against the baseline
x <- setup("corticosteroid_hong2007")
tt <- c(0.25, 1, 2, 4, 6, 8, 12, 16, 24, 32)
a <- sa_profile(x$m, x$p, dose_rec(0, 1, 30000), tt)
b <- mrg("corticosteroid_hong2007", x$truth[1:10], data.frame(time = 0, cmt = 1, amt = 30000), tt)
report("corticosteroid methylprednisolone", a$Y1, b$CMPL)
report("corticosteroid cortisol", a$Y2, b$CORT)
report("corticosteroid lymphocytes", a$Y3, b$LYM)
report("corticosteroid proliferation", a$Y4, b$WBLP)
a <- sa_profile(x$m, x$p, no_dose, tt)
report("corticosteroid placebo cortisol (series)", a$Y2, 80 + 53.13 * cos(2 * pi * tt / 24) - 14.23 * sin(2 * pi * tt / 24))

# Antibody: the lowest and highest dose of the design
x <- setup("mab_mpbpk_tmdd")
tt <- c(1.5, 4, 8, 24, 72, 168, 336, 504, 672, 840, 1008)
for (mgkg in c(0.5, 20)) {
  nmol <- mgkg * 70 * 1e6 / 148000
  a <- sa_profile(x$m, x$p, dose_rec(0, 1, nmol, nmol / 1.5), tt)
  b <- mrg("mab_mpbpk_tmdd", with(x$truth, list(CLp_kg = CLP * 1000 / 70, SIGMA1 = SIGMA1, SIGMA2 = SIGMA2, KSS = KSS,
                                                KSYN = KSYN)),
           data.frame(time = 0, cmt = "CENT", amt = nmol, rate = nmol / 1.5), tt)
  report(sprintf("antibody %.1f mg/kg", mgkg), a$Y1, b$CP_UGML)
}

# T-DM1 in mice: vehicle, 1 and 15 mg/kg
x <- setup("adc_tdm1")
tt <- c(0.02, 1, 3, 7, 10, 14, 21, 28)
mouse <- with(x$truth, list(BW = 0.025, CLADC = CLADC, CLDADC = 0.118, V1ADC = V1ADC, V2ADC = 0.0948, CLDRUG = 11.29,
                            CLDDRUG = 155, V1DRUG = 3.30, V2DRUG = 2.01, KDEC = KDEC, GLIN = 0, DTEXP = DTEXP,
                            KKILL = KKILL, TV0 = TV0 * 1e-6, AG = 1660))
for (mgkg in c(0, 1, 15)) {
  nmol <- mgkg * 0.025 * 1e6 / 148500
  a <- sa_profile(x$m, x$p, if (mgkg > 0) rbind(dose_rec(0, 1, nmol), dose_rec(0, 3, nmol)) else no_dose, tt)
  b <- mrg("adc_tdm1_singh2017", mouse, if (mgkg > 0) data.frame(time = 0, cmt = c("X1TT", "X1ADC"), amt = nmol), tt)
  report(sprintf("T-DM1 %g mg/kg tumour volume", mgkg), a$Y3, b$TV_MM3)
  if (mgkg > 0) {
    report(sprintf("T-DM1 %g mg/kg total antibody", mgkg), a$Y1, b$TT_UGML)
    report(sprintf("T-DM1 %g mg/kg conjugate", mgkg), a$Y2, b$ADC_UGML)
  }
}

# CAR-T: one year
x <- setup("cart_stein2019")
tt <- c(1, 3, 5, 7, 9, 11, 14, 21, 28, 60, 180, 365)
a <- sa_profile(x$m, x$p, no_dose, tt)
b <- mrg("cart_tisagenlecleucel_stein2019", with(x$truth, list(FOLDX = FOLDX, TMAX = TPEAK, CMAX = CPEAK, ALPHA = KALPHA,
                                                               BETA = KBETA, FB = FPERS)), NULL, tt)
report("CAR-T transgene", exp(a$Y1), b$CART)

# PROTAC: the six arms of one experiment, from below the optimum to far
# above it (hook effect)
x <- setup("protac_kcat")
tt <- c(0.5, 1, 2, 4, 8, 16, 24)
conc <- c(1, 10, 100, 1000, 10000, 100000)
a <- sa_profile(x$m, x$p, dose_rec(0, seq_along(conc), conc), tt)
for (k in seq_along(conc)) {
  b <- mrg("protac_degrader_kcat2023", with(x$truth, list(CFIX = conc[k], KCAT = KCAT, THALFP = THALFP)), NULL, tt)
  report(sprintf("PROTAC %g nM target", conc[k]), a[[paste0("Y", k)]], 100 * b$PR)
}

cat(sprintf("\nlargest relative difference: %.2e (%s)\n", max(worst), names(which.max(worst))))
if (max(worst) > 1e-3) stop("a control stream disagrees with its mrgsolve reference")
cat("all control streams agree with mrgsolve\n")
