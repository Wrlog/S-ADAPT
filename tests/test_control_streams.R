# The model files, checked the way SADAPT-TRAN's pre-processor checks them,
# and the datasets against the models they belong to. Base R only; nothing
# is compiled.
#
# Run from the repository root: Rscript tests/test_control_streams.R

for (f in c("ctl.R", "build.R", "simulate.R")) source(file.path("engine", f))

fails <- 0
check <- function(what, ok) {
  cat(sprintf("  %-66s %s\n", what, if (isTRUE(ok)) "ok" else "FAILED"))
  if (!isTRUE(ok)) fails <<- fails + 1
}
refused <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

cat("Translation to double precision Fortran\n")
check("2 -> 2.0D0, 0.7 -> 0.7D0, 1E-6 -> 1.0D-6, .5 -> 0.5D0",
      identical(vapply(c("2", "0.7", "1E-6", ".5"), as_double, ""), c(`2` = "2.0D0", `0.7` = "0.7D0", `1E-6` = "1.0D-6",
                                                                      `.5` = "0.5D0")))
check("LOG is limited and its argument made double",
      ctl_line_fortran("Kout = LOG(2)/Tout12", "t") == "Kout = SALOG(2.0D0)/Tout12")
check("subscripts stay integer, comparisons become double",
      ctl_line_fortran("if (X(1).LT.0) X(1) = 0", "t") == "if (X(1).LT.0.0D0) X(1) = 0.0D0")
check("a literal integer power is refused", refused(ctl_line_fortran("A = B**2", "t")))

cat("\nThe pre-processor's checks\n")
tmpl <- read_ctl("models/idr1_template/final_model.ctl")
par <- read_settings("models/idr1_template/parameter_settings.csv")$par
broken <- function(block, line) {
  bad <- tmpl
  bad$blocks[[block]] <- c(bad$blocks[[block]], line)
  refused(check_ctl(bad, par))
}
check("a name used before it is defined is refused", broken("DIFFEQ_DIF", "XP(2) = XP(2) - KLOSS*X(2)"))
check("assigning to a parameter is refused", broken("DIFFEQ_DIF", "CL = 2"))
check("a state in the time-independent block is refused", broken("OUTPUT_GLB", "Kin = X(2)*Kout"))
check("a residual-error parameter in the equations is refused", broken("DIFFEQ_DIF", "XP(1) = XP(1) - SDin"))
check("a line over 60 characters is refused", broken("OUTPUT_EQN", paste0("Y(2) = X(2)", strrep(" + 0*X(2)", 6))))

cat("\nThe template is the published example (Bulitta et al. 2011, Fig. 4 and Table I)\n")
check("differential equations", identical(trimws(tmpl$blocks$DIFFEQ_DIF), c(
  "DC1 = X(1)/V1", "IF (DC1.GT.0) THEN", "INH  = Imax*DC1**HILL / (DC1**HILL + IC50**HILL)", "ELSE", "INH  = 0", "ENDIF",
  "XP(1) = R(1) -CL*DC1", "XP(2) = Kin * (1 - INH) - Kout * X(2)")))
check("parameter names and initial means", identical(par$PNAME, c("CL", "V1", "Tout12", "Imax", "BASE", "IC50", "HILL", "SDin",
                                                                    "SDsl", "PDin", "PDsl")) &&
        identical(par$PMEAN, c(5, 30, 10, 0.7, 15, 10, 1, 1, 1, 1, 1)))
check("Imax is logistic between 0 and 1; HILL shrinks to variance 0.01",
      par$PTRANSF[4] == "O" && par$PLOW[4] == 0 && par$PHIGH[4] == 1 && par$VARNIT[7] == 0.01 && par$VARIT[7] == 40)

expected <- list(idr1_template = c(2, 2), corticosteroid_hong2007 = c(4, 4), mab_mpbpk_tmdd = c(6, 1), adc_tdm1 = c(14, 3),
                 cart_stein2019 = c(1, 1), protac_kcat = c(12, 6))
for (name in names(expected)) {
  cat("\n", name, "\n", sep = "")
  ctl <- read_ctl(file.path("models", name, "final_model.ctl"))
  set <- read_settings(file.path("models", name, "parameter_settings.csv"))
  dims <- tryCatch(check_ctl(ctl, set$par), error = function(e) conditionMessage(e))
  check("control stream passes the checks", is.list(dims))
  if (!is.list(dims)) { cat("   ", dims, "\n"); next }
  check(sprintf("%d differential equations, %d outputs", expected[[name]][1], expected[[name]][2]),
        identical(c(dims$neq, dims$nout), as.integer(expected[[name]])))
  check("no comments in the control stream", !any(grepl("^\\s*[Cc*!]\\s", unlist(ctl$blocks))))
  check("settings name the model, its data and the run commands",
        set$settings$NFILE == name && set$settings$DATAFILE == sprintf("../../data/%s.csv", name) &&
          identical(set$commands, c("piteraten", "?finish.txt", "stop")))
  design <- read_design(name)
  check("design gives a value for every parameter", setequal(names(design$truth), set$par$PNAME))
  check("design labels every output", length(design$outputs) == dims$nout)

  d <- utils::read.csv(data_paths(name)$data)
  obs <- d$EVID == 0
  check("dataset columns are ID, TIME, DV, AMT, RATE, EVID, MDV, CMT",
        identical(names(d), c("ID", "TIME", "DV", "AMT", "RATE", "EVID", "MDV", "CMT")))
  check("observations name an output, doses a compartment",
        all(d$CMT[obs] %in% seq_len(dims$nout)) && all(d$CMT[!obs] %in% seq_len(dims$neq)))
  check("MDV is 1 on dose rows only", all(d$MDV == as.integer(!obs)))
  check("each subject's rows are in time order, samples before a dose at the same time",
        all(vapply(split(d, d$ID), function(s) !is.unsorted(s$TIME) && !is.unsorted(order(s$TIME, s$EVID)), TRUE)))
  check("an infusion appears as R() in the equations",
        !any(d$RATE > 0) || all(vapply(unique(d$CMT[d$RATE > 0]), function(k)
          any(grepl(sprintf("R(%d)", k), ctl$blocks$DIFFEQ_DIF, fixed = TRUE)), TRUE)))
  cov <- utils::read.csv(data_paths(name)$cov)
  check("covariate file has one row per subject", identical(cov$ID, unique(d$ID)))
}

cat("\n")
if (fails) stop(fails, " check(s) failed")
cat("all checks passed\n")
