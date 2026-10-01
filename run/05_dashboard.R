# Build the results dashboard: a static page, dashboard/index.html, and the
# data it shows, gathered here from the model files, the datasets and
# results/reference/. Prediction curves between the sampling times are
# solved with the reference engine from each subject's estimated
# parameters. Output goes to site/, which GitHub Pages serves.
#
# Run from the repository root: Rscript run/05_dashboard.R

for (f in c("ctl.R", "build.R", "simulate.R")) source(file.path("engine", f))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed")

MODELS <- c("idr1_template", "corticosteroid_hong2007", "mab_mpbpk_tmdd", "adc_tdm1", "cart_stein2019", "protac_kcat")

# What the page says about each model: a short name for the tab, a summary,
# the source, the time unit, and how to word a subject's group
ABOUT <- list(
  idr1_template = list(
    tab = "Template", time = "h", source = "Bulitta et al., AAPS J 2011;13:201 (Fig. 4, Table I)",
    summary = paste("One-compartment PK with a 1-hour infusion, driving indirect response model I. The control stream is",
                    "the worked example published with SADAPT-TRAN."),
    group = function(g) sprintf("%g mg", g)),
  corticosteroid_hong2007 = list(
    tab = "Corticosteroid", time = "h", source = "Hong, Mager, Blum & Jusko, Pharm Res 2007;24:1088",
    summary = paste("Methylprednisolone suppresses circadian cortisol secretion; cortisol and drug together hold lymphocytes",
                    "out of blood; ex vivo lymphocyte proliferation follows through a transit delay. The original",
                    "analysis was done in S-ADAPT."),
    group = function(g) if (g == 0) "placebo" else sprintf("%g mg", g)),
  mab_mpbpk_tmdd = list(
    tab = "Antibody", time = "h", source = "Cao & Jusko, J Pharmacokinet Pharmacodyn 2014;41:375",
    summary = paste("Trastuzumab in a minimal PBPK model (plasma, tight and leaky interstitial fluid, lymph) with",
                    "target-mediated disposition in the interstitial fluid."),
    group = function(g) sprintf("%g mg/kg", g)),
  adc_tdm1 = list(
    tab = "T-DM1", time = "day", source = "Singh & Shah, AAPS J 2017;19:1054",
    summary = paste("Antibody-drug conjugate in KPL-4 xenograft mice: plasma PK of total antibody and conjugate, tumour",
                    "uptake, HER2 binding, internalisation, DM1 release, and tumour killing."),
    group = function(g) if (g == 0) "vehicle" else sprintf("%g mg/kg", g)),
  cart_stein2019 = list(
    tab = "CAR-T", time = "day", source = "Stein et al., CPT Pharmacometrics Syst Pharmacol 2019;8:285",
    summary = "Tisagenlecleucel transgene in blood: expansion to a peak, fast contraction, slow persistence.",
    group = function(g) ""),
  protac_kcat = list(
    tab = "PROTAC", time = "h", source = "Pharmaceutics 2023;15:195",
    summary = paste("Target protein degradation in cells: ternary-complex engagement with a hook effect, and catalytic",
                    "degradation on top of target turnover. One subject is one experiment at six concentrations."),
    group = function(g) "")
)

r4 <- function(x) signif(x, 4)

#' Times for a smooth curve: an even grid, extra points early where the
#' kinetics are fast, and every sampling time
curve_times <- function(obs_times) {
  tmax <- max(obs_times)
  early <- tmax * 10^seq(-3, -0.5, length.out = 12)
  sort(unique(r4(c(0, obs_times, seq(0, tmax, length.out = 48), early))))
}

model_data <- function(name) {
  model <- load_model(name)
  design <- read_design(name)
  about <- ABOUT[[name]]
  dir <- file.path("results", "reference", name)
  est <- utils::read.csv(file.path(dir, "estimates.csv"))
  hist <- utils::read.csv(file.path(dir, "iterations.csv"), check.names = FALSE)
  post <- utils::read.csv(file.path(dir, "posthoc.csv"))
  pred <- utils::read.csv(file.path(dir, "predictions.csv"))
  cov <- utils::read.csv(data_paths(name)$cov)
  subjects <- read_dataset(name)
  ty <- par_types(model)
  means <- model$par$PNAME[ty$p]
  full <- function(p) {
    out <- numeric(nrow(model$par))
    out[ty$p] <- p
    out[ty$v] <- est$ESTIMATE[est$KIND == "residual"]
    out
  }
  typical <- full(est$ESTIMATE[est$KIND == "mean"])

  # an output is drawn on a log axis when it is positive and spans more
  # than two decades
  logged <- vapply(seq_len(model$nout), function(k) {
    z <- pred[pred$CMT == k, ]
    nrow(z) > 0 && all(pmin(z$DV, z$IPRED, z$PRED) > 0) && max(z$DV) / min(z$DV) > 100
  }, TRUE)

  subj <- lapply(seq_along(subjects), function(i) {
    s <- subjects[[i]]
    obs <- s[s$EVID == 0, ]
    pr <- pred[pred$ID == s$ID[1], ]
    times <- curve_times(obs$TIME)
    doses <- s[s$EVID == 1, c("TIME", "EVID", "CMT", "AMT", "RATE")]
    ip <- sa_profile(model, full(unlist(post[i, means])), doses, times)
    pp <- sa_profile(model, typical, doses, times)
    used <- sort(unique(obs$CMT))
    list(id = s$ID[1], group = about$group(cov$GROUP[i]), t = ip$TIME,
         outputs = lapply(used, function(k) {
           o <- pr[pr$CMT == k, ]
           list(k = k, ipred = r4(ip[[paste0("Y", k)]]), pred = r4(pp[[paste0("Y", k)]]),
                obs = list(t = o$TIME, dv = r4(o$DV), ipred = r4(o$IPRED), pred = r4(o$PRED)))
         }))
  })

  rep_file <- file.path(dir, "replicates.csv")
  replicates <- if (file.exists(rep_file)) {
    z <- utils::read.csv(rep_file)
    z <- z[z$KIND == "mean", ]
    list(n = length(unique(z$REPLICATE)),
         rows = unname(lapply(split(z, z$PARAMETER), function(d) list(parameter = d$PARAMETER[1], ratio = r4(d$ESTIMATE / d$SAMPLE)))))
  }
  text <- function(file) paste(readLines(file.path("models", name, file), warn = FALSE), collapse = "\n")
  list(name = name, tab = about$tab, title = design$title, summary = about$summary, source = about$source,
       time = about$time, outputs = lapply(seq_along(design$outputs), function(k) list(label = design$outputs[k],
                                                                                      log = logged[k])),
       nsubjects = length(subjects), nobs = nrow(pred), niter = nrow(hist), objective = round(mean(utils::tail(hist$OBJ, 10)), 1),
       burn = max(model$par$VARBURN, na.rm = TRUE),
       estimates = lapply(seq_len(nrow(est)), function(i) as.list(est[i, ])),
       history = c(list(ITER = hist$ITER, OBJ = r4(hist$OBJ)), lapply(stats::setNames(means, means), function(k) r4(hist[[k]]))),
       means = means, subjects = subj, replicates = replicates,
       ctl = text("final_model.ctl"), settings = text("parameter_settings.csv"))
}

data <- list(generated = format(Sys.Date()), models = lapply(MODELS, model_data))
dir.create("site", showWarnings = FALSE)
json <- jsonlite::toJSON(data, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
writeLines(c("window.DASHBOARD = ", json, ";"), file.path("site", "data.js"), useBytes = TRUE)
invisible(file.copy(file.path("dashboard", "index.html"), file.path("site", "index.html"), overwrite = TRUE))
cat(sprintf("site/ written: %d models, %d subjects, data.js %.0f kB\n", length(MODELS),
            sum(vapply(data$models, function(m) m$nsubjects, 0)), file.size(file.path("site", "data.js")) / 1024))
