# ============================================================================
# Simulating a population dataset from a model and a design
# ============================================================================
# models/<name>/design.R defines `design`: the population the data come from
# (truth, bsv), the random seed, and subjects(), which returns one entry per
# subject with its records and a group label. The dataset is written in the
# NONMEM layout S-ADAPT reads, with the true individual parameters beside it.

#' Dose and observation records for one subject. An observation taken at a
#' dosing time is a pre-dose sample, so it sorts ahead of the dose.
dose_rec <- function(time, cmt, amt, rate = 0) {
  data.frame(TIME = time, EVID = 1L, CMT = as.integer(cmt), AMT = amt, RATE = rate)
}
obs_rec <- function(times, cmt) {
  data.frame(TIME = times, EVID = 0L, CMT = as.integer(cmt), AMT = 0, RATE = 0)
}
subject_rec <- function(...) {
  rec <- do.call(rbind, Filter(Negate(is.null), list(...)))
  rec <- rec[order(rec$TIME, rec$EVID, rec$CMT), ]
  rownames(rec) <- NULL
  rec
}

# ---- parameter scales ------------------------------------------------------
# PTRANSF: L log-normal, N normal, O logistic between PLOW and PHIGH

to_phi <- function(par, p) {
  tr <- rep_len(par$PTRANSF, length(p))
  u <- (p - par$PLOW) / (par$PHIGH - par$PLOW)
  ifelse(tr == "L", log(p), ifelse(tr == "O", log(u / (1 - u)), p))
}

#' Rows of `phi` (transformed scale) back to the natural scale
from_phi <- function(par, phi) {
  if (is.null(dim(phi))) phi <- matrix(phi, nrow = 1)
  out <- phi
  for (j in seq_len(ncol(phi))) {
    out[, j] <- switch(par$PTRANSF[j], L = exp(phi[, j]),
                       O = par$PLOW[j] + (par$PHIGH[j] - par$PLOW[j]) / (1 + exp(-phi[, j])), phi[, j])
  }
  out
}

#' Split the parameter table into random-effect (P) and variance (V) rows
par_types <- function(model) list(p = which(model$par$PTYPE == "P"), v = which(model$par$PTYPE == "V"))

# ---- simulation ------------------------------------------------------------

read_design <- function(name) {
  env <- new.env()
  sys.source(file.path(sa_root(), "models", name, "design.R"), envir = env)
  d <- env$design
  for (k in c("seed", "truth", "bsv", "subjects")) if (is.null(d[[k]])) stop("design.R: `", k, "` is missing")
  d
}

#' Draw the subjects, solve the model for each and add residual error
simulate_dataset <- function(model, design) {
  set.seed(design$seed)
  par <- model$par
  ty <- par_types(model)
  pn <- par$PNAME
  if (!setequal(names(design$truth), pn)) stop("design truth must name every parameter")
  if (!setequal(names(design$bsv), pn[ty$p])) stop("design bsv must name every P parameter")
  truth <- design$truth[pn]
  mu <- to_phi(par[ty$p, ], truth[ty$p])
  sd <- sqrt(design$bsv[pn[ty$p]])
  subj <- design$subjects()
  data <- vector("list", length(subj))
  indiv <- matrix(NA_real_, length(subj), length(pn), dimnames = list(NULL, pn))
  for (i in seq_along(subj)) {
    rec <- subj[[i]]$rec
    p <- truth
    p[ty$p] <- from_phi(par[ty$p, ], mu + stats::rnorm(length(mu), 0, sd))
    sol <- sa_solve(model, p, rec)
    if (sol$failed) stop("the model could not be solved for subject ", i)
    obs <- rec$EVID == 0
    dv <- rep(0, nrow(rec))
    y <- sol$y[obs, 1]
    s <- sqrt(sol$v[obs, 1])
    draw <- y + s * stats::rnorm(sum(obs))
    # an assay cannot report a negative amount: draw those again
    positive <- isTRUE(design$positive)
    while (positive && any(bad <- draw <= 0)) draw[bad] <- y[bad] + s[bad] * stats::rnorm(sum(bad))
    dv[obs] <- draw
    data[[i]] <- data.frame(ID = i, TIME = rec$TIME, DV = signif(dv, 5), AMT = signif(rec$AMT, 6),
                            RATE = signif(rec$RATE, 6), EVID = rec$EVID, MDV = as.integer(!obs), CMT = rec$CMT)
    indiv[i, ] <- p
  }
  groups <- vapply(subj, function(s) as.numeric(s$group), 0)
  list(data = do.call(rbind, data), cov = data.frame(ID = seq_along(subj), GROUP = groups),
       truth = data.frame(ID = seq_along(subj), GROUP = groups, signif(indiv[, ty$p, drop = FALSE], 6)))
}

data_paths <- function(name) {
  d <- file.path(sa_root(), "data")
  list(data = file.path(d, paste0(name, ".csv")), cov = file.path(d, paste0(name, "_cov.csv")),
       truth = file.path(d, paste0(name, "_truth.csv")))
}

write_dataset <- function(name, sim) {
  f <- data_paths(name)
  for (k in names(f)) utils::write.csv(sim[[k]], f[[k]], row.names = FALSE, quote = FALSE)
  invisible(f)
}

#' The dataset as a list of per-subject record tables, in file order
read_dataset <- function(name) {
  d <- utils::read.csv(data_paths(name)$data)
  need <- c("ID", "TIME", "DV", "AMT", "RATE", "EVID", "MDV", "CMT")
  if (!identical(names(d), need)) stop("dataset columns must be ", paste(need, collapse = ","))
  split(d, d$ID)
}
