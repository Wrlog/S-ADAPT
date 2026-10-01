# ============================================================================
# Reference MC-PEM: importance-sampling expectation maximisation
# ============================================================================
# The algorithm S-ADAPT runs for PMETHOD 4 (Bauer & Guzy, in Advanced Methods
# of Pharmacokinetic and Pharmacodynamic Systems Analysis vol. 3, 2004),
# written out plainly so that every model in this repository can be fitted
# before S-ADAPT is installed, and S-ADAPT's answer checked afterwards.
#
# Each subject's parameters phi_i are normal on the transformed scale,
# phi_i ~ N(mu, Omega). One iteration:
#   E step  for every subject, draw phi from a proposal built around that
#           subject's conditional distribution, weight each draw by
#           likelihood x prior / proposal, and form the conditional mean
#           m_i and covariance C_i;
#   M step  mu = mean(m_i), Omega = mean((m_i - mu)(m_i - mu)' + C_i),
#           and the residual variance parameters maximise the weighted
#           log-likelihood of the draws.
# The first proposals come from each subject's posterior mode and its
# curvature (the MAP assistance of PMETHOD 4). The parameter settings file
# drives it as it drives S-ADAPT: PMEAN and PCOV start the search, PBLOCK
# shapes Omega, VARINI/VARBURN hold the variances up for the first
# iterations and VARNIT/VARIT shrink them, PBOUNDL/PBOUNDH bound the draws.
#
# This is not S-ADAPT. It leaves out S-ADAPT's sampling refinements and its
# standard errors, and its objective function is on its own scale.

ldmvn <- function(X, mean, R) {
  z <- backsolve(R, t(X) - mean, transpose = TRUE)
  -0.5 * colSums(z^2) - sum(log(diag(R))) - 0.5 * length(mean) * log(2 * pi)
}
rmvn <- function(n, mean, R) sweep(matrix(stats::rnorm(n * length(mean)), n) %*% R, 2, mean, "+")
safe_chol <- function(S) {
  S <- (S + t(S)) / 2
  for (ridge in c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
    R <- tryCatch(chol(S + diag(ridge * max(diag(S), 1e-12), nrow(S))), error = function(e) NULL)
    if (!is.null(R)) return(R)
  }
  chol(diag(diag(S) + 1e-8, nrow(S)))
}
log_mean_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(mean(exp(x - m)))
}

#' Log-likelihood of one subject's observations for each row of `phi`,
#' with the predictions that produced it
subject_loglik <- function(model, subj, phi, sigma, ty) {
  par <- model$par
  P <- matrix(0, nrow(phi), nrow(par))
  P[, ty$p] <- from_phi(par[ty$p, ], phi)
  P[, ty$v] <- rep(sigma, each = nrow(phi))
  lo <- par$PBOUNDL[ty$p]
  hi <- par$PBOUNDH[ty$p]
  inside <- rep(TRUE, nrow(phi))
  for (j in seq_along(ty$p)) {
    if (!is.na(lo[j])) inside <- inside & P[, ty$p[j]] >= lo[j]
    if (!is.na(hi[j])) inside <- inside & P[, ty$p[j]] <= hi[j]
  }
  sol <- sa_solve(model, P, subj)
  obs <- subj$EVID == 0
  y <- sol$y[obs, , drop = FALSE]
  v <- sol$v[obs, , drop = FALSE]
  ll <- colSums(matrix(stats::dnorm(subj$DV[obs], y, sqrt(v), log = TRUE), nrow(y)))
  ll[!inside | sol$failed | !is.finite(ll)] <- -Inf
  list(ll = ll, y = y)
}

#' The posterior mode of one subject and the curvature there
subject_map <- function(model, subj, mu, Rprior, sigma, ty) {
  nlp <- function(phi) {
    v <- -(subject_loglik(model, subj, matrix(phi, 1), sigma, ty)$ll + ldmvn(matrix(phi, 1), mu, Rprior))
    if (is.finite(v)) v else 1e10
  }
  # a simplex search to get near the mode from a poor start, then a
  # gradient search to settle on it
  fit <- list(par = mu, value = nlp(mu))
  for (method in c("Nelder-Mead", "BFGS")) {
    ctrl <- if (method == "BFGS") list(maxit = 100, reltol = 1e-8) else list(maxit = 200 * length(mu), reltol = 1e-6)
    try <- tryCatch(stats::optim(fit$par, nlp, method = method, control = ctrl), error = function(e) NULL)
    if (!is.null(try) && try$value < fit$value) fit <- try
  }
  H <- tryCatch(stats::optimHess(fit$par, nlp), error = function(e) NULL)
  C <- if (is.null(H)) NULL else tryCatch(solve((H + t(H)) / 2), error = function(e) NULL)
  if (is.null(C) || any(!is.finite(C)) || any(diag(C) <= 0)) C <- crossprod(Rprior) / 4
  list(mean = fit$par, cov = C)
}

#' Omega for this iteration: the block structure, then the variance burn
#' and shrink schedule of the parameter settings file
shape_omega <- function(Omega, par, iter) {
  Omega <- Omega * outer(par$PBLOCK, par$PBLOCK, "==")
  for (j in seq_len(nrow(par))) {
    burn <- if (is.na(par$VARBURN[j])) 0 else par$VARBURN[j]
    held <- NA
    if (iter <= burn) {
      held <- par$VARINI[j]
    } else if (!is.na(par$VARNIT[j])) {
      step <- min(1, (iter - burn) / max(par$VARIT[j], 1))
      held <- par$VARINI[j] + step * (par$VARNIT[j] - par$VARINI[j])
    }
    if (!is.na(held)) {
      Omega[j, ] <- 0
      Omega[, j] <- 0
      Omega[j, j] <- held
    }
  }
  Omega
}

#' Fit a population model. `subjects` is the list read_dataset() returns.
#' Returns the estimates (averaged over the last `navg` iterations), the
#' iteration history, and each subject's conditional mean and predictions.
mcpem <- function(model, subjects, niter = 60, nsamp = 300, navg = 10, seed = 1, prior_share = 0.15, inflate = 2,
                  nresid = 50, verbose = TRUE) {
  set.seed(seed)
  ty <- par_types(model)
  par <- model$par[ty$p, ]
  np <- length(ty$p)
  ns <- length(subjects)
  mu <- to_phi(par, par$PMEAN)
  Omega <- shape_omega(diag(par$PCOV, np), par, 1)
  sigma <- model$par$PMEAN[ty$v]
  cond <- vector("list", ns)
  hist <- matrix(NA_real_, niter, 2 + 2 * np + length(ty$v) + 1,
                 dimnames = list(NULL, c("ITER", "OBJ", par$PNAME, paste0(par$PNAME, "~", par$PNAME),
                                         model$par$PNAME[ty$v], "ESS")))
  keep <- vector("list", niter)
  for (iter in seq_len(niter)) {
    Rprior <- safe_chol(Omega)
    if (iter == 1) for (i in seq_len(ns)) cond[[i]] <- subject_map(model, subjects[[i]], mu, Rprior, sigma, ty)
    m <- matrix(0, ns, np)
    S <- matrix(0, np, np)
    obj <- 0
    ess <- numeric(ns)
    draws <- vector("list", ns)
    for (i in seq_len(ns)) {
      # proposal: mostly the subject's own conditional distribution, widened,
      # with a share of the population distribution so no weight can run away
      Rq <- safe_chol(inflate * cond[[i]]$cov)
      from_prior <- stats::runif(nsamp) < prior_share
      phi <- rmvn(nsamp, cond[[i]]$mean, Rq)
      if (any(from_prior)) phi[from_prior, ] <- rmvn(sum(from_prior), mu, Rprior)
      lprior <- ldmvn(phi, mu, Rprior)
      lq <- log(prior_share * exp(lprior) + (1 - prior_share) * exp(ldmvn(phi, cond[[i]]$mean, Rq)))
      lik <- subject_loglik(model, subjects[[i]], phi, sigma, ty)
      lw <- lik$ll + lprior - lq
      if (!any(is.finite(lw))) stop("subject ", names(subjects)[i], ": no draw has a finite likelihood")
      obj <- obj - 2 * log_mean_exp(lw)
      w <- exp(lw - max(lw))
      w <- w / sum(w)
      ess[i] <- 1 / sum(w^2)
      mi <- colSums(w * phi)
      d <- sweep(phi, 2, mi)
      Ci <- crossprod(d * sqrt(w))
      m[i, ] <- mi
      S <- S + Ci
      # a conditional covariance from too few effective draws is not trusted
      trust <- min(1, ess[i] / (5 * np))
      cond[[i]] <- list(mean = mi, cov = trust * Ci + (1 - trust) * cond[[i]]$cov + diag(1e-8, np))
      # a resample of the draws, in proportion to their weights, carries the
      # predictions to the residual-error update
      use <- sample.int(nsamp, nresid, replace = TRUE, prob = w)
      draws[[i]] <- list(w = rep(1 / nresid, nresid), y = lik$y[, use, drop = FALSE])
    }
    # M step
    mu <- colMeans(m)
    Omega <- shape_omega((crossprod(sweep(m, 2, mu)) + S) / ns, par, iter + 1)
    sigma <- update_sigma(model, subjects, draws, sigma, ty)
    tv <- from_phi(par, mu)[1, ]
    hist[iter, ] <- c(iter, obj, tv, diag(Omega), sigma, mean(ess))
    keep[[iter]] <- list(mu = mu, Omega = Omega, sigma = sigma)
    if (verbose) cat(sprintf("  iter %3d  obj %12.3f  mean ESS %6.1f  %s\n", iter, obj, mean(ess),
                             paste(sprintf("%s=%.4g", par$PNAME, tv), collapse = " ")))
  }
  last <- keep[seq(max(1, niter - navg + 1), niter)]
  mu <- Reduce(`+`, lapply(last, `[[`, "mu")) / length(last)
  Omega <- Reduce(`+`, lapply(last, `[[`, "Omega")) / length(last)
  sigma <- Reduce(`+`, lapply(last, `[[`, "sigma")) / length(last)
  post <- t(vapply(cond, function(z) z$mean, numeric(np)))
  if (np == 1) post <- matrix(post, ncol = 1)
  list(model = model$name, par = par$PNAME, mu = mu, typical = setNames(from_phi(par, mu)[1, ], par$PNAME),
       Omega = `dimnames<-`(Omega, list(par$PNAME, par$PNAME)), sigma = setNames(sigma, model$par$PNAME[ty$v]),
       objective = mean(hist[seq(max(1, niter - navg + 1), niter), "OBJ"]), history = as.data.frame(hist),
       posthoc = `colnames<-`(from_phi(par, post), par$PNAME), settings = list(niter = niter, nsamp = nsamp, navg = navg,
                                                                               seed = seed))
}

#' The residual variance parameters that maximise the weighted
#' log-likelihood of the draws kept from the E step
update_sigma <- function(model, subjects, draws, sigma, ty) {
  if (!length(ty$v)) return(sigma)
  dv <- cmt <- yy <- ww <- vector("list", length(subjects))
  for (i in seq_along(subjects)) {
    obs <- subjects[[i]]$EVID == 0
    k <- length(draws[[i]]$w)
    dv[[i]] <- rep(subjects[[i]]$DV[obs], k)
    cmt[[i]] <- rep(subjects[[i]]$CMT[obs], k)
    yy[[i]] <- as.vector(draws[[i]]$y)
    ww[[i]] <- rep(draws[[i]]$w, each = sum(obs))
  }
  dv <- unlist(dv); cmt <- unlist(cmt); yy <- unlist(yy); ww <- unlist(ww)
  p <- model$par$PMEAN
  nq <- function(ls) {
    p[ty$v] <- exp(ls)
    v <- sa_variance(model, p, cmt, yy)
    if (any(!is.finite(v)) || any(v <= 0)) return(1e300)
    sum(ww * (log(v) + (dv - yy)^2 / v))
  }
  start <- log(pmax(sigma, 1e-8))
  fit <- if (length(start) == 1) {
    stats::optimize(nq, start + c(-3, 3))$minimum
  } else {
    stats::optim(start, nq, method = "Nelder-Mead", control = list(maxit = 300, reltol = 1e-7))$par
  }
  exp(fit)
}

#' Population and individual predictions for every observation
mcpem_predictions <- function(model, subjects, fit) {
  ty <- par_types(model)
  full <- function(p) {
    out <- numeric(nrow(model$par))
    out[ty$p] <- p
    out[ty$v] <- fit$sigma
    out
  }
  do.call(rbind, lapply(seq_along(subjects), function(i) {
    s <- subjects[[i]]
    obs <- s$EVID == 0
    pred <- sa_solve(model, full(fit$typical), s)
    ipred <- sa_solve(model, full(fit$posthoc[i, ]), s)
    data.frame(ID = s$ID[obs], TIME = s$TIME[obs], CMT = s$CMT[obs], DV = s$DV[obs], PRED = pred$y[obs, 1],
               IPRED = ipred$y[obs, 1], IWRES = (s$DV[obs] - ipred$y[obs, 1]) / sqrt(ipred$v[obs, 1]))
  }))
}
