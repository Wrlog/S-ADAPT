# ============================================================================
# Compile a SADAPT-TRAN model for the reference engine, and solve it
# ============================================================================
# load_model() reads models/<name>/, translates the control stream to Fortran
# (engine/ctl.R), compiles it with the driver and LSODA through R CMD SHLIB,
# and loads the library. Compiled files go to a per-user cache, outside the
# repository; set SADAPT_BUILD_DIR to put them somewhere else.

sa_root <- function() {
  root <- getOption("sadapt.root", getwd())
  if (!file.exists(file.path(root, "engine", "driver.f"))) stop("run from the repository root")
  root
}

sa_build_dir <- function() {
  d <- Sys.getenv("SADAPT_BUILD_DIR", tools::R_user_dir("sadapt-reference", "cache"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

sa_shlib <- function(dir, out, files) {
  writeLines("PKG_FFLAGS = -std=legacy -w -ffixed-line-length-132", file.path(dir, "Makevars"))
  owd <- setwd(dir)
  on.exit(setwd(owd))
  log <- suppressWarnings(system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "-o", out, files),
                                  stdout = TRUE, stderr = TRUE))
  if (!file.exists(out)) stop("compilation failed in ", dir, ":\n", paste(log, collapse = "\n"))
  invisible(log)
}

#' LSODA (ODEPACK, public domain) is compiled once and its object files reused
sa_odepack <- function() {
  dir <- file.path(sa_build_dir(), "odepack")
  src <- c("opkdmain.f", "opkda1.f", "opkda2.f")
  obj <- file.path(dir, sub("\\.f$", ".o", src))
  if (!all(file.exists(obj))) {
    dir.create(dir, showWarnings = FALSE)
    file.copy(file.path(sa_root(), "engine", "odepack", src), dir, overwrite = TRUE)
    sa_shlib(dir, paste0("odepack", .Platform$dynlib.ext), src)
  }
  obj
}

#' Read, check, compile and load one model
load_model <- function(name) {
  dir <- file.path(sa_root(), "models", name)
  ctl <- read_ctl(file.path(dir, "final_model.ctl"))
  set <- read_settings(file.path(dir, "parameter_settings.csv"))
  dims <- check_ctl(ctl, set$par)
  code <- ctl_fortran(ctl, set$par)
  driver <- readLines(file.path(sa_root(), "engine", "driver.f"))
  key <- substr(sa_hash(c(code, driver)), 1, 10)
  bdir <- file.path(sa_build_dir(), name)
  dir.create(bdir, showWarnings = FALSE)
  lib <- file.path(bdir, paste0(name, "_", key, .Platform$dynlib.ext))
  if (!file.exists(lib)) {
    obj <- sa_odepack()
    writeLines(code, file.path(bdir, "model.f"))
    writeLines(driver, file.path(bdir, "driver.f"))
    unlink(file.path(bdir, c("model.o", "driver.o")))
    sa_shlib(bdir, basename(lib), c("model.f", "driver.f", file.path("..", "odepack", basename(obj))))
  }
  dll <- dyn.load(lib)
  tol <- function(k, default) 10^-as.numeric(if (is.null(set$settings[[k]])) default else set$settings[[k]])
  list(name = name, project = ctl$project, ctl = ctl, par = set$par, settings = set$settings, commands = set$commands,
       neq = dims$neq, nout = dims$nout, lib = lib, dll = dll[["name"]], fortran = code,
       rtol = tol("RTOL", 6), atol = tol("ATOL", 8))
}

# a content hash without extra packages
sa_hash <- function(lines) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(lines, f)
  unname(tools::md5sum(f))
}

#' Solve one subject's records for each row of P (parameters on their natural
#' scale, columns in the order of the parameter table). Returns predictions
#' and residual variances, records by rows of P; NA where the solver failed.
sa_solve <- function(model, P, rec) {
  if (is.null(dim(P))) P <- matrix(P, nrow = 1)
  if (ncol(P) != nrow(model$par)) stop("P must have one column per parameter")
  if (is.unsorted(rec$TIME)) stop("records must be in time order")
  n <- nrow(rec)
  k <- nrow(P)
  out <- .Fortran("sasim", as.integer(ncol(P)), as.integer(k), as.double(t(P)), as.integer(model$neq), as.integer(n),
                  as.double(rec$TIME), as.integer(rec$EVID), as.integer(rec$CMT), as.double(rec$AMT),
                  as.double(rec$RATE), as.double(model$rtol), as.double(model$atol),
                  y = double(n * k), v = double(n * k), ierr = integer(k), PACKAGE = model$dll)
  y <- matrix(out$y, n, k)
  v <- matrix(out$v, n, k)
  bad <- out$ierr != 0 | colSums(!is.finite(y)) > 0 | colSums(!is.finite(v)) > 0
  y[, bad] <- NA
  v[, bad] <- NA
  list(y = y, v = v, failed = bad)
}

#' Residual variance of predictions `y` of outputs `cmt` under parameters `p`
sa_variance <- function(model, p, cmt, y) {
  .Fortran("savarb", as.integer(length(p)), as.double(p), as.integer(length(y)), as.integer(cmt), as.double(y),
           v = double(length(y)), PACKAGE = model$dll)$v
}

#' Every output on a time grid, for one parameter vector: a data frame with
#' TIME and one column per output. `doses` holds the EVID 1 records.
sa_profile <- function(model, p, doses, times) {
  obs <- expand.grid(TIME = times, CMT = seq_len(model$nout))
  obs <- data.frame(TIME = obs$TIME, EVID = 0L, CMT = obs$CMT, AMT = 0, RATE = 0)
  rec <- rbind(doses[, c("TIME", "EVID", "CMT", "AMT", "RATE")], obs)
  rec <- rec[order(rec$TIME, -rec$EVID), ]
  sol <- sa_solve(model, p, rec)
  keep <- rec$EVID == 0
  wide <- tapply(sol$y[keep, 1], list(rec$TIME[keep], rec$CMT[keep]), function(z) z[1])
  out <- data.frame(TIME = as.numeric(rownames(wide)), unname(wide))
  names(out)[-1] <- paste0("Y", seq_len(model$nout))
  out
}
