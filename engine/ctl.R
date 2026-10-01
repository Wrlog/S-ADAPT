# ============================================================================
# SADAPT-TRAN model files: reading, checking, and translating to Fortran
# ============================================================================
# A model is two files (Bulitta et al., AAPS J 2011;13:201, Fig. 4 and
# Table I): final_model.ctl, the control stream of simplified Fortran in
# $ blocks, and parameter_settings.csv, the parameters and run settings.
#
# SADAPT-TRAN turns these into the Fortran subroutines S-ADAPT compiles.
# This file does the same job for the reference engine, so that the model
# code S-ADAPT will see is compiled and run here first: numbers go to double
# precision, EXP and LOG are limited as SADAPT-TRAN limits them, parameters
# and the globals of $OUTPUT_GLB are shared through common blocks, and every
# name on the right of an equals sign has to be defined before it is used.

CTL_BLOCKS <- c("DIFFEQ_DIF", "OUTPUT_GLB", "OUTPUT_ICS", "OUTPUT_EQN", "VARMOD_EQN", "POPMOD_EQN")
CTL_ARRAYS <- c("X", "XP", "Y", "V", "R", "B")
CTL_KEYWORDS <- c("IF", "THEN", "ELSE", "ELSEIF", "ENDIF", "END")
CTL_FUNCTIONS <- c("EXP", "DEXP", "LOG", "DLOG", "SQRT", "DSQRT", "ABS", "DABS", "MIN", "MAX", "DMIN1", "DMAX1",
                   "SIN", "DSIN", "COS", "DCOS", "TANH", "DTANH", "LOG10", "DLOG10")
# Lines are kept short because SADAPT-TRAN writes them into fixed-form
# Fortran (72 columns) after indenting and expanding the numbers.
CTL_MAX_LINE <- 60

#' Read a control stream into its title and blocks (lines of code, blanks dropped)
read_ctl <- function(file) {
  lines <- sub("\\s+$", "", readLines(file, warn = FALSE))
  head <- grep("^\\$", lines)
  if (!length(head)) stop(file, ": no $ blocks found")
  out <- list(file = file, project = "", blocks = setNames(vector("list", length(CTL_BLOCKS)), CTL_BLOCKS))
  for (i in seq_along(head)) {
    name <- toupper(sub("^\\$(\\S+).*$", "\\1", lines[head[i]]))
    last <- if (i < length(head)) head[i + 1] - 1 else length(lines)
    body <- lines[seq_len(last - head[i]) + head[i]]
    body <- body[nzchar(trimws(body))]
    if (name == "PROJECT") {
      out$project <- trimws(sub("^\\$\\S+", "", lines[head[i]]))
    } else if (name %in% CTL_BLOCKS) {
      out$blocks[[name]] <- body
    } else {
      stop(file, ": unknown block $", name)
    }
  }
  for (b in CTL_BLOCKS) if (is.null(out$blocks[[b]])) out$blocks[[b]] <- character(0)
  out
}

#' Read parameter_settings.csv: the parameter table, then the COVARIATES,
#' RUN_SETTINGS and RUN_COMMANDS sections (Table I of the paper)
read_settings <- function(file) {
  raw <- utils::read.csv(file, header = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                         na.strings = character(0), strip.white = TRUE, fill = TRUE, comment.char = "")
  cols <- sub("^#", "", unname(unlist(raw[1, ])))
  cols <- cols[nzchar(cols)]
  want <- c("PNAME", "PMEAN", "PCOV", "PTYPE", "PBLOCK", "PTRANSF", "PLOW", "PHIGH", "VARINI", "VARBURN", "VARNIT",
            "VARIT", "PBOUNDL", "PBOUNDH")
  if (!identical(cols, want)) stop(file, ": header must be #", paste(want, collapse = ","))
  key <- raw[[1]]
  marks <- sapply(c("COVARIATES", "RUN_SETTINGS", "RUN_COMMANDS"), function(k) match(k, key))
  if (anyNA(marks)) stop(file, ": needs COVARIATES, RUN_SETTINGS and RUN_COMMANDS rows")
  par <- raw[seq(2, marks[["COVARIATES"]] - 1), seq_along(want)]
  names(par) <- want
  num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("", "."), NA, x)))
  for (k in setdiff(want, c("PNAME", "PTYPE", "PTRANSF"))) par[[k]] <- num(par[[k]])
  rownames(par) <- NULL
  if (anyDuplicated(toupper(par$PNAME))) stop(file, ": duplicated parameter name")
  if (!all(par$PTYPE %in% c("P", "V"))) stop(file, ": PTYPE must be P or V")
  if (!all(par$PTRANSF[par$PTYPE == "P"] %in% c("L", "N", "O"))) stop(file, ": PTRANSF must be L, N or O")
  is_o <- par$PTYPE == "P" & par$PTRANSF == "O"
  par$PLOW[is_o & is.na(par$PLOW)] <- 0
  par$PHIGH[is_o & is.na(par$PHIGH)] <- 1
  between <- function(a, b) if (b - a > 1) raw[seq(a + 1, b - 1), , drop = FALSE] else raw[0, , drop = FALSE]
  set <- between(marks[["RUN_SETTINGS"]], marks[["RUN_COMMANDS"]])
  cmd <- between(marks[["RUN_COMMANDS"]], nrow(raw) + 1)
  list(file = file, par = par,
       covariates = between(marks[["COVARIATES"]], marks[["RUN_SETTINGS"]])[[1]],
       settings = setNames(as.list(set[[2]]), set[[1]]),
       commands = cmd[[1]][nzchar(cmd[[1]])])
}

# ---- tokens ----------------------------------------------------------------

CTL_TOKEN <- paste0("\\.[A-Za-z]+\\.",                                         # .GT. .AND. ...
                    "|[A-Za-z][A-Za-z0-9_]*",                                  # names
                    "|[0-9]+(?:\\.(?![A-Za-z]+\\.)[0-9]*)?(?:[EeDd][+-]?[0-9]+)?", # 1  1.5  1E-6  (not the 1. of 1.GT.)
                    "|\\.[0-9]+(?:[EeDd][+-]?[0-9]+)?",                        # .5
                    "|\\*\\*|[-+*/(),=]|\\s+")

ctl_tokens <- function(line, where) {
  tok <- regmatches(line, gregexpr(CTL_TOKEN, line, perl = TRUE))[[1]]
  if (paste(tok, collapse = "") != line) stop(where, ": cannot read `", line, "`")
  tok
}

is_name <- function(tok) grepl("^[A-Za-z]", tok)
is_number <- function(tok) grepl("^[0-9.]", tok) & !grepl("^\\.[A-Za-z]", tok)

#' A number as a double precision constant: 2 -> 2.0D0, 0.7 -> 0.7D0, 1E-6 -> 1.0D-6
as_double <- function(tok) {
  m <- regmatches(tok, regexec("^([0-9.]+)(?:[EeDd]([+-]?[0-9]+))?$", tok))[[1]]
  mant <- m[2]
  if (!grepl(".", mant, fixed = TRUE)) mant <- paste0(mant, ".0")
  if (grepl("\\.$", mant)) mant <- paste0(mant, "0")
  if (grepl("^\\.", mant)) mant <- paste0("0", mant)
  paste0(mant, "D", if (nzchar(m[3])) m[3] else "0")
}

#' One line of control-stream code as double precision Fortran. Subscripts
#' of the S-ADAPT arrays stay integer.
ctl_line_fortran <- function(line, where) {
  tok <- ctl_tokens(line, where)
  solid <- which(!grepl("^\\s+$", tok))
  index <- logical(0)                      # one entry per open bracket: is it a subscript?
  for (k in seq_along(solid)) {
    i <- solid[k]
    t <- tok[i]
    prev <- if (k > 1) tok[solid[k - 1]] else ""
    if (t == "(") {
      index <- c(index, toupper(prev) %in% CTL_ARRAYS)
    } else if (t == ")") {
      if (!length(index)) stop(where, ": unbalanced brackets in `", line, "`")
      index <- index[-length(index)]
    } else if (is_number(t)) {
      if (prev == "**") stop(where, ": write a product or a named exponent, not `**", t, "`, in `", line, "`")
      if (!(length(index) && index[length(index)])) tok[i] <- as_double(t)
    } else if (is_name(t)) {
      if (toupper(t) %in% c("EXP", "DEXP")) tok[i] <- "SAEXP"
      if (toupper(t) %in% c("LOG", "DLOG")) tok[i] <- "SALOG"
    }
  }
  if (length(index)) stop(where, ": unbalanced brackets in `", line, "`")
  paste(tok, collapse = "")
}

#' The names a line assigns to and the names it reads
ctl_line_names <- function(line, where) {
  tok <- ctl_tokens(line, where)
  tok <- tok[!grepl("^\\s+$", tok)]
  depth <- cumsum(tok == "(") - cumsum(tok == ")")
  eq <- which(tok == "=" & depth == 0)
  names_in <- function(t) toupper(t[is_name(t)])
  if (!length(eq)) return(list(sets = character(0), reads = names_in(tok)))
  # left of the equals sign: after the closing bracket of a one-line IF, if any
  start <- 1
  if (toupper(tok[1]) == "IF") start <- which(depth == 0 & tok == ")")[1] + 1
  lhs <- tok[start:(eq[1] - 1)]
  rest <- c(tok[seq_len(start - 1)], if (length(lhs) > 1) lhs[-1], tok[-seq_len(eq[1])])
  list(sets = toupper(lhs[1]), reads = names_in(rest))
}

# ---- checks ----------------------------------------------------------------

#' Check a model the way SADAPT-TRAN's pre-processor does, and return its
#' dimensions: number of differential equations, outputs, and the globals.
check_ctl <- function(ctl, par) {
  pn <- toupper(par$PNAME)
  clash <- intersect(pn, c(CTL_ARRAYS, CTL_KEYWORDS, CTL_FUNCTIONS, "T", "P"))
  if (length(clash)) stop(ctl$file, ": parameter name is reserved: ", paste(clash, collapse = ", "))
  # what each block may read besides parameters, globals and its own locals
  may_read <- list(OUTPUT_GLB = character(0), OUTPUT_ICS = "X", DIFFEQ_DIF = c("X", "XP", "R", "T"),
                   OUTPUT_EQN = c("X", "Y", "T"), VARMOD_EQN = c("X", "Y", "V", "T"), POPMOD_EQN = character(0))
  may_set <- list(OUTPUT_GLB = character(0), OUTPUT_ICS = "X", DIFFEQ_DIF = "XP", OUTPUT_EQN = c("X", "Y"),
                  VARMOD_EQN = "V", POPMOD_EQN = character(0))
  globals <- character(0)
  top <- c(XP = 0L, Y = 0L, V = 0L)
  for (b in c("OUTPUT_GLB", "OUTPUT_ICS", "DIFFEQ_DIF", "OUTPUT_EQN", "VARMOD_EQN")) {
    locals <- character(0)
    for (n in seq_along(ctl$blocks[[b]])) {
      line <- ctl$blocks[[b]][n]
      where <- sprintf("%s $%s line %d", basename(dirname(ctl$file)), b, n)
      if (nchar(line) > CTL_MAX_LINE) stop(where, ": longer than ", CTL_MAX_LINE, " characters")
      if (grepl("[!&;]", line)) stop(where, ": comments, continuation lines and `;` are not used here")
      nm <- ctl_line_names(line, where)
      known <- c(pn, globals, locals, may_read[[b]], CTL_KEYWORDS, CTL_FUNCTIONS)
      # a variance parameter belongs to the residual error model only
      if (b != "VARMOD_EQN") known <- setdiff(known, pn[par$PTYPE == "V"])
      unknown <- setdiff(nm$reads, known)
      if (length(unknown)) stop(where, ": `", paste(unknown, collapse = "`, `"), "` is used before it is defined")
      if (length(nm$sets)) {
        if (nm$sets %in% pn) stop(where, ": assigns to the parameter ", nm$sets)
        if (nm$sets %in% CTL_ARRAYS) {
          if (!nm$sets %in% may_set[[b]]) stop(where, ": ", nm$sets, " cannot be set in $", b)
        } else if (b == "OUTPUT_GLB") {
          globals <- union(globals, nm$sets)
        } else {
          if (nm$sets %in% globals) stop(where, ": assigns to the global ", nm$sets, " outside $OUTPUT_GLB")
          locals <- union(locals, nm$sets)
        }
      }
      for (a in names(top)) {
        idx <- regmatches(line, gregexpr(sprintf("(?i)(?<![A-Za-z0-9_])%s\\(\\s*[0-9]+\\s*\\)", a), line, perl = TRUE))[[1]]
        if (length(idx)) top[[a]] <- max(top[[a]], as.integer(gsub("[^0-9]", "", idx)))
      }
    }
  }
  if (top[["Y"]] < 1) stop(ctl$file, ": $OUTPUT_EQN defines no output")
  if (top[["V"]] != top[["Y"]]) stop(ctl$file, ": $VARMOD_EQN must give one variance per output")
  list(neq = max(top[["XP"]], 1L), nout = top[["Y"]], globals = globals)
}

# ---- Fortran ---------------------------------------------------------------

#' The model as fixed-form Fortran: the routines engine/driver.f calls
ctl_fortran <- function(ctl, par) {
  dims <- check_ctl(ctl, par)
  pn <- toupper(par$PNAME)
  body <- function(b) {
    lines <- ctl$blocks[[b]]
    where <- sprintf("%s $%s", basename(dirname(ctl$file)), b)
    if (!length(lines)) return(character(0))
    paste0("      ", vapply(lines, function(l) ctl_line_fortran(trimws(l), where), ""))
  }
  shared <- c("      IMPLICIT DOUBLE PRECISION (A-Z)",
              paste0("      COMMON /SAPARS/ ", pn),
              if (length(dims$globals)) paste0("      COMMON /SAGLBS/ ", dims$globals),
              "      COMMON /SAINP/ R(50)")
  routine <- function(head, decl, code) c(head, shared, decl, code, "      RETURN", "      END", "")
  c(sprintf("C     %s", ctl$project),
    "C     Generated by engine/ctl.R from final_model.ctl. Do not edit.", "",
    routine("      SUBROUTINE SAPAR(P)", "      DIMENSION P(*)", sprintf("      %s = P(%d)", pn, seq_along(pn))),
    routine("      SUBROUTINE SAGLB", NULL, body("OUTPUT_GLB")),
    routine("      SUBROUTINE SAICS(X)", "      DIMENSION X(*)", body("OUTPUT_ICS")),
    routine("      SUBROUTINE DIFFEQ(T,X,XP)", "      DIMENSION X(*), XP(*)", body("DIFFEQ_DIF")),
    routine("      SUBROUTINE OUTPUT(Y,T,X)", "      DIMENSION Y(*), X(*)", body("OUTPUT_EQN")),
    routine("      SUBROUTINE VARMOD(V,T,X,Y)", "      DIMENSION V(*), X(*), Y(*)", body("VARMOD_EQN")))
}
