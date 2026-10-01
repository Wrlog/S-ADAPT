# Running the models in S-ADAPT

The model files in `models/` are written for S-ADAPT with the SADAPT-TRAN
pre- and post-processor. Neither program is in this repository: S-ADAPT's
licence does not allow it to be passed on. This page says how to get them,
how to run a model, and which points to check on the first run.

Nothing here has been run in S-ADAPT yet. The model files follow the syntax
published for SADAPT-TRAN, and the model code in them compiles and gives
the right answers in the reference engine, but the first S-ADAPT run is
still a first run. The checklist at the end lists every assumption that
could not be confirmed from the published material.

## What to install

1. **S-ADAPT.** Fill in the form at
   <https://bmsr.usc.edu/downloads/s-adapt/s-adapt-download-form/>
   (Biomedical Simulations Resource, University of Southern California).
   The terms include not distributing S-ADAPT, or the ADAPT modules in it,
   without the permission of the BMSR. S-ADAPT is written and maintained
   by Robert J. Bauer.
2. **A Fortran compiler.** S-ADAPT compiles each model. The BMSR page lists
   Compaq Visual Fortran 6.x, Intel Fortran 9.x to 11.x and g77, on Windows
   NT to XP and on Linux. Whether it builds with a current Intel compiler or
   with gfortran on Windows 11 is the first thing to find out; ask the BMSR
   or Dr. Bauer which compiler they recommend now.
3. **SADAPT-TRAN**, the Perl scripts of Bulitta, Bingölbali, Shin and
   Landersdorfer (*AAPS J* 2011;13:201-211), distributed through the BMSR.
   It needs Perl, and R for the diagnostic plots.

## Running one model

A model is a folder with two files, as in the paper:

- `final_model.ctl`: the control stream. `$DIFFEQ_DIF` holds the
  differential equations, `$OUTPUT_GLB` the constants and derived
  quantities, `$OUTPUT_ICS` the initial conditions, `$OUTPUT_EQN` the
  outputs and `$VARMOD_EQN` the residual variance of each output.
- `parameter_settings.csv`: one row per parameter (initial mean and
  variance, type, block of the variance matrix, transformation, variance
  burn settings, bounds), then the run settings and the run commands.

In the model's folder, type:

```
sago
```

SADAPT-TRAN translates the two files to Fortran, compiles the model, writes
`run.txt` and runs it in S-ADAPT (loading the data, setting the initial
estimates, estimating with `piteraten`, then the standard errors, the
post-hoc step and the exports of `finish.txt`), and finally runs its
evaluation and plotting scripts, `saev` and `sapl`. The results are csv
files named `<NFILE>_<VERS>_mean.csv`, `_iter.csv`, `_par.csv`,
`_ipred.csv`, `_varc.csv`, `_var.csv`, `_pred.csv` and `_hoc.csv`, and pdf
files of diagnostic plots.

The data are in `data/<name>.csv`, which the settings file points to as
`../../data/<name>.csv`.

## What to compare

`results/reference/<name>/estimates.csv` holds, for the same data and the
same initial estimates, the values the data were simulated from and the
estimates of the reference MC-PEM in this repository. S-ADAPT's population
means (`_mean.csv`) and variances (`_var.csv`) should land close to both.
The two will not agree to the last digit: both are Monte Carlo methods, and
S-ADAPT's sampling is more refined. S-ADAPT also gives standard errors
(`_varc.csv`), which the reference does not.

Put S-ADAPT's exported csv files in `results/sadapt/<name>/`.

## First-run checklist

These follow from the published examples but are not stated in them. Each
is quick to check on the first run, and SADAPT-TRAN's error messages will
point at most of them.

**The dataset**

- [ ] *Column names.* The data files use `ID, TIME, DV, AMT, RATE, EVID,
      MDV, CMT`. S-ADAPT reads data in the NONMEM layout; Dr. Bauer's
      examples name the observation and dose columns `CONC` and `DOSE`.
      Rename the columns if S-ADAPT asks for other names.
- [ ] *CMT.* On a dose row it is the compartment dosed; on an observation
      row it is the number of the output, as in Dr. Bauer's two-output
      example.
- [ ] *Bolus and infusion.* A dose row with `RATE = 0` is taken as a bolus
      into `X(CMT)`. A row with `RATE > 0` is an infusion that the model
      reads as `R(CMT)`; only the template and the antibody model use one,
      and both have `R(1)` in their first equation.
- [ ] *A sample at a dosing time* is a pre-dose sample, and its row comes
      before the dose row.
- [ ] *No dose rows.* The CAR-T data and the placebo subjects of the
      corticosteroid data have observations only.
- [ ] *Line endings.* The files in this repository end lines with LF. If
      S-ADAPT or SADAPT-TRAN on Windows misreads a file, convert it to
      CRLF.
- [ ] *The covariate file.* The published example loads a covariate file
      (`COVFILE`) even with no covariates in the model, and does not show
      its layout. `data/<name>_cov.csv` holds `ID` and `GROUP` (the dose or
      dose group). If SADAPT-TRAN wants another layout, start from the
      covariate file in its own examples.

**The control stream**

- [ ] *Time in the equations.* The corticosteroid model uses `T` in
      `$DIFFEQ_DIF` for the circadian input.
- [ ] *An empty `$OUTPUT_ICS` block*, in the CAR-T model, where every state
      starts at zero.
- [ ] *A block IF in `$OUTPUT_EQN`*, in the CAR-T model. The published
      example has a block IF in `$DIFFEQ_DIF` and one-line IFs in
      `$OUTPUT_EQN`.
- [ ] *Reserved names.* SADAPT-TRAN refuses variable names that S-ADAPT
      uses itself. The names here avoid the obvious ones; rename any it
      objects to.

**The settings**

- [ ] *ATOL and RTOL* are taken to be the exponents of the tolerances
      (`8` for 1e-8), as the published values of 5 suggest.
- [ ] *Parallel run.* The settings keep the `BEOFILE` and `beosetup` rows of
      the published example, which set up S-ADAPT's parallel run on one
      computer.
- [ ] *NPOPITER.* The number of iterations in each settings file is what
      the reference MC-PEM needed on that dataset. Look at SADAPT-TRAN's
      convergence plots before trusting any run, and raise it if a mean is
      still moving.
