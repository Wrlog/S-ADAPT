$PROB
Tisagenlecleucel (CAR-T) cellular kinetics

// Stein AM, Grupp SA, Levine JE, et al. Tisagenlecleucel model-based
// cellular kinetic analysis of chimeric antigen receptor-T cells.
// CPT Pharmacometrics Syst Pharmacol 2019;8:285-295 (open access).
//
// CAR-T cells are a living drug: after infusion they expand exponentially
// (rate RHO) until TMAX, then contract in two phases - a fast one (ALPHA)
// as activated effector cells die, and a very slow one (BETA) carried by a
// small persisting fraction FB (memory-like cells). The published model is
// the analytical function
//   f(t) = R0 exp(RHO t)                                 t < TMAX
//   f(t) = A exp(-ALPHA (t-TMAX)) + B exp(-BETA (t-TMAX)) t >= TMAX
// with R0 = CMAX / FOLDX, A = (1 - FB) CMAX, B = FB CMAX; here it is written
// as the equivalent differential equations: two populations that expand
// together and then contract at their own rates. RHO follows from the fold
// expansion over TMAX. The states are the natural logs of the two levels
// (the paper also fitted in log space), so each rate is piecewise constant.
//
// Comedication (the paper's MLXTRAN model, Supplementary Material): from the
// first dose of tocilizumab (TTOCI) or corticosteroids (TSTER) before TMAX,
// the expansion rate is multiplied by FTOCI or FSTER. Use a time beyond TMAX
// (e.g. 99999) for a patient who did not receive the drug. The estimates
// (1.2 and 1.0) mean neither slowed expansion.
//
// Parameters: population estimates, Table 1 of the paper. Covariate effects
// on CMAX and the log-normal between-patient variability are applied by the
// app, as in the paper (Equations 4-6 of the Supplementary Material).
//
// Units: time in days; transgene level in DNA copies per ug genomic DNA
// (qPCR in peripheral blood).

$PARAM @annotated
FOLDX : 3900   : Fold expansion from t = 0 to TMAX
TMAX  : 9.3    : Time of maximal expansion (days)
CMAX  : 24000  : Maximal transgene level (copies/ug)
ALPHA : 0.16   : Rate of the fast contraction phase (1/day)
BETA  : 0.0032 : Rate of the persistence phase (1/day)
FB    : 0.0079 : Fraction of cells in the persistence phase
FTOCI : 1.2    : Expansion-rate factor after tocilizumab
FSTER : 1      : Expansion-rate factor after corticosteroids
TTOCI : 99999  : Day of the first tocilizumab dose
TSTER : 99999  : Day of the first corticosteroid dose

$CMT @annotated
LEFF : Log of the contracting (effector) CAR-T transgene (copies/ug)
LPER : Log of the persisting CAR-T transgene (copies/ug)

$MAIN
double R0 = CMAX / FOLDX;
double RHO = log(FOLDX) / TMAX;
LEFF_0 = log((1 - FB) * R0);
LPER_0 = log(FB * R0);

$ODE
double FT = SOLVERTIME >= TTOCI ? FTOCI : 1.0;
double FS = SOLVERTIME >= TSTER ? FSTER : 1.0;
double GROW = RHO * FT * FS;
dxdt_LEFF = SOLVERTIME < TMAX ? GROW : -ALPHA;
dxdt_LPER = SOLVERTIME < TMAX ? GROW : -BETA;

$TABLE
double EFF = exp(LEFF);
double PER = exp(LPER);
double CART = EFF + PER;

$CAPTURE CART EFF PER
