$PROB
PROTAC targeted protein degradation: the kcat model

// Mechanistic pharmacodynamic framework for proteolysis targeting chimeras:
//   A Mechanistic Pharmacodynamic Modeling Framework for the Assessment and
//   Optimization of Proteolysis Targeting Chimeras (PROTACs).
//   Pharmaceutics 2023;15:195 (open access; PMC9865105).
//
// A PROTAC bridges the target protein (P) and an E3 ubiquitin ligase (E).
// Only the ternary complex P-PROTAC-E leads to ubiquitination and proteasomal
// degradation, and the PROTAC is released to act again (catalytic). At high
// concentrations binary complexes (P-PROTAC, PROTAC-E) outcompete the ternary
// one: the "hook effect".
//
// Pharmacodynamics, as published:
//   target engagement TE = TC/P (rapid equilibrium, TC << E0), Equation (3)
//     TE = a E0 / (a E0 + KD,P + KD,E + KD,P KD,E / C + C)
//   target turnover with PROTAC-catalysed degradation, Equation (8)
//     d(P/P0)/dt = kdeg,P - kdeg,P (P/P0) - kcat TE (P/P0)
//   degradation D = 1 - P/P0 (7); occupancy-driven inhibition I = C/(KD,P + C)
//   (19); total target modulation TM = D + I - D I (6); downstream response
//   PD (20), normalised to 1 at baseline.
// Default parameters: BTK degrader Cpd. D of Zorba et al. (PNAS 2018) in
// Ramos cells (Tables S1-S2 of the paper) and kcat = 4.6 1/h, the value the
// authors fitted to the lead compound and used to predict the whole series.
//
// Pharmacokinetics: the paper is an in vitro framework (constant C). For an
// in vivo illustration this file adds a one-compartment oral PK model with
// ILLUSTRATIVE parameters (not from the paper) that drives the same PD
// through the unbound concentration. Set CFIX >= 0 for the in vitro setting
// (constant unbound concentration, nM).
//
// Units: time h; concentrations nM; dose mg.

$PARAM @annotated
KDP    : 71    : Binary dissociation constant, PROTAC-target (nM)
KDE    : 2500  : Binary dissociation constant, PROTAC-E3 ligase (nM)
COOP   : 0.86  : Cooperativity factor alpha (-)
E0     : 203   : Total E3 ligase concentration (nM)
THALFP : 16    : Target protein half-life (h)
KCAT   : 4.6   : Catalytic degradation rate constant (1/h)
INHIB  : 1     : Does binding also inhibit the target? (1 yes, 0 no)
PDMIN  : 0     : Target-independent residual response (-)
P50    : 0.5   : Target level giving half the response (-)
NPD    : 1     : Shape factor of the response (-)
CFIX   : -1    : In vitro, constant unbound concentration (nM); < 0 uses PK
KA     : 0.8   : Illustrative PK, absorption rate constant (1/h)
CL     : 12    : Illustrative PK, apparent clearance CL/F (L/h)
VD     : 150   : Illustrative PK, apparent volume V/F (L)
FU     : 0.05  : Illustrative PK, fraction unbound in plasma (-)
MW     : 950   : Molecular weight (g/mol)

$INIT @annotated
GUT  : 0 : PROTAC in the gut (mg)
CENT : 0 : PROTAC in the body (mg)
PR   : 1 : Target protein relative to baseline, P/P0 (-)

$ODE
double CU = CFIX >= 0 ? CFIX : FU * (CENT / VD) / MW * 1e6;
double TE = COOP * E0 * CU / (CU * (COOP * E0 + KDP + KDE + CU) + KDP * KDE);
double KDEG = log(2.0) / THALFP;
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - CL / VD * CENT;
dxdt_PR = KDEG - KDEG * PR - KCAT * TE * PR;

$TABLE
double CUNB = CFIX >= 0 ? CFIX : FU * (CENT / VD) / MW * 1e6;
double CPLASMA = (CENT / VD) / MW * 1e6;
double TENG = COOP * E0 * CUNB / (CUNB * (COOP * E0 + KDP + KDE + CUNB) + KDP * KDE);
double DEG = 1 - PR;
double INH = INHIB * CUNB / (KDP + CUNB);
double TM = DEG + INH - DEG * INH;
double RES = pow(1 - TM, NPD);
double PDR = PDMIN + (1 - PDMIN) * RES * (1 - pow(P50, NPD)) / (pow(P50, NPD) + RES * (1 - 2 * pow(P50, NPD)));

$CAPTURE CUNB CPLASMA TENG DEG INH TM PDR
