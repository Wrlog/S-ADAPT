$PROB
Methylprednisolone: cortisol suppression, lymphocyte trafficking and ex vivo
whole blood lymphocyte proliferation

// Hong Y, Mager DE, Blum RA, Jusko WJ. Population pharmacokinetic/
//   pharmacodynamic modeling of systemic corticosteroid inhibition of whole
//   blood lymphocytes: modeling interoccasion pharmacodynamic variability.
//   Pharm Res 2007;24:1088-1097 (Eqs. 1-3 and 6-8; Tables I-IV,
//   methylprednisolone).
//
// Written independently of models/corticosteroid_hong2007/final_model.ctl,
// as the check on it.
//
// The paper describes each subject's cortisol baseline with a three-harmonic
// Fourier series whose coefficients it does not report. Here the baseline is
// one harmonic, with ILLUSTRATIVE values: mean 80 ng/mL, amplitude 55 ng/mL,
// peak one hour before dosing.
//
// Units: time h; methylprednisolone in ug and ng/mL; cortisol ng/mL;
// lymphocytes cells/uL; proliferation in percent of the pre-dose value.

$PARAM @annotated
CLMPL  : 22.8   : Methylprednisolone clearance (L/h)
VMPL   : 78.4   : Methylprednisolone volume (L)
KCOUT  : 0.300  : Cortisol elimination rate (1/h)
IC50MC : 1.59   : Methylprednisolone IC50, cortisol secretion (ng/mL)
KINL   : 1125   : Lymphocyte zero-order return rate (cells/uL/h)
KBE    : 0.283  : Lymphocyte trafficking rate (1/h)
IC50CL : 91.1   : Cortisol IC50, lymphocyte trafficking (ng/mL)
IC50ML : 19.6   : Methylprednisolone IC50, lymphocyte trafficking (ng/mL)
KT     : 0.217  : Transit rate to the active concentration (1/h)
IC50MW : 0.984  : Methylprednisolone IC50, ex vivo proliferation (ng/mL)
FA0    : 80     : Cortisol baseline, mean (ng/mL; illustrative)
FA1    : 53.13  : Cortisol baseline, cosine coefficient (ng/mL; illustrative)
FB1    : -14.23 : Cortisol baseline, sine coefficient (ng/mL; illustrative)
DILF   : 0.05   : Dilution of blood in the ex vivo assay

$CMT @annotated
MPL  : Methylprednisolone in the body (ug)
CORT : Endogenous cortisol (ng/mL)
LYM  : Blood lymphocytes (cells/uL)
CA   : Active methylprednisolone concentration (ng/mL)

$MAIN
double W = 2.0 * 3.14159265358979 / 24.0;
double CEN0 = FA0 + FA1;
double LYM0 = KINL / KBE * (1.0 - CEN0 / (IC50CL + CEN0));
CORT_0 = CEN0;
LYM_0 = LYM0;

$ODE
double CP = MPL > 0 ? MPL / VMPL : 0.0;
double KIN = KCOUT * FA0 + (KCOUT * FA1 + FB1 * W) * cos(W * SOLVERTIME) + (KCOUT * FB1 - FA1 * W) * sin(W * SOLVERTIME);
double JOINT = CP * (IC50CL / IC50ML) + CORT;
dxdt_MPL = -CLMPL * CP;
dxdt_CORT = KIN * (1.0 - CP / (IC50MC + CP)) - KCOUT * CORT;
dxdt_LYM = KINL * (1.0 - JOINT / (IC50CL + JOINT)) - KBE * LYM;
dxdt_CA = KT * (CP - CA);

$TABLE
double CMPL = MPL / VMPL;
double WBLP = 100.0 * LYM / LYM0 * (1.0 - CA * DILF / (IC50MW + CA * DILF));

$CAPTURE CMPL WBLP
