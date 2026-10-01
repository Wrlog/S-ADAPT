$PROB
One-compartment PK with zero-order infusion and indirect response model I

// Written independently of models/idr1_template/final_model.ctl, as the
// check on it. The model is the example of Bulitta et al., AAPS J
// 2011;13:201 (Fig. 4).

$PARAM @annotated
CL     : 5   : Clearance (L/h)
V1     : 30  : Volume of distribution (L)
TOUT12 : 10  : Half-life of response loss (h)
IMAX   : 0.7 : Maximum inhibition
BASE   : 15  : Baseline response
IC50   : 10  : Concentration giving half of IMAX (mg/L)
HILL   : 1   : Hill coefficient

$CMT @annotated
CENT : Drug in the body (mg)
RESP : Response

$MAIN
double KOUT = log(2.0) / TOUT12;
double KIN = BASE * KOUT;
RESP_0 = BASE;

$ODE
double C = CENT / V1;
double INH = C > 0 ? IMAX * pow(C, HILL) / (pow(C, HILL) + pow(IC50, HILL)) : 0.0;
dxdt_CENT = -CL * C;
dxdt_RESP = KIN * (1 - INH) - KOUT * RESP;

$TABLE
double CP = CENT / V1;

$CAPTURE CP
