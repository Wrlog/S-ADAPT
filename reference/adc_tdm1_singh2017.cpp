$PROB
Antibody-drug conjugate T-DM1: systemic PK, tumour disposition, payload
release and tumour growth inhibition (mouse and patient)

// Singh AP, Shah DK. Application of a PK-PD modeling and simulation-based
//   strategy for clinical translation of antibody-drug conjugates: a case
//   study with trastuzumab emtansine (T-DM1). AAPS J 2017;19:1054-1070
//   (Eqs. 1-21, Tables I and III).
// Tumour disposition and cellular parameters: Singh AP, Maass KF, Betts AM,
//   et al. Evolution of antibody-drug conjugate tumor disposition model to
//   predict preclinical tumor pharmacokinetics of trastuzumab-emtansine
//   (T-DM1). AAPS J 2016;18:861-875.
//
// Plasma: total trastuzumab (TT, two compartments) and conjugated T-DM1 (two
// compartments, with nonspecific deconjugation Kdec from the central
// compartment); the released DM1 catabolites have their own two-compartment
// model and receive DAR payloads per deconjugated or cleared ADC. The average
// DAR falls with Kdec.
// Tumour (a sphere of volume TV): ADC and DM1 enter across the vessels
// (permeability, Krogh cylinder) and the tumour surface (diffusion), free
// ADC binds HER2, the complex is internalised and degraded in lysosomes,
// releasing DAR DM1 inside the cell, where DM1 binds tubulin or diffuses
// out. Intracellular DM1 (free + tubulin-bound) kills tumour cells with a
// first-order constant KKILL.
// Growth: exponential in mice (GLIN = 0); in patients the Haddish-Berhane
// form, exponential when small, linear when large, capped at VMAX (Eq. 19).
//
// Transcription notes. Four printed equations disagree with the model they
// describe and are used here in their mass-balanced form: Eq. 3's
// deconjugation acts on the central compartment (X1, not X2), as Eq. 16's
// DM1 input requires; Eqs. 4 and 9 use CLD/V2 (as Eq. 2 does), not CL/V2;
// Eq. 8 keeps the deconjugation term of Eq. 3; and the tumour exchange of
// DM1 in Eq. 16 is scaled by TV/V1 like the ADC exchange in Eq. 8. In Eq. 13
// the efflux term is KOUT (zero in Table I). TV is in litres: that is the
// unit in which Eqs. 8 and 16 balance, and in which the clinical growth
// parameters of Table III give KGLIN = ln2/DTLIN in L/day. Rate constants of
// Table I given per hour are converted to per day.
//
// Units: time day; antibody amounts nmol, concentrations nM; volumes L
// (tumour TV too); doses: mg/kg x BW x 1e6 / MW nmol into X1TT and X1ADC.
// Defaults: patient, HER2 3+ (1660 nM antigen), translated human PK.

$PARAM @annotated
BW      : 70       : Body weight (kg)
MW      : 148500   : T-DM1 molecular weight (g/mol)
DAR0    : 3.5      : Initial average drug:antibody ratio
CLADC   : 0.0043   : ADC (and total antibody) clearance (L/day/kg)
CLDADC  : 0.014    : ADC distributional clearance (L/day/kg)
V1ADC   : 0.034    : ADC central volume (L/kg)
V2ADC   : 0.04     : ADC peripheral volume (L/kg)
KDEC    : 0.241    : Systemic nonspecific deconjugation (1/day)
CLDRUG  : 2.23     : DM1 catabolite clearance (L/day/kg)
CLDDRUG : 1.0      : DM1 distributional clearance (L/day/kg)
V1DRUG  : 0.034    : DM1 central volume (L/kg)
V2DRUG  : 5.0      : DM1 peripheral volume (L/kg)
KONADC  : 8.88     : ADC-HER2 association (1/nM/day; 0.37 1/nM/h)
KOFFADC : 2.328    : ADC-HER2 dissociation (1/day; 0.097 1/h)
KINT    : 2.16     : Internalisation of HER2-ADC (1/day; 0.09 1/h)
KDEG    : 0.72     : Lysosomal degradation of ADC (1/day; 0.03 1/h)
KDECT   : 0.5352   : Deconjugation in tumour interstitium (1/day; 0.0223 1/h)
KONTUB  : 0.72     : DM1-tubulin association (1/nM/day; 0.03 1/nM/h)
KOFFTUB : 0.684    : DM1-tubulin dissociation (1/day; 0.0285 1/h)
TUB     : 65       : Intracellular tubulin (nM)
KDIFF   : 2.208    : DM1 diffusion across the cell membrane (1/day; 0.092 1/h)
KOUT    : 0        : Active DM1 efflux (1/day)
AG      : 1660     : Total HER2 in tumour (nM): 3+ 1660, 2+ 830, 1+ 166
RCAP    : 8.0      : Tumour capillary radius (um)
RKROGH  : 75.0     : Half-distance between capillaries (um)
PADC    : 334      : ADC vascular permeability (um/day)
PDRUG   : 21000    : DM1 vascular permeability (um/day)
DADC    : 0.022    : ADC diffusivity in tumour (cm2/day)
DDRUG   : 0.25     : DM1 diffusivity in tumour (cm2/day)
EPSADC  : 0.24     : Tumour void fraction for ADC
EPSDRUG : 0.44     : Tumour void fraction for DM1
TV0     : 0.002432 : Initial tumour volume (L)
GLIN    : 1        : 1 = exponential-linear-plateau growth (patients), 0 = exponential (mice)
DTEXP   : 25       : Doubling time, exponential growth (day)
DTLIN   : 621      : Linear growth: KGLIN = ln2/DTLIN (L/day)
PSI     : 20       : Switch between exponential and linear growth
VMAX    : 0.5238   : Maximum tumour volume (L)
KKILL   : 1.8e-4   : Killing per nM of intracellular DM1 (1/day/nM)

$CMT @annotated
X1TT    : Total trastuzumab, central (nmol)
X2TT    : Total trastuzumab, peripheral (nmol)
X1ADC   : Conjugated T-DM1, central (nmol)
X2ADC   : Conjugated T-DM1, peripheral (nmol)
C1DRUG  : DM1 catabolites, central (nM)
C2DRUG  : DM1 catabolites, peripheral (nM)
DAR     : Average drug:antibody ratio
ADCF    : Free ADC in tumour (nM of tumour)
ADCB    : HER2-bound ADC in tumour (nM)
ADCE    : ADC in endosomes/lysosomes (nM)
DRUGFC  : Free DM1 in tumour cells (nM)
DRUGBC  : Tubulin-bound DM1 in tumour cells (nM)
DRUGFX  : DM1 in tumour interstitium (nM)
TV      : Tumour volume (L)

$MAIN
double V1A = V1ADC * BW;
double V2A = V2ADC * BW;
double CLA = CLADC * BW;
double CLDA = CLDADC * BW;
double V1D = V1DRUG * BW;
double V2D = V2DRUG * BW;
double CLD1 = CLDRUG * BW;
double CLDD = CLDDRUG * BW;
double KGEX = log(2.0) / DTEXP;
double KGLIN = log(2.0) / DTLIN;
double VASA = 2.0 * PADC * RCAP / (RKROGH * RKROGH);        // 1/day
double VASD = 2.0 * PDRUG * RCAP / (RKROGH * RKROGH);
DAR_0 = DAR0;
TV_0 = TV0;

$ODE
// tumour radius (cm) from its volume (L): Eq. 21 with TV in mm3
double TVMM3 = fmax(TV, 1e-12) * 1e6;
double RTUM = pow(3.0 * TVMM3 / (4.0 * 3.14159265358979), 1.0 / 3.0) / 10.0;
double SURA = 6.0 * DADC / (RTUM * RTUM);
double SURD = 6.0 * DDRUG / (RTUM * RTUM);
double CPADC = X1ADC / V1A;
double XADC = (VASA + SURA) * (CPADC * EPSADC - ADCF);       // per tumour volume
double XDRUG = (VASD + SURD) * (C1DRUG * EPSDRUG - DRUGFX);
double GROWTH = GLIN > 0.5 ? KGEX * (1.0 - TV / VMAX) / pow(1.0 + pow(KGEX / KGLIN * TV, PSI), 1.0 / PSI) : KGEX;

dxdt_X1TT = -CLA / V1A * X1TT - CLDA / V1A * X1TT + CLDA / V2A * X2TT;
dxdt_X2TT = CLDA / V1A * X1TT - CLDA / V2A * X2TT;
dxdt_X1ADC = -CLA / V1A * X1ADC - CLDA / V1A * X1ADC + CLDA / V2A * X2ADC - KDEC * X1ADC - XADC * TV;
dxdt_X2ADC = CLDA / V1A * X1ADC - CLDA / V2A * X2ADC;
dxdt_C1DRUG = -CLD1 / V1D * C1DRUG - CLDD / V1D * C1DRUG + CLDD / V1D * C2DRUG + X1ADC * DAR * KDEC / V1D
              + CLA * DAR * X1ADC / (V1A * V1D) - XDRUG * TV / V1D;
dxdt_C2DRUG = CLDD / V2D * C1DRUG - CLDD / V2D * C2DRUG;
dxdt_DAR = -KDEC * DAR;
dxdt_ADCF = XADC - KONADC * ADCF * (AG - ADCB) / EPSADC + KOFFADC * ADCB - KDECT * ADCF;
dxdt_ADCB = KONADC * ADCF * (AG - ADCB) / EPSADC - (KOFFADC + KINT + KDECT) * ADCB;
dxdt_ADCE = KINT * ADCB - KDEG * ADCE;
dxdt_DRUGFC = KDEG * DAR * ADCE - KONTUB * DRUGFC * (TUB - DRUGBC) + KOFFTUB * DRUGBC - KOUT * DRUGFC
              + KDIFF * (DRUGFX - DRUGFC);
dxdt_DRUGBC = KONTUB * DRUGFC * (TUB - DRUGBC) - KOFFTUB * DRUGBC;
dxdt_DRUGFX = XDRUG + KOUT * DRUGFC + KDECT * DAR * (ADCF + ADCB) - KDIFF * (DRUGFX - DRUGFC);
dxdt_TV = (GROWTH - KKILL * (DRUGFC + DRUGBC)) * TV;

$TABLE
double TT_UGML = X1TT / V1A * MW / 1e6;                     // total trastuzumab (ug/mL)
double ADC_UGML = X1ADC / V1A * MW / 1e6;                   // T-DM1 (ug/mL)
double DM1_NGML = C1DRUG * 737.5 / 1e3;                     // DM1 catabolites (ng/mL; DM1 737.5 g/mol)
double DM1_TUMOUR = DRUGFC + DRUGBC;                        // intracellular DM1 (nM)
double TV_MM3 = fmax(TV, 0.0) * 1e6;
double DIAM_MM = 2.0 * pow(3.0 * TV_MM3 / (4.0 * 3.14159265358979), 1.0 / 3.0);
double RO_HER2 = AG > 0 ? ADCB / AG : 0.0;

$CAPTURE TT_UGML ADC_UGML DM1_NGML DM1_TUMOUR TV_MM3 DIAM_MM RO_HER2
