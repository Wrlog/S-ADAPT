$PROJECT One-compartment PK with zero-order infusion and indirect response model I

$DIFFEQ_DIF

DC1 = X(1)/V1

IF (DC1.GT.0) THEN
   INH  = Imax*DC1**HILL / (DC1**HILL + IC50**HILL)
ELSE
   INH  = 0
ENDIF

XP(1) = R(1) -CL*DC1
XP(2) = Kin * (1 - INH) - Kout * X(2)


$OUTPUT_GLB

Kout = LOG(2)/Tout12
Kin  = BASE * KOUT


$OUTPUT_ICS

X(2) = X(2) + BASE


$OUTPUT_EQN

if (X(1).LT.0) X(1) = 0
if (X(2).LT.0) X(2) = 0

Y(1) = X(1)/V1
Y(2) = X(2)


$VARMOD_EQN

V(1) = ( SDin + SDsl*Y(1) ) * ( SDin + SDsl*Y(1) )
V(2) = ( PDin + PDsl*Y(2) ) * ( PDin + PDsl*Y(2) )


$POPMOD_EQN
