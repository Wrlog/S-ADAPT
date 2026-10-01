$PROJECT PROTAC targeted protein degradation in cells: the kcat model, six concentrations per experiment

$DIFFEQ_DIF

TE1 = AE0*X(1)/(X(1)*(AE0 + KDSUM + X(1)) + KDPROD)
TE2 = AE0*X(2)/(X(2)*(AE0 + KDSUM + X(2)) + KDPROD)
TE3 = AE0*X(3)/(X(3)*(AE0 + KDSUM + X(3)) + KDPROD)
TE4 = AE0*X(4)/(X(4)*(AE0 + KDSUM + X(4)) + KDPROD)
TE5 = AE0*X(5)/(X(5)*(AE0 + KDSUM + X(5)) + KDPROD)
TE6 = AE0*X(6)/(X(6)*(AE0 + KDSUM + X(6)) + KDPROD)

XP(1) = 0
XP(2) = 0
XP(3) = 0
XP(4) = 0
XP(5) = 0
XP(6) = 0
XP(7) = KDEGP - KDEGP*X(7) - KCAT*TE1*X(7)
XP(8) = KDEGP - KDEGP*X(8) - KCAT*TE2*X(8)
XP(9) = KDEGP - KDEGP*X(9) - KCAT*TE3*X(9)
XP(10) = KDEGP - KDEGP*X(10) - KCAT*TE4*X(10)
XP(11) = KDEGP - KDEGP*X(11) - KCAT*TE5*X(11)
XP(12) = KDEGP - KDEGP*X(12) - KCAT*TE6*X(12)


$OUTPUT_GLB

KDP = 71
KDE = 2500
COOP = 0.86
EZERO = 203
AE0 = COOP*EZERO
KDSUM = KDP + KDE
KDPROD = KDP*KDE
KDEGP = LOG(2)/THALFP


$OUTPUT_ICS

X(7) = X(7) + 1
X(8) = X(8) + 1
X(9) = X(9) + 1
X(10) = X(10) + 1
X(11) = X(11) + 1
X(12) = X(12) + 1


$OUTPUT_EQN

Y(1) = 100*X(7)
Y(2) = 100*X(8)
Y(3) = 100*X(9)
Y(4) = 100*X(10)
Y(5) = 100*X(11)
Y(6) = 100*X(12)


$VARMOD_EQN

V(1) = (SDIN + SDSL*Y(1))*(SDIN + SDSL*Y(1))
V(2) = (SDIN + SDSL*Y(2))*(SDIN + SDSL*Y(2))
V(3) = (SDIN + SDSL*Y(3))*(SDIN + SDSL*Y(3))
V(4) = (SDIN + SDSL*Y(4))*(SDIN + SDSL*Y(4))
V(5) = (SDIN + SDSL*Y(5))*(SDIN + SDSL*Y(5))
V(6) = (SDIN + SDSL*Y(6))*(SDIN + SDSL*Y(6))


$POPMOD_EQN
