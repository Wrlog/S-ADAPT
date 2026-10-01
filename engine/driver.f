C     ==================================================================
C     Reference driver for a model written as a SADAPT-TRAN control
C     stream. engine/ctl.R turns the control stream into the routines
C     SAPAR, SAGLB, SAICS, DIFFEQ, OUTPUT and VARMOD; this file
C     integrates them over one subject's records with LSODA, the solver
C     ADAPT and S-ADAPT use, for a batch of parameter vectors.
C
C     Records follow the NONMEM layout S-ADAPT reads:
C       EVID 1, RATE = 0   bolus of AMT into compartment CMT
C       EVID 1, RATE > 0   zero-order input R(CMT) for AMT/RATE
C       EVID 0             observation of output Y(CMT)
C     Records must be in time order. Integration starts at time zero,
C     or at the first record if that is earlier.
C     ==================================================================
      SUBROUTINE SASIM(NPAR, NSAMP, P, NEQ, NREC, TIME, EVID, CMT,
     &                 AMT, RATE, RTOL, ATOL, YPRED, VPRED, IERR)
      IMPLICIT NONE
      INTEGER NPAR, NSAMP, NEQ, NREC
      INTEGER EVID(NREC), CMT(NREC), IERR(NSAMP)
      DOUBLE PRECISION P(NPAR,NSAMP), TIME(NREC), AMT(NREC)
      DOUBLE PRECISION RATE(NREC), RTOL, ATOL
      DOUBLE PRECISION YPRED(NREC,NSAMP), VPRED(NREC,NSAMP)
      INTEGER MAXEQ, MAXOUT, MAXINF, LRW, LIW
      PARAMETER (MAXEQ=60, MAXOUT=30, MAXINF=200)
      PARAMETER (LRW=22+9*MAXEQ+MAXEQ*MAXEQ, LIW=20+MAXEQ)
      DOUBLE PRECISION X(MAXEQ), XC(MAXEQ), Y(MAXOUT), V(MAXOUT)
      DOUBLE PRECISION RWORK(LRW), T, TNEXT, TE
      DOUBLE PRECISION TEND(MAXINF), RINF(MAXINF)
      INTEGER IWORK(LIW), ICINF(MAXINF)
      INTEGER I, J, K, M, NINF, ISTATE, IFAIL
      DOUBLE PRECISION R
      COMMON /SAINP/ R(50)

      CALL XSETF(0)
      DO 900 K = 1, NSAMP
        IERR(K) = 0
        IFAIL = 0
        CALL SAPAR(P(1,K))
        CALL SAGLB
        DO 10 I = 1, MAXEQ
          X(I) = 0.0D0
   10   CONTINUE
        DO 20 I = 1, 50
          R(I) = 0.0D0
   20   CONTINUE
        CALL SAICS(X)
        T = DMIN1(0.0D0, TIME(1))
        NINF = 0
        ISTATE = 1
        DO 800 J = 1, NREC
          YPRED(J,K) = 0.0D0
          VPRED(J,K) = 0.0D0
          IF (IFAIL.NE.0) GOTO 800
          TNEXT = TIME(J)
C         stop at the end of every infusion that finishes on the way
  100     CONTINUE
          M = 0
          TE = TNEXT
          DO 110 I = 1, NINF
            IF (TEND(I).LE.TE) THEN
              M = I
              TE = TEND(I)
            ENDIF
  110     CONTINUE
          IF (M.GT.0) THEN
            IF (TE.GT.T) CALL SASTEP(NEQ, X, T, TE, RTOL, ATOL,
     &          ISTATE, RWORK, LRW, IWORK, LIW, IFAIL)
            IF (IFAIL.NE.0) GOTO 800
            R(ICINF(M)) = R(ICINF(M)) - RINF(M)
            TEND(M) = TEND(NINF)
            RINF(M) = RINF(NINF)
            ICINF(M) = ICINF(NINF)
            NINF = NINF - 1
            ISTATE = 1
            GOTO 100
          ENDIF
          IF (TNEXT.GT.T) CALL SASTEP(NEQ, X, T, TNEXT, RTOL, ATOL,
     &        ISTATE, RWORK, LRW, IWORK, LIW, IFAIL)
          IF (IFAIL.NE.0) GOTO 800
          IF (EVID(J).EQ.1) THEN
            IF (RATE(J).GT.0.0D0) THEN
              IF (NINF.GE.MAXINF) THEN
                IFAIL = 2
                GOTO 800
              ENDIF
              NINF = NINF + 1
              R(CMT(J)) = R(CMT(J)) + RATE(J)
              TEND(NINF) = T + AMT(J)/RATE(J)
              RINF(NINF) = RATE(J)
              ICINF(NINF) = CMT(J)
            ELSE
              X(CMT(J)) = X(CMT(J)) + AMT(J)
            ENDIF
            ISTATE = 1
          ELSE
C           the output code may clamp the states it is given, so it
C           works on a copy and the integration carries on undisturbed
            DO 200 I = 1, NEQ
              XC(I) = X(I)
  200       CONTINUE
            CALL OUTPUT(Y, T, XC)
            CALL VARMOD(V, T, XC, Y)
            YPRED(J,K) = Y(CMT(J))
            VPRED(J,K) = V(CMT(J))
          ENDIF
  800   CONTINUE
        IERR(K) = IFAIL
  900 CONTINUE
      RETURN
      END

C     One LSODA call from T to TOUT. IFAIL = 1 if the solver gave up.
      SUBROUTINE SASTEP(NEQ, X, T, TOUT, RTOL, ATOL, ISTATE, RWORK,
     &                  LRW, IWORK, LIW, IFAIL)
      IMPLICIT NONE
      INTEGER NEQ, ISTATE, LRW, LIW, IWORK(LIW), IFAIL, I
      DOUBLE PRECISION X(*), T, TOUT, RTOL, ATOL, RWORK(LRW)
      EXTERNAL SAFEX, SAJAC
      IF (ISTATE.EQ.1) THEN
        DO 10 I = 5, 10
          RWORK(I) = 0.0D0
          IWORK(I) = 0
   10   CONTINUE
        IWORK(6) = 100000
      ENDIF
      CALL DLSODA(SAFEX, NEQ, X, T, TOUT, 1, RTOL, ATOL, 1, ISTATE,
     &            1, RWORK, LRW, IWORK, LIW, SAJAC, 2)
      IF (ISTATE.LT.0) IFAIL = 1
      RETURN
      END

      SUBROUTINE SAFEX(NEQ, T, X, XP)
      IMPLICIT NONE
      INTEGER NEQ
      DOUBLE PRECISION T, X(*), XP(*)
      CALL DIFFEQ(T, X, XP)
      RETURN
      END

      SUBROUTINE SAJAC(NEQ, T, X, ML, MU, PD, NROWPD)
      IMPLICIT NONE
      INTEGER NEQ, ML, MU, NROWPD
      DOUBLE PRECISION T, X(*), PD(NROWPD,*)
      RETURN
      END

C     Residual variances for predictions already in hand, so that the
C     variance parameters can be updated without solving the model
C     again. Each output's variance may depend on that output only.
      SUBROUTINE SAVARB(NPAR, P, N, CMT, YP, VP)
      IMPLICIT NONE
      INTEGER NPAR, N, CMT(N), I, J
      DOUBLE PRECISION P(NPAR), YP(N), VP(N)
      DOUBLE PRECISION X(60), Y(30), V(30)
      CALL SAPAR(P)
      CALL SAGLB
      DO 10 J = 1, 60
        X(J) = 0.0D0
   10 CONTINUE
      DO 30 I = 1, N
        DO 20 J = 1, 30
          Y(J) = 0.0D0
   20   CONTINUE
        Y(CMT(I)) = YP(I)
        CALL VARMOD(V, 0.0D0, X, Y)
        VP(I) = V(CMT(I))
   30 CONTINUE
      RETURN
      END

C     SADAPT-TRAN limits the argument of EXP to 40 and of LOG to 1E-15
C     in the code it generates; the translated model calls these.
      DOUBLE PRECISION FUNCTION SAEXP(A)
      IMPLICIT NONE
      DOUBLE PRECISION A
      SAEXP = DEXP(DMIN1(A, 40.0D0))
      RETURN
      END

      DOUBLE PRECISION FUNCTION SALOG(A)
      IMPLICIT NONE
      DOUBLE PRECISION A
      SALOG = DLOG(DMAX1(A, 1.0D-15))
      RETURN
      END
