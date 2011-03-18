* SUBROUTINE PNINT3                ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* EXTRAPOLATION OR INTERPOLATION FOR LINE SEARCH WITHOUT DIRECTIONAL
* DERIVATIVES.
*
* PARAMETERS :
*  RI  RO  INITIAL VALUE OF THE STEPSIZE PARAMETER.
*  RI  RL  LOWER VALUE OF THE STEPSIZE PARAMETER.
*  RI  RU  UPPER VALUE OF THE STEPSIZE PARAMETER.
*  RI  RI  INNER VALUE OF THE STEPSIZE PARAMETER.
*  RI  FO  VALUE OF THE OBJECTIVE FUNCTION FOR R=RO.
*  RI  FL  VALUE OF THE OBJECTIVE FUNCTION FOR R=RL.
*  RI  FU  VALUE OF THE OBJECTIVE FUNCTION FOR R=RU.
*  RI  FI  VALUE OF THE OBJECTIVE FUNCTION FOR R=RI.
*  RO  PO  INITIAL VALUE OF THE DIRECTIONAL DERIVATIVE.
*  RO  R  VALUE OF THE STEPSIZE PARAMETER OBTAINED.
*  II  MODE  MODE OF LINE SEARCH.
*  II  MTYP  METHOD SELECTION. MTYP=1-BISECTION. MTYP=2-TWO POINT
*         QUADRATIC INTERPOLATION. MTYP=2-THREE POINT QUADRATIC
*         INTERPOLATION.
*  IO  MERR  ERROR INDICATOR. MERR=0 FOR NORMAL RETURN.
*
* METHOD :
* EXTRAPOLATION OR INTERPOLATION WITH STANDARD MODEL FUNCTIONS.
*
      SUBROUTINE PNINT3(RO,RL,RU,RI,FO,FL,FU,FI,PO,R,MODE,MTYP,MERR)
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,THREE,C1L,C1U,C2L,C2U,C3L
      PARAMETER        (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,
     +                 THREE=3.0D0,C1L=1.1D0,C1U=1.0D3,C2L=1.0D-2,
     +                 C2U=0.9D0,C3L=1.0D-1)
      DOUBLE PRECISION FI,FL,FO,FU,PO,R,RI,RL,RO,RU
      INTEGER          MERR,MODE,MTYP
      DOUBLE PRECISION AI,AL,AU,DEN,DIS
      INTEGER          NTYP
      LOGICAL          L1,L2
      MERR = 0
      IF (MODE.LE.0) RETURN
      IF (PO.GE.ZERO) THEN
          MERR = 2
          RETURN

      ELSE IF (RU.LE.RL) THEN
          MERR = 3
          RETURN
      END IF
      L1 = RL .LE. RO
      L2 = RI .LE. RL
      DO 10 NTYP = MTYP,1,-1
          IF (NTYP.EQ.1) THEN
*
*     BISECTION
*
              IF (MODE.EQ.1) THEN
                  R = TWO*RU
                  RETURN
              ELSE IF (RI-RL.LE.RU-RI) THEN
                  R = HALF* (RI+RU)
                  RETURN
              ELSE
                  R = HALF* (RL+RI)
                  RETURN
              END IF
          ELSE IF (NTYP.EQ.MTYP .AND. L1) THEN
              IF (.NOT.L2) AI = (FI-FO)/ (RI*PO)
              AU = (FU-FO)/ (RU*PO)
          END IF
          IF (L1 .AND. (NTYP.EQ.2.OR.L2)) THEN
*
*     TWO POINT QUADRATIC EXTRAPOLATION OR INTERPOLATION
*
              IF (AU.GE.ONE) GO TO 10
              R = HALF*RU/ (ONE-AU)
          ELSE IF (.NOT.L1 .OR. .NOT.L2 .AND. NTYP.EQ.3) THEN
*
*     THREE POINT QUADRATIC EXTRAPOLATION OR INTERPOLATION
*
              AL = (FI-FL)/ (RI-RL)
              AU = (FU-FI)/ (RU-RI)
              DEN = AU - AL
              IF (DEN.LE.ZERO) GO TO 10
              R = RI - HALF* (AU* (RI-RL)+AL* (RU-RI))/DEN
          ELSE IF (L1 .AND. .NOT.L2 .AND. NTYP.EQ.4) THEN
*
*     THREE POINT CUBIC EXTRAPOLATION OR INTERPOLATION
*
              DIS = (AI-ONE)* (RU/RI)
              DEN = (AU-ONE)* (RI/RU) - DIS
              DIS = AU + AI - DEN - TWO* (ONE+DIS)
              DIS = DEN*DEN - THREE*DIS
              IF (DIS.LT.ZERO) GO TO 10
              DEN = DEN + SQRT(DIS)
              IF (DEN.EQ.ZERO) GO TO 10
              R = (RU-RI)/DEN
          ELSE
              GO TO 10
          END IF
          IF (MODE.EQ.1 .AND. R.GT.RU) THEN
*
*     EXTRAPOLATION ACCEPTED
*
              R = MAX(R,C1L*RU)
              R = MIN(R,C1U*RU)
              RETURN
          ELSE IF (MODE.EQ.2 .AND. R.GT.RL .AND. R.LT.RU) THEN
*
*     INTERPOLATION ACCEPTED
*
              IF (RI.EQ.ZERO .AND. NTYP.NE.4) THEN
                  R = MAX(R,RL+C2L* (RU-RL))
              ELSE
                  R = MAX(R,RL+C3L* (RU-RL))
              END IF
              R = MIN(R,RL+C2U* (RU-RL))
              IF (R.EQ.RI) GO TO 10
              RETURN
          END IF
   10 CONTINUE
      END
