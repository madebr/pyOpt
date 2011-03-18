* SUBROUTINE PLMINS             ALL SYSTEMS                   91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DETERMINATION OF THE NEW ACTIVE SIMPLE BOUND.
*
* PARAMETERS :
*  II  NF DECLARED NUMBER OF VARIABLES.
*  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RI  XO(NF)  SAVED VECTOR OF VARIABLES.
*  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  RI  S(NF)  DIRECTION VECTOR.
*  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
*         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
*  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
*  IO  KNEW  SIGNUM OF THE NEW NORMAL.
*  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
*  RA  PAR  AUXILIARY VARIABLE.
*
      SUBROUTINE PLMINS(NF,IX,XO,XL,XU,S,KBF,INEW,KNEW,EPS9,PAR)
      DOUBLE PRECISION EPS9,PAR
      INTEGER          INEW,KBF,KNEW,NF
      DOUBLE PRECISION S(*),XL(*),XO(*),XU(*)
      INTEGER          IX(*)
      DOUBLE PRECISION POM,TEMP
      INTEGER          I
      IF (KBF.GT.0) THEN
          DO 10 I = 1,NF
              IF (IX(I).GT.0) THEN
                  TEMP = 1.0D0
                  IF (IX(I).EQ.1 .OR. IX(I).GE.3) THEN
                      POM = XO(I) + S(I)*TEMP - XL(I)
                      IF (POM.LT.MIN(PAR,-EPS9*MAX(ABS(XL(I)),
     +                    TEMP))) THEN
                          INEW = -I
                          KNEW = 1
                          PAR = POM
                      END IF
                  END IF
                  IF (IX(I).EQ.2 .OR. IX(I).GE.3) THEN
                      POM = XU(I) - S(I)*TEMP - XO(I)
                      IF (POM.LT.MIN(PAR,-EPS9*MAX(ABS(XU(I)),
     +                    TEMP))) THEN
                          INEW = -I
                          KNEW = -1
                          PAR = POM
                      END IF
                  END IF
              END IF
   10     CONTINUE
      END IF
      RETURN
      END
