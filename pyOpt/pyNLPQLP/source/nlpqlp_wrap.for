      SUBROUTINE NLPQLP_WRAP (L,M,ME,MMAX,N,NMAX,MNN2,
     /                 X,F,G,DF,DG,U,XL,XU,C,D,
     /                 ACC,ACCQP,STPMIN,MAXFUN,MAXIT,MAXNM,RHOB,MODE,
     /                 IFAIL,IPRINT,IOUT,IFILE,
     /                 WA,LWA,KWA,LKWA,ACTIVE,LACTIV,LQL,
     /                 NFUN,NGRD,NLFUNC,NLGRAD)
C
      IMPLICIT         NONE
      INTEGER          NMAX, MMAX, MNN2, LWA, LKWA, LACTIV,
     /                 KWA(LKWA), N, ME, M, L, MAXIT, MAXFUN, 
     /                 IPRINT,MAXNM,IOUT,MODE,IFAIL,NFUN,NGRD
      DOUBLE PRECISION X(NMAX,L), F(L), G(MMAX,L), DF(NMAX),
     /                 DG(MMAX,NMAX), U(MNN2), XL(NMAX), XU(NMAX),
     /                 C(NMAX,NMAX), D(NMAX), WA(LWA), ACC, ACCQP, 
     /                 STPMIN, RHOB
      LOGICAL          ACTIVE(LACTIV), LQL
      CHARACTER*(*)    IFILE
      EXTERNAL         QL
      EXTERNAL         NLFUNC,NLGRAD
C
C  OPEN WRITE FILE
C
      IF (IPRINT.NE.0) THEN
        OPEN(UNIT=IOUT,FILE=IFILE,STATUS='UNKNOWN')
      ENDIF
C
C  EVALUATE OBJECTIVE AND CONSTRAINTS
C
   30 CONTINUE
      CALL NLFUNC(L,NMAX,MMAX,X,LACTIV,ACTIVE,F,G)
      NFUN=NFUN+1
      IF (IFAIL.EQ.-1) GOTO 50
C
C  EVALUATE GRADIENTS OF OBJECTIVE AND CONSTRAINTS
C
   40 CONTINUE
      CALL NLGRAD(L,NMAX,MMAX,X,LACTIV,ACTIVE,F,G,DF,DG)
      NGRD=NGRD+1
C
C  CALL NLPQLP
C
   50 CONTINUE
      CALL NLPQLP (     L,      M,     ME,   MMAX,      N, 
     /               NMAX,   MNN2,      X,      F,      G,  
     /                 DF,     DG,      U,     XL,     XU,
     /                  C,      D,    ACC,  ACCQP, STPMIN, 
     /             MAXFUN,  MAXIT,  MAXNM,   RHOB, IPRINT,
     /               MODE,   IOUT,  IFAIL,     WA,    LWA,
     /                KWA,   LKWA, ACTIVE, LACTIV,    LQL, 
     /                 QL)
C
C  EVALUATIONS
C
	  IF (IFAIL.EQ.-1) GOTO 30
	  IF (IFAIL.EQ.-2) GOTO 40
C
      RETURN
      END
