* SUBROUTINE MXDPRB                ALL SYSTEMS                89/12/01
* PORTABILITY : ALL SYSTEMS
* 89/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* SOLUTION OF A SYSTEM OF LINEAR EQUATIONS WITH A DENSE SYMMETRIC
* POSITIVE DEFINITE MATRIX A USING THE FACTORIZATION A=TRANS(R)*R.
*
* PARAMETERS :
*  II  N ORDER OF THE MATRIX A.
*  RI  A(N*(N+1)/2) FACTORIZATION A=TRANS(R)*R.
*  RU  X(N)  ON INPUT THE RIGHT HAND SIDE OF A SYSTEM OF LINEAR
*         EQUATIONS. ON OUTPUT THE SOLUTION OF A SYSTEM OF LINEAR
*         EQUATIONS.
*  II  JOB  OPTION. IF JOB=0 THEN X:=A**(-1)*X. IF JOB>0 THEN
*         X:=TRANS(R)**(-1)*X. IF JOB<0 THEN X:=R**(-1)*X.
*
* METHOD :
* BACK SUBSTITUTION
*
      SUBROUTINE MXDPRB(N,A,X,JOB)
      INTEGER          JOB,N
      DOUBLE PRECISION A(*),X(*)
      INTEGER          I,II,IJ,J
      IF (JOB.GE.0) THEN
*
*     PHASE 1 : X:=TRANS(R)**(-1)*X
*
          IJ = 0
          DO 20 I = 1,N
              DO 10 J = 1,I - 1
                  IJ = IJ + 1
                  X(I) = X(I) - A(IJ)*X(J)
   10         CONTINUE
              IJ = IJ + 1
              X(I) = X(I)/A(IJ)
   20     CONTINUE
      END IF
      IF (JOB.LE.0) THEN
*
*     PHASE 2 : X:=R**(-1)*X
*
          II = N* (N+1)/2
          DO 40 I = N,1,-1
              IJ = II
              DO 30 J = I + 1,N
                  IJ = IJ + J - 1
                  X(I) = X(I) - A(IJ)*X(J)
   30         CONTINUE
              X(I) = X(I)/A(II)
              II = II - I
   40     CONTINUE
      END IF
      RETURN
      END
