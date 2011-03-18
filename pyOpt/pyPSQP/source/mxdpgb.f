* SUBROUTINE MXDPGB                ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* SOLUTION OF A SYSTEM OF LINEAR EQUATIONS WITH A DENSE SYMMETRIC
* POSITIVE DEFINITE MATRIX A+E USING THE FACTORIZATION A+E=L*D*TRANS(L)
* OBTAINED BY THE SUBROUTINE MXDPGF.
*
* PARAMETERS :
*  II  N ORDER OF THE MATRIX A.
*  RI  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
*         SUBROUTINE MXDPGF.
*  RU  X(N)  ON INPUT THE RIGHT HAND SIDE OF A SYSTEM OF LINEAR
*         EQUATIONS. ON OUTPUT THE SOLUTION OF A SYSTEM OF LINEAR
*         EQUATIONS.
*  II  JOB  OPTION. IF JOB=0 THEN X:=(A+E)**(-1)*X. IF JOB>0 THEN
*         X:=L**(-1)*X. IF JOB<0 THEN X:=TRANS(L)**(-1)*X.
*
* METHOD :
* BACK SUBSTITUTION
*
      SUBROUTINE MXDPGB(N,A,X,JOB)
      INTEGER JOB,N
      DOUBLE PRECISION A(*),X(*)
      INTEGER I,II,IJ,J
      IF (JOB.GE.0) THEN
*
*     PHASE 1 : X:=L**(-1)*X
*
          IJ = 0
          DO 20 I = 1,N
              DO 10 J = 1,I - 1
                  IJ = IJ + 1
                  X(I) = X(I) - A(IJ)*X(J)
   10         CONTINUE
              IJ = IJ + 1
   20     CONTINUE
      END IF
      IF (JOB.EQ.0) THEN
*
*     PHASE 2 : X:=D**(-1)*X
*
          II = 0
          DO 30 I = 1,N
              II = II + I
              X(I) = X(I)/A(II)
   30     CONTINUE
      END IF
      IF (JOB.LE.0) THEN
*
*     PHASE 3 : X:=TRANS(L)**(-1)*X
*
          II = N* (N-1)/2
          DO 50 I = N - 1,1,-1
              IJ = II
              DO 40 J = I + 1,N
                  IJ = IJ + J - 1
                  X(I) = X(I) - A(IJ)*X(J)
   40         CONTINUE
              II = II - I
   50     CONTINUE
      END IF
      RETURN
      END
