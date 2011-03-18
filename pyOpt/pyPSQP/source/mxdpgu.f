* SUBROUTINE MXDPGU                ALL SYSTEMS                89/12/01
* PORTABILITY : ALL SYSTEMS
* 89/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* CORRECTION OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX A+E IN THE
* FACTORED FORM A+E=L*D*TRANS(L) OBTAINED BY THE SUBROUTINE MXDPGF.
* THE CORRECTION IS DEFINED AS A+E:=A+E+ALF*X*TRANS(X) WHERE ALF IS A
* GIVEN SCALING FACTOR AND X IS A GIVEN VECTOR.
*
* PARAMETERS :
*  II  N ORDER OF THE MATRIX A.
*  RU  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
*         SUBROUTINE MXDPGF.
*  RI  ALF  SCALING FACTOR IN THE CORRECTION TERM.
*  RI  X(N)  VECTOR IN THE CORRECTION TERM.
*  RA  Y(N) AUXILIARY VECTOR.
*
* METHOD :
* P.E.GILL, W.MURRAY, M.SAUNDERS: METHODS FOR COMPUTING AND MODIFYING
* THE LDV FACTORS OF A MATRIX, MATH. OF COMP. 29 (1974) PP. 1051-1077.
*
      SUBROUTINE MXDPGU(N,A,ALF,X,Y)
      DOUBLE PRECISION ZERO,ONE,FOUR,CON
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,FOUR=4.0D0,CON=1.0D-8)
      DOUBLE PRECISION ALF,ALFR
      INTEGER N
      DOUBLE PRECISION A(*),X(*),Y(*)
      DOUBLE PRECISION B,D,P,R,T,TO
      INTEGER I,II,IJ,J
      IF (ALF.GE.ZERO) THEN
*
*     FORWARD CORRECTION IN CASE WHEN THE SCALING FACTOR IS NONNEGATIVE
*
          ALFR = SQRT(ALF)
          CALL MXVSCL(N,ALFR,X,Y)
          TO = ONE
          II = 0
          DO 30 I = 1,N
              II = II + I
              D = A(II)
              P = Y(I)
              T = TO + P*P/D
              R = TO/T
              A(II) = D/R
              B = P/ (D*T)
              IF (A(II).LE.FOUR*D) THEN
*
*     AN EASY FORMULA FOR LIMITED DIAGONAL ELEMENT
*
                  IJ = II
                  DO 10 J = I + 1,N
                      IJ = IJ + J - 1
                      D = A(IJ)
                      Y(J) = Y(J) - P*D
                      A(IJ) = D + B*Y(J)
   10             CONTINUE
              ELSE
*
*     A MORE COMPLICATE BUT NUMERICALLY STABLE FORMULA FOR UNLIMITED
*     DIAGONAL ELEMENT
*
                  IJ = II
                  DO 20 J = I + 1,N
                      IJ = IJ + J - 1
                      D = A(IJ)
                      A(IJ) = R*D + B*Y(J)
                      Y(J) = Y(J) - P*D
   20             CONTINUE
              END IF
              TO = T
   30     CONTINUE
      ELSE
*
*     BACKWARD CORRECTION IN CASE WHEN THE SCALING FACTOR IS NEGATIVE
*
          ALFR = SQRT(-ALF)
          CALL MXVSCL(N,ALFR,X,Y)
          TO = ONE
          IJ = 0
          DO 50 I = 1,N
              D = Y(I)
              DO 40 J = 1,I - 1
                  IJ = IJ + 1
                  D = D - A(IJ)*Y(J)
   40         CONTINUE
              Y(I) = D
              IJ = IJ + 1
              TO = TO - D*D/A(IJ)
   50     CONTINUE
          IF (TO.LE.ZERO) TO = CON
          II = N* (N+1)/2
          DO 70 I = N,1,-1
              D = A(II)
              P = Y(I)
              T = TO + P*P/D
              A(II) = D*TO/T
              B = -P/ (D*TO)
              TO = T
              IJ = II
              DO 60 J = I + 1,N
                  IJ = IJ + J - 1
                  D = A(IJ)
                  A(IJ) = D + B*Y(J)
                  Y(J) = Y(J) + P*D
   60         CONTINUE
              II = II - I
   70     CONTINUE
      END IF
      RETURN
      END
