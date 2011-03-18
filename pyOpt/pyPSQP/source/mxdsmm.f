* SUBROUTINE MXDSMM                ALL SYSTEMS                89/12/01
* PORTABILITY : ALL SYSTEMS
* 89/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* MULTIPLICATION OF A DENSE SYMMETRIC MATRIX A BY A VECTOR X.
*
* PARAMETERS :
*  II  N  ORDER OF THE MATRIX A.
*  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
*  RI  X(N)  INPUT VECTOR.
*  RO  Y(N)  OUTPUT VECTOR EQUAL TO  A*X.
*
      SUBROUTINE MXDSMM(N,A,X,Y)
      INTEGER          N
      DOUBLE PRECISION A(*),X(*),Y(*)
      DOUBLE PRECISION TEMP
      INTEGER          I,J,K,L
      K = 0
      DO 30 I = 1,N
          TEMP = 0.0D0
          L = K
          DO 10 J = 1,I
              L = L + 1
              TEMP = TEMP + A(L)*X(J)
   10     CONTINUE
          DO 20 J = I + 1,N
              L = L + J - 1
              TEMP = TEMP + A(L)*X(J)
   20     CONTINUE
          Y(I) = TEMP
          K = K + I
   30 CONTINUE
      RETURN
      END
