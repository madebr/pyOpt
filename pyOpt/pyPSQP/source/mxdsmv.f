* SUBROUTINE MXDSMV                ALL SYSTEMS                91/12/01
* PORTABILITY : ALL SYSTEMS
* 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* K-TH ROW OF A DENSE SYMMETRIC MATRIX A IS COPIED TO THE VECTOR X.
*
* PARAMETERS :
*  II  N  ORDER OF THE MATRIX A.
*  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
*  RO  X(N)  OUTPUT VECTOR.
*  II  K  INDEX OF COPIED ROW.
*
      SUBROUTINE MXDSMV(N,A,X,K)
      INTEGER          K,N
      DOUBLE PRECISION A(*),X(*)
      INTEGER          I,L
      L = K* (K-1)/2
      DO 10 I = 1,N
          IF (I.LE.K) THEN
              L = L + 1
          ELSE
              L = L + I - 1
          END IF
          X(I) = A(L)
   10 CONTINUE
      RETURN
      END
