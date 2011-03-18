* SUBROUTINE MXDSMI                ALL SYSTEMS                88/12/01
* PORTABILITY : ALL SYSTEMS
* 88/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
* DENSE SYMMETRIC MATRIX A IS SET TO THE UNIT MATRIX WITH THE SAME
* ORDER.
*
* PARAMETERS :
*  II  N  ORDER OF THE MATRIX A.
*  RO  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM
*         WHICH IS SET TO THE UNIT MATRIX (I.E. A:=I).
*
      SUBROUTINE MXDSMI(N,A)
      INTEGER          N
      DOUBLE PRECISION A(*)
      INTEGER          I,M
      M = N* (N+1)/2
      DO 10 I = 1,M
          A(I) = 0.0D0
   10 CONTINUE
      M = 0
      DO 20 I = 1,N
          M = M + I
          A(M) = 1.0D0
   20 CONTINUE
      RETURN
      END
