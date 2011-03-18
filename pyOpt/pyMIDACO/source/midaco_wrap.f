CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   
C   MIDACO 0.3 beta version
C   -----------------------
C
C   MIDACO solves the general mixed integer non-linear program (MINLP):
C
C             minimize    f(x,y)              (x continuous, y integer)
C
C             subject to  g(j)  =  0          ( j = 1, ..., me)
C                         g(j) >=  0          ( j = me + 1, ..., m)
C                         xl   <=  x  <=  xu  
C                         yl   <=  x  <=  yu 
C
C   MIDACO is a stochastic global optimization solver based on a mixed
C   integer extension of the ant colony optimization heuristic described
C   in [1] and [2]. Constraints are handled by the oracle penalty method 
C   approach described in [3].
C   The optimization problem is considered as a black-box modell and no
C   derivative information is required. Only the objective and constraint
C   function values must be delivered. Hence MIDACO is well applicable
C   on non-convex, discontinuous and general non-smooth problems.   
C   MIDACO also excepts purely continuous (NLP) and purely combinatorial 
C   (COP) optimization problems. 
C   
C   With regard to the upmost user friendliness MIDACO sets all its algo-
C   rithmic parameters by itself (like kernel size, population size, 
C   oracles,...). This selftuning comes with the price of an increased 
C   amount of function evaluations. Nevertheless the MIDACO kernel (based 
C   on a FORTRAN 77 reverse communication interface) provides fast 
C   internal calculation times.
C
C   MIDACO 0.3 and higher features an option to parallelize the objective
C   functions calls. 
C
C
C   References:
C
C       [1] Schlueter M., Egea J.A., Banga J.R.: "Extended ant colony
C           optimization for nonconvex mixed integer nonlinear programming"
C           Computers and Operations Research, Vol. 36(7), pp 2217-2229, 
C           2009.
C
C       [2] Schlueter M., Egea J.A., Antelo L.T., Alonso A.A., Banga J.R.: 
C           "An extended ant colony optimization algorithm for integrated 
C           process and control system design", Industrial & Engineering
C           Chemistry, Vol. 48 (14), pp 6723–6738, 2009.
C
C       [3] Schlueter M., Gerdts M.: "The oracle penalty method",
C           Journal of Global Optimization,accepted (2009)
C
C   Author: Martin Schlueter           
C           Theoretical and Computational Optimization Group,
C           School of Mathematics, University of Birmingham, UK.
C
C   Email:  info@midaco-solver.com
C   URL:    www.midaco-solver.com
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine midaco_wrap(N,NINT,M,ME,XL,XU,X,F,G,
     &  ACC,ISEED,IFAIL,QSTART,ISTOP,IPRINT,IOUT,IFILE,
     &  EVAL,MAXEVAL,MAXTIME,
     &  LIW,IW,LRW,RW,
     &  OBJFUN)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C       This parameter declares the amount of ants (iterates) which can 
C       be processed at once by MIDACO within one reverse communication 
C       step. Hence, this parameter defines the amount of objective  
C       function evaluation within one reverse communication step.  
C       Therefore this is the crucial parameter for parallelization.
        INTEGER L
        PARAMETER ( L = 1 )  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
C       Dimensions of test problem        
        INTEGER N,NINT,M,ME
C       Bounds and optimization variable   
        DOUBLE PRECISION XL(N),XU(N),X(N)
C       Objective and constraint function values
        EXTERNAL OBJFUN
        DOUBLE PRECISION F,G(M)  
c       Dummy fields for parallel ants (iterates)        
        DOUBLE PRECISION XXX(N*L),FFF(L),GGG(M*L+1)              
C       MIDACO parameter
        INTEGER IFAIL,ISEED,ISTOP,QSTART
        DOUBLE PRECISION ACC
C       MIDACO workspace        
        INTEGER LIW,LRW      
        INTEGER IW(LIW)
        DOUBLE PRECISION RW(LRW)
C       Parameter for reverse communication loop
        INTEGER EVAL,MAXEVAL,I,J
        DOUBLE PRECISION TIME,STARTTIME,MAXTIME 
C       Parameter needed for print        
        INTEGER ENDRW,ENDIW,KMAX,FIRST
        INTEGER IPRINT,IOUT
        CHARACTER*(*) IFILE   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  
C       Copy user defined starting point in dummy fields       
        DO I = 1,N   
          XXX(I) = X(I)             
        ENDDO
        
C       Parameter needed for print         
        KMAX  = 2*N+2  
        ENDRW = 2*N+M+(N+5)*KMAX+8 
        ENDIW = 31+L+N+N 
        FIRST = 1
                          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         

C       Initialize time ticker
        CALL CPU_TIME(STARTTIME)  ! Get starting time
                  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF (IPRINT.EQ.0)THEN
          WRITE(*,1) ! Print headline
      ENDIF
      IF (IPRINT.GT.0)THEN
           OPEN(UNIT=IOUT,FILE=IFILE(1:LEN_TRIM(IFILE)),
     &        STATUS='UNKNOWN')
          WRITE(IOUT,1) ! Print headline
      ENDIF
      
    1 FORMAT(/,' MIDACO 0.3 BETA VERSION   (www.midaco-solver.com)',/,
     &   /, ' [     EVAL,   TIME]        OBJECTIVE FUNCTION VALUE',
     &       '         RESIDUAL VALUE',/,
     &       ' ---------------------------------------------------',
     &       '-----------------------')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
C
C     Call MIDACO solver by reverse communication
C                   
  100 CONTINUE !~~~~~~~~~~~~~~~~~~~~~start~of~reverse~communication~loop
  
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
C     !!!!!!!!!!!!!!!!!! This loop can be parallelized !!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO J = 1,L    
          DO I =1,N
              X(I) = XXX((J-1)*N+I)   ! Get X out of XXX  
          ENDDO
          
          CALL OBJFUN(N,M,X,F,G)      ! Evaluate objective function
                                      ! F & constraints G for iterate X
                                                   
          EVAL = EVAL + 1             ! Count evaluations          
          FFF(J) = F                  ! Store F in FFF
          DO I =1,M
              GGG((J-1)*M+I) = G(I)   ! Store G in GGG
          ENDDO
      ENDDO
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

C     Dummy for MIDACO
      GGG(23*L+1) = 0.0D0
      
C     Stopping criteria      
      CALL CPU_TIME(TIME)                          ! Get the time	   
      IF(EVAL.GE.MAXEVAL)           ISTOP = 1      ! MAXEVAL criteria	
      IF(TIME-STARTTIME.GE.MAXTIME) ISTOP = 1      ! MAXTIME criteria	  
      
      CALL MIDACO(L,N,NINT,M,ME,XXX,FFF,GGG,XL,XU,ACC,
     &            IFAIL,ISTOP,QSTART,ISEED,RW,LRW,IW,LIW)
          
C     Check if MIDACO returns ERROR-MESSAGE
      IF(IFAIL.GE.10)THEN
          print*,'***ERROR-MESSAGE***   IFAIL = ',IFAIL
          pause
      ENDIF       
      
C     Print current best solution
      IF(IFAIL.EQ.-3.AND.FIRST.EQ.1.OR.EVAL.LE.L)THEN
          IF (IPRINT.EQ.0)THEN
             WRITE(*,2) EVAL,TIME-STARTTIME,RW(2+N),RW(2+N+1+M)
          ENDIF
          IF (IPRINT.GT.0)THEN
             WRITE(IOUT,2) EVAL,TIME-STARTTIME,RW(2+N),RW(2+N+1+M)
          ENDIF
      ENDIF
      IF(IFAIL.EQ.-300.OR.IFAIL.GE.0)THEN
       IF(FIRST.EQ.1.AND.IFAIL.LT.0)THEN
        IF(IPRINT.EQ.0)THEN
         WRITE(*,*) 'Print frequency is slowed down from now on...'
        ENDIF
        IF(IPRINT.GT.0)THEN
         WRITE(IOUT,*) 'Print frequency is slowed down from now on...'
        ENDIF
        FIRST = 0
       ENDIF
       IF(IPRINT.EQ.0)THEN
        WRITE(*,2) EVAL,TIME-STARTTIME,RW(ENDRW+1+N),RW(ENDRW+1+N+1+M) 
       ENDIF
       IF(IPRINT.GT.0)THEN
        WRITE(IOUT,2) EVAL,TIME-STARTTIME,RW(ENDRW+1+N),
     &    RW(ENDRW+1+N+1+M) 
       ENDIF
      ENDIF
    2 FORMAT(' [',I9,',',F7.2,']        F:',F22.8,'         RES:',D10.3)
       
      IF(IFAIL.LT.0) GOTO 100 !~~~~~~~~end~of~reverse~communication~loop
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Get solution out of dummy fields
      F = FFF(1)
      DO I = 1,n
          X(I) = XXX(I)
      ENDDO
      DO I = 1,M
          G(I) = GGG(I)
      ENDDO            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
C     
C     Independent test of solution X given back by MIDACO     
C
      CALL OBJFUN(N,M,X,F,G)      
C      
C     Print solution  
C   
      IF (IPRINT.EQ.0)THEN
         WRITE(*,3) F
      ENDIF
      IF (IPRINT.GT.0)THEN
         WRITE(IOUT,3) F
      ENDIF
      DO I = 1,ME
          IF (IPRINT.EQ.0)THEN
             WRITE(*,4) I,G(I)
          ENDIF
          IF (IPRINT.GT.0)THEN
             WRITE(IOUT,4) I,G(I)
          ENDIF
      ENDDO   
      IF (IPRINT.EQ.0)THEN
         WRITE(*,*)         
      ENDIF
      IF (IPRINT.GT.0)THEN
         WRITE(IOUT,*)
      ENDIF  
      DO I = ME+1 , M
        IF (IPRINT.EQ.0)THEN
          WRITE(*,5) I,G(I)
        ENDIF
        IF (IPRINT.GT.0)THEN
          WRITE(IOUT,5) I,G(I)
        ENDIF
      ENDDO
      IF (IPRINT.EQ.0)THEN
         WRITE(*,*)         
      ENDIF
      IF (IPRINT.GT.0)THEN
         WRITE(IOUT,*)
      ENDIF
      DO I = 1,N
        IF (IPRINT.EQ.0)THEN
          WRITE(*,6) I,X(I)
        ENDIF
        IF (IPRINT.GT.0)THEN
          WRITE(IOUT,6) I,X(I)
        ENDIF
      ENDDO 
      IF (IPRINT.EQ.0)THEN
         WRITE(*,*)         
      ENDIF
      IF (IPRINT.GT.0)THEN
         WRITE(IOUT,*)
      ENDIF   
      
    3 FORMAT(/,/,/,' Independent test of MIDACO solution',/,
     &             ' -----------------------------------',/,
     &             ' Objective function value      F(X) = ',F22.8,/) 
    4 FORMAT(' Equality   constraint   0  = G(',I2,') = ',F22.8)     
    5 FORMAT(' Inequality constraint   0 <= G(',I2,') = ',F22.8)   
    6 FORMAT(' Solution variable            X(',I2,') = ',F22.8)    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RETURN
      END
