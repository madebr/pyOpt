CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  MIDACO Wrapper for pyOpt
C  ------------------------
C
C  Authors: 
C    - Martin Schlueter           
C      Theoretical and Computational Optimization Group,
C      School of Mathematics, University of Birmingham, UK.
C    - Ruben E. Perez
C      Department of Mechanical and Aerospace Engineering,
C      Royal Military College of Canada, CANADA.
C
C   MIDACO (Mixed Integer Distributed Ant Colony Optimization)
C   Email: info@midaco-solver.com
C   URL: www.midaco-solver.com
C
C   pyOpt (PYthon OPTimization framework)
C   Email: Ruben E. Perez Ruben.Perez@rmc.ca
C   URL: www.pyopt.org
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine midaco_wrap(L,N,NINT,M,ME,X,XL,XU,F,G,
     &  ACC,PARAM,MAXEVAL,MAXTIME,IFAIL,EVAL,
     &  IPRINT,PRINTEVAL,PRINTBEST,IOUT1,IOUT2,IFILE1,IFILE2,
     &  LICENSE_KEY,
     &  LIW,IW,LRW,RW,OBJFUN)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
      IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     Dimensions of optimization problem        
      INTEGER L,N,NINT,M,ME
C     Lower and upper bounds ('XL' and 'XU') and optimization variable 'X'   
      DOUBLE PRECISION X(L*N),XL(N),XU(N)
C     Objective 'F(X)' and constraints 'G(X)' 
      DOUBLE PRECISION F(L),G(L*M)
C     MIDACO flags, parameters and counters
      INTEGER MAXEVAL,IFAIL,ISTOP,EVAL
      INTEGER IPRINT,PRINTEVAL,PRINTBEST,IOUT1,IOUT2
      DOUBLE PRECISION ACC,PARAM(7),MAXTIME
      CHARACTER*(*) IFILE1,IFILE2
C     License LICENSE_KEY
      CHARACTER*60 LICENSE_KEY
C     MIDACO workspace        
      INTEGER LIW,LRW      
      INTEGER IW(LIW)
      DOUBLE PRECISION RW(LRW)
C     Objective and constraint function
      EXTERNAL OBJFUN
C     Other counters
      INTEGER I
      DOUBLE PRECISION TIME,STARTTIME    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Stopflag, MUST be 0 on first call 
      ISTOP    = 0              
      
C     Get starting time
      CALL CPU_TIME(STARTTIME)

C     Print headline and basic information
      IF (IPRINT.GE.0)THEN
        CALL MIDACOPRINT(1,IPRINT,PRINTEVAL,PRINTBEST,IOUT1,IOUT2,
     &    IFILE1,IFILE2,IFAIL,EVAL,F,G,X,XL,N,NINT,M,ME,TIME,STARTTIME,
     &    RW,LRW,ACC,MAXEVAL,MAXTIME,PARAM,LICENSE_KEY)
      ENDIF

C     Call MIDACO solver by reverse communication       
  100 CONTINUE

C     Evaluate objective function
      CALL OBJFUN(L,N,M,X,F,G)                   ! Evaluate F(X) and G(X)
      EVAL = EVAL + L                            ! Count evaluation

C     Stopping criteria      
      CALL CPU_TIME(TIME)                        ! Get current time	  
      IF(EVAL.GE.MAXEVAL)            ISTOP = 1   ! MAXEVAL criteria	
      IF(TIME-STARTTIME.GE.MAXTIME)  ISTOP = 1   ! MAXTIME criteria	  

C     MIDACO call within the rev.com.loop
      CALL MIDACO(L,N,NINT,M,ME,X,F,G,XL,XU,ACC,
     &     IFAIL,ISTOP,PARAM,RW,LRW,IW,LIW,LICENSE_KEY)

C     Check if MIDACO returns WARNING or ERROR message
      IF(IFAIL.GT.10)THEN
          IF(IFAIL.LT.100)THEN
C              WRITE(*,110) IFAIL
C  110         FORMAT(/,' **** MIDACO-WARNING: IFAIL =',I5,' ****',/) 
              GOTO 100          
          ELSE
C              WRITE(*,120) IFAIL
C  120         FORMAT(/,' **** MIDACO-ERROR: IFAIL =',I5,' ****',/)    
              GOTO 999
          ENDIF
      ENDIF

C     Print best solution after every *PRINTEVAL* and update *BESTX* file
      IF (IPRINT.GE.0)THEN
        CALL MIDACOPRINT(2,IPRINT,PRINTEVAL,PRINTBEST,IOUT1,IOUT2,
     &    IFILE1,IFILE2,IFAIL,EVAL,F,G,X,XL,N,NINT,M,ME,TIME,STARTTIME,
     &    RW,LRW,ACC,MAXEVAL,MAXTIME,PARAM,LICENSE_KEY)
      ENDIF
      
C     Continue rev.com.loop
      IF(ISTOP.EQ.0) GOTO 100

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     Independent check of the MIDACO solution 
  999 CALL OBJFUN(L,N,M,X,F,G)

C     Print the MIDACO solution
      IF (IPRINT.GE.0)THEN
        CALL MIDACOPRINT(3,IPRINT,PRINTEVAL,PRINTBEST,IOUT1,IOUT2,
     &    IFILE1,IFILE2,IFAIL,EVAL,F,G,X,XL,N,NINT,M,ME,TIME,STARTTIME,
     &    RW,LRW,ACC,MAXEVAL,MAXTIME,PARAM,LICENSE_KEY)
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END
