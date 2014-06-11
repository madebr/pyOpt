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
     &  PARAM,MAXEVAL,MAXTIME,IFLAG,NEVAL,
     &  IPRINT,PRINTEVAL,IOUT1,IOUT2,IFILE1,IFILE2,
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
C     MIDACO information and stop flags
      INTEGER IFLAG, ISTOP
C     MIDACO parameter      
      DOUBLE PRECISION PARAM(9)
C     MIDACO workspace        
      INTEGER LIW,LRW      
      INTEGER IW(LIW)
      DOUBLE PRECISION RW(LRW)
C     Parameter for stopping criteria, and printing
      INTEGER MAXTIME,MAXEVAL,NEVAL,PRINTEVAL,IPRINT
      INTEGER IOUT1,IOUT2
      CHARACTER*(*) IFILE1,IFILE2
C     Other counters
      INTEGER I
C     License LICENSE_KEY
      CHARACTER*60 LICENSE_KEY
C     Objective and constraint function
      EXTERNAL OBJFUN
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         
C     Call MIDACO by Reverse Communication
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Stopflag, MUST be 0 on first call 
      ISTOP = 0              
C
C     Print headline and basic information
      CALL MIDACO_PRINT_WRAP(1,PRINTEVAL,IPRINT,IFLAG,ISTOP,
     &  F(1),G,X,XL,XU,N,NINT,M,ME,RW,LRW,MAXEVAL,MAXTIME,PARAM,
     &  1,LICENSE_KEY,IOUT1,IOUT2,IFILE1,IFILE2,NEVAL)
C
C     Call MIDACO solver by reverse communication       
      DO WHILE(ISTOP.EQ.0)

C       Evaluate Objective F(X) and constraints G(X)
        CALL OBJFUN(L,N,M,X,F,G)
C
C       MIDACO call
        CALL MIDACO(L,N,NINT,M,ME,X,F,G,XL,XU,IFLAG,
     &              ISTOP,PARAM,RW,LRW,IW,LIW,LICENSE_KEY)  
C
C        Print best solution after every PRINTEVAL evaluation
         CALL MIDACO_PRINT_WRAP(2,PRINTEVAL,IPRINT,IFLAG,ISTOP,
     &     F(1),G,X,XL,XU,N,NINT,M,ME,RW,LRW,MAXEVAL,MAXTIME,PARAM,
     &     1,LICENSE_KEY,IOUT1,IOUT2,IFILE1,IFILE2,NEVAL)
C
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     Independent check of the MIDACO solution 
      CALL OBJFUN(L,N,M,X,F,G)
C
C     Print the MIDACO solution
      CALL MIDACO_PRINT_WRAP(3,PRINTEVAL,IPRINT,IFLAG,ISTOP,
     &  F(1),G,X,XL,XU,N,NINT,M,ME,RW,LRW,MAXEVAL,MAXTIME,PARAM,
     &  1,LICENSE_KEY,IOUT1,IOUT2,IFILE1,IFILE2,NEVAL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      RETURN
      END





CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine handles all printing commands for MIDACO.
C     Note that this subroutine is called independently from MIDACO and
C     MIDACO itself does not include any print commands (due to 
C     compiler portability and robustness).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MIDACO_PRINT_WRAP(C,PRINTEVAL,SAVE2FILE,IFLAG,ISTOP,
     &           F,G,X,XL,XU,N,NI,M,ME,RW,LRW,MAXEVAL,MAXTIME,PARAM,P,
     &           KEY,IOUT1,IOUT2,IFILE1,IFILE2,NEVAL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC              
      IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER C,PRINTEVAL,SAVE2FILE,IFLAG,ISTOP,EVAL,N,NI,M,ME,LRW
      INTEGER TIC,Q,I,IOUT1,IOUT2,MAXEVAL,MAXTIME,KMAX,UPDATE
      INTEGER KF,KG,KRES,KX,WF,WG,WRES,WX,KBEST,WBEST,P,NEVAL
      DOUBLE PRECISION TNOW,TSTART,TMAX,RW(LRW),BESTF,BESTR,F,G(M),X(N)
      DOUBLE PRECISION XL(N),XU(N),ACC,PARAM(9),DUMMY_F,DUMMY_VIO
C     Increase size, if problems with N > 1000 are solved      
      DOUBLE PRECISION BESTX(1000),BESTG(1000)
      CHARACTER*60 KEY
      CHARACTER*(*) IFILE1,IFILE2
      DATA KF,KG,KRES,KX,WF,WG,WRES,WX /0,0,0,0,0,0,0,0/      
      DATA KBEST,WBEST,TIC,Q,EVAL /0,0,0,0,0/    
      DATA BESTF,BESTR,DUMMY_F,DUMMY_VIO /0.0D0,0.0D0,0.0D0,0.0D0/
      DATA TNOW,TSTART,TMAX,ACC /0.0D0,0.0D0,0.0D0,0.0D0/
      IF(C.EQ.2)THEN                  
          CALL GET_TIME(TNOW)
          TNOW = TNOW - TSTART 
          EVAL = EVAL + P                  
          IF(IFLAG.GE.10)THEN
              IF(SAVE2FILE.EQ.0)THEN          
                  CALL WARNINGS_AND_ERRORS( IFLAG, 0 )              
              ENDIF
              IF(SAVE2FILE.GE.1)THEN
                  CALL WARNINGS_AND_ERRORS( IFLAG, IOUT1 )
                  CALL WARNINGS_AND_ERRORS( IFLAG, IOUT2 )                            
              ENDIF 
              IF(IFLAG.GE.100) RETURN     
          ENDIF  
          IF(PRINTEVAL.GE.1)THEN          
              TIC = TIC + P
              IF(TIC.GE.PRINTEVAL.OR.EVAL.EQ.P.OR.IFLAG.GE.1)THEN
                  IF(EVAL.GT.P) TIC = 0      
                  IF(RW(KRES).EQ.RW(WRES))THEN            
                    KBEST=KF
                    WBEST=WF
                  ELSE                 
                    KBEST=KRES
                    WBEST=WRES
                  ENDIF                 
                  IF(RW(WBEST).LT.RW(KBEST).OR.
     &              (IFLAG.GE.1.OR.IFLAG.EQ.-300))THEN                     
                      BESTF = RW(WF)
                      BESTR = RW(WRES)
                      DO I=1,M
                        BESTG(I) = RW(WG+I)                                       
                      ENDDO
                      DO I=1,N
                        BESTX(I) = RW(WX+I)
                      ENDDO
                  ELSE                  
                      BESTF = RW(KF)
                      BESTR = RW(KRES)
                      DO I=1,M
                        BESTG(I) = RW(KG+I)                    
                      ENDDO      
                      DO I=1,N
                        BESTX(I) = RW(KX+I)
                      ENDDO                      
                  ENDIF
                  IF(SAVE2FILE.EQ.0)THEN                                     
                    CALL PRINT_LINE(EVAL,TNOW,BESTF,BESTR,0)
                  ENDIF             
                  IF(SAVE2FILE.GE.1)THEN                    
                    CALL PRINT_LINE(EVAL,TNOW,BESTF,BESTR,IOUT1)
                  ENDIF               
                  IF(SAVE2FILE.GE.1)THEN   
                    UPDATE = 0                 
                    IF( (BESTR.LT.DUMMY_VIO) .OR.
     &                  (BESTR.EQ.DUMMY_VIO.AND.BESTF.LT.DUMMY_F) )THEN 
                       DUMMY_F   = BESTF
                       DUMMY_VIO = BESTR
                       UPDATE    = 1 
                    ENDIF                
                    
                    IF(UPDATE.EQ.1)THEN
                     WRITE(IOUT2,31)
   31                FORMAT(/,/,'            CURRENT BEST SOLUTION') 
                     CALL PRINT_SOLUTION( N, M, ME, BESTX, BESTG, BESTF, 
     &                     BESTR, XL, XU, ACC, EVAL, TNOW, IFLAG, IOUT2)
                     CALL FORCE_OUTPUT( IOUT1 )
                     CALL FORCE_OUTPUT( IOUT2 )
                     CONTINUE     
                    ENDIF
                  ENDIF                       
              ENDIF  
          ENDIF    
          IF(ISTOP.EQ.0)THEN
              IF(TNOW.GE.TMAX     ) IFLAG = -999
              IF(EVAL.GE.MAXEVAL-1)THEN            
                IF(MAXEVAL.LE.99999999) IFLAG = -999
              ENDIF
          ENDIF       
          RETURN          
      ENDIF      
      IF(C.EQ.1)THEN             
          IFLAG = 0
          ISTOP = 0      
          TMAX = DBLE(MAXTIME)           
          CALL GET_TIME(TSTART)        
          EVAL = 0                
          if(param(1).le.0.0D0)then
              acc = 1.0D-3
          else
              acc = param(1)
          endif              
          KMAX = 100
          Q    = 2*N+M+(N+5)*KMAX+8                            
          KX   = 1           
          KF   = 2+N
          KG   = 2+N
          KRES = 2+N+1+M 
          WX   = Q
          WF   = Q+1+N
          WG   = Q+1+N
          WRES = Q+1+N+1+M          
          IF(SAVE2FILE.GE.1)THEN                     
              OPEN(IOUT1,FILE=IFILE1,STATUS='UNKNOWN')
              OPEN(IOUT2,FILE=IFILE2,STATUS='UNKNOWN')              
          ENDIF     
          BESTF = 1.0D+32
          BESTR = 1.0D+32 
          DUMMY_F   = 1.0D+32
          DUMMY_VIO = 1.0D+32      
          TIC = 0             
          IF(PRINTEVAL.GE.1)THEN
              IF(SAVE2FILE.EQ.0)THEN
                CALL PRINT_HEAD( N, NI, M, ME, PARAM, MAXEVAL, MAXTIME,
     &                         PRINTEVAL, SAVE2FILE, KEY, 0)
              ENDIF 
              IF(SAVE2FILE.GE.1)THEN
                CALL PRINT_HEAD( N, NI, M, ME, PARAM, MAXEVAL, MAXTIME, 
     &                        PRINTEVAL, SAVE2FILE, KEY, IOUT1)
              ENDIF            
          ENDIF      
          IF(SAVE2FILE.GE.1) THEN
             WRITE(IOUT2,3)
    3        FORMAT('MIDACO - SOLUTION',/,
     &           '-----------------',/,
     &'This file saves the current best solution X found by MIDACO.',/,
     &'This file is updated after every PRINTEVAL function evaluation,',
     &/,'if X has been improved.',/,/)              
             CALL FORCE_OUTPUT( IOUT1 )
             CALL FORCE_OUTPUT( IOUT2 )
           ENDIF                                           
      ENDIF
      IF(C.EQ.3.AND.PRINTEVAL.GE.1)THEN  
          IF(SAVE2FILE.EQ.0)THEN
             CALL PRINT_FINAL( IFLAG,TNOW,TMAX,EVAL,MAXEVAL,
     &                      N,M,ME,X,G,F,XL,XU,RW,ACC,WRES,PARAM, 
     &                      0 ) 
          ENDIF
          IF(SAVE2FILE.GE.1)THEN
             CALL PRINT_FINAL( IFLAG,TNOW,TMAX,EVAL,MAXEVAL,
     &                         N,M,ME,X,G,F,XL,XU,RW,ACC,WRES,PARAM, 
     &                         IOUT1 )
             CALL PRINT_FINAL( IFLAG,TNOW,TMAX,EVAL,MAXEVAL,
     &                         N,M,ME,X,G,F,XL,XU,RW,ACC,WRES,PARAM, 
     &                         IOUT2 )  
             CALL FORCE_OUTPUT( IOUT1 )
             CALL FORCE_OUTPUT( IOUT2 )
          ENDIF
          NEVAL = EVAL
      ENDIF     
      END      
