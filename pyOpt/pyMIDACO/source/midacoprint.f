CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine handles all printing commands for MIDACO used in
C     the example templates that are distributed along with MIDACO.
C     Note that this subroutine is called independently from MIDACO and
C     MIDACO itself does not include any print commands (due to 
C     compiler portability and robustness). DO NOT CHANGE THIS FILE!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MIDACOPRINT(C,IPRINT,PRINTEVAL,PRINTBEST,IOUT1,IOUT2,
     &    IFILE1,IFILE2,IFAIL,EVAL,F,G,X,XL,N,NINT,M,ME,TIME,STARTTIME,
     &    RW,LRW,ACC,MAXEVAL,MAXTIME,PARAM,LICENSE_KEY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC              
      IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER C,PRINTEVAL,IPRINT,IFAIL,EVAL,N,NINT,M,ME,LRW
      INTEGER TIC,TOC,R,I,IOUT1,IOUT2,PRINTBEST,MAXEVAL
      INTEGER KF,KRES,KX,WF,WRES,WX,KFRES,WFRES
      DOUBLE PRECISION TIME,STARTTIME,RW(LRW),BESTF,BESTR,F(1),G(M),X(N)
      DOUBLE PRECISION XL(N),Z(N),ACC,MAXTIME,PARAM(7)
      CHARACTER*60 LICENSE_KEY
      CHARACTER*(*) IFILE1,IFILE2
      DATA KFRES,WFRES /0,0/
      DATA TIC,TOC,R /0,0,0/
      DATA KF,KRES,KX,WF,WRES,WX /0,0,0,0,0,0/
      DATA BESTF,BESTR /0.0D0,0.0D0/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    1 FORMAT(/,' MIDACO 3.0   (www.midaco-solver.com)',/,
     &         ' ------------------------------------',/,/,
     &         ' LICENSE-KEY:  ',A60,/,/,           
     &         ' ------------------------------------',/,
     &         ' | N',I7,'   | ACC',F16.9,' |',/,
     &         ' | NINT',I4,'   | MAXEVAL',I12,' |',/,     
     &         ' | M',I7,'   | MAXTIME',F12.1,' |',/,     
     &         ' | ME',I6,'   | PRINTEVAL',I10,' |',/,
     &         ' |----------------------------------|')        
  101 FORMAT(' | PARAMETER:  All by default (0)   |')            
  102 FORMAT(' | PARAM(1)',F10.1,' (RANDOM-SEED) |',/,  
     &       ' | PARAM(2)',F10.1,' (QSTART)      |',/,
     &       ' | PARAM(3)',F10.1,' (AUTOSTOP)    |',/,     
     &       ' | PARAM(4)',F10.1,' (ORACLE)      |',/,     
     &       ' | PARAM(5)',F10.1,' (ANTS)        |',/,    
     &       ' | PARAM(6)',F10.1,' (KERNEL)      |',/,     
     &       ' | PARAM(7)',F10.1,' (CHARACTER)   |')     
  103 FORMAT(' ------------------------------------',/,/,  
     &' [     E',
     &'VAL,    TIME]        OBJECTIVE FUNCTION VALUE         RESIDUAL ',
     &'VALUE',/,' ----------------------------------------------------',
     &'-----------------------') 
  104 FORMAT(' ------------------------------------',/,/,  
     &' [     E',
     &'VAL,    TIME]        OBJECTIVE FUNCTION VALUE         RESIDUAL ',
     &'VALUE  |   SOLUTION VECTOR',/,' -------------------------------',
     &'----------------------------------------------|----------------',
     &'--') 
    2 FORMAT(' [',I9,',',F8.2,']        F(X):',F19.8,'         RES:',
     &D10.3)
   20 FORMAT(' [',I9,',',F8.2,']        F(X):',F19.8,'         RES:',
     &D10.3,'  |   X:',1000D16.8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
    3 FORMAT('MIDACO - Current Best Solution',/,
     &         '------------------------------',/,/,
     &         'This file saves the current best solution X every:',I7,
     &         ' evaluation (if X is improved)',/,/,
     &         'Hint: Use these solutions by copy and paste for later',
     &         ' refinements (QSTART=PARAM(2))')
   31 FORMAT(/,/,'EVAL:',I10,/,'TIME:',F10.2,/,'IFAIL:',I9,/,
     &'-----------------------------------')
   32 FORMAT('F(X)    = ',D25.18,/,'RES(X)  = ',F25.12,/,
     &'-----------------------------------')
   33 FORMAT('X(',I4,') = ',D25.18)      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
    4 FORMAT(/,/,' INDEPENDENT TEST OF SOLUTION GIVEN BY MIDACO',/,
     &       ' --------------------------------------------')
   41 FORMAT(' [EVAL:',I10,', TIME:',F8.2,', IFAIL:',I4,']',/,
     &       ' --------------------------------------------')
   42 FORMAT(' F(X) = ',F37.18)
  142 FORMAT(' --------------------------------------------')
   43 FORMAT(' G(',I4,') = ',F15.9,'  (equality constr)')
   44 FORMAT(' G(',I4,') = ',F15.9,'  (in-equal constr)')
  431 FORMAT(' G(',I4,') = ',F15.9,'  (equality constr) <- infeasible')
  441 FORMAT(' G(',I4,') = ',F15.9,'  (in-equal constr) <- infeasible')
   
  145 FORMAT(' --------------------------------------------',
     &                              '        ---BOUNDS-PROFILER---')  
  400 FORMAT(' X(',I4,') = ',F34.18,'     !  XL___________________') 
  401 FORMAT(' X(',I4,') = ',F34.18,'     !  x____________________')
  402 FORMAT(' X(',I4,') = ',F34.18,'     !  _x___________________')
  403 FORMAT(' X(',I4,') = ',F34.18,'     !  __x__________________')  
  404 FORMAT(' X(',I4,') = ',F34.18,'     !  ___x_________________')
  405 FORMAT(' X(',I4,') = ',F34.18,'     !  ____x________________')
  406 FORMAT(' X(',I4,') = ',F34.18,'     !  _____x_______________')
  407 FORMAT(' X(',I4,') = ',F34.18,'     !  ______x______________')
  408 FORMAT(' X(',I4,') = ',F34.18,'     !  _______x_____________')
  409 FORMAT(' X(',I4,') = ',F34.18,'     !  ________x____________')
  410 FORMAT(' X(',I4,') = ',F34.18,'     !  _________x___________')
  411 FORMAT(' X(',I4,') = ',F34.18,'     !  __________x__________')
  412 FORMAT(' X(',I4,') = ',F34.18,'     !  ___________x_________')
  413 FORMAT(' X(',I4,') = ',F34.18,'     !  ____________x________')
  414 FORMAT(' X(',I4,') = ',F34.18,'     !  _____________x_______')
  415 FORMAT(' X(',I4,') = ',F34.18,'     !  ______________x______')
  416 FORMAT(' X(',I4,') = ',F34.18,'     !  _______________x_____')
  417 FORMAT(' X(',I4,') = ',F34.18,'     !  ________________x____')
  418 FORMAT(' X(',I4,') = ',F34.18,'     !  _________________x___')
  419 FORMAT(' X(',I4,') = ',F34.18,'     !  __________________x__')
  420 FORMAT(' X(',I4,') = ',F34.18,'     !  ___________________x_')
  421 FORMAT(' X(',I4,') = ',F34.18,'     !  ____________________x')  
  422 FORMAT(' X(',I4,') = ',F34.18,'     !  ___________________XU')  
  490 FORMAT(' X(',I4,') = ',F34.18,'     !  WARNING: XL = XU     ')  
  491 FORMAT(' X(',I4,') = ',F34.18,' ***ERROR*** (X > XU)        ')  
  492 FORMAT(' X(',I4,') = ',F34.18,' ***ERROR*** (X < XL)        ')  
  493 FORMAT(' X(',I4,') = ',F34.18,' ***ERROR*** (XL > XU)       ')  
   47 FORMAT(/,' ')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(C.EQ.1)THEN      
          R    = 2*N+M+(N+5)*(2*N+10)+8          
          KF   = 2+N
          KRES = 2+N+1+M
          KX   = 1                    
          WF   = R+1+N
          WRES = R+1+N+1+M
          WX   = R
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Open files for saving 
      IF(C.EQ.1)THEN
          IF(IPRINT.GE.1)THEN
              OPEN(IOUT1,FILE=IFILE1,STATUS='UNKNOWN')
              IF(PRINTBEST.GE.1)THEN
                  OPEN(IOUT2,FILE=IFILE2,STATUS='UNKNOWN')
                  BESTF = 1.0D+32
                  BESTR = 1.0D+32
              ENDIF  
          ENDIF
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Print Headlines
      IF(C.EQ.1)THEN
          IF(IPRINT.EQ.0)THEN
              IF(PRINTEVAL.GE.1) WRITE(*,1) LICENSE_KEY,N,ACC,
     &                     NINT,MAXEVAL,M,MAXTIME,ME,PRINTEVAL
              DO I = 1,7
                  IF(PARAM(I).NE.0.0D0) GOTO 66
              ENDDO
              IF(PRINTEVAL.GE.1) WRITE(*,101)     
              GOTO 67
   66         IF(PRINTEVAL.GE.1) WRITE(*,102) PARAM(1),PARAM(2),
     &              PARAM(3),PARAM(4),PARAM(5),PARAM(6),PARAM(7)
   67         CONTINUE
              IF(PRINTEVAL.GE.1) WRITE(*,103)
          ENDIF 
          IF (IPRINT.GT.0)THEN
              IF(PRINTEVAL.GE.1) WRITE(IOUT1,1) LICENSE_KEY,N,ACC,
     &                         NINT,MAXEVAL,M,MAXTIME,ME,PRINTEVAL
              DO I = 1,7
                  IF(PARAM(I).NE.0.0D0) GOTO 661
              ENDDO
              IF(PRINTEVAL.GE.1) WRITE(IOUT1,101)     
              GOTO 671
  661         IF(PRINTEVAL.GE.1) WRITE(IOUT1,102) PARAM(1),PARAM(2),
     &                  PARAM(3),PARAM(4),PARAM(5),PARAM(6),PARAM(7)
  671         CONTINUE
              IF(PRINTEVAL.GE.1) WRITE(IOUT1,104)
              IF(PRINTBEST.GE.1) WRITE(IOUT2,3) PRINTBEST
          ENDIF
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     Print best solution after every *PRINTEVAL* evaluation
      IF(C.EQ.2.AND.PRINTEVAL.GE.1)THEN
        TIC = TIC + 1
        IF(TIC.GE.PRINTEVAL.OR.EVAL.EQ.1.OR.IFAIL.GE.1)THEN
            IF(EVAL.NE.1)  TIC = 0              
            IF(RW(KRES).EQ.RW(WRES))THEN
              KFRES=KF
              WFRES=WF
            ELSE
              KFRES=KRES
              WFRES=WRES
            ENDIF
            IF(RW(KFRES).GT.RW(WFRES)
     &            .OR.(IFAIL.GE.1.OR.IFAIL.EQ.-300))THEN
              IF(IPRINT.EQ.0)THEN
                  WRITE(*,2) EVAL,TIME-STARTTIME,RW(WF),RW(WRES)
              ENDIF
              IF(IPRINT.GT.0)THEN
                DO I = 1,N
                  Z(I) = RW(WX+I)          
                ENDDO
                WRITE(IOUT1,20) EVAL,TIME-STARTTIME,RW(WF),RW(WRES),Z
              ENDIF                     
            ELSE
              IF(IPRINT.EQ.0)THEN
                  WRITE(*,2) EVAL,TIME-STARTTIME,RW(KF),RW(KRES)
              ENDIF
              IF(IPRINT.GT.0)THEN
                DO I = 1,N
                  Z(I) = RW(KX+I)          
                ENDDO
                WRITE(IOUT1,20) EVAL,TIME-STARTTIME,RW(KF),RW(KRES),Z
              ENDIF
            ENDIF
c           FLUSH() is a non-standard command to flush output units.
c           Comment below flush commands in case of compilation problems.               
c            FLUSH(IOUT1)
c            FLUSH(IOUT2)
        ENDIF
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Save the current best solution X to file
      IF(C.EQ.2.AND.PRINTBEST.GE.1)THEN
          TOC = TOC + 1
          IF(TOC.GE.PRINTBEST.OR.EVAL.EQ.1.OR.IFAIL.GE.1)THEN
              IF(EVAL.NE.1) TOC = 0              
              IF(IFAIL.GE.1.OR.IFAIL.EQ.-300)THEN
                          BESTF = RW(WF)
                          BESTR = RW(WRES)
                          IF(IPRINT.GT.0)THEN
                              WRITE(IOUT2,31) EVAL,TIME-STARTTIME,IFAIL
                              WRITE(IOUT2,32) BESTF,BESTR  
                              DO I = 1,N
                                WRITE(IOUT2,33) I,RW(WX+I)
                              ENDDO
                              GOTO 555
                          ENDIF
              ENDIF
              IF(RW(KRES).EQ.RW(WRES))THEN
                  IF(RW(KF).LE.RW(WF))THEN
                      IF(RW(KRES).LE.BESTR.AND.RW(KF).LT.BESTF)THEN
                          BESTF = RW(KF)
                          BESTR = RW(KRES)
                          IF(IPRINT.GT.0)THEN
                              WRITE(IOUT2,31) EVAL,TIME-STARTTIME,IFAIL
                              WRITE(IOUT2,32) BESTF,BESTR  
                              DO I = 1,N
                                  WRITE(IOUT2,33) I,RW(KX+I)
                              ENDDO
                          ENDIF
                      ENDIF
                  ELSE
                      IF(RW(WRES).LE.BESTR.AND.RW(WF).LT.BESTF)THEN
                          BESTF = RW(WF)
                          BESTR = RW(WRES)
                          IF(IPRINT.GT.0)THEN
                              WRITE(IOUT2,31) EVAL,TIME-STARTTIME,IFAIL
                              WRITE(IOUT2,32) BESTF,BESTR  
                              DO I = 1,N
                                  WRITE(IOUT2,33) I,RW(WX+I)
                              ENDDO
                          ENDIF
                      ENDIF   
                  ENDIF
              ELSE
                  IF(RW(KRES).LT.RW(WRES))THEN
                      IF(RW(KRES).LT.BESTR.OR.
     &                  (RW(KRES).EQ.BESTR.AND.RW(KF).LT.BESTF))THEN
                          BESTF = RW(KF)
                          BESTR = RW(KRES)
                          IF(IPRINT.GT.0)THEN
                              WRITE(IOUT2,31) EVAL,TIME-STARTTIME,IFAIL
                              WRITE(IOUT2,32) BESTF,BESTR  
                              DO I = 1,N
                                  WRITE(IOUT2,33) I,RW(KX+I)
                              ENDDO
                          ENDIF
                      ENDIF   
                  ELSE
                      IF(RW(WRES).LT.BESTR.OR.
     &                  (RW(WRES).EQ.BESTR.AND.RW(WF).LT.BESTF))THEN
                          BESTF = RW(WF)
                          BESTR = RW(WRES)
                          IF(IPRINT.GT.0)THEN
                              WRITE(IOUT2,31) EVAL,TIME-STARTTIME,IFAIL
                              WRITE(IOUT2,32) BESTF,BESTR  
                              DO I = 1,N
                                  WRITE(IOUT2,33) I,RW(WX+I)
                              ENDDO
                          ENDIF
                      ENDIF      
                  ENDIF
             ENDIF  
c            FLUSH() is a non-standard command to flush output units.
c            Comment below flush commands in case of compilation problems.               
c             FLUSH(IOUT1)
c             FLUSH(IOUT2)      
         ENDIF
      ENDIF   
  555 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     Print information on independet solution check
      IF(C.EQ.3.AND.PRINTEVAL.GT.0)THEN
          IF(IPRINT.EQ.0)THEN
              WRITE(*,4)
              WRITE(*,41) EVAL,TIME-STARTTIME,IFAIL
              WRITE(*,42) F  
              IF(M.GT.0) WRITE(*,142)          
                  DO I = 1,ME
                      IF(DABS(G(I)).LE.ACC)THEN
                          WRITE(*,43) I,G(I)
                      ELSE
                          WRITE(*,431) I,G(I)                  
                      ENDIF
                  ENDDO
                  DO I = ME+1,M
                      IF(G(I).GE.-ACC)THEN
                          WRITE(*,44) I,G(I)
                      ELSE
                          WRITE(*,441) I,G(I)                  
                      ENDIF
                  ENDDO              
                  WRITE(*,145)
                  DO I = 1,N
                     IF(XL(I).EQ. 0.0D0) WRITE(*,400) I,X(I)
                     IF(XL(I).EQ. 1.0D0) WRITE(*,401) I,X(I)
                     IF(XL(I).EQ. 2.0D0) WRITE(*,402) I,X(I)
                     IF(XL(I).EQ. 3.0D0) WRITE(*,403) I,X(I)
                     IF(XL(I).EQ. 4.0D0) WRITE(*,404) I,X(I)
                     IF(XL(I).EQ. 5.0D0) WRITE(*,405) I,X(I)
                     IF(XL(I).EQ. 6.0D0) WRITE(*,406) I,X(I)
                     IF(XL(I).EQ. 7.0D0) WRITE(*,407) I,X(I)
                     IF(XL(I).EQ. 8.0D0) WRITE(*,408) I,X(I)
                     IF(XL(I).EQ. 9.0D0) WRITE(*,409) I,X(I)
                     IF(XL(I).EQ.10.0D0) WRITE(*,410) I,X(I)
                     IF(XL(I).EQ.11.0D0) WRITE(*,411) I,X(I)
                     IF(XL(I).EQ.12.0D0) WRITE(*,412) I,X(I)
                     IF(XL(I).EQ.13.0D0) WRITE(*,413) I,X(I)
                     IF(XL(I).EQ.14.0D0) WRITE(*,414) I,X(I)
                     IF(XL(I).EQ.15.0D0) WRITE(*,415) I,X(I)
                     IF(XL(I).EQ.16.0D0) WRITE(*,416) I,X(I)
                     IF(XL(I).EQ.17.0D0) WRITE(*,417) I,X(I)
                     IF(XL(I).EQ.18.0D0) WRITE(*,418) I,X(I)
                     IF(XL(I).EQ.19.0D0) WRITE(*,419) I,X(I)
                     IF(XL(I).EQ.20.0D0) WRITE(*,420) I,X(I)    
                     IF(XL(I).EQ.21.0D0) WRITE(*,421) I,X(I)
                     IF(XL(I).EQ.22.0D0) WRITE(*,422) I,X(I)  
                     IF(XL(I).EQ.90.0D0) WRITE(*,490) I,X(I)
                     IF(XL(I).EQ.91.0D0) WRITE(*,491) I,X(I)
                     IF(XL(I).EQ.92.0D0) WRITE(*,492) I,X(I)
                     IF(XL(I).EQ.93.0D0) WRITE(*,493) I,X(I)         
                  ENDDO
              WRITE(*,47)
          ENDIF
          IF(IPRINT.GT.0)THEN
              WRITE(IOUT1,4)
              WRITE(IOUT1,41) EVAL,TIME-STARTTIME,IFAIL
              WRITE(IOUT1,42) F  
              IF(M.GT.0) WRITE(IOUT1,142)              
              DO I = 1,ME
                  IF(DABS(G(I)).LE.ACC)THEN
                      WRITE(IOUT1,43) I,G(I)
                  ELSE
                      WRITE(IOUT1,431) I,G(I)                  
                  ENDIF
              ENDDO
              DO I = ME+1,M
                  IF(G(I).GE.-ACC)THEN
                      WRITE(IOUT1,44) I,G(I)
                  ELSE
                      WRITE(IOUT1,441) I,G(I)                  
                  ENDIF
              ENDDO
              WRITE(IOUT1,145)   
              DO I = 1,N
                 IF(XL(I).EQ. 0.0D0) WRITE(IOUT1,400) I,X(I)
                 IF(XL(I).EQ. 1.0D0) WRITE(IOUT1,401) I,X(I)
                 IF(XL(I).EQ. 2.0D0) WRITE(IOUT1,402) I,X(I)
                 IF(XL(I).EQ. 3.0D0) WRITE(IOUT1,403) I,X(I)
                 IF(XL(I).EQ. 4.0D0) WRITE(IOUT1,404) I,X(I)
                 IF(XL(I).EQ. 5.0D0) WRITE(IOUT1,405) I,X(I)
                 IF(XL(I).EQ. 6.0D0) WRITE(IOUT1,406) I,X(I)
                 IF(XL(I).EQ. 7.0D0) WRITE(IOUT1,407) I,X(I)
                 IF(XL(I).EQ. 8.0D0) WRITE(IOUT1,408) I,X(I)
                 IF(XL(I).EQ. 9.0D0) WRITE(IOUT1,409) I,X(I)
                 IF(XL(I).EQ.10.0D0) WRITE(IOUT1,410) I,X(I)
                 IF(XL(I).EQ.11.0D0) WRITE(IOUT1,411) I,X(I)
                 IF(XL(I).EQ.12.0D0) WRITE(IOUT1,412) I,X(I)
                 IF(XL(I).EQ.13.0D0) WRITE(IOUT1,413) I,X(I)
                 IF(XL(I).EQ.14.0D0) WRITE(IOUT1,414) I,X(I)
                 IF(XL(I).EQ.15.0D0) WRITE(IOUT1,415) I,X(I)
                 IF(XL(I).EQ.16.0D0) WRITE(IOUT1,416) I,X(I)
                 IF(XL(I).EQ.17.0D0) WRITE(IOUT1,417) I,X(I)
                 IF(XL(I).EQ.18.0D0) WRITE(IOUT1,418) I,X(I)
                 IF(XL(I).EQ.19.0D0) WRITE(IOUT1,419) I,X(I)
                 IF(XL(I).EQ.20.0D0) WRITE(IOUT1,420) I,X(I)    
                 IF(XL(I).EQ.21.0D0) WRITE(IOUT1,421) I,X(I)
                 IF(XL(I).EQ.22.0D0) WRITE(IOUT1,422) I,X(I) 
                 IF(XL(I).EQ.90.0D0) WRITE(IOUT1,490) I,X(I)
                 IF(XL(I).EQ.91.0D0) WRITE(IOUT1,491) I,X(I)
                 IF(XL(I).EQ.92.0D0) WRITE(IOUT1,492) I,X(I)
                 IF(XL(I).EQ.93.0D0) WRITE(IOUT1,493) I,X(I)           
              ENDDO
              WRITE(IOUT1,47)           
          ENDIF
      ENDIF 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RETURN
      END
