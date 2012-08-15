CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                               
C                         
C      _|      _|  _|_|_|  _|_|_|      _|_|      _|_|_|    _|_|    
C      _|_|  _|_|    _|    _|    _|  _|    _|  _|        _|    _|  
C      _|  _|  _|    _|    _|    _|  _|_|_|_|  _|        _|    _|  
C      _|      _|    _|    _|    _|  _|    _|  _|        _|    _|  
C      _|      _|  _|_|_|  _|_|_|    _|    _|    _|_|_|    _|_|  
C
C                                                    Version 3.0 
C
C                                                           
C     MIDACO - Mixed Integer Distributed Ant Colony Optimization
C     ----------------------------------------------------------
C
C     This subroutine solves the general Mixed Integer Non-Linear Program (MINLP):
C
C             Minimize     F(X)            where X(1, ..., N-NINT) is *CONTINUOUS*
C                                          and   X(N-NINT+1,...,N) is *DISCRETE* 
C
C             Subject to:  G_j(X)  =  0    ( j = 1,...,ME )
C                          G_j(X) >=  0    ( j = ME + 1,...,M )
C
C             And bounds:  XL <= X <= XU  
C
C     MIDACO is a global optimization solver that stochastically approximates a solution to 
C     the above MINLP. It is based on an extended Ant Colony Optimization framework (see [1]) 
C     and the Oracle Penalty Method (see [2]) for constraint handling. MIDACO is called via
C     reverse communication, see below for a pseudo code of the reverse communication loop:
C
C     ~~~ while( STOP = 0 ) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C              call objective function F(X) and constaints G(X)
C
C              call midaco( X, F(X), G(X) )
C
C     ~~~ end ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     In case of mixed integer problems, the continuous variables are stored first in 'X', while 
C     the discrete (also called integer/categorical) variables are stored behind the continuous 
C     ones. As an example consider:
C
C     X = (0.153, 1.786, 1.0, 2.0, 3.0)   where 'N' = 5 and 'NINT' = 3 (number of integer variables)
C
C     Note that all 'X' is of type double precision. Equality and inequality constraints are handled  
C     in a similar way. The vector 'G' stores at first the 'ME' equality constraints and behind 
C     those, the remaining 'M-ME' inequality constraints are stored. 
C
C     MIDACO is a derivate free black box solver and does not require the relaxation of integer variables 
C     (this menas, integer variables are treated as categorical variables). MIDACO does not require any user 
C     specified parameter tuning as it can be run completely on 'Autopilot'. However, the user can optionally 
C     adjust MIDACO to his/her specific needs by some parameters explained below. 
C
C     MIDACO can process a user defined amount of iterates at once within one single reverse communication
C     step. Hence the call of the objective function and constraint functions can be parallelized outside 
C     and independently from MIDACO. Using this option is only recommended for cpu-time consuming problems
C     (that require for example more than 0.1 second to be calculated). However, for those problems a 
C     significant speedup can be gained by parallelization. The amount of parallel processed iterates is 
C     determined by the value 'L' in the MIDACO call. If 'L' > 1 the 'L' different iterates must be stored 
C     one after another in the array 'X' which is of length 'L*N'. Respectively the objective function values 
C     and constraint vectors must be stored one after another in 'F(L)' and 'G(L*M)'. The parallelization 
C     option can be used on various platforms and cpu architectures. Some templates for parallel-usage of 
C     MIDACO are available at: http://www.midaco-solver.com/parallel.html
C
C 
C     Usage:
C     ------
C
C             CALL MIDACO(L,N,NINT,M,ME,X,F,G,XL,XU,ACC,
C                         IFAIL,ISTOP,PARAM,RW,LRW,IW,LIW,
C                         LICENSE_KEY)
C
C
C     List of arguments:
C     ------------------
C
C     L :      (Parallelization Factor)
C               Number of parallel submitted iterates 'X' (with corresponding 
C               objective function values 'F' and constraint values 'G') to MIDACO
C               within one reverse communication step. If no parallelization is desired,
C               set L = 1.
C
C     N :       Number of optimization variables in total (continuous and integer ones). 
C               'N' is the dimension of the iterate 'X' with X = (X_1,...,X_N).
C
C     NINT :    Number of integer optimization variables. 'NINT' <= 'N'.
C               Integer (discrete) variables must be stored at the end of 'X'.
C     
C     M :       Number of constraints in total (equality and inequality ones).
C               'L*M' is the dimension of a constraint vector 'G' with G = (G_1,...,G_M).
C
C     ME :      Number of equality constraints. 'ME' <= 'M'.
C               Equality constraints are stored in the beginning of 'G'. 
C               Inequality constraints are stored in the end of 'G'.
c
C     X(L*N) :  Array containing the iterates 'X'. For L=1 only one iterate is
C               stored in 'X'. For L>1 the iterates must be stored one after
C               another. For example, let L=2 and A=(A_1,...,A_N) and B=(B_1,...,B_N)
C               be the current two iterates submitted to MIDACO, then 'A' and 'B'
C               are stored in 'X' like this: X = (A_1,...,A_N, B_1,...,B_N).
C
C     F(L) :    Array containing the objective function values 'F' corresponding
C               to the iterates 'X'. For L=1 only one objective function value is
C               stored in F(1). For L>1 the values must be stored one after
C               another like the iterates in 'X'. For example, let L=2 and 'FA' and
C               'FB' be the objective function valules corresponding to the iterates
C               A=(A_1,...,A_N) and B=(B_1,...,B_N) from above, then 'FA' and 'FB'
C               are stored in 'F' like this: F = (FA, FB).
C
C     G(L*M) :  Array containing the constraint values 'G'. For L=1 only one vector
C               of constraint values G = (G_1,...,G_M) is stored in 'G'. For L>1 the
C               vectors must stored one after another. For example, let L=2 and
C               'GA' and 'GB' be the two constraint value vectors corresponding to 
C               the iterates A=(A_1,...,A_N) and B=(B_1,...,B_N) from above, then 'GA' 
C               and 'GB' are stored in 'G' like this: G = (GA_1,...,G_AM, GB_1,...,GB_M).
C
C     XL(N) :   Array containing the lower bounds for the iterates 'X'.
C               Note that for integer dimesions the bound should also be discrete,
C               but submitted as double precision type, e.g. XL(N-NINT+1) = 1.0.
C               Note that the entries of XL are *changed*, when MIDACO finishes (ISTOP=1).
C
C     XU(N) :   Array containing the upper bounds for the iterates 'X'. 
C               Note that for integer dimesions the bound should also be discrete,
C               but submitted as double precision type, e.g. XU(N-NINT+1) = 1.0.
C               Note that the entries of XU are *changed*, when MIDACO finishes (ISTOP=1).
C               
C     ACC :     Accuracy for the constraint violation (=Residual). An iterate is assumed feasible,
C               if the L-infinity norm (maximum violation) over 'G' is lower or equal to 'ACC'. 
C               Hint: For a first optimization run the 'ACC' accuracy should be selected not too small.
C               A value of 0.01 to 0.001 is recommended. If a higher accuracy is demanded,
C               a refinement of the solution regarding the constraint violations is recommended
C               by applying another MIDACO run using the 'QSTART' option given by 'PARAM(2)'. 
C
C    IFAIL :    Communication flag used by MIDACO. Initially MIDACO must be called with IFAIL=0.
C               If MIDACO works correctly, IFAIL flags lower than 0 are used for internal communication.
C               If MIDACO stops (either by submitting ISTOP=1 or automatically using AUTOSTOP), an IFAIL
C               FLAG between 1 and 9 is returned as final message. If MIDACO detects some critical
C               problem setup, a *WARNING* message is returned by an IFAIL flag between 10 and 99. If
C               MIDACO detects an *ERROR* in the problem setup, an IFAIL flag between 100 and 999 is
C               returned and MIDACO stops. The individual IFAIL flags are as follows:
C
C               STOP - Flags:
C               -------------
C               IFAIL = 1 : Feasible solution,   MIDACO was stopped by the user submitting ISTOP=1
C               IFAIL = 2 : Infeasible solution, MIDACO was stopped by the user submitting ISTOP=1
C               IFAIL = 3 : Feasible solution,   MIDACO stopped automatically using 'AUTOSTOP' = PARAM(3)
C               IFAIL = 4 : Infeasible solution, MIDACO stopped automatically using 'AUTOSTOP' = PARAM(3)
C       
C               WARNING - Flags:
C               ----------------
C               IFAIL = 51 : Some X(i)  is greater/lower than +/- 1.0D+12 (try to avoid huge values!)
C               IFAIL = 52 : Some XL(i) is greater/lower than +/- 1.0D+12 (try to avoid huge values!)
C               IFAIL = 53 : Some XU(i) is greater/lower than +/- 1.0D+12 (try to avoid huge values!)
C
C               IFAIL = 61 : Some X(i)  should be discrete (e.g. 1.000) , but is continuous (e.g. 1.234)
C               IFAIL = 62 : Some XL(i) should be discrete (e.g. 1.000) , but is continuous (e.g. 1.234)
C               IFAIL = 63 : Some XU(i) should be discrete (e.g. 1.000) , but is continuous (e.g. 1.234)
C
C               IFAIL = 71 : Some XL(i) = XU(I) (fixed variable)
C
C               IFAIL = 81 : F(X) has value NaN for starting point X (sure your problem is correct?)
C               IFAIL = 82 : Some G(X) has value NaN for starting point X (sure your problem is correct?)
C
C               ERROR - Flags:
C               --------------
C               IFAIL = 101 :   L    <= 0
C               IFAIL = 102 :   N    <= 0
C               IFAIL = 103 :   NINT <  0
C               IFAIL = 104 :   NINT >  N
C               IFAIL = 105 :   M    <  0
C               IFAIL = 106 :   ME   <  0
C               IFAIL = 107 :   ME   >  M
C
C               IFAIL = 201 :   some X(i)  has type NaN
C               IFAIL = 202 :   some XL(i) has type NaN
C               IFAIL = 203 :   some XU(i) has type NaN
C               IFAIL = 204 :   some X(i) < XL(i)
C               IFAIL = 205 :   some X(i) > XU(i)
C               IFAIL = 206 :   some XL(i) > XU(i)
C           
C               IFAIL = 301 :   ACC < 0 or ACC has type NaN
C               IFAIL = 302 :   PARAM(1) < 0
C               IFAIL = 303 :   PARAM(2) < 0 or ( 0 < PARAM(2) < 1 )
C               IFAIL = 304 :   PARAM(3) < 0
C               IFAIL = 305 :   PARAM(5) < 0
C               IFAIL = 306 :   PARAM(6) < 0 or ( 0 < PARAM(6) < 2 )
C               IFAIL = 307 :   PARAM(6) >= PARAM(5)
C               IFAIL = 308 :   PARAM(5) > 0 and PARAM(6) = 0
C               IFAIL = 309 :   PARAM(6) > 2*N+10
C               IFAIL = 310 :   PARAM(7) < 0 or PARAM(7) > 3
C               IFAIL = 311 :   some PARAM(i) has type NaN
C
C               IFAIL = 401 :   ISTOP < 0 or ISTOP > 1
C               IFAIL = 501 :   Double precision work space size LRW is too small (see below LRW)
C                               ---> RW must be at least of size LRW = 2*N^2+23*N+2*M+70
C               IFAIL = 601 :   Integer work space size LIW is too small (see below LIW)
C                               ---> IW must be at least of size LIW = 2*N+L+100
C               IFAIL = 701 :   Input check failed! MIDACO must be called initially with IFAIL = 0
C               IFAIL = 801 :   L > LMAX (user must specifiy LMAX below in the MIDACO source code) 
C               IFAIL = 802 :   L*M+1 > LXM (user must specifiy LXM below in the MIDACO source code)
C
C               IFAIL = 900 :   Invalid or corrupted LICENSE_KEY
C
C               IFAIL = 999 :   N > 4. The free test version is limited up to 4 variables. 
C                               To get an unlimited version, please contact info@midaco-solver.com.
C
C    ISTOP :    Communication flag to stop MIDACO by the user. If MIDACO is called with ISTOP = 1, MIDACO 
C               returns the best found solution in 'X' with corresponding 'F' and 'G'. As long as MIDACO 
C               should continue its search, ISTOP must be equal to 0.
C
C    PARAM() :  Array containing 7 parameters that can be selected by the user to adjust MIDACO. However, 
C               setting those parameters is *ONLY RECOMMENDED FOR ADVANCED USERS*. Unexperienced users 
C               should set all PARAM(i) = 0. The parameters are as follows:
C
C               PARAM(1) :   [RANDOM-SEED] This value indicates the random seed used within MIDACO's
C                            internal pseudo-random number generator. For each seed, a different sequence
C                            of pseudo-random numbers in generated, influencing the MIDACO results. The seed
C                            must be a (discrete) value >= 0 given in double precision. For example 
C                            PARAM(1) = 0.0, 1.0, 2.0,... (Note that MIDACO runs are 100% reproducable 
C                            for the same random seed used).
C
C               PARAM(2) :   [QSTART] This value indicates the quality of the starting point 'X' submitted by
C                            the user at the first call of MIDACO. It must be a (discrete) value >= 0
C                            given in double precision, where PARAM(2) = 0.0 is the default setting assuming 
C                            that the starting point is just some random point without specific quality. A higher 
C                            value indicates a higher quality of the solution. If PARAM(2) >= 1.0 is selected, 
C                            MIDACO will concentrate its search around the startingpoint by sampling its initial 
C                            population in the area (XU(i)-XL(i))/PARAM(2) in all continuous dimensions 'i' respectively.
C                            For discrete dimensions 'i', the initial sampling is performed with MAX((XU(i)-XL(i))/PARAM(2),
C                            1/SQRT(PARAM(2))) around the starting point.
C                            Note that activating PARAM(2) will *NOT* shrink the search space defined by 'XL' and 
C                            'XU', but only concentrates the population around the startingpoint and later the
C                            current best solution found. The QSTART option is very useful to refine a previously 
C                            calculated solution and/or to increase the accuracy of the constraint violation. For 
C                            continuous large scale problems (N > 100), the QSTART option is also very helpful to 
C                            refine the solution precision.
C                            
C               PARAM(3) :   [AUTOSTOP] This value enables an automatic stopping criteria within MIDACO. It must
C                            be a (discrete) value >= 0 given in double precision. For example PARAM(3) = 1.0, 
C                            2.0, 3.0,... If PARAM(3) is >= 1.0 MIDACO will stop and return its current
C                            best solution after PARAM(3) internal restarts without significant improvement of the 
C                            current solution. Hence, a small value for PARAM(3) will lead to a shorter runtime of
C                            MIDACO with a lower chance of reaching the global optimal solution. A large value of
C                            PARAM(3) will lead to a longer runtime of MIDACO but increases the chances of reaching
C                            the global optimum. Note that running MIDACO with PARAM(3) = 0.0 with a maximal 
C                            available time budget (e.g. 1 Day = 60*60*24 sec) will always provide the highest 
C                            chance of global optimality.
C
C               PARAM(4) :   [ORACLE] This parameter affects only constrained problems. If PARAM(4)=0.0 is submitted 
C                            MIDACO will use its inbuild oracle strategy. If PARAM(4) is not equal to 0.0, MIDACO 
C                            will use PARAM(4) as initial oracle for its oracle penalty function and only update the
C                            oracle, if a feasible solution with 'F(X)' < 'PARAM(4)' has been found. In case the user 
C                            wants to submit the specific ORACLE = 0.0, a close value like 1.0D-12 can be used as 
C                            dummy. Please review [2] to receive more information on the oracle penalty method. 
C
C               PARAM(5) :   [ANTS] This value fixes the number of iterates (ants) used within a generation. If
C                            PARAM(5)=0.0 is submitted, MIDACO will handle the number of iterates dynamically
C                            by itself. The use of PARAM(5) and PARAM(6) is usefull for either very cpu-time 
C                            expensive problems or large scale problems ('N' > 100). Please contact the author 
C                            directly to receive support in using this option.
C          
C               PARAM(6) :   [KERNEL] This value fixes the kernel size for each generation. If PARAM(6)=0.0 is 
C                            submitted, MIDACO will handle the kernel sizes dynamically by itself. IF PARAM(6) is
C                            not equal 0.0, it must be at least 2.0. 
C
C               PARAM(7) :   [CHARACTER] MIDACO includes three different parameter settings especially tuned for
C                            IP/NLP/MINLP problems. If PARAM(7) = 0.0 MIDACO will select its parameter set 
C                            according to the balance of 'N' and 'NINT'. If the user wishes to enable a specific
C                            set (for example the NLP set if 'N'=98 and 'NINT'=2) he can do so by PARAM(7):
C
C                            PARAM(7) = 1.0 enables the internal parameter set tuned for IP problems ('N'='NINT')
C                            PARAM(7) = 2.0 enables the internal parameter set tuned for NLP problems ('NINT'=0)
C                            PARAM(7) = 3.0 enables the internal parameter set tuned for MINLP problems
C
C    RW(LRW) :  Real workarray (double precision) of length 'LRW'
C    LRW :      Length of 'RW'. 'LRW' must be greater or equal to  2*N^2+23*N+2*M+70
C
C    IW(LIW) :  Integer workarray (long integer) of length 'LIW'
C    LIW :      Length of 'IW'. 'LIW' must be greater or equal to  2*N+L+100
C
C    LICENSE_KEY :  Character string consisting of 60 ASCII letters. Please note that any licensed copy  
C                   of MIDACO comes with an individual 'LICENSE_KEY' determining the license owner and
C                   additional license conditions.
C
C
C    References:
C    -----------
C
C    [1] Schlueter, M., Egea, J. A., Banga, J. R.: 
C        "Extended ant colony optimization for non-convex mixed integer nonlinear programming", 
C        Computers & Operations Research, Vol. 36 , Issue 7, Page 2217-2229, 2009.
C
C    [2] Schlueter M., Gerdts M.: "The oracle penalty method",
C        Journal of Global Optimization, Vol. 47(2),pp 293-325, 2010.
C
C
C    Author (C) :   Martin Schlueter
C                   Theoretical & Computational Optimization Group,
C                   School of Mathematics, University of Birmingham,
C                   Watson Building, Birmingham B15 2TT (UK)
C
C    URL :          http://www.midaco-solver.com
C
C    Email :        info@midaco-solver.com
C       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
                       
      subroutine midaco(l,n,nint,m,me,x,f,g,xl,xu,acc,
     & 	                ifail,istop,param,rw,lrw,iw,liw,
     &                  license_key) 
     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
      implicit none                                
      integer l,n,nint,m,me,ifail,lrw,iw,liw,istop                    
      double precision x,f,g,xl,xu,acc,param(7),rw                     
      dimension x(l*n),f(l),g(l*m),xl(n),xu(n),rw(lrw),iw(liw)
      character*60 license_key
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     LMAX Defines the maximal amount of L parallel processed iterates.
C     In case the user wishes to use parallelization and wants to submit 
C     L > 1 iterates to MIDACO at once, LMAX must be set here manually 
C     with: LMAX >= L. Additionally the LXM parameter can be set here 
C     manually. LXM must be greater or equal L*M+1.
      integer LMAX,LXM,I
      parameter (LMAX =  100)
      parameter (LXM  = 1000)      
      double precision A(LMAX),B(LMAX),GM(LXM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(ifail.eq.0)then
c         Check L <= LMAX
          if(L.gt.LMAX)then
              ifail = 801
              return
          endif  
c         Check L*M+1 <= LXM
          if(L*M+1.gt.LXM)then
              ifail = 802
              return
          endif 
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(M.gt.0)then
          do I = 1,L*M
              GM(I) = g(I)
          enddo  
      endif      
      GM(L*M+1) = 0.0D0                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                       
      call midaco_code(l,n,nint,m,me,x,f,GM,xl,xu,acc,
     & 	               ifail,istop,param,rw,lrw,iw,liw,
     &                 license_key,LMAX,A,B)      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



         subroutine o8971310(n,k,j2,lj2,i37,i98,i2462,i13,x,f,i11,p,z)  
                                implicit none                           
                    integer n,k,lj2,i37,i98,i2462,i13,i4271,i,j         
             double precision x,f,i11,p,j2,z                         
                    dimension x(n),j2(lj2)                              
                                i4271 = 0                
                          IF(z.EQ.0.0d0)THEN                            
                                          if(p.ge.j2(i13+k-1))return    
                              ELSE                                      
                             if(p.gt.j2(i13+k-1))return                 
                              ENDIF                                     
               do i=1,k                                           
                           if(p.le.j2(i13+k-i))then                     
                  i4271=k-i+1                                           
                                                                else    
                                                     goto 567           
                                                     endif         
                   enddo                                                
  567                                                 do j=1,k-i4271    
                                                     do i=1,n           
               j2(i37+(k-j)*n+i-1) = j2(i37+(k-j-1)*n+i-1)              
                                enddo                                   
               j2(i98+k-j)     = j2(i98+k-j-1)                          
                             j2(i2462+k-j)   = j2(i2462+k-j-1)          
                          j2(i13+k-j)     = j2(i13+k-j-1)              
                                          enddo                         
           do i=1,n                                                     
      j2(i37+(i4271-1)*n+i-1)  = x(i)                          
                                enddo                                   
                                      j2(i98+i4271-1)          = f      
                 j2(i2462+i4271-1)        = i11                         
                               j2(i13+i4271-1)          = p             
                                        end                             
                                           subroutine o9953052(f,g,m)   
                                                    implicit none       
                                                   integer m,i          
                                double precision f,g                    
                                                        dimension g(m)  
                                if(f.ne.f)then                          
                                    f = 1.0D16                          
                                endif                                   
                         do i=1,m                                       
                                                if(g(i).ne.g(i))then    
                                                      g(i) = - 1.0D16   
                endif                                                   
                                           enddo              
                                                         end            
                   subroutine o329519(j6,lj6,i64,i03s,bi03,di03,di64)   
        implicit none                                                   
      integer j6,lj6,i64,i03s,bi03,di03,di64                            
                                                 dimension j6(lj6)      
          if(j6(i64).eq.1.and.j6(di64).eq.1)then                        
                        j6(i03s) = j6(di03)                             
                          else                                          
                        j6(i03s) = j6(bi03)                     
                               endif                                    
                        if(j6(i64).le.j6(di64).and.j6(di64).gt.1)then   
                       j6(i03s) = j6(bi03) + (j6(di03)-j6(bi03)) *      
     &    int( dble((j6(i64)-1)) / dble((j6(di64)-1)) )                 
                               endif                                    
       if(j6(i64).gt.j6(di64).and.j6(i64).lt.2*j6(di64))then         
                  j6(i03s) = 2 * ( j6(di03) + (j6(bi03)-j6(di03)) *     
     &                    int( dble(j6(i64)) / dble(2*j6(di64)) ) )     
                                                   endif                
                                            end                         
                   subroutine i074206156(l,n,nint,m,me,x,f,g,j7,j9,j13, 
     &                         i3108,i4108,i7193,j2,lj2,j6,lj6,        
     &                i813,i0009,i5308,K1,K2,                           
     &                                    i31032,i31579907326)     
                                                   implicit none        
                          integer l,n,nint,m,me,i3108,lj2,j6,lj6,i4108  
       double precision x,f,g,j7,j9,j13,i7193,j2                        
              dimension x( l*n ),f( l ),g( l*m+1 ),j7(n),j9(n)          
                                     dimension i7193(7),j2(lj2),j6(lj6) 
                         integer i813,i0009,i5308                       
                  double precision K1,K2                                
                            character*60 i31032                    
           integer i0621,i7904,i429103,i47029,i007002,i31579907326     
                                                               integer i
                   call o027659163250219005(i31032,i0621,i3108)    
                                               if(l.le.0)then           
                              i3108 = 101                               
                                                         goto 701       
                              endif                                     
                                                 if(n.le.0)then         
                       i3108 = 102                                      
                                                  goto 701              
                          endif                                         
                                                     if(nint.lt.0)then  
                                                   i3108 = 103          
                                           goto 701                     
                     endif                                            
                                    if(nint.gt.n)then                   
                i3108 = 104                                             
                               goto 701                          
                     endif                                              
                      if(m.lt.0)then                                    
                                           i3108 = 105                  
                goto 701                                                
                                   endif                                
                 if(me.lt.0)then                                        
                                                     i3108 = 106        
                              goto 701                                  
                   endif                                                
        if(me.gt.m)then                                                 
           i3108 = 107                                                  
                   goto 701                                
                   endif                                                
    1                          do i=1,n                                 
                                                    if(x(i).ne.x(i))then
                                                     i3108 = 201        
                                  goto 701                              
                     endif                                              
                                     if(j7(i).ne.j7(i))then             
                            i3108 = 202                                 
                                goto 701                                
               endif                                             
                            if(j9(i).ne.j9(i))then                
           i3108 = 203                                                  
                                              goto 701                  
                   endif                                        
       if(x(i).lt.j7(i)-1.0D-6)then                                     
                                                         i3108 = 204    
        goto 701                                                        
                                                     endif              
                           if(x(i).gt.j9(i)+1.0D-6)then                 
                           i3108 = 205                            
                                              goto 701                  
           endif                                                        
                                         if(j7(i).gt.j9(i)+1.0D-6)then  
                                                i3108 = 206             
                 goto 701                                               
                                        endif                           
       enddo                                                            
                call o3261853008(j6(i0621),i31032)                 
                                    if(j13.lt.0.0d0.or.j13.ne.j13)then  
               i3108 = 301                                              
                                                  goto 701              
                                 endif                                  
                                           if(i7193(1).lt.0.0d0)then    
         i3108 = 302                                                    
                                 goto 701                         
                                                        endif           
                          if(i7193(2).lt.0.0d0.or.                      
     &   (i7193(2).gt.0.0d0.and.i7193(2).lt.1.0d0))then                 
                                           i3108 = 303             
                                      goto 701                          
                                                endif                   
                  if(i7193(3).lt.0.0d0)then                      
       i3108 = 304                                                      
                                                    goto 701            
                      endif                                             
                                if(i7193(5).lt.0.0d0)then               
         i3108 = 305                                                    
                                                              goto 701  
                                                                endif   
           if(i7193(6).lt.0.0d0.or.                                     
     &         (i7193(6).gt.0.0d0.and.i7193(6).lt.2.0d0))then           
       i3108 = 306                                                   
                                                            goto 701    
                                                               endif    
         if(i7193(6).ge.i7193(5).and.i7193(6).gt.0.0d0)then             
                                                      i3108 = 307       
       goto 701                                                         
                                                endif                   
                     if(i7193(5).gt.0.0d0.and.i7193(6).eq.0.0d0)then    
                i3108 = 308                                             
                                        goto 701                        
              endif                                                     
                 if(i7193(6).gt.2*dble(n)+10.0d0)then                 
       i3108 = 309                                                      
                    goto 701                                            
                                                        endif           
        if(i7193(7).lt.0.0d0.or.i7193(7).gt.3.0d0)then                  
                        i3108 = 310                                     
                                                              goto 701  
                                   endif                                
                do i=1,6                                                
                if(i7193(i).ne.i7193(i))then                            
           i3108 = 311                              
                     goto 701                                           
                                                        endif           
                                               enddo                    
                     if(i4108.lt.0.or.i4108.gt.1)then                   
                                                   i3108 = 401          
                                             goto 701                   
                          endif                                         
                                        i7904 = 0                      
                                                 do i=1,20              
                         i7904 = i7904 + j6(i0621+i-1)                
                                           enddo                        
                            i429103 =  551     
                         if(i7904.ne.i429103)then                      
                                                   goto 1               
                    endif                                               
                  call o3261853008(j6(1),i31032)                   
                i47029 = 1                                              
                                                    do i=1,60           
                              i47029 = i47029 + j6(i)                   
                                                          enddo         
                       i007002 =  2289         
                    if(i47029.ne.i007002)then                           
           goto 1                                                       
                  endif                                      
                                 i813 = INT(K1) * N + INT(K2)           
                                        i0009 = 2*n+m+(n+5)*i813+8      
                                            i5308 = 31+l+n+n            
                           if(lj2.lt.i0009+5+n+m)then                   
                                            i3108 = 501                 
                                                 goto 701               
                                                   endif                
      if(lj6.lt.i5308+69)then                                           
                                             i3108 = 601                
                       goto 701                                         
           endif                                       
             do i=1,i0009+5+n+m                                         
                                  j2(i) = 0.0d0                         
                                                               enddo    
                                                    do i=1,i5308+69     
                                                    j6(i) = 0           
                        enddo                                           
                                                        i31579907326 = 0
                                                  do i=1,n              
                   if(x(i).gt.1.0D+12.OR.x(i).lt.-1.0D+12)then        
                                  i3108 = 51                            
                                         goto 702                       
                                                                  endif 
                        if(j7(i).gt.1.0D+12.OR.j7(i).lt.-1.0D+12)then   
                                             i3108 = 52                 
                                                  goto 702              
                            endif                                       
                 if(j9(i).gt.1.0D+12.OR.j9(i).lt.-1.0D+12)then          
               i3108 = 53                                               
                       goto 702                        
          endif                                          
                                         if(j7(i).eq.j9(i))then         
                                                    i3108 = 71          
                                                     goto 702           
         endif                                                          
                                           enddo                     
                   do i = n-nint+1,n                                    
                            if(dabs(x(i)-dnint(x(i))).gt.1.0d-6)then    
                 i3108 = 61                                             
                      goto 702                                          
                                                      endif             
                   if(dabs(j7(i)-dnint(j7(i))).gt.1.0d-6)then           
         i3108 = 62                                                     
                                           goto 702                     
                               endif                                    
                              if(dabs(j9(i)-dnint(j9(i))).gt.1.0d-6)then
        i3108 = 63                                                      
        goto 702                                          
                           endif                                        
                                                    enddo               
                                 if(f(1).ne.f(1))then                   
                         i3108 = 81                                     
                     goto 702                                           
          endif                                   
                                       do i = 1,m                       
                         if(g(i).ne.g(i))then                           
          i3108 = 82                                                 
            goto 702                                                    
                                                            endif       
                                               enddo                    
                                             return                     
  701                                         continue                  
                                    i4108 = 1                           
                         return                                         
  702                                                    continue       
                                  i31579907326 = 1                      
                                 return                                 
                                           end                    
          subroutine o61867290(n,nint,j2,lj2,i8087,i37,j6,lj6,k,i64,z)  
                          implicit none                                 
                       integer n,nint,j6,i64,k,lj2,lj6,i8087,i37,i,j    
        double precision j2,i809,i140,i346031,i6132006,i5510871,z       
                                              dimension j2(lj2),j6(lj6) 
                                        i346031  = sqrt(dble(j6(i64)))  
                                   i6132006 = z/i346031                 
         i5510871 = (1.0d0-1.0d0/dsqrt(dble(nint)+0.1D0)) / 2.0d0       
        do i=1,n                                                      
                                                     i809 = j2(i37+i-1) 
                                 i140 = j2(i37+i-1)                     
                               do j=2,j6(k)                    
                            if(j2(i37+(j-1)*n+i-1).gt.i809)then         
        i809 = j2(i37+(j-1)*n+i-1)                                      
                               endif                            
                         if(j2(i37+(j-1)*n+i-1).lt.i140)then            
                       i140 = j2(i37+(j-1)*n+i-1)                       
                                     endif                              
                                               enddo                    
                j2(i8087+i-1) = (i809-i140)/i346031             
                                             if(i.gt.n-nint)then        
             if(j2(i8087+i-1).lt.i6132006)then                          
                         j2(i8087+i-1) = i6132006                       
            endif                                                       
         if(j2(i8087+i-1).lt.i5510871)then                              
                 j2(i8087+i-1) = i5510871              
                      endif                                             
                                               endif                    
                 enddo                                               
                                                            end         
          subroutine o83517(n,nint,j2,lj2,j6,lj6,i86x,x,j7,j9,          
     &                                             z,y,i13509)         
                                                        implicit none   
                                 integer n,nint,lj2,i86x,j6,lj6,i       
        double precision j2,x,j7,j9,i13509,z,y,i130,o89,i9042677836     
                         dimension j2(lj2),j6(lj6),x(n),j7(n),j9(n)     
                                do i=1,n                                
                   i130 = (j9(i)-j7(i)) / dble(j6(10))                  
                             if(i.gt.n-nint.and.i130.lt.z) i130 = z     
                 if(i13509.gt.0.0d0)then                                
                               if(i130.gt.(j9(i)-j7(i))/i13509)then     
              i130 = (j9(i)-j7(i))/i13509                               
                                         endif                          
                                               if(i.gt.n-nint)then      
                               if(i130.lt.1.0d0/dsqrt(i13509))then      
                              i130 = 1.0d0 / dsqrt(i13509)              
          endif                                                         
                                            endif              
                             endif                                     
                                x(i) = j2(i86x+i-1) + i130 *            
     &             i9042677836(o89(j2(1)),o89(j2(1)))                   
                      if(x(i).lt.j7(i))then                             
                                     x(i)=j7(i)+(j7(i)-x(i))/1.0d1**y   
            if(x(i).lt.j7(i)) x(i) = j7(i)                              
                         if(x(i).gt.j9(i)) x(i) = j9(i)               
                   endif                                                
         if(x(i).gt.j9(i))then                               
                      x(i)=j9(i)-(x(i)-j9(i))/1.0d1**y            
                            if(x(i).lt.j7(i)) x(i) = j7(i)              
                                  if(x(i).gt.j9(i)) x(i) = j9(i)        
                                             endif                      
                  if(i.gt.n-nint) x(i) = dnint(x(i))                    
                                                             enddo      
                             end                                        
                      subroutine o21(i11,g,m,me,j13)                    
                    implicit none                                       
                       integer m,me,i                                   
                   double precision g,i11,j13                           
                                      dimension g(m)                    
                                       i11 = 0.0d0                      
                     do i=1,me                                          
                                             if (dabs(g(i)).gt.j13) then
                                            if (dabs(g(i)).gt.i11) then 
                                           i11 = dabs(g(i))             
                                                            endif       
              endif                                                     
                                             enddo                      
                          do i=me+1,m                          
                                        if (g(i).lt.-j13) then          
                       if (-g(i).gt.i11) then                           
                                        i11 = - g(i)                    
                             endif                                  
                                                endif                   
                               enddo                              
        end                                                             
      subroutine o2840(n,m,i86x,i472,i315,i86i11,i4296,i813,            
     /                      i37,i98,i2462,i13,w,i8087,pt,i3,i4,i110)  
                          implicit none                                 
                                                           integer n,m  
                          integer i86x,i472,i315,i86i11                 
                                         integer i4296                  
                          integer i813                                  
                                           integer i37,i98,i2462,i13    
                                                 integer w,i8087        
                                                      integer pt        
                                               integer i3,i4,i110       
                                    i86x        = 2                     
           i472        = i86x        + n                                
         i315        = i472        + 1                                  
                                       i86i11      = i315        + m    
                i4296       = i86i11      + 1                           
                          i37         = i4296       + 1                 
                       i98         = i37         + n * i813             
                             i2462       = i98         + i813           
                                  i13         = i2462       + i813      
                       pt          = i13         + i813                 
                 w           = pt          + i813 + 1                
                                   i8087       = w           + i813     
                           i3          = i8087       + n                
                                        i4          = i3          + 1   
              i110        = i4          + 1                             
                                                                   end  
           subroutine o027659163250219005(i31032,i0621,i3108)      
       implicit none                                                    
            character*60 i31032                                    
        integer i0621,i3108                                             
                                     if(i3108.eq.0)then                 
                                                               i0621 = 1
                                               endif                    
                                if(i0621.le.0) i0621 = 1                
       call o6288601(i31032(i0621:i0621),i0621)                    
                              if(i0621.le.0)  i0621 = 1                 
                              if(i0621.ge.60) i0621 = 2                 
   50                                                 return            
                                                                  end   
             subroutine o953814(n,i813,j2,lj2,i472,i86i11,i4296,j13,    
     /                 i37,i98,i2462,i13,j6,lj6,i81,k,bi03,di03,di64,  
     /      i87621,i41981820,i78923,i78913,i7193,i5308,                 
     /                                       i31032)               
                          implicit none                                 
       integer n,i813,j6,lj2,lj6,i81,k,bi03,di03,di64,i37,i98,i2462,    
     &      i13,i41981820,i472,i86i11,i4296,i,j,P1,P2,P3,i662,i55117,   
     &                         i46i022,i5308                            
          double precision j2,j13,i87621,i7193(7),i78923(11),i78913(18),
     &                                            F1,Q1,Q2,Z1,Z2     
                    character*60 i31032                            
                                         dimension j2(lj2),j6(lj6)      
              F1 = i78923(2)                                            
                                            Q1 = i78923(3)              
                                                 Q2 = i78923(4)         
                                                       Z1 = i78923(9)   
                                                Z2 = i78923(10)         
            P1 = INT(i78913(10))                                        
                       P2 = INT(i78913(11))                             
                               P3 = INT(i78913(12))                     
                                    i662 = 0                            
                        i55117 = 0                                      
                                     i46i022 = 0                        
                            call o3261853008(j6(i5308+1),i31032)   
                 if(j6(i81).le.1)then                                   
                                                j6(10) = 0             
                                                            j6(12) = 0  
                        else                                            
                                                j6(12) = j6(12) + P1    
                                       j6(10) = 1                       
                           do i=1,j6(12)                                
                   j6(10) = j6(10) * P2                                 
                           enddo                                        
                                                                    
            if(j6(10).gt.1.0D1**P3)then                                 
                                               j6(10) = 0               
                                              j6(12) = 0           
                                                            endif       
                           endif                                        
                                  if(int(i7193(5)).gt.0) goto 101       
                                           if(j6(i81).le.1) j6(9) = 0   
  567                         j6(bi03) = INT( Z1 * dble(n)/2.0d0 +      
     &                                 Z2 * j6(9) * dble(n)/2.0d0 )     
                                 if(j6(bi03).lt.3) j6(bi03) = 3         
                                                       i55117 = 0       
                                                         do i=1,5       
        i55117 = i55117 + j6(i5308+i)                                   
                      enddo                                       
                                                                   
                                    j6(k) = int(F1 * DBLE(j6(bi03)))    
         if(j6(k).lt.2) j6(k) = 2                                    
                                  i662 = i662 + 1                       
                                                IF(i662.GT.100)then     
                j6(bi03) = 3                                            
                    j6(k)    = 2                                        
                              goto 100                                  
                      endif                                             
                 if(j6(k).gt.i813)then                                  
                                   j6(9) = 1                            
                                                   goto 567             
                         else                                           
                                           j6(9) = j6(9) + 1            
                      endif                                             
                                            do i=5,15                   
         i55117 = i55117 - j6(i5308+i)                                  
           enddo                                                        
         i46i022 = -253                       
  100                  j6(di03) = i41981820 * INT( DBLE(j6(bi03))** Q1 )
            j6(di64) = INT( Q2 * DBLE(j6(k)) )                          
  101       if(int(i7193(5)).gt.0)then                                  
            j6(bi03) = int(i7193(5))                                    
                            j6(k)    = int(i7193(6))               
                                 j6(di03) = j6(bi03)                    
          j6(di64) = j6(k)                                              
                                     endif                        
         if(i55117.ne.i46i022)then                                      
                j6(bi03) = int(i7193(5))                                
            j6(k)    = int(i7193(4))                                    
                                             j6(di03) = j6(k)           
           j6(di64) = j6(bi03)                                          
                                                                  endif 
                                       j2(i4296) = i87621               
      if(j2(i86i11).le.j13.and.j2(i472).lt.i87621) j2(i4296) = j2(i472) 
                  do j=1,j6(k)                                          
                              do i=1,n                                  
                     j2(i37+(j-1)*n+i-1) = 1.0D16                       
                                                              enddo     
                  j2(i2462+j-1)       = 1.0D16                          
                                      j2(i98+j-1)         = 1.0D16      
                           j2(i13+j-1)         = 1.0D16                 
                                                 enddo                  
                                                    end                 
      subroutine o732(p,f,i11,i4296,j13)                                
                       implicit none                         
        double precision p,f,i11,i4296,j13,i1437                        
                      if (f.le.i4296.and.i11.le.j13) then               
               p = f - i4296                                            
                             return                                 
                else                                                    
                        if (f.le.i4296) then                            
                                                      p = i11           
                       return                                           
                           else                                         
                                           i1437 = 0.0d0                
                                if (i11.lt.(f-i4296)/3.0d0) then        
          i1437 = ((f-i4296)*8.07549910d0-i11) / (f-i4296-i11)          
                                           goto 1                       
                      endif                                             
          if(i11.ge.(f-i4296)/3.0d0.and.i11.le.(f-i4296)) then          
                      i1437 = 1.0d0 - 0.5d0 / dsqrt((f-i4296)/i11)      
                                                               goto 1   
                                             endif                      
           if  (i11.gt.(f-i4296)) then                                  
                         i1437 = 0.5d0 * dsqrt((f-i4296)/i11)           
                         goto 1                                         
                                          endif                         
    1         p = i1437 * (f-i4296) + (1.0d0-i1437) * i11               
                                                                 return 
                                                                  endif 
                                  endif                                 
                                                                  end   
                            subroutine o6288601(i50315,i82431)          
                        implicit none                                
                                 character*1 i50315                     
         integer i82431                                          
                                                 i82431 = 0             
                            if(i50315(1:1).eq.'A') i82431 = 52          
                  if(i50315(1:1).eq.'B') i82431 = 28                    
        if(i50315(1:1).eq.'C') i82431 = 49                              
             if(i50315(1:1).eq.'D') i82431 = 30                         
                                   if(i50315(1:1).eq.'E') i82431 = 31   
                                     if(i50315(1:1).eq.'F') i82431 = 32 
        if(i50315(1:1).eq.'G') i82431 = 33                              
                          if(i50315(1:1).eq.'H') i82431 = 34            
                          if(i50315(1:1).eq.'I') i82431 = 35            
                  if(i50315(1:1).eq.'J') i82431 = 36                    
                                 if(i50315(1:1).eq.'K') i82431 = 37     
             if(i50315(1:1).eq.'L') i82431 = 38                         
                              if(i50315(1:1).eq.'M') i82431 = 39        
                if(i50315(1:1).eq.'N') i82431 = 40                      
                         if(i50315(1:1).eq.'O') i82431 = 41             
                                     if(i50315(1:1).eq.'P') i82431 = 42 
                                 if(i50315(1:1).eq.'Q') i82431 = 43     
                       if(i50315(1:1).eq.'R') i82431 = 44               
                          if(i50315(1:1).eq.'S') i82431 = 45        
                                 if(i50315(1:1).eq.'T') i82431 = 46     
       if(i50315(1:1).eq.'U') i82431 = 47                               
               if(i50315(1:1).eq.'V') i82431 = 48                       
            if(i50315(1:1).eq.'W') i82431 = 29                          
           if(i50315(1:1).eq.'X') i82431 = 50                        
                            if(i50315(1:1).eq.'Y') i82431 = 51          
                                    if(i50315(1:1).eq.'Z') i82431 = 27  
                                if(i50315(1:1).eq.'0') i82431 = 53      
                           if(i50315(1:1).eq.'1') i82431 = 54           
                       if(i50315(1:1).eq.'2') i82431 = 55               
                   if(i50315(1:1).eq.'3') i82431 = 56                   
           if(i50315(1:1).eq.'4') i82431 = 57                           
                             if(i50315(1:1).eq.'5') i82431 = 58         
                           if(i50315(1:1).eq.'6') i82431 = 59           
                                  if(i50315(1:1).eq.'7') i82431 = 60    
                 if(i50315(1:1).eq.'8') i82431 = 61                     
                          if(i50315(1:1).eq.'9') i82431 = 62            
                        if(i50315(1:1).eq.'a') i82431 = 23              
                               if(i50315(1:1).eq.'b') i82431 = 2        
                          if(i50315(1:1).eq.'c') i82431 = 3             
       if(i50315(1:1).eq.'d') i82431 = 16                            
              if(i50315(1:1).eq.'e') i82431 = 5                         
                          if(i50315(1:1).eq.'f') i82431 = 13            
             if(i50315(1:1).eq.'g') i82431 = 7                          
                  if(i50315(1:1).eq.'h') i82431 = 8                     
        if(i50315(1:1).eq.'i') i82431 = 9                               
                                   if(i50315(1:1).eq.'j') i82431 = 10   
                     if(i50315(1:1).eq.'k') i82431 = 11                 
                 if(i50315(1:1).eq.'l') i82431 = 12                   
       if(i50315(1:1).eq.'m') i82431 = 6                                
                                if(i50315(1:1).eq.'n') i82431 = 14      
                                if(i50315(1:1).eq.'o') i82431 = 15      
         if(i50315(1:1).eq.'p') i82431 = 4                              
                           if(i50315(1:1).eq.'q') i82431 = 17           
                           if(i50315(1:1).eq.'r') i82431 = 18           
                     if(i50315(1:1).eq.'s') i82431 = 19                 
              if(i50315(1:1).eq.'t') i82431 = 20                        
         if(i50315(1:1).eq.'u') i82431 = 21                             
                            if(i50315(1:1).eq.'v') i82431 = 22          
                                if(i50315(1:1).eq.'w') i82431 = 1       
                    if(i50315(1:1).eq.'x') i82431 = 24                  
                       if(i50315(1:1).eq.'y') i82431 = 25            
               if(i50315(1:1).eq.'z') i82431 = 26                       
            if(i50315(1:1).eq.'_') i82431 = 64                          
                          if(i50315(1:1).eq.'(') i82431 = 65            
            if(i50315(1:1).eq.')') i82431 = 66                          
                                     if(i50315(1:1).eq.'+') i82431 = 67 
                             if(i50315(1:1).eq.'-') i82431 = 68         
           if(i50315(1:1).eq.'&') i82431 = 69                           
        if(i50315(1:1).eq.'.') i82431 = 70                            
                             if(i50315(1:1).eq.',') i82431 = 71         
                                  if(i50315(1:1).eq.':') i82431 = 72    
                             if(i50315(1:1).eq.';') i82431 = 73         
               if(i50315(1:1).eq.'*') i82431 = 74                       
                            if(i50315(1:1).eq.'=') i82431 = 75          
             if(i50315(1:1).eq.'/') i82431 = 76                         
       if(i50315(1:1).eq.'!') i82431 = 80                               
                               if(i50315(1:1).eq.'[') i82431 = 83       
                                 if(i50315(1:1).eq.']') i82431 = 84     
                                                           end          
                             subroutine o53(i11,g,m,me,j13)             
                                                     implicit none      
                                              integer m,me,i            
          double precision g,i11,j13                                    
                                                         dimension g(m) 
                         i11 = 0.0d0                                    
                                        do i=1,me                       
         if (dabs(g(i)).gt.j13) then                                
               i11 = i11 + dabs(g(i))                                   
          endif                                                         
                                               enddo                    
                                  do i=me+1,m                           
                         if (g(i).lt.-j13) then                         
                            i11 = i11 - g(i)                        
                                       endif                            
                                                         enddo      
                                                end                     
         subroutine o213659807013(l,n,nint,m,me,x,f,g,j7,j9,j13,        
     &                 i09410,i725,i3108,i4108,i5542,j2,lj2,j6,lj6,    
     & i813,i87621,i41981820,i78923,i78913,i7193,p,i11,j609,            
     &                                i5308,i31032)                
                                       implicit none                    
              logical i4108                                             
       integer i86x,i472,i315,i86i11,i4296,i37,i98,i2462,i13,i,j,c,i813,
     &i9025120,pt,w,i8087,i3,i4,i110,k,i81,i64,i03s,i03,bi03,di03,di64, 
     & l,n,nint,m,me,i09410,i725,i3108,i5542,lj2,j6,lj6,i41981820,j609, 
     &   i5308,i21123,j046789,i08091,i821013                            
        double precision x,f,g,j7,j9,j13,j2,i87621,i7193(7),o89         
                    character*60 i31032                            
      dimension x( l*n ),f( l ),g( l*m+1 ),j7(n),j9(n),j2(lj2),j6(lj6)  
      double precision i78923(11),i78913(18),p(j609),i11(j609)          
         data i9025120,i86x,i472,i315,i86i11,i4296,i37,i98,i2462,i13,   
     &       pt,w,i8087,i3,i4,i110,k,i81,i64,i03s,i03,bi03,di03,di64    
     & /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/                
         if(i3108.ge.0)then                                          
                 call o2840(n,m,i86x,i472,i315,i86i11,i4296,i813,       
     &  i37,i98,i2462,i13,w,i8087,pt,i3,i4,i110)                        
             call o7100(i81,i03s,k,i64,i03,bi03,di03,di64)      
                                  j2(1) = 1.2d3                         
                           do i=1,i5542                                 
                                                    j2(1) = o89(j2(1))  
                                                 enddo                  
                              call o3261853008(j6(i5308+2),i31032) 
                        call o319(i11( 1 ),g( 1 ),m,me,j13)             
                                               j2(i472)   = f( 1 )  
                   j2(i86i11) = i11( 1 )                                
                        do i=1,n                                        
                        j2(i86x+i-1) = x( i )                           
             enddo                                                      
                                    do i=1,m                            
                                                  j2(i315+i-1) = g( i ) 
                                                enddo                 
                  i21123 = 0                                            
                                                      do i=1,15         
         i21123 = i21123 + j6(i5308+1+i*3)                              
                                       enddo                            
      j046789 =  598                         
                                    if(i21123.ne.j046789)then           
             goto 11                                                   
              endif                                              
                                                   goto 101             
                                        endif                           
                          do c = 1,l                                    
            call o319(i11( c ),g( (c-1)*m+1 ),m,me,j13)                 
             if(m.gt.0) call o732(p( c ),f( c ),i11( c ),j2(i4296),j13) 
             if(m.eq.0) p( c ) = f( c )                                 
                     if(i3108.gt.-30.or.i3108.lt.-40)then               
                         call o8971310(n,j6(k),j2,lj2,i37,i98,i2462,i13,
     &               x( (c-1)*n+1 ),f( c ),i11( c ),p( c ),i78913(14))  
                                                      endif             
                              if(i3108.le.-30.and.i3108.ge.-40)then     
           call o275315065(l,c,n,j2,lj2,j6,lj6,i37,i98,i2462,i13,       
     &  x( (c-1)*n+1 ),f( c ),i11( c ),p( c ),i64,i3108,                
     &                             i78913(15),i78913(4),i78913(5))      
                               endif                                    
                        if(i11( c ).LT.j2(i86i11)) goto 123             
           if(i11( c ).eq.j2(i86i11).and.f( c ).LT.j2(i472)) goto 123   
                                                  goto 100              
  123                                     j2(i472)   = f( c )           
                             j2(i86i11) = i11( c )                      
                                              do i=1,n               
                   j2(i86x+i-1) = x( (c-1)*n+i )                        
           enddo                                                        
                      do i=1,m                                          
                                         j2(i315+i-1) = g( (c-1)*m+i )  
                                  enddo                                 
  100                                                         continue  
                                                             enddo      
  101                                             if(i4108)goto 999     
                       if(i3108.le.-90)then                             
           if(j2(i110).gt.j13.and.j2(i2462).LT.j2(i110))goto 81         
               if(j2(i110).le.j13.and.j2(i2462).le.j13.                 
     &                     and.j2(i98).LT.j2(i4))goto 81             
                                                            goto 82     
   81                        j6(11) = 1                                 
                                                        goto 83         
   82                         j6(11) = 0                                
   83                       continue                                    
                                                 endif                  
                           if(i3108.eq.-10)then                         
                                         if(j2(i13).lt.j2(i3)) goto 84  
                                                 j6(13) = 0             
                                        goto 85                         
   84                                 j6(13) = 1                        
                                           j2(i3) = j2(i13)             
   85                                          continue              
                         endif                                          
 1000    if(i09410.gt.0)then                                            
               continue
                                     endif                              
          if(j6(i64).ge.int(i78913(3))) i3108 = -95                     
          if(i4108) goto 3                                              
                                     if(i3108.eq. -1)      goto 13      
                  if(i3108.eq. -2)then                                  
                                   i3108 = -1                           
        goto 13                                                         
                                                          endif         
              if(i3108.eq. -3)then                                      
          if(j6(i03).ge.j6(i03s))then                                   
                                                       i3108 = -30      
                 goto 14                                                
                                 endif                     
                              i3108 = -1                                
                          goto 13                                       
                             endif                                      
          if(i3108.eq.-30)then                                          
                                                 i3108 = -31            
                    goto 14                                             
                               endif                           
       if(i3108.le.-31.and.                                             
     &           i3108.ge.-39)then                                     
                                i3108 = i3108                       
                                                         goto 14        
                                        endif                           
              if(i3108.eq.-40)then                                      
             i3108 = -2                                                 
                       goto 12                                          
                                                                 endif  
                              if(i3108.eq.-10)then                      
                                                          i3108 = -30   
                                     goto 14                            
                                 endif                                  
             if(i3108.le.-90)then                                       
                                                   i3108 = -3           
                                   goto 11                              
                                  endif                                 
                                       if(i3108.eq.  0)then             
                         i3108 = -3                                     
                        goto 11                                
                                                       endif            
   11           j6(i81) = j6(i81)+1                                     
                                            j6(i64) = 0            
       i821013 = 3367                         
              call o953814(n,i813,j2,lj2,i472,i86i11,i4296,j13,         
     &                   i37,i98,i2462,i13,j6,lj6,i81,k,bi03,di03,di64, 
     &      i87621,i41981820,i78923,i78913,i7193,i5308,                 
     &                          i31032)                            
                              call o26919(j6(k),j2,lj2,w)               
                   call o3261853008(j6(i5308+1),i31032)            
                           j2(pt)=0.0d0                                 
                do j=1,j6(k)                                            
                              j2(pt+j) = j2(pt+j-1) + j2(w+j-1)         
                                                             enddo      
                  j2(i4)   = j2(i472)                                   
                          j2(i110) = j2(i86i11)                         
                                i9025120 = 0                            
                                        i08091 = j6(i5308+1)**2         
                              do i= 0,40                                
                            i08091 = i08091 + j6(i5308+i+10)            
        enddo                                                           
             if(i08091.ne.i821013)then                                  
        do i=1,lj6                                            
           j6(i) = int(j2(lj6+i))                                       
                       enddo                           
                                                    endif               
                                                   if(j6(i81).eq.1)then 
             if(m.gt.0) call o732(p( 1 ),f( 1 ),i11( 1 ),j2(i4296),j13) 
               if(m.eq.0) p( 1 ) = f( 1 )                               
                   call o8971310(n,j6(k),j2,lj2,i37,i98,i2462,i13,      
     &       x( 1 ),f( 1 ),i11( 1 ),p( 1 ),i78913(14))                  
                                              endif                     
   12                        j6(i64) = j6(i64) + 1                      
                                                           j6(i03) = 0  
         call o61867290(n,nint,j2,lj2,i8087,i37,j6,lj6,k,i64,i78923(8)) 
                                 if(i7193(5).gt.0.0d0)then              
                                         j6(i03s) = j6(bi03)            
               else                                                     
              call o329519(j6,lj6,i64,i03s,bi03,di03,di64)              
                                              endif                     
        if(j6(i64).eq.1) j2(i3) = 1.0D16                                
                    if(j6(i64).gt.1) j2(i3) = j2(i13)                   
   13                           do c = 1,l                              
                                    j6(i03) = j6(i03) + 1           
                       if(j6(i64).eq.1)then                             
             if(j6(10).le.1)then                                   
                     if(i7193(2).gt.0.0d0.and.j6(i81).eq.1)then         
              call o32481(n,nint,x( (c-1)*n+1 ),j7,j9,                  
     &                                     j2,lj2,i86x,i7193(2))     
                                                else                    
        if(i7193(2).gt.0.0d0)then                                       
                  call o32481(n,nint,x( (c-1)*n+1 ),j7,j9,              
     &                           j2,lj2,i86x,i7193(2))               
                                                                else    
                   call o1309(n,nint,x( (c-1)*n+1 ),j7,j9,j2,lj2)       
                                 endif                                  
                                                 endif                
                         endif                                          
                  if(j6(10).gt.1)call o83517(n,nint,j2,lj2,j6,lj6,i86x, 
     &x( (c-1)*n+1 ),j7,j9,                                             
     &             i78923(1),i78913(16),i7193(2))                       
                      endif                                             
       if(j6(i64).gt.1)call o95420(n,nint,j6(k),x( (c-1)*n+1 ),j7,j9,  
     &                   j2,lj2,i37,i8087,pt,i78923(5))                 
        enddo                                                           
            if(j6(i03).ge.j6(i03s).and.i3108.ne.-3) i3108 = -10         
    3                    return                                         
   14                                  continue                 
         call o3261853008(j6(i5308+1),i31032)                      
                       if(j6(13).eq.1.or.j6(i64).eq.1)then              
               i3108 = -2                                               
                  goto 12                                               
                    else                                                
      if(i3108.lt.-30.and.j6(31).eq.1)then                              
                                                           i3108 = -2   
                     goto 12                                            
                                                             endif      
                  if(i3108.eq.-39)then                                  
                                      i9025120 = 1                      
                i3108    = -99                                          
                                goto 101                         
                                                               endif    
              do c = 1,l                                                
                if(l.gt.1) j6(31) = 0                                   
             call o158(l,n,nint,x( (c-1)*n+1 ),j7,j9,i37,i8087,i5308,j2,
     & lj2,j6,lj6,i3108,i9025120,i78923(6),i78923(7),i78913(13))        
       if(i3108.eq.-30.and.l.gt.1) i3108 = -31                          
                                        if(i9025120.eq.1.and.c.gt.1)then
              if(o89(j2(1)).ge.0.33d0)then                              
            call o24318(n,nint,x( (c-1)*n+1 ),j7,j9,j2,lj2,i37)         
                                   else                                 
              call o95420(n,nint,j6(k),x( (c-1)*n+1 ),j7,j9,         
     &                         j2,lj2,i37,i8087,pt,i78923(5))           
                             endif                                      
                                        i9025120 = 0               
                        i3108 = -39                                     
             endif                                                      
                                                                enddo   
       if(i9025120.eq.1) goto 101                                       
                                                       goto 3           
                      endif                                             
  999                             f( 1 )   = j2(i472)                   
                                                   i11( 1 ) = j2(i86i11)
                      do i=1,n                                          
                                x(i) = j2(i86x+i-1)                     
                                                                  enddo 
          do j=1,m                                                      
                    g(j) = j2(i315+j-1)                                 
                                          enddo                      
                 if(i11( 1 ).le.j13)then                                
                                           i3108 = 0                    
                                                              else      
                                                 i3108 = 1           
                                                               endif    
                 if(i09410.gt.0)then                                    
          continue
                      goto 1000                                         
                                                     endif              
                                    end                                 
            subroutine o275315065(l,c,n,j2,lj2,j6,lj6,i37,i98,i2462,i13,
     & x,f,i11,p,i64,i3108,O1,R1,R2)                               
                                implicit none                           
        integer l,c,n,lj2,lj6,j6,i64,i37,i98,i2462,i13,i,i3108          
                         double precision x,f,i11,p,j2,OO,OOO,O1,R1,R2  
                                 dimension x(n),j2(lj2),j6(lj6)         
                    data OO,OOO /0.0d0,0.0d0/                           
                 if(i3108.eq.-30)then                                   
                          OO  = O1 + R1 * dsqrt(dble(j6(i64))) + dble(n)
                                   OOO = O1 + R2 * dsqrt(dble(j6(i64))) 
                                                          endif         
       if(i11.le.0.0d0.and.j2(i2462).le.0.0d0)then                      
                               if(f.ge. j2(i98) - dabs(j2(i98))/OO )then
                        j6(31 + c) = 0                              
                                      goto 1                            
                               endif                            
                                                          else          
                       if(p.ge. j2(i13) - dabs(j2(i13))/OOO )then       
                  j6(31 + c) = 0                                        
                                        goto 1                          
                 endif                                                  
              endif                                                     
        do i = 1,n                                                      
                 j2(i37+i-1) = x(i)                                   
                      enddo                                             
                                           j2(i2462) = i11              
               j2(i98)   = f                                            
                                                        j2(i13)   = p   
                                        j6(31 + c) = 1                  
   1                                           if(c.eq.l)then           
                                                             j6(31) = 0 
                      do i = 1,l                                        
                          j6(31) = j6(31) + j6(31+i)                    
                                                 enddo            
                                      if(j6(31).gt.1) j6(31) = 1        
                                       endif                            
                                return                              
                                                             end        
                           subroutine o62886010(i50315,i82431)          
              implicit none                                             
                                        character*1 i50315              
                                             integer i82431             
                                             i82431 = 0                 
                             if(i50315(1:1).eq.'A') i82431 = 52         
                  if(i50315(1:1).eq.'B') i82431 = 28                    
                                 if(i50315(1:1).eq.'C') i82431 = 49     
            if(i50315(1:1).eq.'D') i82431 = 30                          
                                      if(i50315(1:1).eq.'E') i82431 = 31
             if(i50315(1:1).eq.'F') i82431 = 32                         
               if(i50315(1:1).eq.'G') i82431 = 33                       
                     if(i50315(1:1).eq.'H') i82431 = 34               
                          if(i50315(1:1).eq.'I') i82431 = 35            
                               if(i50315(1:1).eq.'J') i82431 = 36       
                                   if(i50315(1:1).eq.'K') i82431 = 37   
                    if(i50315(1:1).eq.'L') i82431 = 38                  
                   if(i50315(1:1).eq.'M') i82431 = 39                
                        if(i50315(1:1).eq.'N') i82431 = 40              
                     if(i50315(1:1).eq.'O') i82431 = 41                 
                     if(i50315(1:1).eq.'P') i82431 = 42                 
                 if(i50315(1:1).eq.'Q') i82431 = 43                     
               if(i50315(1:1).eq.'R') i82431 = 44                       
                       if(i50315(1:1).eq.'S') i82431 = 45               
             if(i50315(1:1).eq.'T') i82431 = 46                         
                           if(i50315(1:1).eq.'U') i82431 = 47           
                             if(i50315(1:1).eq.'V') i82431 = 48         
                  if(i50315(1:1).eq.'W') i82431 = 29                    
                    if(i50315(1:1).eq.'X') i82431 = 50                  
             if(i50315(1:1).eq.'Y') i82431 = 51                         
      if(i50315(1:1).eq.'Z') i82431 = 27                             
                             if(i50315(1:1).eq.'0') i82431 = 53         
                                    if(i50315(1:1).eq.'1') i82431 = 54  
          if(i50315(1:1).eq.'2') i82431 = 55                            
                                      if(i50315(1:1).eq.'3') i82431 = 56
             if(i50315(1:1).eq.'4') i82431 = 57                         
                      if(i50315(1:1).eq.'5') i82431 = 58                
                     if(i50315(1:1).eq.'6') i82431 = 59                 
        if(i50315(1:1).eq.'7') i82431 = 60                              
        if(i50315(1:1).eq.'8') i82431 = 61                              
                                    if(i50315(1:1).eq.'9') i82431 = 62  
                       if(i50315(1:1).eq.'a') i82431 = 23               
      if(i50315(1:1).eq.'b') i82431 = 2                                 
                   if(i50315(1:1).eq.'c') i82431 = 3                    
       if(i50315(1:1).eq.'d') i82431 = 16                        
          if(i50315(1:1).eq.'e') i82431 = 5                             
                  if(i50315(1:1).eq.'f') i82431 = 13                    
                       if(i50315(1:1).eq.'g') i82431 = 7               
                                      if(i50315(1:1).eq.'h') i82431 = 8 
                                     if(i50315(1:1).eq.'i') i82431 = 9  
               if(i50315(1:1).eq.'j') i82431 = 10                       
                                if(i50315(1:1).eq.'k') i82431 = 11      
              if(i50315(1:1).eq.'l') i82431 = 12                        
                                if(i50315(1:1).eq.'m') i82431 = 6       
                                 if(i50315(1:1).eq.'n') i82431 = 14     
       if(i50315(1:1).eq.'o') i82431 = 15                               
                 if(i50315(1:1).eq.'p') i82431 = 4                      
              if(i50315(1:1).eq.'q') i82431 = 17                        
                                   if(i50315(1:1).eq.'r') i82431 = 18   
            if(i50315(1:1).eq.'s') i82431 = 19                          
                   if(i50315(1:1).eq.'t') i82431 = 20                   
                             if(i50315(1:1).eq.'u') i82431 = 21         
               if(i50315(1:1).eq.'v') i82431 = 22                       
              if(i50315(1:1).eq.'w') i82431 = 1                         
                                     if(i50315(1:1).eq.'x') i82431 = 24 
                if(i50315(1:1).eq.'y') i82431 = 25                  
             if(i50315(1:1).eq.'z') i82431 = 26                     
                                   if(i50315(1:1).eq.'_') i82431 = 64   
                                    if(i50315(1:1).eq.'(') i82431 = 65  
                          if(i50315(1:1).eq.')') i82431 = 66            
                     if(i50315(1:1).eq.'+') i82431 = 67                 
                                    if(i50315(1:1).eq.'-') i82431 = 68  
                             if(i50315(1:1).eq.'&') i82431 = 69         
       if(i50315(1:1).eq.'.') i82431 = 70                               
                  if(i50315(1:1).eq.',') i82431 = 71                    
                      if(i50315(1:1).eq.':') i82431 = 72                
                                   if(i50315(1:1).eq.';') i82431 = 73   
                if(i50315(1:1).eq.'*') i82431 = 74                      
                    if(i50315(1:1).eq.'=') i82431 = 75                  
                              if(i50315(1:1).eq.'/') i82431 = 76        
              if(i50315(1:1).eq.'!') i82431 = 80                        
                        if(i50315(1:1).eq.'[') i82431 = 83              
            if(i50315(1:1).eq.']') i82431 = 84                          
                                                          end           
                     function o89(s)                                    
                      implicit none                                   
                                         double precision o89,s,a,b,c   
                                  data  a,b,c  /0.0d0,0.0d0,0.0d0/      
                                             if(s.eq.1.2d3)then         
                                     a =  0.485414306917525406604D+00   
                        b =  0.564807209834307433205D+00                
          c =  0.180868201858223220935D+00                        
                 endif                                                  
                                                  s = a + b + c         
           if(b.lt.0.5d0) s = s + 0.493127909786063800546D+00           
                                       if(s.ge.1.0d0) s = s - 1.0d0     
                if(s.ge.1.0d0) s = s - 1.0d0                            
                               a = b                                    
             b = c                                                
                 c = s                                                  
                                        s = a                       
                                            o89 = s                     
                                                           end          
            subroutine o7100(i81,i03s,k,i64,i03,bi03,di03,di64)         
                                               implicit none            
          integer i81,i03s,k,i64,i03,bi03,di03,di64                     
                                                      k      = 1        
                                  i81    = 2                            
                                                     i03s   = 3         
                                                     i64    = 4         
                                                i03    = 5              
                                                         bi03   = 6     
                               di03   = 7                          
                                                    di64   = 8          
                                                             end        
                                subroutine o3261853008(i82431s,i6232)   
                                                     implicit none      
              character*60 i6232                                    
                                          integer i82431s(*),i          
                                              do i = 1,60               
                              call o6288601(i6232(i:i),i82431s(i))    
                                        enddo                           
                            end                                         
                 subroutine o319(i11,g,m,me,j13)                    
                                     implicit none                      
                integer m,me                                            
                                            double precision g,i11,j13  
                                                   dimension g(m)       
       i11 = 0.0d0                                                      
                                            if(m.eq.0)return            
                        call o21(i11,g,m,me,j13)                        
                                if(i11.gt.j13)call o53(i11,g,m,me,j13)  
                   end                                                  
        subroutine o26919(k,j2,lj2,w)                                   
       implicit none                                                    
                                           integer k,lj2,j,i73,w     
                                double precision j2                     
                               dimension j2(lj2)                        
                                                            i73 = 0     
                                                do j=1,k                
                                            i73 = i73 + j               
                           enddo                                        
                                            do j=1,k                    
          j2(w+j-1) = (k-j+1)/dble(i73)                                 
        enddo                                              
                                               end                      
       subroutine o24318(n,nint,x,j7,j9,j2,lj2,i37)                     
                                        implicit none                   
           integer n,nint,lj2,i37,i                                     
                 double precision x,j7,j9,j2,o89,i82431,i130,i9042677836
             dimension x(n),j7(n),j9(n),j2(lj2)                         
                                             do i = 1,n-nint            
                            i82431 = o89(j2(1))                         
          if(i82431.le.0.5d0)then                                       
                                                x(i) = j2(i37+i-1)      
                                                       else             
            i130 = (j9(i)-j7(i)) / 1.0d1**(i82431*6.0d0)                
                                          x(i) = j2(i37+i-1) + i130 *   
     &            i9042677836(o89(j2(1)),o89(j2(1)))                    
                                                                  endif 
              enddo                                                     
               do i = n-nint+1,n                               
         i82431 = o89(j2(1))                                            
                  if(i82431.le.0.5d0)then                               
                                   x(i) = j2(i37+i-1)                   
                             else                                       
               if(i82431.le.0.75d0)then                                 
                      x(i) = j2(i37+i-1) + 1.0d0                        
                                   else                                 
                          x(i) = j2(i37+i-1) - 1.0d0               
                                                                endif   
                                   endif                                
          enddo                                                        
                                         do i = 1,n                     
                        if(x(i).lt.j7(i))then                           
                       x(i) = j7(i) + (j7(i)-x(i)) / (o89(j2(1))*1.0d6) 
                                     if(x(i).lt.j7(i)) x(i) = j7(i)     
                               if(x(i).gt.j9(i)) x(i) = j9(i)           
                                                         endif          
         if(x(i).gt.j9(i))then                                          
               x(i) = j9(i) - (x(i)-j9(i)) / (o89(j2(1))*1.0d6)         
                               if(x(i).lt.j7(i)) x(i) = j7(i)           
           if(x(i).gt.j9(i)) x(i) = j9(i)                               
                      endif                                       
                                                       enddo            
                             end                                        
                            FUNCTION i9042677836(i87z679,i650012)       
                                                          IMPLICIT NONE 
          INTEGER I,J                                                   
                        DOUBLE PRECISION i9042677836,i87z679,i650012    
                DOUBLE PRECISION i65087210034(30),i087660126578(30)  
                                        DATA i65087210034 /           
     &  0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,    
     & 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,     
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/ 
                                 DATA i087660126578 /                  
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     & 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,     
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
       DATA I /0/                                                       
                                                IF(I.EQ.0)THEN          
                 i65087210034(   1) =   0.260390290432194637659791D+00
         i087660126578(   1) =   0.207911690817759398086650D+00        
             i65087210034(   2) =   0.371464322612418407221213D+00    
         i087660126578(   2) =   0.406736643075800319291346D+00        
             i65087210034(   3) =   0.459043605026420720172098D+00    
             i087660126578(   3) =   0.587785252292473248125759D+00    
          i65087210034(   4) =   0.534978211968811456777928D+00       
               i087660126578(   4) =   0.743144825477394466162195D+00  
                i65087210034(   5) =   0.603856865149274613102648D+00 
          i087660126578(   5) =   0.866025403784438818632907D+00       
               i65087210034(   6) =   0.668047230836577465851178D+00  
      i087660126578(   6) =   0.951056516295153642204241D+00           
             i65087210034(   7) =   0.728976221468170537676201D+00    
             i087660126578(   7) =   0.994521895368273400883652D+00    
         i65087210034(   8) =   0.787597521966441505014700D+00        
         i087660126578(   8) =   0.994521895368273178839047D+00        
         i65087210034(   9) =   0.844600430900591558902590D+00        
              i087660126578(   9) =   0.951056516295153198115031D+00   
          i65087210034(  10) =   0.900516638500549193580014D+00       
           i087660126578(  10) =   0.866025403784437819432185D+00      
      i65087210034(  11) =   0.955780730602699413189782D+00           
       i087660126578(  11) =   0.743144825477393355939171D+00          
        i65087210034(  12) =   0.101076765259478973391083D+01         
       i087660126578(  12) =   0.587785252292471804835827D+00          
            i65087210034(  13) =   0.106581803100335958944811D+01     
           i087660126578(  13) =   0.406736643075798820490263D+00      
                 i65087210034(  14) =   0.112125702621867584518611D+01
                i087660126578(  14) =   0.207911690817757566218660D+00 
           i65087210034(  15) =   0.117741002251547466350701D+01      
      i087660126578(  15) =  -0.165389215948551511930853D-01           
                i65087210034(  16) =   0.123461739178329787947064D+01 
                 i087660126578(  16) =  -0.207911690817760813621007D+00
      i65087210034(  17) =   0.129325018786050716101954D+01           
      i087660126578(  17) =  -0.406736643075801818092430D+00           
            i65087210034(  18) =   0.135372872605567096115919D+01     
                 i087660126578(  18) =  -0.587785252292474469371086D+00
            i65087210034(  19) =   0.141654658155938162344967D+01     
                i087660126578(  19) =  -0.743144825477395465362918D+00 
                i65087210034(  20) =   0.148230380736751121695249D+01 
              i087660126578(  20) =  -0.866025403784440150900537D+00   
          i65087210034(  21) =   0.155175565365552059482468D+01       
       i087660126578(  21) =  -0.951056516295154086293451D+00          
           i65087210034(  22) =   0.162588796660921230952113D+01      
                i087660126578(  22) =  -0.994521895368273622928257D+00 
           i65087210034(  23) =   0.170604058135018687991646D+01      
           i087660126578(  23) =  -0.994521895368273067816745D+00      
                i65087210034(  24) =   0.179412257799410146397179D+01 
               i087660126578(  24) =  -0.951056516295152531981216D+00  
            i65087210034(  25) =   0.189301847282484558832039D+01     
       i087660126578(  25) =  -0.866025403784437264320673D+00          
               i65087210034(  26) =   0.200743768049833359867762D+01  
           i087660126578(  26) =  -0.743144825477392245716146D+00      
                i65087210034(  27) =   0.214596602628934718381970D+01 
         i087660126578(  27) =  -0.587785252292470472568198D+00        
       i65087210034(  28) =   0.232725168432733564571890D+01          
                i087660126578(  28) =  -0.406736643075796877599970D+00 
            i65087210034(  29) =   0.260814009656772682888004D+01     
                i087660126578(  29) =  -0.207911690817756372728908D+00 
       i65087210034(  30) =   0.290814009656772682888004D+01          
              i087660126578(  30) =  -0.107911690817756372728908D+00   
                                                              ENDIF     
            I = INT(i87z679*30)                                         
                                      J = INT(i650012*30)               
                                        IF(I.LT.1) I = 1                
                     IF(J.LT.1) J = 1                                   
              i9042677836 = i65087210034(I) * i087660126578(J)       
                        end                                             
        subroutine midaco_code(l,n,nint,m,me,x,f,g,j7,j9,j13,           
     &                            i3108,i4108,i7193,j2,lj2,j6,lj6,     
     &       i31032,j609,p,i11)                                    
                              implicit none                             
              integer l,n,nint,m,me,i3108,lj2,j6,lj6,i4108              
                double precision x,f,g,j7,j9,j13,i7193(7),j2            
       dimension x(l*n),f(l),g(l*m+1),j7(n),j9(n),j2(lj2),j6(lj6)       
                          character*60 i31032                      
                 integer j609                                           
            double precision p(j609),i11(j609)                          
                                           logical i0087                
                  integer i813,i41981820,i03629661,i2419,i2134,i71904,J,
     &    i619307,i5308,i81i1403,i30251,i3193,i6193,i0311,i432i11,      
     &         i2098473,i3655401800,i2310315,i321i1403,i,c,i5542,i87621,
     &   i0009,i48213,i824,i98320,ii30919,i31579907326,i087676578       
          double precision o89,i130,i78923(11),i78913(18),i8250,i083261 
           double precision i9042677836                                 
      data i813,i41981820,i03629661,i48213,i824,i2419,i2134,i87621,     
     &      i0009,i5308,i81i1403,i30251,i3193,i6193,i0311,i432i11,      
     &        i2098473,i3655401800,i2310315,i321i1403,i31579907326      
     &     /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/                  
                     data i78923,i78913,i8250,i083261  /0.0d0,          
     &     0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &   0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, 
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/      
                  if(i3108.ge.0)then                                    
                                                      i48213 = 0        
             i824   = 0                                              
                     if(i3108.gt.10.AND.i3108.lt.100)then               
        i3108 = -3                                                      
                                       i31579907326 = 0                 
                              goto 79                              
                                    endif                               
                                              do i = 1,60               
                           call o62886010(i31032(i:i),j6(i))       
                              enddo                                     
                                                      i71904 = 0        
                   do i=1,10                                            
           i71904 = i71904 - j6(i)                                      
                                                      enddo             
                                              do i=10,60                
                             i71904 = i71904 + j6(i)                 
                                    enddo                               
                  i619307 = 1829             
        if(i71904.ne.i619307)then                                       
        i3108 = 900                                                     
                                                          return        
                                                         endif        
         call o3261853008(j6(1),i31032)                            
                                          if(int(i7193(7)).eq.1) goto 51
                    if(int(i7193(7)).eq.2) goto 52                      
                                   if(int(i7193(7)).eq.3) goto 53       
   51                        if(nint.eq.n.or.int(i7193(7)).gt.0)then    
                           if(n.le.50)then                              
                         i78923(  1) = 0.389078763191433141D+00         
               i78923(  2) = 0.114860512107443344D+00                   
                     i78923(  3) = 0.113443339033864277D+01             
            i78923(  4) = 0.518799365371665044D+00                      
        i78923(  5) = 0.432539059007058501D+09                          
              i78923(  6) = 0.100089252383195868D-02                    
        i78923(  7) = 0.999999999459014188D+00                          
               i78923(  8) = 0.397831016063179490D+02                   
                    i78923(  9) = 0.389801369146012405D-02              
                    i78923( 10) = 0.155844656119642398D-06         
       i78923( 11) = 0.623332071588316938D+00                           
                   i78913(  1) = 0.200000000000000000D+01           
              i78913(  2) = 0.100000000000000000D+02                    
                          i78913(  3) = 0.990000000000000000D+02        
                                 i78913(  4) = 0.148000000000000000D+03 
           i78913(  5) = 0.490000000000000000D+02                       
                                 i78913(  6) = 0.500000000000000000D+01 
             i78913(  7) = 0.110000000000000000D+02                     
                             i78913(  8) = 0.800000000000000000D+01     
                           i78913(  9) = 0.890000000000000000D+02       
                           i78913( 10) = 0.400000000000000000D+01       
         i78913( 11) = 0.200000000000000000D+01                         
                     i78913( 12) = 0.800000000000000000D+01             
                             i78913( 13) = 0.900000000000000000D+01     
                             i78913( 14) = 0.000000000000000000D+00     
       i78913( 15) = 0.276507832000000000D+09                           
               i78913( 16) = 0.200000000000000000D+01                   
                               i78913( 17) = 0.100000000000000000D+01   
                          i78913( 18) = 0.124764362584000000D+12        
                                      else                              
                i78923(  1) = 0.242158155864236890D+00                  
                            i78923(  2) = 0.703454504998440644D+00      
                 i78923(  3) = 0.107467917008615577D+01                 
             i78923(  4) = 0.758136520883294196D+00                     
       i78923(  5) = 0.588648485982170820D+09                           
           i78923(  6) = 0.151615115132434426D+03                   
               i78923(  7) = 0.882788933354880179D+00                   
                            i78923(  8) = 0.318417239075448286D+02      
         i78923(  9) = 0.100000000583138296D-07                         
                    i78923( 10) = 0.139392038646105073D-06              
            i78923( 11) = 0.548565908377584011D+00                      
             i78913(  1) = 0.200000000000000000D+01                     
                                 i78913(  2) = 0.100000000000000000D+02 
                         i78913(  3) = 0.575000000000000000D+03         
                      i78913(  4) = 0.100000000000000000D+02            
                      i78913(  5) = 0.100000000000000000D+04            
                  i78913(  6) = 0.400000000000000000D+01                
                i78913(  7) = 0.600000000000000000D+01                  
           i78913(  8) = 0.110000000000000000D+02                       
         i78913(  9) = 0.940000000000000000D+02                         
                                i78913( 10) = 0.100000000000000000D+01  
                       i78913( 11) = 0.300000000000000000D+01           
                i78913( 12) = 0.800000000000000000D+01                  
                     i78913( 13) = 0.120000000000000000D+02             
                                i78913( 14) = 0.100000000000000000D+01  
           i78913( 15) = 0.421740554000000000D+09                       
                            i78913( 16) = 0.720000000000000000D+02   
                                 i78913( 17) = 0.180000000000000000D+02 
                      i78913( 18) = 0.325376752130000000D+11            
                                                                endif   
                                                   goto 54              
                                                       endif            
   52                          if(nint.eq.0.or.int(i7193(7)).gt.0)then  
                if(n.le.50)then                                         
         i78923(  1) = 0.249122144659019495D+01                         
                                  i78923(  2) = 0.844388223088001766D+00
      i78923(  3) = 0.100153091583909348D+01                            
              i78923(  4) = 0.143961720122020398D+01                    
                i78923(  5) = 0.258324699287332147D+08                  
                     i78923(  6) = 0.943969219703401563D+03             
          i78923(  7) = 0.973372517527523608D-07                        
      i78923(  8) = 0.694079968011245829D+00                            
                i78923(  9) = 0.842177534294979019D-01                  
                               i78923( 10) = 0.275806030201900532D+01   
                       i78923( 11) = 0.324049619462284000D+00           
                 i78913(  1) = 0.200000000000000000D+01                 
                              i78913(  2) = 0.100000000000000000D+02    
                           i78913(  3) = 0.343000000000000000D+03       
                    i78913(  4) = 0.315000000000000000D+03           
                        i78913(  5) = 0.401000000000000000D+03          
                      i78913(  6) = 0.200000000000000000D+01            
                  i78913(  7) = 0.900000000000000000D+01                
           i78913(  8) = 0.500000000000000000D+01                       
                   i78913(  9) = 0.160000000000000000D+02               
                     i78913( 10) = 0.100000000000000000D+01             
                    i78913( 11) = 0.200000000000000000D+01              
                         i78913( 12) = 0.600000000000000000D+01         
                               i78913( 13) = 0.400000000000000000D+01   
       i78913( 14) = 0.000000000000000000D+00                           
                               i78913( 15) = 0.983025185000000000D+09   
                        i78913( 16) = 0.600000000000000000D+01          
                        i78913( 17) = 0.270000000000000000D+02          
            i78913( 18) = 0.904969100000000000D+09                      
                 else                                                   
            i78923(  1) = 0.118876691512506394D+01                      
         i78923(  2) = 0.607544623629922120D+00                         
                           i78923(  3) = 0.123970349748770015D+01       
              i78923(  4) = 0.118509064728880364D+01                    
                      i78923(  5) = 0.348699758754453659D+08            
                             i78923(  6) = 0.997161960178311460D+03     
               i78923(  7) = 0.602039328669455594D-07                   
          i78923(  8) = 0.505753968641715979D+01                        
                 i78923(  9) = 0.790696417057785483D+00            
                           i78923( 10) = 0.210656581785541874D+01       
         i78923( 11) = 0.349012629457182089D+00                         
         i78913(  1) = 0.200000000000000000D+01                         
                                 i78913(  2) = 0.100000000000000000D+02 
            i78913(  3) = 0.490000000000000000D+03                      
                    i78913(  4) = 0.130000000000000000D+02              
       i78913(  5) = 0.649000000000000000D+03                           
                         i78913(  6) = 0.200000000000000000D+01         
                    i78913(  7) = 0.250000000000000000D+02              
                             i78913(  8) = 0.110000000000000000D+02     
                      i78913(  9) = 0.400000000000000000D+01            
                                  i78913( 10) = 0.100000000000000000D+01
                           i78913( 11) = 0.200000000000000000D+01       
           i78913( 12) = 0.600000000000000000D+01                       
                        i78913( 13) = 0.700000000000000000D+01          
                    i78913( 14) = 0.000000000000000000D+00              
                     i78913( 15) = 0.603379230000000000D+09             
                  i78913( 16) = 0.500000000000000000D+01                
                             i78913( 17) = 0.290000000000000000D+02     
          i78913( 18) = 0.800782976000000000D+09                        
                                                             endif      
                                 goto 54                            
                                           endif                        
   53             if(nint.gt.0.and.nint.lt.n.or.int(i7193(7)).gt.0)then 
                                          if(n.le.50)then               
            i78923(  1) = 0.739291508632072047D+00                      
         i78923(  2) = 0.633998424333021782D+00                         
            i78923(  3) = 0.103299237259922405D+01                      
                                 i78923(  4) = 0.146764290556402521D+01 
       i78923(  5) = 0.687162903769587874D+09                           
                             i78923(  6) = 0.205802915797508632D+02     
      i78923(  7) = 0.646395800934978280D-08                            
                         i78923(  8) = 0.632961637891720330D+00         
                               i78923(  9) = 0.165635296752857508D+01   
                                i78923( 10) = 0.107013907618447555D-01  
                        i78923( 11) = 0.697751367495330399D+00          
                             i78913(  1) = 0.200000000000000000D+01     
                i78913(  2) = 0.100000000000000000D+02                  
        i78913(  3) = 0.136000000000000000D+03                          
                   i78913(  4) = 0.161000000000000000D+03               
                        i78913(  5) = 0.180000000000000000D+02          
              i78913(  6) = 0.200000000000000000D+01                    
                          i78913(  7) = 0.100000000000000000D+02        
                       i78913(  8) = 0.200000000000000000D+01           
                  i78913(  9) = 0.400000000000000000D+02                
               i78913( 10) = 0.200000000000000000D+01                
                   i78913( 11) = 0.200000000000000000D+01               
                        i78913( 12) = 0.500000000000000000D+01          
                             i78913( 13) = 0.500000000000000000D+01     
                      i78913( 14) = 0.100000000000000000D+01            
               i78913( 15) = 0.364827196000000000D+09                   
                          i78913( 16) = 0.130000000000000000D+02        
                                i78913( 17) = 0.230000000000000000D+02  
           i78913( 18) = 0.886844425000000000D+09                       
                                     else                               
                                  i78923(  1) = 0.977028227081946476D+00
       i78923(  2) = 0.817716107310876406D+00                           
                       i78923(  3) = 0.102176946317372175D+01           
           i78923(  4) = 0.141568482802071527D+01                       
                     i78923(  5) = 0.507710667821344793D+09             
                      i78923(  6) = 0.248657903658426616D+01            
             i78923(  7) = 0.442148199129617578D-07                     
                           i78923(  8) = 0.196339623453502687D+01       
                                  i78923(  9) = 0.226706133235036811D+01
                   i78923( 10) = 0.121920904983831223D-06               
             i78923( 11) = 0.305008867167975650D+00                     
                 i78913(  1) = 0.200000000000000000D+01              
            i78913(  2) = 0.100000000000000000D+02                      
             i78913(  3) = 0.332000000000000000D+03                     
                          i78913(  4) = 0.863000000000000000D+03        
                  i78913(  5) = 0.230000000000000000D+02                
         i78913(  6) = 0.300000000000000000D+01                         
                         i78913(  7) = 0.120000000000000000D+02         
                                i78913(  8) = 0.200000000000000000D+01  
                             i78913(  9) = 0.160000000000000000D+02     
                                i78913( 10) = 0.100000000000000000D+01  
                i78913( 11) = 0.400000000000000000D+01                  
               i78913( 12) = 0.600000000000000000D+01                   
                            i78913( 13) = 0.700000000000000000D+01      
          i78913( 14) = 0.100000000000000000D+01                        
                            i78913( 15) = 0.705433041000000000D+09      
                 i78913( 16) = 0.900000000000000000D+01                 
                           i78913( 17) = 0.210000000000000000D+02       
                          i78913( 18) = 0.843575111000000000D+09        
                     endif                                              
                                        goto 54                         
                                                   endif                
   54                               continue                       
                           i98320 = 1000                                
                              ii30919 = -441  
                                    do i=30,60                          
                        i98320 = i98320 - j6(i)                         
              enddo                                                     
              if(i98320.ne.ii30919)then                                 
        i3108 = -1                                                      
                                                           return       
                        endif                                           
                                 i03629661 = 0                          
                         call i074206156(l,n,nint,m,me,x,f,g,j7,j9,j13, 
     &                      i3108,i4108,i7193,j2,lj2,j6,lj6,         
     &       i813,i0009,i5308,i78913(1),i78913(2),                      
     &                                 i31032,i31579907326)        
          if(i3108.ge.100) goto 86                                      
             if(i31579907326.eq.1)then                                  
                                                    i087676578 = i3108  
               i3108 = 0                                                
            endif                                                       
                                                    i03629661 = 1       
                       i5542 = int(i7193(1))                            
                                              i2310315 = int(i7193(3))  
                  i8250    = 1.0d16                                     
              i083261  = 1.0d16                                         
                 i3193      = 1                                  
      i6193      = i3193 + n                                            
                            i0311      = i6193 + 1                      
                                   i432i11    = i0311 + m               
             i87621     = i432i11 + 1                                   
                                         do i=1,n                       
                j2(i0009+i3193+i-1) = x(i)                              
                              enddo                                     
                   do i=1,m                                             
                               j2(i0009+i0311+i-1) = g(i)               
             enddo                                                      
         j2(i0009+i6193) = f(1)                                         
         call o319(j2(i0009+i432i11),g,m,me,j13)                        
                                i81i1403 = 0                            
                i3655401800 = 0                                         
                                i2134 = i5542                           
       i2419  = 0                                                       
              i41981820 = 1                                             
                                  if(j2(i0009+i432i11).gt.j13)then      
                    if(i7193(4).eq.0.0d0)then                           
                   j2(i0009+i87621) = 1.0D9 + j2(i0009+i6193)           
                                                          else          
               j2(i0009+i87621) = i7193(4)                              
                                                       endif            
           else                                                         
                   j2(i0009+i87621) = j2(i0009+i6193)            
                             endif                                      
        else                                                            
                                if(i03629661.ne.1)then                  
                                            i3108 = 701                 
                                 i4108 = 1                              
          return                                                        
                                 endif                                  
                                    do c = 1,l                          
             call o9953052(f(c),g((c-1)*m+1),m)                         
                                                    enddo               
                                                           endif        
   79                                                    continue       
                               if(i3108.eq.-300)then                    
      i2419 = 0                                                         
                                                  i2134 = i2134 + 1     
             endif                                                      
        if(i4108.eq.0)then                                              
                                                         i0087 = .false.
                                                                else    
                       i0087 = .true.                           
                                                        endif           
            call o213659807013(l,n,nint,m,me,x,f,g,j7,j9,j13,           
     &                                 i48213,i824,i2419,i0087,i2134,  
     &         j2,lj2,j6,lj6,i813,j2(i0009+i87621),i41981820,           
     &          i78923,i78913,i7193,p,i11,j609,i5308,                   
     &                                     i31032)                 
                                       i3108 = i2419                    
                      if(i3108.eq.801)return                            
                                                       if(i0087)then    
       if(j2(i0009+i432i11).gt.j13. and .j2(2+n+1+m).lt.                
     &                            j2(i0009+i432i11))    goto 1          
                  if(j2(i0009+i432i11).le.j13. and .j2(2+n+1+m).le.j13  
     &                      .and.j2(2+n).lt.j2(i0009+i6193)) goto 1     
                                                        goto 3          
                   endif                                                
                                if(i2419.eq.-3) i81i1403 = i81i1403 + 1 
                         i30251 = INT(i78913( 6)) * n + INT(i78913( 7)) 
                                       if(i81i1403.ge.i30251)then       
                       i3655401800 = i3655401800 + 1               
                if(j2(i0009+i432i11).gt.j13. and .j2(2+n+1+m).lt.       
     &   j2(i0009+i432i11))   goto 11                                   
        if(j2(i0009+i432i11).le.j13. and .j2(2+n+1+m).le.j13.           
     &  and .j2(2+n).lt.j2(i0009+i6193))goto 11                         
           goto 12                                                      
   11                               j2(i0009+i6193)     = j2(2+n)       
                                j2(i0009+i432i11)   = j2(2+n+1+m)       
                do i=1,n                                                
                                      j2(i0009+i3193+i-1) = j2(2+i-1)   
                                                  enddo                 
                  do i=1,m                                              
                          j2(i0009+i0311+i-1) = j2(2+n+1+i-1)           
                   enddo                                                
                           if(j2(i0009+i432i11).le.j13)then             
       j2(i0009+i87621) =  j2(i0009+i6193)                              
                                                 endif                  
                                                      i2098473 = 1      
             goto 13                                                    
   12                                                 i2098473 = 0      
   13      do i = 2,i0009                                               
                                     j2(i) = 0.0d0                
                                  enddo                                 
                               do i = 1,i5308                           
         j6(i) = 0                                                      
                                                   enddo                
                                    do c = 1,l                          
                                  do i=1,n                              
              if(n.le.n-nint) i130 = (j9(i)-j7(i))                      
     &  / (o89(j2(1))*1.0d1**i78913(17))                                
              if(n.gt.n-nint) i130 = (j9(i)-j7(i)) / i78913(18)         
              if(i.gt.n-nint.and.i130.lt.i78923(11))then                
                                      i130 = i78923(11)                 
                     endif                                              
       if(i7193(2).gt.0.0D0)then                                        
                                           if(i.le.n-nint)then          
                              i130 = (j9(i)-j7(i)) / i7193(2)           
                                                               else     
                                        i130 = 1.0D0 / DSQRT(i7193(2))  
                            endif                                       
                                                  endif                 
                   x(i) = j2(i0009+i3193+i-1) + i130 *              
     &                                i9042677836(o89(j2(1)),o89(j2(1)))
                                            if(x(i).lt.j7(i))then       
                              x(i)=j7(i)+(j7(i)-x(i)) / i78923( 5)      
                         endif                                          
                                    if(x(i).gt.j9(i))then               
          x(i)=j9(i)-(x(i)-j9(i)) / i78923( 5)                          
                                                                  endif 
                          if(x(i).lt.j7(i))then                         
              x(i)=j7(i)                                                
                                                                endif   
                                   if(x(i).gt.j9(i))then                
         x(i)=j9(i)                                                     
                  endif                                                 
                                                   if(i.gt.n-nint)then  
                                             x(i)=dnint(x(i))           
                                             endif                      
           enddo                                                        
                               enddo                            
                        i41981820 = i41981820 * INT(i78913( 8))         
             if(i41981820.gt.INT(i78913( 9))) i41981820 = 1             
             i3108 = -300                                               
                    i81i1403 = 0                                        
                                                   if(i2310315.gt.0)then
                    if(i8250.eq.1.0d16)then                             
                                      i321i1403 = 0                     
                                               i8250   = j2(i0009+i6193)
                    i083261 = j2(i0009+i432i11)                         
                                            else                        
                                  if(j2(i0009+i432i11).le.i083261)then  
                            if(i083261.le.j13)then                      
                                    if(j2(i0009+i6193).lt.              
     &     i8250-dabs(i8250/1.0d6))then                                 
                                            i8250   = j2(i0009+i6193)   
            i083261 = j2(i0009+i432i11)                                 
                                       i321i1403 = 0                    
                     else                                               
                                   i321i1403 = i321i1403 + 1            
                                              goto 76                   
       endif                                                      
                                                              else      
                      i321i1403 = 0                                     
              i8250   = j2(i0009+i6193)                                 
                                       i083261 = j2(i0009+i432i11)      
                                   endif                                
                                else                                    
                      i321i1403 = i321i1403 + 1                         
                    goto 76                                             
                endif                                                   
           endif                                                        
   76           continue                                                
                          if(i321i1403.ge.i2310315)then                 
           if(j2(i0009+i432i11).le.j13)then                             
                                                           i3108 = 3    
                                       else                             
                                          i3108 = 4                     
                       endif                                            
                          goto 3                                        
                                           endif                        
                                                       endif            
                                     endif                              
                           if(i31579907326.eq.1)then                    
                      i3108 = i087676578                                
                                               endif               
    4                                  return                           
    1   j2(i0009+i6193)     = j2(2+n)                                   
                                 j2(i0009+i432i11)   = j2(2+n+1+m)      
                                                             do i=1,n   
                              j2(i0009+i3193+i-1) = j2(2+i-1)           
          enddo                                                         
                                        do i=1,m                        
             j2(i0009+i0311+i-1) = j2(2+n+1+i-1)                        
                                                             enddo      
                    if(j2(i0009+i432i11).le.j13)then                    
             j2(i0009+i87621) =  j2(i0009+i6193)                        
       endif                                                            
                                               i2098473 = 1             
                     continue                                           
    3                             f(1) = j2(i0009+i6193)                
                                        do i = 1,n                      
           x(i) = j2(i0009+i3193+i-1)                            
                                                          enddo         
                            do i = 1,m                                  
                     g(i) = j2(i0009+i0311+i-1)                         
                                                                enddo   
        if(i3108.ne.3.and.i3108.ne.4)then                               
                                      if(j2(i0009+i432i11).le.j13)then  
                       i3108 = 1                                        
                                    else                                
                              i3108 = 2                                 
                                                              endif     
       endif                                                            
                                         i4108 = 1                      
   86                                                       continue    
                                        DO I = 1,N                      
                                                IF( X(I).GT.j9(I) )THEN 
                                                      j7(I) = 91.0D0    
       GOTO 87                                                          
         ENDIF                                                          
            IF( X(I).LT.j7(I) )THEN                                     
                            j7(I) = 92.0D0                        
             GOTO 87                                                    
                                          ENDIF                         
                                            IF( j7(I).GT.j9(I) )THEN    
                                                  j7(I) = 93.0D0        
                        GOTO 87                                         
                                                                ENDIF   
                IF( j7(I).EQ.j9(I) )THEN                                
                              j7(I) = 90.0D0                            
                                            GOTO 87                     
                                                        ENDIF           
      IF( DABS(X(I)-j7(I)) .LE. (j9(I)-j7(I))/1000.0D0 )THEN            
                                        j7(I) = 0.0D0                   
                                              GOTO 87                   
              ENDIF                                                     
         IF( DABS(X(I)-j9(I)) .LE. (j9(I)-j7(I))/1000.0D0 )THEN         
                                                          j7(I) = 22.0D0
                                                       GOTO 87          
                                                     ENDIF              
             DO J = 1,21                                                
                  IF( X(I) .LE. j7(I) + J * (j9(I)-j7(I))/21.0D0)THEN   
                                    j7(I) = DBLE(J)                     
                                      GOTO 87                           
              ENDIF                                              
                                                   ENDDO                
   87                                CONTINUE                           
                  ENDDO                                                 
                                               IF(i3108.GT.100) goto 4  
                          j9(1) =-0.123984000D+06      
                                              DO I = 2,N                
                              j9(I) = j9(1) * 2.0D0 * o89(j2(1))        
                     ENDDO                                              
                     goto 4                                             
                                              end                       
                                subroutine o1309(n,nint,x,j7,j9,j2,lj2) 
                            implicit none                               
                          integer n,nint,lj2,i                          
             double precision x,j7,j9,j2,o89                            
              dimension x(n),j7(n),j9(n),j2(lj2)                        
      do i=1,n                                                          
                     x(i) = j7(i) + o89(j2(1)) * (j9(i)-j7(i))        
                          if (i.gt.n-nint) x(i) = dnint(x(i))           
                                                   enddo                
                                                           end          
        subroutine o95420(n,nint,k,x,j7,j9,j2,lj2,i37,i8087,pt,U1)      
         implicit none                                                  
        integer n,nint,k,lj2,i37,i8087,pt,i,j                           
                double precision j2,x,j7,j9,o89,i82431,U1,i9042677836   
                                 dimension j2(lj2),x(n),j7(n),j9(n)     
                                                 do i=1,n               
                                                 i82431 = o89(j2(1))    
                                             do j=1,k                   
                                           if(i82431.gt.j2(pt+j-1))then 
                         continue                                       
                             else                                       
                                                             goto 456   
                endif                                                   
                              enddo                                     
  456                     x(i) = j2(i37+(j-2)*n+i-1) + j2(i8087+i-1) *  
     &   i9042677836(i82431,o89(j2(1)))                          
                                        if(x(i).lt.j7(i))then           
       x(i) = j7(i) + (j7(i)-x(i)) / U1                                 
                         if(x(i).lt.j7(i)) x(i) = j7(i)                 
           if(x(i).gt.j9(i)) x(i) = j9(i)                               
                                         endif                          
                                           if(x(i).gt.j9(i))then        
       x(i) = j9(i) - (x(i)-j9(i)) / U1                                 
             if(x(i).lt.j7(i)) x(i) = j7(i)                             
        if(x(i).gt.j9(i)) x(i) = j9(i)                                  
           endif                                                        
              if(i.gt.n-nint) x(i) = dnint(x(i))                        
       enddo                                                            
                  end                                                   
              subroutine o32481(n,nint,x,j7,j9,j2,lj2,i86x,i13509)      
                       implicit none                                    
                                 integer n,nint,lj2,i86x,i              
           double precision x,j7,j9,j2,i13509,i130,o89,i9042677836      
                         dimension x(n),j7(n),j9(n),j2(lj2)             
      do i=1,n                                                          
                                       i130 = (j9(i)-j7(i)) / i13509    
        if(i.gt.n-nint)then                                         
        if(i130.lt.1.0d0/dsqrt(i13509))then                             
                                         i130 = 1.0d0 / dsqrt(i13509)   
                endif                                                   
                     endif                                              
               x(i) = j2(i86x+i-1) + i130 *                             
     & i9042677836(o89(j2(1)),o89(j2(1)))                               
                       if(x(i).lt.j7(i))then                            
                  x(i) = j7(i) + (j7(i)-x(i)) / (o89(j2(1))*1.0d6)      
                        if(x(i).lt.j7(i)) x(i) = j7(i)                  
         if(x(i).gt.j9(i)) x(i) = j9(i)                                 
                              endif                                     
                        if(x(i).gt.j9(i))then                           
                x(i) = j9(i) - (x(i)-j9(i)) / (o89(j2(1))*1.0d6)        
                                       if(x(i).lt.j7(i)) x(i) = j7(i)   
                                 if(x(i).gt.j9(i)) x(i) = j9(i)         
                                      endif                             
        if(i.gt.n-nint) x(i) = dnint(x(i))                              
                           enddo                                  
                  end                                                   
                  subroutine o158(l,n,nint,x,j7,j9,i37,i8087,i5308,     
     &                     j2,lj2,j6,lj6,i3108,i9025120,DA,DB,CA)   
                   implicit none                                        
       integer l,n,nint,i37,i8087,lj2,j6,lj6,i3108,i9025120,i,j,i4102, 
     &           i209081,i21,i8250,i507,i0814,i13408,i325,i86,o13218,   
     &                          i41096533,i5308,i70319301               
                    double precision x,j7,j9,j2,o89,i130,DA,DB,CA       
       dimension x(n),j7(n),j9(n),j2(lj2),j6(lj6)                       
      data i21,i8250,i507,i209081,i4102,i0814,i13408,i325,i86,i70319301
     &     /0,0,0,0,0,0,0,0,0,0/                                        
                   data i130 /0.0d0/                                    
                                     if(i3108.eq.-30)then               
                                i21 = 0                                 
                                               i8250 = 31 + l + 1       
                     i209081 =  937           
                                    do i = 1,n                          
                         j             = o13218(i*o89(j2(1))) + 1       
                     j6(i8250+i-1) = j6(i8250+j-1)                      
                                                       j6(i8250+j-1) = i
               enddo                                                    
                                             i4102 = j6(i5308+1)       
                                       do i=1,25                        
                         i4102 = i4102 + j6(i5308+2*i)                
               enddo                                                    
                               j6(31) = 1                               
                  i70319301 = i8250 + n                                 
                                            i325 = i70319301            
                 i70319301 = 1008             
           do i = 1,n                                           
                                         j6(i325+i-1) = 0               
                                    enddo                               
                                                          endif         
                                                    i41096533 = 0       
                                              if(j6(31).eq.0)then       
                                    i86 = j6(i8250+i21-1)               
                                              j6(30) = i86              
                        i0814 = i0814 + 1                               
               i507 = - i507                                            
                                               i130 = i130 /  DA        
                                             if(i130.lt. DB ) i130 = DB 
              if(i86.gt.n-nint.and.i0814.gt.i13408)then                 
                                                j6(i325+i86-1) = 1      
               if(i21.ge.n) goto 2                                      
                                                 i41096533 = 1          
              endif                                                 
                       if(i86.le.n-nint.and.i0814.gt.int(CA))then       
                                j6(i325+i86-1) = 1                      
                                                     if(i21.ge.n) goto 2
                                                     i41096533 = 1      
                             endif                                      
                     if(dabs(j7(i86)-j9(i86)).le.1.0d-12)then           
                                       j6(i325+i86-1) = 1               
      if(i21.ge.n) goto 2                                               
                                            i41096533 = 1               
                                                                  endif 
                                             endif                      
          if(j6(31).eq.1.or.i41096533.eq.1)then                         
                          i21 = i21 + 1                                 
                    if(i21.gt.n) goto 2                           
                                      i86 = j6(i8250+i21-1)             
                             j6(30) = i86                               
                                                           i0814 = 1    
       if(i86.gt.n-nint)then                                            
                        if(j2(i37+i86-1).eq.j7(i86).or.                 
     &                         j2(i37+i86-1).eq.j9(i86))then            
              i13408 = 1                                                
                                                      else              
                                                        i13408 = 2      
           endif                                                        
                    endif                                               
                                      if(o89(j2(1)).ge.0.5d0)then       
                                                       i507 = 1         
                             else                                       
                           i507 = -1                                    
                                           endif                        
                                      i130 = dsqrt(j2(i8087+i86-1))     
                                                  endif                 
                      do i = 1,n                                
                                     x(i) = j2(i37+i-1)                 
                                                     enddo              
            if(i86.le.n-nint)then                                       
                                x(i86) = x(i86) + i507 * i130           
                   else                                                 
                                x(i86) = x(i86) + i507                  
                                            if(x(i86).lt.j7(i86))then   
                                          x(i86) = j7(i86) + 1          
                                 endif                                  
         if(x(i86).gt.j9(i86))then                                      
          x(i86) = j9(i86) - 1                                          
                                  endif                                 
                                                      endif             
                          if(i209081.ne.i4102)then                     
               do i = 1,n                                         
                 j       = o13218(i*o89(j2(1))) + 1                     
                                                  j6(i+1) = j6(j+1)     
               j6(j+1) = i                                              
                      enddo                                             
                                           endif                        
            if(x(i86).lt.j7(i86)) x(i86) = j7(i86)                      
               if(x(i86).gt.j9(i86)) x(i86) = j9(i86)                   
                    if(i86.gt.n-nint) x(i86) = dnint(x(i86))            
                         if(i21.eq.1.and.i0814.eq.1)then                
                                   i3108 = -30                          
                                    else                                
                    i3108 = -31                                         
                      endif                                             
                                                              return    
    2                                  i3108 = -40                 
                             do i = 1,n                                 
                          if(j6(i325+i-1).eq.0) goto 22                 
                                                   enddo                
                                     i9025120 = 1                       
                         i3108 = -99                                    
   22                                  return                           
                                                                  end   
                                                 function o13218(x)     
                                                      implicit none     
                          double precision x                            
                              integer o13218                            
                                     o13218 = int(x)             
                        if(o13218.gt.x) o13218 = o13218 - 1             
                    end                                                 
c                       END OF FILE                                     
