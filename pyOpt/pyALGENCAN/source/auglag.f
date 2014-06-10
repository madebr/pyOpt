C     ******************************************************************
C     ******************************************************************

      subroutine auglag(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,
     +efacc,eoacc,f,c,csupn,snorm,nl,nlpsupn,fu,cu,csupnu,fub,csupnub,
     +fb,csupnb,nlpsupnb,ncsupn,rsupn,outiter,totiter,nwcalls,nwtotit,
     +msqcalls,msqtotit,innfail,alinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical innfail
      integer alinfo,inform,m,msqcalls,msqtotit,n,nwcalls,nwtotit,
     +        outiter,totiter
      double precision csupn,csupnb,csupnu,csupnub,efacc,eoacc,epsfeas,
     +        epsopt,f,fb,fu,fub,ncsupn,nlpsupn,nlpsupnb,rsupn,snorm

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision c(m),cu(m),l(n),lambda(m),nl(n),u(n),x(n)

C     Solves the nonlinear programming problem
C
C     min f(x)
C
C     subject to
C
C             c_j(x)  = 0, j in E,
C             c_j(x) <= 0, j in I,
C             l <= x <= u,
C
C     where E is the set of indices of the equality constraints, I is
C     the set of indices of the inequality constraints, and there are
C     n variables and m constraints, using the method of multipliers
C     described in
C
C     R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
C     "On Augmented Lagrangian methods with general lower-level
C     constraints", SIAM Journal on Optimization 18, pp. 1286-1309, 2007.

C     alinfo:
C
C     0: Feasibility, optimality and complementarity satisfied
C     2: Too large penalty parameter. Infeasible problem?
C     3: Maximum number of algencan iterations reached

      include "dim.par"
      include "machconst.inc"
      include "algconst.par"
      include "algparam.inc"
      include "outtyp.inc"
      include "scaling.inc"

C     LOCAL SCALARS
      character innstp
      logical scaletmp
      integer cifacnt,geninfo,i,iter,maxit,msqiter,nwinfo,nwiter,
     +        outitnp,xind
      double precision al,csupna,epsfeas12,epsfeas14,epsopt12,epsopt14,
     +        epsopk,fa,nlpi,nlpsupna,rhoa,rhob,rhoc,rhoini,snormb,
     +        snormprev,xsupn

C     LOCAL ARRAYS
      double precision lambar(mmax),lambdaa(mmax),lambdab(mmax),
     +        rho(mmax),xa(nmax),xb(nmax),nc(nmax)

C     DATA BLOCKS
      character genstp(0:9)
      data genstp(0) /'C'/
      data genstp(1) /'M'/
      data genstp(2) /' '/
      data genstp(3) /'P'/
      data genstp(4) /' '/
      data genstp(5) /' '/
      data genstp(6) /'U'/
      data genstp(7) /' '/
      data genstp(8) /' '/
      data genstp(9) /' '/

C     ==================================================================
C     Print initial information
C     ==================================================================

      if ( iprintout .ge. 1 ) then
C          write(* ,1000) n,m
          write(10,1000) n,m

          if ( iprintout .ge. 4 .and. nprint .ne. 0 ) then
C              write(* ,1010) nprint,(l(i),i=1,nprint)
C              write(* ,1020) nprint,(u(i),i=1,nprint)

              write(10,1010) nprint,(l(i),i=1,nprint)
              write(10,1020) nprint,(u(i),i=1,nprint)
          end if

      end if

C     ==================================================================
C     Initialization
C     ==================================================================

C     Counters

      outiter  = 0
      outitnp  = 0

      totiter  = 0

      nwcalls  = 0
      nwtotit  = 0

      msqcalls = 0
      msqtotit = 0

      cifacnt  = 0
      innfail  = .false.

C     Constants

      epsopt12  = sqrt( epsopt    )
      epsopt14  = sqrt( epsopt12  )
      epsfeas12 = sqrt( epsfeas   )
      epsfeas14 = sqrt( epsfeas12 )

C     Project initial point

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Compute current point norm

      xsupn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind  = i
          end if
      end do

C     Compute objective function and constraints

      call ssetp(n,x)

      call sevalobjc(n,x,f,fu,m,c,cu,inform)
      if ( inform .lt. 0 ) return

C     Compute safeguarded Lagrange multipliers

      do i = 1,m
          lambar(i) = max( lammin, min( lambda(i), lammax ) )
      end do

C     Compute complementarity and feasibility violations

      snormprev = bignum

      csupn  = 0.0d0
      csupnu = 0.0d0
      snorm  = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              csupn  = max( csupn , abs( c(i)  ) )
              csupnu = max( csupnu, abs( cu(i) ) )
              snorm  = max( snorm , abs( c(i)  ) )
          else
              csupn  = max( csupn , c(i)  )
              csupnu = max( csupnu, cu(i) )
              snorm  = max( snorm , abs( min( - c(i), lambda(i) ) ) )
          end if
      end do

C     Compute continuous projected Lagrangian gradient norm

      call sevalnl(n,x,m,lambda,equatn,linear,nl,inform)
      if ( inform .lt. 0 ) return

      nlpsupn = 0.0d0
      do i = 1,n
          nlpi = x(i) - nl(i)
          if ( l(i) .le. nlpi .and. nlpi .le. u(i) ) then
              nlpi = - nl(i)
          else
              nlpi = max( l(i), min( nlpi, u(i) ) ) - x(i)
          end if
          nlpsupn = max( nlpsupn, abs( nlpi ) )
      end do

C     Compute squared-infeasibility gradient norm

      call stafeas(n,x,l,u,m,c,lambda,equatn,nc,ncsupn,inform)
      if ( inform .lt. 0 ) return

C     Save best solution

      fub      = fu
      csupnub  = csupnu

      fb       = f
      csupnb   = csupn
      snormb   = snorm
      nlpsupnb = nlpsupn

      do i = 1,n
          xb(i) = x(i)
      end do

      do i = 1,m
          lambdab(i) = lambda(i)
      end do

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Print information of this iteration
C     ==================================================================

      if ( iprintout .eq. 1 ) then

          if ( outiter .eq. 0 ) then
              innstp = ' '
          else
              innstp = genstp(geninfo)
          end if

          if ( .not. scale ) then
              if ( mod(outiter,10) .eq. 0 ) then
C                  write(* ,1030)
                  write(10,1030)
              end if

              if ( outiter .eq. 0 ) then
C                  write(* ,1041) outiter,      f,csupn,snorm,nlpsupn,
C     +                           xsupn,ncsupn,totiter,innstp,nwcalls,
C     +                           nwtotit
                  write(10,1041) outiter,      f,csupn,snorm,nlpsupn,
     +                           xsupn,ncsupn,totiter,innstp,nwcalls,
     +                           nwtotit
              else
C                  write(* ,1040) outiter,rsupn,f,csupn,snorm,nlpsupn,
C     +                           xsupn,ncsupn,totiter,innstp,nwcalls,
C     +                           nwtotit
                  write(10,1040) outiter,rsupn,f,csupn,snorm,nlpsupn,
     +                           xsupn,ncsupn,totiter,innstp,nwcalls,
     +                           nwtotit
              end if
          else
              if ( mod(outiter,10) .eq. 0 ) then
C                  write(* ,1050)
                  write(10,1050)
              end if

              if ( outiter .eq. 0 ) then
C                  write(* ,1061) outiter,      fu,csupnu,f,csupn,snorm,
C     +                           nlpsupn,ncsupn,totiter,innstp,nwcalls,
C     +                           nwtotit
                  write(10,1061) outiter,      fu,csupnu,f,csupn,snorm,
     +                           nlpsupn,ncsupn,totiter,innstp,nwcalls,
     +                           nwtotit
              else
C                  write(* ,1060) outiter,rsupn,fu,csupnu,f,csupn,snorm,
C     +                           nlpsupn,ncsupn,totiter,innstp,nwcalls,
C     +                           nwtotit
                  write(10,1060) outiter,rsupn,fu,csupnu,f,csupn,snorm,
     +                           nlpsupn,ncsupn,totiter,innstp,nwcalls,
     +                           nwtotit
              end if
          end if

      else if ( iprintout .ge. 2 ) then
          if ( outiter .eq. 0 ) then
C              write(* ,1071) outiter
              write(10,1071) outiter
          else
C              write(* ,1070) outiter,rsupn
              write(10,1070) outiter,rsupn
          end if

          if ( iprintout .ge. 3 ) then
C              write(* ,1080) totiter,nwcalls,nwtotit
              write(10,1080) totiter,nwcalls,nwtotit
          end if

          if ( .not. scale ) then
C              write(* ,1090) f,csupn
              write(10,1090) f,csupn
          else
C              write(* ,1100) f,fu,csupn,csupnu
              write(10,1100) f,fu,csupn,csupnu
          end if

C          write(* ,1110) snorm,nlpsupn,ncsupn,xind,xsupn
          write(10,1110) snorm,nlpsupn,ncsupn,xind,xsupn

          if ( iprintout .ge. 4 ) then
              if ( nprint .ne. 0 ) then
C                  write(* ,1120) nprint,(x(i),i=1,nprint)
                  write(10,1120) nprint,(x(i),i=1,nprint)
              end if

              if ( mprint .ne. 0 ) then
C                  write(* ,1130) mprint,(lambda(i),i=1,mprint)
                  write(10,1130) mprint,(lambda(i),i=1,mprint)

                  if ( .not. scale ) then
C                      write(* ,1150) mprint,(c(i),i=1,mprint)
                      write(10,1150) mprint,(c(i),i=1,mprint)
                  else
C                      write(* ,1160) mprint,(cu(i),i=1,mprint)
                      write(10,1160) mprint,(cu(i),i=1,mprint)
C                      write(* ,1170) mprint,(c(i),i=1,mprint)
                      write(10,1170) mprint,(c(i),i=1,mprint)
                  end if
              end if
          end if

      end if

C     Save intermediate data for crash report

      if ( iprintctl(6) ) then
          open(20,file='algencan-tabline.out')
          write(20,1400) fu,csupnu,f,csupn,nlpsupn,fub,csupnub,fb,
     +                   csupnb,nlpsupnb,ncsupn,rsupn,inform,9,innfail,
     +                   n,m,outiter,totiter,0,nwcalls,nwtotit,msqcalls,
     +                   msqtotit,999.9d0
          close(20)
      end if

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     Test feasibility, optimality and complementarity

      if ( max( snorm, csupnu ) .le. epsfeas .and.
     +     nlpsupn .le. epsopt ) then
          alinfo = 0

          if ( iprintout .ge. 1 ) then
C              write(*, 1300)
              write(10,1300)
          end if

          return
      end if

C     Test whether we are at a stationary of the infeasibility 
C     infeasible point

c      if ( csupn .gt. epsfeas .and. ncsupn .le. epsopt ** 1.5d0  ) then

c          alinfo = 1

c          if ( iprintout .ge. 1 ) then
c              write(*, 1310)
c              write(10,1310)
c          end if

c          return
c      end if

C     Test whether the penalty parameter is too large

      if ( outiter .gt. 0 .and. rsupn .gt. rhomax ) then

          alinfo = 2

          if ( iprintout .ge. 1 ) then
C              write(*, 1320)
              write(10,1320)
          end if

          return
      end if

C     Test whether the number of iterations is exhausted

      if ( outiter .ge. maxoutit ) then
          alinfo = 3

          if ( iprintout .ge. 1 ) then
C              write(*, 1330)
              write(10,1330)
          end if

          return
      end if

C     ==================================================================
C     Near the solution, try to solve the KKT system by Newton's method
C     ==================================================================

      if ( seconde .and. .not. skipacc .and. 
     +     .not. ( efacc .lt. 0.0d0 .and. eoacc .lt. 0.0d0 ) .and. 
     +     ( nwcalls .gt. 0 .or.
     +     ( snorm .le. efacc     .and. nlpsupn .le. eoacc    ) .or.
     +     ( snorm .le. epsfeas12 .and. nlpsupn .le. epsopt12 ) .or.
     +     ( snorm .le. epsfeas14 .and. nlpsupn .le. epsopt14 .and.
     +       outiter .gt. 0 .and. geninfo .ne. 0 ) ) ) then

          nwcalls = nwcalls + 1

          do i = 1,n
              xa(i) = x(i)
          end do

          do i = 1,m
              lambdaa(i) = lambar(i)
          end do

          scaletmp = scale

          if ( scale ) then
              do i = 1,m
                  lambdaa(i) = lambdaa(i) * sc(i) / sf
              end do

              scale = .false.
          end if

          call newtonkkt(n,xa,l,u,m,lambdaa,equatn,linear,epsfeas,
     +    epsopt,fa,csupna,nlpsupna,nwiter,msqiter,nwinfo,inform)

          nwtotit = nwtotit  + nwiter

          if ( nwinfo .ge. 6 ) then
              skipacc = .true.
          end if

          if ( msqiter .gt. 0 ) then
              msqcalls = msqcalls + 1
              msqtotit = msqtotit + msqiter
          end if

          if ( inform .lt. 0 ) return

C         Save best solution

          if ( ( csupnub .gt. epsfeas .and. csupna .le. csupnub ) .or.
     +         ( csupnub .le. epsfeas .and. csupna .le. epsfeas .and.
     +           fa .le. fub ) ) then

              fub      = fa
              csupnub  = csupna

              fb       = fa
              csupnb   = csupna
              snormb   = 0.0d0
              nlpsupnb = nlpsupna

              do i = 1,n
                  xb(i) = xa(i)
              end do

              do i = 1,m
                  lambdab(i) = lambdaa(i)
              end do

          end if

C         Test feasibility, optimality and complementarity

          if ( csupna .le. epsfeas .and. nlpsupna .le. epsopt ) then

              fu      = fa
              csupnu  = csupna

              f       = fa
              csupn   = csupna
              snorm   = 0.0d0
              nlpsupn = nlpsupna

              do i = 1,n
                  x(i) = xa(i)
              end do

              do i = 1,m
                  lambda(i) = lambdaa(i)
              end do

              alinfo = 0

              if ( iprintout .ge. 1 ) then
C                  write(*, 1300)
                  write(10,1300)
              end if

              return

          end if

          scale = scaletmp

      end if

C     ==================================================================
C     Iteration
C     ==================================================================

      outiter = outiter + 1

C     ==================================================================
C     Set penalty parameter
C     ==================================================================

C     Set (or update) penalty parameters and compute maximum rho

      if ( outiter .eq. 1 .or. outiter .eq. 2 ) then

          call comprhoini(f,m,c,equatn,rhoini)

          rhoini = max( rhoinimin, min( rhoini, rhoinimax ) )

          do i = 1,m
              rho(i) = rhoini
          end do

          if ( iprintout .ge. 5 ) then
C              write(*, 1200)
              write(10,1200)
          end if

          rhoa = rhoinimin
          rhob = rhoinimax
          rhoc = rhoinimin

      else if ( max( snorm, csupnu ) .gt. epsfeas ) then

          if ( snorm .gt. rhofrac * snormprev ) then

              do i = 1,m
                  rho(i) = max( rhomult * rho(i), rhoc )
              end do

              if ( iprintout .ge. 5 ) then
C                  write(*, 1220) rhomult
                  write(10,1220) rhomult
              end if

          else

              if ( iprintout .ge. 5 ) then
C                  write(*, 1210)
                  write(10,1210)
              end if

          end if

      else if ( rhorestart ) then
 
          if ( cifacnt .ge. 2 ) then

              rhoa = min( rhoa * rhomult, 1.0d0 )
              rhob = max( rhob / rhomult, 1.0d0 )
              rhoc = rhoc * rhomult

              call comprhoini(f,m,c,equatn,rhoini)

              rhoini = max( rhoa, min( rhoini, rhob ) )

              do i = 1,m
                  rho(i) = min( rhoini, rho(i) )
              end do

              if ( iprintout .ge. 5 ) then
C                  write(*, 1230)
                  write(10,1230)
              end if

          else

              if ( iprintout .ge. 5 ) then
C                  write(*, 1210)
                  write(10,1210)
              end if

          end if

      end if

C     Compute maximum rho

      rsupn = 0.0d0
      do i = 1,m
          rsupn = max( rsupn, rho(i) )
      end do

C     ==================================================================
C     Solve the augmented Lagrangian subproblem
C     ==================================================================

C     Set optimality requeriment for the subproblem

      if ( outiter .eq. 1 ) then
          epsopk = sqrt( epsopt )
      else if ( snorm .le. epsfeas12 .and. nlpsupn .le. epsopt12 ) then
          epsopk = min( rhofrac * nlpsupn, 0.1d0 * epsopk )
          epsopk = max( epsopk, epsopt )
      end if

      if ( outiter .eq. 1 ) then
          maxit =   10
      else
          maxit = 1000
      end if

C     Call the inner-solver

      call gencan(n,x,l,u,m,lambar,equatn,linear,rho,epsfeas,epsopk,
     +maxit,iter,al,nl,nlpsupn,csupn,csupnu,geninfo,inform)

      totiter = totiter + iter

      if ( outiter .eq. 1 ) then
          geninfo = 0
      end if

      if ( geninfo .ne. 0 ) then
          innfail = .true.
      end if

      if ( inform .lt. 0 ) return

C     ==================================================================
C     Prepare for the next iteration
C     ==================================================================

C     Compute current point norm and its difference with the previous one

      xsupn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind  = i
          end if
      end do

C     Compute objective function and constraints

      call ssetp(n,x)

      call sevalobjc(n,x,f,fu,m,c,cu,inform)
      if ( inform .lt. 0 ) return

C     Compute feasibility violation

      csupn  = 0.0d0
      csupnu = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              csupn  = max( csupn , abs( c(i)  ) )
              csupnu = max( csupnu, abs( cu(i) ) )
          else
              csupn  = max( csupn , c(i)  )
              csupnu = max( csupnu, cu(i) )
          end if
      end do

C     Update Lagrange multipliers approximation

      do i = 1,m
          call evaldpdy(c(i),rho(i),lambar(i),equatn(i),lambda(i))
          lambar(i) = max( lammin, min( lambda(i), lammax ) )
      end do

C     Compute complementarity violation

      snormprev = snorm

      snorm = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              snorm = max( snorm, abs( c(i) ) )
          else
              snorm = max( snorm, abs( min( - c(i), lambda(i) ) ) )
          end if
      end do

C     Track consecutve failures of the inner solver at feasible points

      if ( max( snorm, csupnu ) .le. epsfeas .and. geninfo .ne. 0 ) then
          cifacnt = cifacnt + 1
      else
          cifacnt = 0
      end if
 
C     Compute continuous projected Lagrangian gradient norm

      call sevalnl(n,x,m,lambda,equatn,linear,nl,inform)
      if ( inform .lt. 0 ) return

      nlpsupn = 0.0d0
      do i = 1,n
          nlpi = x(i) - nl(i)
          if ( l(i) .le. nlpi .and. nlpi .le. u(i) ) then
              nlpi = - nl(i)
          else
              nlpi = max( l(i), min( nlpi, u(i) ) ) - x(i)
          end if
          nlpsupn = max( nlpsupn, abs( nlpi ) )
      end do

C     Compute squared-infeasibility gradient norm

      call stafeas(n,x,l,u,m,c,lambda,equatn,nc,ncsupn,inform)
      if ( inform .lt. 0 ) return

C     Save best solution

      if ( ( csupnub .gt. epsfeas .and. csupnu .le. csupnub ) .or.
     +     ( csupnub .le. epsfeas .and. csupnu .le. epsfeas .and.
     +       f .lt. fb ) ) then

          fub      = fu
          csupnub  = csupnu

          fb       = f
          csupnb   = csupn
          snormb   = snorm
          nlpsupnb = nlpsupn

          do i = 1,n
              xb(i) = x(i)
          end do

          do i = 1,m
              lambdab(i) = lambda(i)
          end do

      end if

C     ==================================================================
C     Iterate
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 1000 format(/,1X,'Entry to ALGENCAN.',
     +       /,1X,'Number of variables  : ',I7,
     +       /,1X,'Number of constraints: ',I7)
 1010 format(/,1X,'Lower bounds  (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))
 1020 format(/,1X,'Upper bounds  (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))

 1030 format(/,'out',1X,'penalt',2X,'objective',1X,'infeas',1X,'infeas',
     +         1X,'norm',3X,'norm',3X,'|Grad|',2X,'inner',1X,'Newton',
     +       /,'ite',9X,'function',2X,'ibilty',1X,'+compl',1X,'graLag',
     +         1X,'point',2X,'infeas',2X,'totit',1X,'forKKT')

 1040 format(I3,1X,1P,D6.0,1X,1P,D10.3,5(1X,1P,D6.0),1X,I5,A1,1X,I2,1X,
     +       I3)
 1041 format(I3,7X,        1X,1P,D10.3,5(1X,1P,D6.0),1X,I5,A1,1X,I2,1X,
     +       I3)

 1050 format(/,'out',1X,'penalt',2X,'objective',1X,'infeas',2X,'scaled',
     +         4X,'scaled',1X,'infeas',1X,'norm',3X,'|Grad|',1X,'inner',
     +         1X,'Newton',/,'ite',9X,'function',2X,'ibilty',2X,
     +         'obj-funct',1X,'infeas',1X,'+compl',1X,'graLag',1X,
     +         'infeas',1X,'totit',1X,'forKKT')

 1060 format(I3,1X,1P,D6.0,1X,1P,D10.3,1X,1P,D6.0,1X,1P,D10.3,4(1X,1P,
     +       D6.0),I5,A1,1X,I2,1X,I3)
 1061 format(I3,7X,        1X,1P,D10.3,1X,1P,D6.0,1X,1P,D10.3,4(1X,1P,
     +       D6.0),I5,A1,1X,I2,1X,I3)

 1070 format(/,1X,'ALGENCAN OUTER ITERATION                         = ',
     +         I7,
     +       /,1X,'Penalty parameter                                = ',
     +         1P,D7.1)

 1071 format(/,1X,'ALGENCAN OUTER ITERATION                         = ',
     +         I7)

 1080 format(/,1X,'Up-to-now total number of iterations             = ',
     +         I7,
     +       /,1X,'Up-to-now acceleration trials                    = ',
     +         I7,
     +       /,1X,'Up-to-now total number of Newton iterations      = ',
     +         I7)

 1090 format(/,1X,'Functional value                                 = ',
     +             1P,D24.16,
     +       /,1X,'Sup-norm of infeasibility                        = ',
     +             17X,1P,D7.1)

 1100 format(/,1X,'Functional value             (scaled = ',1P,D8.1,
     +            ') = ',1P,D24.16,
     +       /,1X,'Sup-norm of infeasibility    (scaled = ',1P,D8.1,
     +            ') = ',17X,1P,D7.1)

 1110 format(  1X,'Sup-norm of scaled infeasibility-complementarity = ',
     +             17X,1P,D7.1,
     +       /,1X,'Sup-norm of scaled Lagrangian projected gradient = ',
     +             17X,1P,D7.1,
     +       /,1X,'Sup-norm of scaled infeasibilities gradient      = ',
     +             17X,1P,D7.1,
     +       /,1X,'Sup-norm of x (attained at x_{', I7,'})          = ',
     +             17X,1P,D7.1)

 1120 format(/,1X,'Current point (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))
 1130 format(/,1X,'Current Lagrange multipliers (first ',I7,
     +         1X,'components):',
     +       /,6(1X,1P,D11.4))
 1150 format(/,1X,'Infeasibility (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))
 1160 format(/,1X,'Unscaled infeasibility (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))
 1170 format(/,1X,'Scaled infeasibility (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))

 1200 format(/,1X,'Penalty parameter was computed from scratch as it ',
     +            'is usual in the first two',/,1X,'iterations.')
 1210 format(/,1X,'The desired infeasibility improvement was achieved.',
     +       /,1X,'The penalty parameter was not modified.')
 1220 format(/,1X,'The desired infeasibility improvement was not ',
     +            'achieved.',
     +       /,1X,'The penalty parameter was increased by multiplying ',
     +            'it by rhomult = ',1PD11.4,'.')
 1230 format(/,1X,'Penalty parameter was computed from scratch ',
     +            'because its value seems to be preventing ',
     +            'convergence of the inner solver.')

 1300 format(/,1X,'Flag of ALGENCAN: Solution was found.')
 1310 format(/,1X,'Flag of ALGENCAN: It seems that a stationary-of-',
     +            'the-infeasibility probably',/,1X,'infeasible point ',
     +            'was found. Whether the final iterate is a solution ',
     +            'or not',/,1X,'requires further analysis.')
 1320 format(/,1X,'Flag of ALGENCAN: The penalty parameter is too ',
     +            'large.',/,1X,'The problem may be infeasibile or ',
     +            'bad scaled. Further analysis is required.')
 1330 format(/,1X,'Flag of ALGENCAN: Maximum of iterations reached. ',
     +            'The feasibility-',/,1X,'complementarity and ',
     +            'optimality tolerances could not be achieved. ',
     +            'Whether the',/,1X,'final iterate is a solution or ',
     +            'not requires further analysis.')

 1400 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D7.1,1X,1P,D7.1,1X,I3,1X,I1,1X,L1,1X,I6,1X,I6,1X,I3,
     +       1X,I7,1X,I7,1X,I2,1X,I7,1X,I7,1X,I7,0P,F8.2)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine comprhoini(f,m,c,equatn,rhoini)

C     SCALAR ARGUMENTS
      integer m
      double precision f,rhoini

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision c(m)

C     Consider the Augmented Lagrangian function
C
C         al = f(x) + \sum_i P(lambda_i,c_i(x),rho),
C
C     where
C
C         P(lambda,y,rho) = y ( lambda + 0.5 rho y ),
C
C     If c_i(x) is an equality constraint or lambda + rho y > 0, and
C
C         P(lambda,y,rho) = - 0.5 lambda^2 / rho,
C
C     otherwise.
C
C     Assuming that lambda_i = 0 for all i, it is clear that
C
C         P(lambda_i,c_i(x),rho) = 0.5 rho c_i(x)^2
C
C     and that the value of  rho that balances f(x) and
C
C     \sum_i P(lambda_i,c_i(x),rho) is given by
C
C     rho = f(x) / ( 0.5 \sum_i c_i(x)^2 ).

      include "machconst.inc"

C     LOCAL SCALARS
      integer i
      double precision sumc

      sumc = 0.0d0
      do i = 1,m
          if ( equatn(i) .or. c(i) .gt. 0.0d0 ) then
              sumc = sumc + 0.5d0 * c(i) ** 2
          end if
      end do

      rhoini = 10.0d0 * max( 1.0d0, abs( f ) ) / max( 1.0d0, sumc )

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalfeas(n,x,m,equatn,csupn,csupnu,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision csupn,csupnu
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision x(n)

      include "dim.par"
      include "graddat.inc"
      include "scaling.inc"

C     LOCAL SCALARS
      integer i

C     Subroutine sevalnal is always called by GENCAN with the same point
C     before calling this subroutine. Then, constraints are computed and
C     saved.

C     Compute infeasibility

      csupn = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              csupn = max( csupn, abs( c(i) ) )
          else
              csupn = max( csupn, c(i) )
          end if
      end do

      if ( scale ) then
          csupnu = 0.0d0
          do i = 1,m
              if ( equatn(i) ) then
                  csupnu = max( csupnu, abs( c(i) / sc(i) ) )
              else
                  csupnu = max( csupnu, c(i) / sc(i) )
              end if
          end do

      else
          csupnu = csupn
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine stafeas(n,x,l,u,m,r,lambda,equatn,nc,ncsupn,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision ncsupn
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision l(n),lambda(m),nc(n),r(m),u(n),x(n)

      include "dim.par"
      include "algparam.inc"
      include "graddat.inc"

C     LOCAL SCALARS
      integer i,j,nrnnz
      double precision dum,nci

C     LOCAL ARRAYS
      integer nrvar(nmax)
      double precision nrval(nmax),p(mmax)

      if ( gjacpcoded ) then

          do j = 1,m
              if ( equatn(j) .or. r(j) .gt. 0.0d0 ) then
                  p(j) = r(j)
              else
                  p(j) = 0.0d0
              end if
          end do

          call sevalgjacp(n,x,dum,m,p,nc,'t',gotj,inform)
          if ( inform .lt. 0 ) return

      else if ( gjaccoded ) then

          do i = 1,n
              nc(i) = 0.0d0
          end do

          do j = 1,m
              if ( equatn(j) .or. r(j) .gt. 0.0d0 ) then
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      nc(jcvar(i)) = nc(jcvar(i)) + r(j) * jcval(i)
                  end do
              end if
          end do

      else ! if ( gcoded .and. jaccoded ) then
c     This else is the choice when gcoded and jaccoded and also when
c     first derivatives were not provided by the user. In such case,
c     sevaljac below will ultimately call subrotuines to compute 
c     derivatives by finite differences.
          do i = 1,n
              nc(i) = 0.0d0
          end do

          do j = 1,m
              if ( equatn(j) .or. r(j) .gt. 0.0d0 ) then
                  if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then
                      do i = jcsta(j),jcsta(j) + jclen(j) - 1
                          nc(jcvar(i)) = nc(jcvar(i)) + r(j) * jcval(i)
                      end do
                  else
                      call sevaljac(n,x,j,nrvar,nrval,nrnnz,inform)
                      if ( inform .lt. 0 ) return
                      do i = 1,nrnnz
                          nc(nrvar(i)) = nc(nrvar(i)) + r(j) * nrval(i)
                      end do
                  end if
              end if
          end do

      end if

      ncsupn = 0.0d0
      do i = 1,n
          nci = x(i) - nc(i)
          if ( l(i) .le. nci .and. nci .le. u(i) ) then
              nci = - nc(i)
          else
              nci = max( l(i), min( nci, u(i) ) ) - x(i)
          end if

          ncsupn = max( ncsupn, abs( nci ) )
      end do

      end
