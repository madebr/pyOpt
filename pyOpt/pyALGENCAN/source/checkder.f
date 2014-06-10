C     *****************************************************************
C     *****************************************************************

      subroutine checkd(n,l,u,m,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision l(n),u(n)

C     This subrotutine checks the user supplied first and second
C     derivatives subroutines (evalg, evalh, evaljac and evalhc) for
C     computing the objective function gradient and Hessian and the
C     constraints gradients and Hessians, respectively.

      include "dim.par"
      include "algparam.inc"

C     LOCAL SCALARS
      character answer
      integer i,j
      double precision drand,seed,smalll,smallu

C     LOCAL ARRAYS
      double precision x(nmax)

C     EXTERNAL FUNCTIONS
      external drand

C     SET A RANDOM POINT

      seed = 123456.0d0
      do i = 1,n
          smalll = max( l(i), - 10.0d0 )
          smallu = min( u(i),   10.0d0 )
          if ( .not. smalll .lt. smallu ) then
              smalll = l(i)
              smallu = u(i)
          end if
          x(i) = smalll + ( smallu - smalll ) * drand(seed)
      end do

C      write(* ,100)
      write(10,100)

      do i = 1,n
C          write(* ,110) i,x(i)
          write(10,110) i,x(i)
      end do

C     CHECK OBJECTIVE FUNCTION GRADIENT

      if ( .not. gcoded ) then
C          write(* ,120) 'evalg'
          write(10,120) 'evalg'

          go to 1000
      end if

C      write(* ,130) 'evalg'
      write(10,130) 'evalg'

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
          go to 1000

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkg(n,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK JACOBIAN OF CONSTRAINTS

 1000 continue

      if ( .not. jaccoded ) then
C          write(* ,120) 'evaljac'
          write(10,120) 'evaljac'

          go to 1020
      end if

C      write(* ,130) 'evaljac'
      write(10,130) 'evaljac'

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
          go to 1020

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then

          j = 1

 1010     if ( j .le. m ) then

C              write(* ,220) j
              write(10,220) j

              read(*,*) answer

              if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
                  return

              else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
                  go to 1020

              else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
                  call checkjac(n,x,j,inform)
                  if ( inform .lt. 0 ) return
              end if

              j = j + 1

              go to 1010

          end if

      end if

C     CHECK HESSIAN OF THE OBJECTIVE FUNCTION

 1020 continue

      if ( .not. hcoded ) then
C          write(* ,120) 'evalh'
          write(10,120) 'evalh'

          go to 1030
      end if

C      write(* ,130) 'evalh'
      write(10,130) 'evalh'

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
          go to 1030

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkh(n,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK HESSIANS OF THE CONSTRAINTS

 1030 continue

      if ( .not. hccoded ) then
C          write(* ,120) 'evalhc'
          write(10,120) 'evalhc'

          go to 1050
      end if

C      write(* ,130) 'evalhc'
      write(10,130) 'evalhc'

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
          go to 1050

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then

          j = 1

 1040     if ( j .le. m ) then

C              write(* ,230) j
              write(10,230) j

              read(*,*) answer

              if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
                  return

              else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
                  go to 1050

              else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
                  call checkhc(n,x,j,inform)
                  if ( inform .lt. 0 ) return
              end if

              j = j + 1

              go to 1040

          end if

      end if

C     CHECK GRADIENT OF OBJECTIVE FUNCTION PLUS JACOBIAN OF CONSTRAINTS

 1050 continue

      if ( .not. gjaccoded ) then
C          write(* ,120) 'evalgjac'
          write(10,120) 'evalgjac'

          go to 1060
      end if

C      write(* ,130) 'evalgjac'
      write(10,130) 'evalgjac'

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
          go to 1060

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkgjac(n,m,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK GRADIENT OF OBJECTIVE FUNCTION PLUS PRODUCT OF 
C     JACOBIAN OF CONSTRAINTS TIMES A GIVEN VECTOR

 1060 continue

      if ( .not. gjacpcoded ) then
C          write(* ,120) 'evalgjacp'
          write(10,120) 'evalgjacp'

          go to 1070
      end if

C      write(* ,130) 'evalgjacp'
      write(10,130) 'evalgjacp'

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
          go to 1070

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
C          write(* ,*) 'Test of evalgjacp not implemented yet!'
          write(10,*) 'Test of evalgjacp not implemented yet!'
      end if

C     CHECK HESSIAN OF THE LAGRANGIAN

 1070 continue

      if ( .not. hlcoded ) then
C          write(* ,120) 'evalhl'
          write(10,120) 'evalhl'

          go to 1080
      end if

C      write(* ,130) 'evalhl'
      write(10,130) 'evalhl'

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
          go to 1080

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
C          write(* ,*) 'Test of evalhl not implemented yet!'
          write(10,*) 'Test of evalhl not implemented yet!'
      end if

C     CHECK HESSIAN OF THE LAGRANGIAN TIMES A VECTOR

 1080 continue

      if ( .not. hlpcoded ) then
C          write(* ,120) 'evalhlp'
          write(10,120) 'evalhlp'

          go to 1090
      end if

C      write(* ,130) 'evalhlp'
      write(10,130) 'evalhlp'

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
          go to 1090

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
C          write(* ,*) 'Test of evalhlp not implemented yet!'
          write(10,*) 'Test of evalhlp not implemented yet!'
      end if

C     NON-EXECUTABLE STATEMENTS

 1090 continue

 100  format(/,1X,'Derivatives will be tested at the random point: ')
 110  format(  1X,'x(',I6,') = ',1P,D15.8)
 120  format(/,1X,'Skipping checking of uncoded ',A9,' subroutine.')
 130  format(/,1X,'Would you like to ckeck subroutine ',A9,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip checking this subroutine): ')

 200  format(/,1X,'Check gradient of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 210  format(/,1X,'Check Hessian matrix of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')

 220  format(/,1X,'Check gradient of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip gradients of constraints): ')
 230  format(/,1X,'Check Hessian matrix of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip Hessians of constraints): ')

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkg(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalg for
C     computing the gradient of the objective function using central
C     finite differences with two different discretization steps.

      include "dim.par"
      include "machconst.inc"

C     LOCAL SCALARS
      integer i
      double precision fminus,fplus,gdiff1,gdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      double precision g(nmax)

      call vsetp(n,x)
      call vevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

C      write(* ,100)
      write(10,100)

      maxerr = 0.0d0

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vsetp(n,x)
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vsetp(n,x)
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vsetp(n,x)
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vsetp(n,x)
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )

C          write(* ,110) i,g(i),gdiff1,gdiff2,tmp
          write(10,110) i,g(i),gdiff1,gdiff2,tmp

          maxerr = max( maxerr, tmp )

      end do

C      write(* ,120) maxerr
      write(10,120) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of the objective function.',
     +       /,1X,'Index',13X,'evalg',2X,'Central diff (two different ',
     +            'steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkh(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalh for
C     computing the Hessian of the objective function using central
C     finite differences with two different discretization steps.

      include "dim.par"
      include "machconst.inc"

C     LOCAL SCALARS
      logical nullcol
      integer i,j,hnnz
      double precision elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2)
      double precision g(nsmax),gplus1(nsmax),gplus2(nsmax),
     +        H(nsmax,nsmax),hval(nsmax**2),maxcoe(nsmax)

C     Check viability of the test

      if ( n .gt. nsmax ) then
C          write(*, 100) nsmax,nsmax
          write(10,100) nsmax,nsmax

          return
      end if

C     Compute the gradient of the objective function at x

      call vsetp(n,x)
      call vevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

C     Compute the Hessian of the objective function at x and save in a
C     dense matrix

      call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,hnnz
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

C     Test column by column

C      write(* ,200)
      write(10,200)

      maxerr = 0.0d0

      do j = 1,n

          tmp  = x(j)

          step1 = macheps12 * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call vsetp(n,x)
          call vevalg(n,x,gplus1,inform)
          if ( inform .lt. 0 ) return

          step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call vsetp(n,x)
          call vevalg(n,x,gplus2,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp

C          write(* ,210) j
          write(10,210) j

          maxcoe(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )
              if ( elem   .ne. 0.0d0 .or.
     +             hdiff1 .ne. 0.0d0 .or.
     +             hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
C                      write(* ,220)
                      write(10,220)
                  end if
C                  write(* ,230) i,elem,hdiff1,hdiff2,tmp
                  write(10,230) i,elem,hdiff1,hdiff2,tmp
              end if
              maxcoe(j) = max( maxcoe(j), tmp )
          end do

          maxerr = max( maxerr, maxcoe(j) )

          if ( nullcol ) then
C              write(* ,240)
              write(10,240)
          else
C              write(* ,250) maxcoe(j)
              write(10,250) maxcoe(j)
          end if

      end do

C      write(* ,*)
      write(10,*)

      do j = 1,n
C          write(* ,260) j,maxcoe(j)
          write(10,260) j,maxcoe(j)
      end do

C      write(* ,270) maxerr
      write(10,270) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Subroutine CHECKH uses dense matrices up to ',
     +            'dimension ',I6,' times ',I6,'. The Hessian ',
     +            'checking will be skipped.')

 200  format(/,1X,'Hessian matrix of the objective function column by ',
     +            'column.')
 210  format(/,1X,'Column:  ',I6)
 220  format(/,1X,'Index',13X,'evalh',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 230  format(  1X,I5,4(3X,1P,D15.8))
 240  format(  1X,'All the elements of this column are null.')
 250  format(  1X,'Maximum absolute error = ',1P,D15.8)
 260  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
 270  format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkjac(n,x,ind,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evaljac for
C     computing the gradients of the constraints using central finite
C     differences with two different discretization steps.

      include "dim.par"
      include "machconst.inc"

C     LOCAL SCALARS
      logical nullcol
      integer i,jcnnz
      double precision cminus,cplus,jacdiff1,jacdiff2,maxerr,step1,
     +        step2,tmp

C     LOCAL ARRAYS
      integer jcvar(nmax)
      double precision g(nmax),jcval(nmax)

C     COMPUTE THE GRADIENT OF THE CONSTRAINT AND SAVE IT INTO A DENSE
C     VECTOR

      call vsetp(n,x)
      call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,jcnnz
          g(jcvar(i)) = g(jcvar(i)) + jcval(i)
      end do

C     COMPARE WITH CENTRAL FINITE DIFFERENCES

C      write(* ,100) ind
      write(10,100) ind

      maxerr = 0.0d0

      nullcol = .true.

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vsetp(n,x)
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vsetp(n,x)
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          jacdiff1 = ( cplus - cminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vsetp(n,x)
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vsetp(n,x)
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          jacdiff2 = ( cplus - cminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - jacdiff1 ), abs( g(i) - jacdiff2 ) )

          if ( g(i)     .ne. 0.0d0 .or.
     +         jacdiff1 .ne. 0.0d0 .or.
     +         jacdiff2 .ne. 0.0d0 ) then
              if ( nullcol ) then
                  nullcol = .false.
C                  write(* ,110)
                  write(10,110)
              end if
C              write(* ,120) i,g(i),jacdiff1,jacdiff2,tmp
              write(10,120) i,g(i),jacdiff1,jacdiff2,tmp
          end if

          maxerr = max( maxerr, tmp )
      end do

      if ( nullcol ) then
C          write(* ,130)
          write(10,130)
      else
C          write(* ,140) maxerr
          write(10,140) maxerr
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of constraints ',I5,'.')
 110  format(/,1X,'Index',11X,'evaljac',2X,'Central diff (two ',
     +            'different steps)',4X,'Absolute error')
 120  format(  1X,I5,4(3X,1P,D15.8))
 130  format(  1X,'All the elements of this gradient are null.')
 140  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************
      subroutine checkhc(n,x,ind,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalhc for
C     computing the Hessians of the constraints using finite
C     differences.

      include "dim.par"
      include "machconst.inc"

C     LOCAL SCALARS
      logical nullcol
      integer i,j,hnnz,jcnnz
      double precision elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2),jcvar(nsmax)
      double precision g(nsmax),gplus1(nsmax),gplus2(nsmax),
     +        H(nsmax,nsmax),hval(nsmax**2),jcval(nsmax),maxcoe(nsmax)

C     Check viability of the test

      if ( n .gt. nsmax ) then
C          write(*, 100) nsmax,nsmax
          write(10,100) nsmax,nsmax

          return
      end if

C     Compute the gradient of constraint ind at x and save it in a
C     dense vector

      call vsetp(n,x)
      call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,jcnnz
          g(jcvar(i)) = g(jcvar(i)) + jcval(i)
      end do

C     Compute the Hessian of constraint ind at x and save it in a
C     dense matrix

      call vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,hnnz
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

C      write(* ,200) ind
      write(10,200) ind

      maxerr = 0.0d0

      do j = 1,n

          tmp  = x(j)

C         Compute the gradient of constraint ind at xplus1 and
C         save in a dense vector

          step1 = macheps12 * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call vsetp(n,x)
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              gplus1(i) = 0.0d0
          end do

          do i = 1,jcnnz
              gplus1(jcvar(i)) = jcval(i)
          end do

C         Compute the gradient of constraint ind at xplus2 and
C         save in a dense vector

          step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call vsetp(n,x)
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              gplus2(i) = 0.0d0
          end do

          do i = 1,jcnnz
              gplus2(jcvar(i)) = jcval(i)
          end do

          x(j) = tmp

C          write(* ,210) j
          write(10,210) j

          maxcoe(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )

              if ( elem   .ne. 0.0d0 .or.
     +             hdiff1 .ne. 0.0d0 .or.
     +             hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
C                      write(* ,220)
                      write(10,220)
                  end if
C                  write(* ,230) i,elem,hdiff1,hdiff2,tmp
                  write(10,230) i,elem,hdiff1,hdiff2,tmp
              end if

              maxcoe(j) = max( maxcoe(j), tmp )
          end do

          maxerr = max( maxerr, maxcoe(j) )

          if ( nullcol ) then
C              write(* ,240)
              write(10,240)
          else
C              write(* ,250) maxcoe(j)
              write(10,250) maxcoe(j)
          end if

      end do

C      write(* ,*)
      write(10,*)

      do j = 1,n
C          write(* ,260) j,maxcoe(j)
          write(10,260) j,maxcoe(j)
      end do

C      write(* ,270) maxerr
      write(10,270) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Subroutine CHECKHC uses dense matrices up to ',
     +            'dimension ',I6,' times ',I6,'. The Hessian ',
     +            'checking will be skipped.')

 200  format(/,1X,'Hessian matrix of constraint ',I5,' column by ',
     +            'column.')
 210  format(/,1X,'Column:  ',I6)
 220  format(/,1X,'Index',12X,'evalhc',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 230  format(  1X,I5,4(3X,1P,D15.8))
 240  format(  1X,'All the elements of this column are null.')
 250  format(  1X,'Maximum absolute error = ',1P,D15.8)
 260  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
 270  format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkgjac(n,m,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalgjac for
C     computing the gradient of the objective function plus the Jacobian
C     of the constraints. It uses central finite differences with two 
C     different discretization steps.

      include "dim.par"
      include "machconst.inc"

C     LOCAL SCALARS
      logical nullcol
      integer i,ind,jcnnz
      double precision fminus,fplus,gdiff1,gdiff2,jacdiff1,jacdiff2,
     +        maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer jcfun(jcnnzmax),jcvar(jcnnzmax)
      double precision cminus(mmax),cplus(mmax),g(nmax),jcval(jcnnzmax)


C     COMPUTE GRADIENT OF OBJECTIVE FUNCTION AND JACOBIAN OF 
C     CONSTRAINTS AT THE CURRENT POINT

      call vsetp(n,x)
      call vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

C     CHECK GRADIENT OF THE OBJECTIVE FUNCTION

C      write(* ,100)
      write(10,100)

      maxerr = 0.0d0

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vsetp(n,x)
          call vevalfc(n,x,fplus,m,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vsetp(n,x)
          call vevalfc(n,x,fminus,m,cminus,inform)
          if ( inform .lt. 0 ) return

          gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vsetp(n,x)
          call vevalfc(n,x,fplus,m,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vsetp(n,x)
          call vevalfc(n,x,fminus,m,cminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )

C          write(* ,110) i,g(i),gdiff1,gdiff2,tmp
          write(10,110) i,g(i),gdiff1,gdiff2,tmp

          maxerr = max( maxerr, tmp )

      end do

C      write(* ,120) maxerr
      write(10,120) maxerr

C     CHECK JACOBIAN OF CONSTRAINTS

      do ind = 1,m

C         COPY GRADIENT OF IND-TH CONSTRAINT INTO A DENSE VECTOR

          do i = 1,n
              g(i) = 0.0d0
          end do

          do i = 1,jcnnz
              if ( jcfun(i) .eq. ind ) then
                  g(jcvar(i)) = g(jcvar(i)) + jcval(i)
              end if
          end do

C         COMPARE WITH CENTRAL FINITE DIFFERENCES

C          write(* ,200) ind
          write(10,200) ind

          maxerr = 0.0d0

          nullcol = .true.

          do i = 1,n
              tmp  = x(i)

              step1 = macheps13 * max( abs( tmp ), 1.0d0 )

              x(i) = tmp + step1
              call vsetp(n,x)
              call vevalfc(n,x,fplus,m,cplus,inform)
              if ( inform .lt. 0 ) return

              x(i) = tmp - step1
              call vsetp(n,x)
              call vevalfc(n,x,fminus,m,cminus,inform)
              if ( inform .lt. 0 ) return

              jacdiff1 = ( cplus(ind) - cminus(ind) ) / ( 2.0d0 * step1)

              step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

              x(i) = tmp + step2
              call vsetp(n,x)
              call vevalfc(n,x,fplus,m,cplus,inform)
              if ( inform .lt. 0 ) return

              x(i) = tmp - step2
              call vsetp(n,x)
              call vevalfc(n,x,fminus,m,cminus,inform)
              if ( inform .lt. 0 ) return

              x(i) = tmp

              jacdiff2 = ( cplus(ind) - cminus(ind) ) / ( 2.0d0 * step2)

              tmp = min( abs(g(i) - jacdiff1), abs(g(i) - jacdiff2) )

              if ( g(i)     .ne. 0.0d0 .or.
     +             jacdiff1 .ne. 0.0d0 .or.
     +             jacdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
C                      write(* ,210)
                      write(10,210)
                  end if
C                  write(* ,220) i,g(i),jacdiff1,jacdiff2,tmp
                  write(10,220) i,g(i),jacdiff1,jacdiff2,tmp
              end if

              maxerr = max( maxerr, tmp )
          end do

          if ( nullcol ) then
C              write(* ,230)
              write(10,230)
          else
C              write(* ,120) maxerr
              write(10,120) maxerr
          end if

      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of the objective function.',
     +       /,1X,'Index',13X,'evalgjac',1X,'Central diff (two ',
     +            'different steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)
 200  format(/,1X,'Gradient vector of constraints ',I5,'.')
 210  format(/,1X,'Index',11X,'evalgjac',1X,'Central diff (two ',
     +            'different steps)',4X,'Absolute error')
 220  format(  1X,I5,4(3X,1P,D15.8))
 230  format(  1X,'All the elements of this gradient are null.')

      end

