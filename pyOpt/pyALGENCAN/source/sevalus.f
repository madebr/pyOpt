C     SCALE OBJECTIVE FUNCTION AND CONSTRAINTS

C     ******************************************************************
C     ******************************************************************

      subroutine sinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(11),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

      include "dim.par"
      include "outtyp.inc"
      include "scaling.inc"
      include "algparam.inc"
      include "machconst.inc"
      include "probdata.inc"
      include "graddat.inc"

C     LOCAL SCALARS
      integer i,fun,j,jcnnz,neq
      double precision scmin
      double precision dum(1)

C     LOCAL ARRAYS
c     integer jcfun(jcnnzmax),jcvar(jcnnzmax)
c     double precision g(nmax),jcval(jcnnzmax),p(mmax),q(nmax)
      integer jcfun(jcnnzmax)
      double precision p(mmax),q(nmax)

      neq = 0
      do j = 1,m
          if ( equatn(j) ) neq = neq + 1
      end do

      nbds = 0
      do i = 1,n
          if ( l(i) .gt. - 1.0d+20 ) nbds = nbds + 1
          if ( u(i) .lt.   1.0d+20 ) nbds = nbds + 1
      end do

      if ( iprintctl(2) ) then
C          write(* ,100) n,neq,m-neq,nbds
          write(10,100) n,neq,m-neq,nbds
      end if

      call tinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

C     Write classification line of final model

      if ( iprintctl(6) ) then
          open(50,file='class-tabline.out')
          write(50,400) n,neq,m-neq,nbds
          close(50)
      end if

C     Scaling

      usf = 1.0d0
      do j = 1,m
          usc(j) = 1.0d0
      end do

      if ( scale ) then

          if ( m .eq. 0 ) then
              sf = 1.0d0
              if ( iprintctl(2) ) then
C                  write(* ,200) sf
                  write(10,200) sf
              end if

              return
          end if

          call tsetp(n,x)

          if ( gjacpcoded ) then

C             Scaling the constraints may be avoided if computing all
C             individual gradients of constraints (even only once) is
C             a very expensive task. In this case, the code below may
C             be replaced by:
C
C             do j = 1,m
C                 sc(j) = 1.0d0
C             end do
C
C             and a single call to tevalgjacp with p=0 and work equal
C             to 'J' or 'T' may be preserved to compute the gradient
C             of the objective function g.

C             Scale constraints

              do j = 1,m
                  p(j) = 0.0d0
              end do

              do j = 1,m
                  p(j) = 1.0d0
                  if ( j .eq. 1 ) then
                      call tevalgjacp(n,x,g,m,p,q,'T',gotj,inform)
                      if ( inform .lt. 0 ) return
                  else
                      call tevalgjacp(n,x,dum,m,p,q,'t',gotj,inform)
                      if ( inform .lt. 0 ) return
                  end if
                  p(j) = 0.0d0

                  sc(j) = 1.0d0
                  do i = 1,n
                      sc(j) = max( sc(j), abs( q(i) ) )
                  end do
              end do

              do j = 1,m
                  sc(j) = 1.0d0 / sc(j)
              end do

          else if ( gjaccoded ) then

              call tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
              if ( inform .lt. 0 ) return

C             Scale constraints

              do j = 1,m
                  sc(j) = 1.0d0
              end do

              do i = 1,jcnnz
                  fun = jcfun(i)
                  sc(fun) = max( sc(fun), abs( jcval(i) ) )
              end do

              do j = 1,m
                  sc(j) = 1.0d0 / sc(j)
              end do

          else

              call tevalg(n,x,g,inform)
              if ( inform .lt. 0 ) return

C             Scale constraints

              do j = 1,m
                  call tevaljac(n,x,j,jcvar,jcval,jcnnz,inform)
                  if ( inform .lt. 0 ) return

                  sc(j) = 1.0d0
                  do i = 1,jcnnz
                      sc(j) = max( sc(j), abs( jcval(i) ) )
                  end do
                  sc(j) = 1.0d0 / sc(j)
              end do

          end if

C         Scale objective function

          sf = 1.0d0
          do i = 1,n
              sf = max( sf, abs( g(i) ) )
          end do
          sf = 1.0d0 / sf

C         Report scaling factors

          scmin = bignum
          do j = 1,m
              scmin = min( scmin, sc(j) )
          end do

          if ( iprintctl(2) ) then
C              write(* ,300) sf,scmin
              write(10,300) sf,scmin
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Number of variables               : ',I7,
     +       /,1X,'Number of equality constraints    : ',I7,
     +       /,1X,'Number of inequality constraints  : ',I7,
     +       /,1X,'Number of bound constraints       : ',I7)

 200  format(/,1X,'Objective function scale factor   : ',1P,D7.1,
     +       /,1X,'The scaling feature was mainly developed for ',
     +            'constrained problems. For',/,1X,'unconstrained and ',
     +            'bound-constrained problem, please, set the ',
     +            'optimality',/,1X,'tolerance (related to the ',
     +            'sup-norm of the projected gradient of the',/,1X,
     +            'objective function) with a convenient value.')

 300  format(/,1X,'Objective function scale factor   : ',1P,D7.1,
     +       /,1X,'Smallest constraints scale factor : ',1P,D7.1)

 400  format(  1X,I6,1X,I6,1X,I6,1X,I6)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

      include "dim.par"
      include "scaling.inc"

C     LOCAL SCALARS
      integer i

      if ( scale ) then
          do i = 1,m
              lambda(i) = lambda(i) * sc(i) / sf
          end do

          scale = .false.
      end if

      call tendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalobjc(n,x,f,fu,m,c,cu,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision f,fu
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision c(m),cu(m),x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer j

      if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then

          if ( ignoref ) then
              fu = 0.0d0
          else
              call tevalf(n,x,fu,inform)
              if ( inform .lt. 0 ) return
          end if

          do j = 1,m
              call tevalc(n,x,j,cu(j),inform)
              if ( inform .lt. 0 ) return
          end do

      else ! if ( fccoded ) then
          call tevalfc(n,x,fu,m,cu,inform)
          if ( inform .lt. 0 ) return

          if ( ignoref ) fu = 0.0d0
      end if

      if ( scale ) then
          f = fu * sf

          do j = 1,m
              c(j) = cu(j) * sc(j)
          end do

      else
          f = fu

          do j = 1,m
              c(j) = cu(j)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

      if ( ignoref ) then
          f = 0.0d0
          return
      end if

      call tevalf(n,x,f,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) f = f * sf

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer i

      if ( ignoref ) then
          do i = 1,n
              g(i) = 0.0d0
          end do
          return
      end if

      call tevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,n
              g(i) = g(i) * sf
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer i

      if ( ignoref ) then
          hnnz = 0
          return
      end if

      call tevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,hnnz
              hval(i) = hval(i) * sf
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

      include "dim.par"
      include "scaling.inc"

      call tevalc(n,x,ind,c,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) c = c * sc(ind)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

      include "dim.par"
      include "scaling.inc"

C     LOCAL SCALARS
      integer i

      call tevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,jcnnz
              jcval(i) = jcval(i) * sc(ind)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

      include "dim.par"
      include "scaling.inc"

C     LOCAL SCALARS
      integer i

      call tevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,hnnz
              hval(i) = hval(i) * sc(ind)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhl(n,x,m,lambda,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      double precision tsf

      if ( scale ) then
          if ( ignoref ) then
              tsf = 0.0d0
          else
              tsf = sf
          end if

          call tevalhl(n,x,m,lambda,tsf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( ignoref ) then
              tsf = 0.0d0
          else
              tsf = usf
          end if

          call tevalhl(n,x,m,lambda,tsf,usc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhlp(n,x,m,lambda,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      double precision tsf

      if ( scale ) then
          if ( ignoref ) then
              tsf = 0.0d0
          else
              tsf = sf
          end if

          call tevalhlp(n,x,m,lambda,tsf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

      else
          if ( ignoref ) then
              tsf = 0.0d0
          else
              tsf = usf
          end if

          call tevalhlp(n,x,m,lambda,tsf,usc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer j

      call tevalfc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

      if ( ignoref ) f = 0.0d0

      if ( scale ) then
          f = f * sf

          do j = 1,m
              c(j) = c(j) * sc(j)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer i

      call tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      if ( ignoref ) then
          do i = 1,n
              g(i) = 0.0d0
          end do
      end if

      if ( scale ) then
          do i = 1,n
              g(i) = g(i) * sf
          end do

          do i = 1,jcnnz
              jcval(i) = jcval(i) * sc(jcfun(i))
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalgjacp(n,x,g,m,p,q,work,gotj,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gotj
      integer inform,m,n
      character work

C     ARRAY ARGUMENTS
      double precision g(n),p(m),q(n),x(n)

      include "dim.par"
      include "scaling.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer i,j

      if ( scale ) then
          if ( work .eq. 't' .or. work .eq. 'T' ) then
              do j = 1,m
                  p(j) = p(j) * sc(j)
              end do
          end if
      end if

      call tevalgjacp(n,x,g,m,p,q,work,gotj,inform)
      if ( inform .lt. 0 ) return

      if ( ignoref ) then
          if ( work .eq. 'J' .or. work .eq. 'T' ) then
              do i = 1,n
                  g(i) = 0.0d0
              end do
          end if
      end if

      if ( scale ) then
          if ( work .eq. 'j' .or. work .eq. 'J' ) then
              do j = 1,m
                  p(j) = p(j) * sc(j)
              end do
          end if

          if ( work .eq. 'J' .or. work .eq. 'T' ) then
              do i = 1,n
                  g(i) = g(i) * sf
              end do
          end if
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ssetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     EXTERNAL SUBROUTINES
      external tsetp

      call tsetp(n,x)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sunsetp()

      implicit none

C     EXTERNAL SUBROUTINES
      external tunsetp

      call tunsetp()

      end
