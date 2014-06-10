C     SAFE CALL TO THE USER PROVIDED SUBROUTINES

C     ******************************************************************
C     ******************************************************************

      subroutine vinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(11),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

      include "dim.par"

C     LOCAL SCALARS
      integer i

C     INTRINSEC FUNCTIONS
      intrinsic max,min

C     EXTERNAL SUBROUTINES
      external checkd

C     Avoid huge bounds

      do i = 1,n
          l(i) = max( l(i), - 1.0d+20 )
          u(i) = min( u(i),   1.0d+20 )
      end do

C     Project initial guess

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Check derivatives

      if ( checkder ) then
          call checkd(n,l,u,m,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

      include "outtyp.inc"

C     LOCAL SCALARS
      integer i

      if ( iprintctl(2) ) then

C         Save solution

C          open(20,file='solution.txt')

C         Point

          write(10,100)
          do i = 1,n
              write(10,200) i,x(i)
          end do

C         Lagrange multipliers

          if ( m .gt. 0 ) then
              write(10,300)
              do i = 1,m
                  write(10,200) i,lambda(i)
              end do
          end if

          close(20)

      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,'FINAL POINT:',//,2X,'INDEX',16X,'X(INDEX)')
 200  format(I7,1P,D24.16)
 300  format(/,'FINAL ESTIMATION OF THE LAGRANGE MULTIPLIERS: ',
     +       //,2X,'INDEX',11X,'LAMBDA(INDEX)')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer ercode,flag

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalf

      call evalf(n,x,f,flag)

      efcnt = efcnt + 1
      fcnt  = fcnt  + 1

      if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
C              write(* ,100)
              write(10,100)
          end if

          if ( safemode ) then
              inform = - 80
              call reperr(inform)
              return
          end if
      end if

      if ( .not. IsANumber(f) ) then
          if ( iprintctl(3) ) then
C              write(* ,200)
C              write(* ,300) f
              write(10,200)
              write(10,300) f
          end if

          if ( safemode ) then
              inform = - 80
              call reperr(inform)
              return
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALF WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALF WARNING: The objective function value ',
     +            'computed by the user-supplied subroutine EVALF is ',
     +            '+Inf, -Inf or NaN.')

 300  format(/,1X,'Value: ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalg,ivevalg

      if ( gcoded ) then
          call evalg(n,x,g,flag)

          egcnt = egcnt + 1

          if ( flag .ne. 0 ) then
              if ( iprintctl(3) ) then
C                  write(* ,100)
                  write(10,100)
              end if

              if ( safemode ) then
                  inform = - 81
                  call reperr(inform)
                  return
              end if
          end if

          do i = 1,n
              if ( .not. IsANumber(g(i)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,200)
C                      write(* ,300) n,i,g(i)
                      write(10,200)
                      write(10,300) n,i,g(i)
                  end if

                  if ( safemode ) then
                      inform = - 81
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ivevalg(n,x,g,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALG WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALG WARNING: There is an element whose value is ',
     +            '+Inf, -Inf or NaN in the gradient of the objective ',
     +            'function computed by the user-supplied subroutine ',
     +            'EVALG.')

 300  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ivevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

      include "machconst.inc"

C     LOCAL SCALARS
      integer j
      double precision fminus,fplus,step,tmp

C     INTRINSIC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external vsetp,vevalf

      do j = 1,n
          tmp  = x(j)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(j) = tmp + step
          call vsetp(n,x)
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp - step
          call vsetp(n,x)
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          g(j) = ( fplus - fminus ) / ( 2.0d0 * step )

          x(j) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalh,reperr

      call evalh(n,x,hlin,hcol,hval,hnnz,flag)

      ehcnt = ehcnt + 1

      if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
C              write(* ,100)
              write(10,100)
          end if

          if ( safemode ) then
              inform = - 82
              call reperr(inform)
              return
          end if
      end if

      do i = 1,hnnz
          if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +         hcol(i) .lt. 1 .or. hcol(i) .gt. n .or.
     +         hcol(i) .gt. hlin(i) ) then

              if ( iprintctl(3) ) then
C                  write(* ,200)
C                  write(* ,400) n,i,hlin(i),hcol(i),hval(i)
                  write(10,200)
                  write(10,400) n,i,hlin(i),hcol(i),hval(i)
              end if

              hlin(i) = 1
              hcol(i) = 1
              hval(i) = 0.0d0
          end if

          if ( .not. IsANumber(hval(i)) ) then
              if ( iprintctl(3) ) then
C                  write(* ,300)
C                  write(* ,400) n,i,hlin(i),hcol(i),hval(i)
                  write(10,300)
                  write(10,400) n,i,hlin(i),hcol(i),hval(i)
              end if

              if ( safemode ) then
                  inform = - 82
                  call reperr(inform)
                  return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALH WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALH WARNING: There is an element out of range, ',
     +            'or in the upper triangle, of the Hessian of the ',
     +            'objetive function computed by the user-supplied ',
     +            'subroutine EVALH. It will be ignored.')

 300  format(/,1X,'VEVALH WARNING: There is an element whose value is ',
     +            '+Inf, -Inf or NaN in the Hessian of the objetive ',
     +            'function computed by the user-supplied subroutine ',
     +            'EVALH.')

 400  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalc,reperr

      call evalc(n,x,ind,c,flag)

      eccnt(ind) = eccnt(ind) + 1

      if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
C              write(* ,100)
              write(10,100)
          end if

          if ( safemode ) then
              inform = - 83
              call reperr(inform)
              return
          end if
      end if

      if ( .not. IsANumber(c) ) then
          if ( iprintctl(3) ) then
C              write(* ,200) ind
C              write(* ,300) c
              write(10,200) ind
              write(10,300) c
          end if

          if ( safemode ) then
              inform = - 83
              call reperr(inform)
              return
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALC WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALC WARNING: The value of constraint ',I16,' ',
     +            'computed by the user-supplied subroutine EVALC is ',
     +            '+Inf, -Inf or NaN.')

 300  format(/,1X,'Value: ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evaljac,ivevaljac,reperr

      if ( jaccoded ) then
          call evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

          ejccnt(ind) = ejccnt(ind) + 1

          if ( flag .ne. 0 ) then
              if ( iprintctl(3) ) then
C                  write(* ,100)
                  write(10,100)
              end if

              if ( safemode ) then
                  inform = - 84
                  call reperr(inform)
                  return
              end if
          end if

          do i = 1,jcnnz
              if ( jcvar(i) .lt. 1 .or. jcvar(i) .gt. n ) then

                  if ( iprintctl(3) ) then
C                      write(* ,200) ind
C                      write(* ,400) n,i,jcvar(i),jcval(i)
                      write(10,200) ind
                      write(10,400) n,i,jcvar(i),jcval(i)
                  end if

                  jcvar(i) = 1
                  jcval(i) = 0.0d0
              end if

              if ( .not. IsANumber(jcval(i)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,300) ind
C                      write(* ,400) n,i,jcvar(i),jcval(i)
                      write(10,300) ind
                      write(10,400) n,i,jcvar(i),jcval(i)
                  end if

                  if ( safemode ) then
                      inform = - 84
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ivevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALJAC WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALJAC WARNING: There is an element out of ',
     +            'range in the gradient of constraint ',I16,' ',
     +            'computed by the user-supplied subroutine EVALJAC. ',
     +            'It will be ignored.')

 300  format(/,1X,'VEVALJAC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the gradient of constraint ',
     +            I16,'computed by the user-supplied subroutine ',
     +            'EVALJAC.')

 400  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Variable : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ivevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision jcval(n),x(n)

      include "dim.par"
      include "machconst.inc"

C     LOCAL SCALARS
      integer j
      double precision cminus,cplus,step,tmp

C     INTRINSEC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external vsetp,vevalc

      jcnnz = 0

      do j = 1,n
          tmp  = x(j)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(j) = tmp + step
          call vsetp(n,x)
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp - step
          call vsetp(n,x)
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          jcvar(jcnnz + 1) = j
          jcval(jcnnz + 1) = ( cplus - cminus ) / ( 2.0d0 * step )

          if ( abs( jcval(jcnnz + 1) ) .gt. 0.0d0 ) then
              jcnnz = jcnnz + 1
          end if

          x(j) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhc,reperr

      call evalhc(n,x,ind,hlin,hcol,hval,hnnz,flag)

      ehccnt(ind) = ehccnt(ind) + 1

      if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
C              write(* ,100)
              write(10,100)
          end if

          if ( safemode ) then
              inform = - 85
              call reperr(inform)
              return
          end if
      end if

      do i = 1,hnnz
          if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +         hcol(i) .lt. 1 .or. hcol(i) .gt. n .or.
     +         hcol(i) .gt. hlin(i) ) then

              if ( iprintctl(3) ) then
C                  write(* ,200) ind
C                  write(* ,400) n,i,hlin(i),hcol(i),hval(i)
                  write(10,200) ind
                  write(10,400) n,i,hlin(i),hcol(i),hval(i)
              end if

              hlin(i) = 1
              hcol(i) = 1
              hval(i) = 0.0d0
          end if

          if ( .not. IsANumber(hval(i)) ) then
              if ( iprintctl(3) ) then
C                  write(* ,300) ind
C                  write(* ,400) n,i,hlin(i),hcol(i),hval(i)
                  write(10,300) ind
                  write(10,400) n,i,hlin(i),hcol(i),hval(i)
              end if

              if ( safemode ) then
                  inform = - 85
                  call reperr(inform)
                  return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHC WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALHC WARNING: There is an element out of range ',
     +            'or in the upper triangle of the Hessian of ',
     +            'constraint ',I16,' computed by the user-supplied ',
     +            'subroutine EVALHC. It will be ignored.')

 300  format(/,1X,'VEVALHC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the Hessian of constraint ',
     +            I16,' computed by the user-supplied subroutine ',
     +            'EVALHC.')

 400  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUCNTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalfc,reperr

      call evalfc(n,x,f,m,c,flag)

      efccnt = efccnt + 1
      fcnt   = fcnt   + 1

      if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
C              write(* ,100)
              write(10,100)
          end if

          if ( safemode ) then
              inform = - 86
              call reperr(inform)
              return
          end if
      end if

      if ( .not. IsANumber(f) ) then
          if ( iprintctl(3) ) then
C              write(* ,200)
C              write(* ,400) f
              write(10,200)
              write(10,400) f
          end if

          if ( safemode ) then
              inform = - 86
              call reperr(inform)
              return
          end if
      end if

      do i = 1,m
          if ( .not. IsANumber(c(i)) ) then
              if ( iprintctl(3) ) then
C                  write(* ,300)
C                  write(* ,500) n,m,i,c(i)
                  write(10,300)
                  write(10,500) n,m,i,c(i)
              end if

              if ( safemode ) then
                  inform = - 86
                  call reperr(inform)
                  return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALFC WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALFC WARNING: The objective function value ',
     +            'computed by the user-supplied',/,1X,'subroutine ',
     +            'EVALFC is +Inf, -Inf or NaN.')

 300  format(/,1X,'VEVALFC WARNING: The value of a constraint ',
     +            'computed by the user-supplied',/,1X,'subroutine ',
     +            'EVALFC is +Inf, -Inf or NaN.')

 400  format(/,1X,'Value: ',1P,D24.16)

 500  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Number of constraints : ',I16,
     +       /,1X,'Constraint            : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUCNTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalgjac,ivevalgjac,reperr

      if ( gjaccoded ) then
          call evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

          egjccnt = egjccnt + 1

          if ( flag .ne. 0 ) then
              if ( iprintctl(3) ) then
C                  write(* ,100)
                  write(10,100)
              end if

              if ( safemode ) then
                  inform = - 87
                  call reperr(inform)
                  return
              end if
          end if

          do i = 1,n
              if ( .not. IsANumber(g(i)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,300)
C                      write(* ,400) n,i,g(i)
                      write(10,300)
                      write(10,400) n,i,g(i)
                  end if

                  if ( safemode ) then
                      inform = - 87
                      call reperr(inform)
                      return
                  end if
              end if
          end do

          do i = 1,jcnnz
              if ( jcfun(i) .lt. 1 .or. jcfun(i) .gt. m .or.
     +             jcvar(i) .lt. 1 .or. jcvar(i) .gt. n ) then

                  if ( iprintctl(3) ) then
C                      write(* ,200)
C                      write(* ,500) n,m,i,jcfun(i),jcvar(i),jcval(i)
                      write(10,200)
                      write(10,500) n,m,i,jcfun(i),jcvar(i),jcval(i)
                  end if

                  jcfun(i) = 1
                  jcvar(i) = 1
                  jcval(i) = 0.0d0
              end if

              if ( .not. IsANumber(jcval(i)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,300)
C                      write(* ,500) n,m,i,jcfun(i),jcvar(i),jcval(i)
                      write(10,300)
                      write(10,500) n,m,i,jcfun(i),jcvar(i),jcval(i)
                  end if

                  if ( safemode ) then
                      inform = - 87
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ivevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALGJAC WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALGJAC WARNING: There is an element out of ',
     +            'range in the gradient of the objective function or ',
     +            'in the Jacobian of the constraints computed by the ',
     +            'user-supplied subroutine EVALGJAC. It will be ',
     +            'ignored.')

 300  format(/,1X,'VEVALGJAC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the',/,1X,'gradient of the ',
     +            'objective function or in the Jacobian of the ',
     +            'constraints',/,1X,'computed by the user-supplied ',
     +            'subroutine EVALGJAC.')

 400  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

 500  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Number of constraints : ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Constraint            : ',I16,
     +       /,1X,'Variable              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ivevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

      include "dim.par"
      include "machconst.inc"

C     LOCAL SCALARS
      integer i,j
      double precision fminus,fplus,step,tmp

C     LOCAL ARRAYS
      double precision cminus(mmax),cplus(mmax)

C     INTRINSIC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external vsetp,vevalfc

      jcnnz = 0

      do i = 1,n
          tmp  = x(i)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(i) = tmp + step
          call vsetp(n,x)
          call vevalfc(n,x,fplus,m,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step
          call vsetp(n,x)
          call vevalfc(n,x,fminus,m,cminus,inform)
          if ( inform .lt. 0 ) return

          do j = 1,m
              jcfun(jcnnz + 1) = j
              jcvar(jcnnz + 1) = i
              jcval(jcnnz + 1) = ( cplus(j) - cminus(j) ) /
     +                           ( 2.0d0 * step )

              if ( abs( jcval(jcnnz + 1) ) .gt. 0.0d0 ) then
                  jcnnz = jcnnz + 1
              end if
          end do

          g(i) = ( fplus - fminus ) / ( 2.0d0 * step )

          x(i) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalgjacp(n,x,g,m,p,q,work,gotj,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gotj
      integer inform,m,n
      character work

C     ARRAY ARGUMENTS
      double precision g(n),p(m),q(n),x(n)

C     The meaning of argument work follows: work = 'j' or 'J' means that 
C     q is an input array argument and that p = Jacobian x q must be 
C     computed, while work = 't' or 'T' means that p is an input array 
C     argument and that q = Jacobian^t x p must be computed. Moreover, a 
C     capital letter (i.e. 'J' or 'T') means that the gradient of the 
C     objective function g must also be computed. A lower letter (i.e. 
C     'j' or 't') means that only the product with the Jacobian is 
C     required. In the later case, input array argument g MUST NOT be 
C     modified nor referenced.

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i,j

C     EXTERNAL FUCNTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalgjacp,reperr

      call evalgjacp(n,x,g,m,p,q,work,gotj,flag)

      egjcpcnt = egjcpcnt + 1

      if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
C              write(* ,100)
              write(10,100)
          end if

          if ( safemode ) then
              inform = - 88
              call reperr(inform)
              return
          end if
      end if

      if ( work .eq. 'J' .or. work .eq. 'T' ) then
          do i = 1,n
              if ( .not. IsANumber(g(i)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,200)
C                      write(* ,400) n,i,g(i)
                      write(10,200)
                      write(10,400) n,i,g(i)
                  end if

                  if ( safemode ) then
                      inform = - 88
                      call reperr(inform)
                      return
                  end if
              end if
          end do
      end if

      if ( work .eq. 'j' .or. work .eq. 'J' ) then
          do j = 1,m
              if ( .not. IsANumber(p(j)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,300)
C                      write(* ,400) m,j,p(j)
                      write(10,300)
                      write(10,400) m,j,p(j)
                  end if

                  if ( safemode ) then
                      inform = - 88
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else ! if ( work .eq. 't' .or. work .eq. 'T' ) then
          do i = 1,n
              if ( .not. IsANumber(q(i)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,300)
C                      write(* ,400) n,i,q(i)
                      write(10,300)
                      write(10,400) n,i,q(i)
                  end if

                  if ( safemode ) then
                      inform = - 88
                      call reperr(inform)
                      return
                  end if
              end if
          end do
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALGJACP WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALGJACP WARNING: There is an element whose ',
     +            'value is +Inf, -Inf or NaN in the gradient of the ',
     +            'objective function computed by the user-supplied ',
     +            'subroutine EVALGJACP.')

 300  format(/,1X,'VEVALGJACP WARNING: There is an element whose ',
     +            'value is +Inf, -Inf or NaN in the product of the ',
     +            'Jacobian (or transpose Jacobian) times a given ',
     +            'vector, computed by the user-supplied subroutine ',
     +            'EVALGJACP.')

 400  format(/,1X,'Dimension of computed array: ',I16,
     +       /,1X,'Position                   : ',I16,
     +       /,1X,'Value                      : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhl,ivevalhl,reperr

      if ( hlcoded ) then

          call evalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,flag)

          ehlcnt = ehlcnt + 1

          if ( flag .ne. 0 ) then
              if ( iprintctl(3) ) then
C                  write(* ,100)
                  write(10,100)
              end if

              if ( safemode ) then
                  inform = - 89
                  call reperr(inform)
                  return
              end if
          end if

          do i = 1,hnnz
              if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +             hcol(i) .lt. 1 .or. hcol(i) .gt. n .or.
     +             hcol(i) .gt. hlin(i) ) then

                  if ( iprintctl(3) ) then
C                      write(* ,200)
C                      write(* ,400) n,i,hlin(i),hcol(i),hval(i)
                      write(10,200)
                      write(10,400) n,i,hlin(i),hcol(i),hval(i)
                  end if

                  hlin(i) = 1
                  hcol(i) = 1
                  hval(i) = 0.0d0
              end if

              if ( .not. IsANumber(hval(i)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,300)
C                      write(* ,400) n,i,hlin(i),hcol(i),hval(i)
                      write(10,300)
                      write(10,400) n,i,hlin(i),hcol(i),hval(i)
                  end if

                  if ( safemode ) then
                      inform = - 89
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else if ( hcoded .and. ( hccoded .or. m .eq. 0 ) ) then
         call ivevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
         if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHL WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALHL WARNING: There is an element out of range, ',
     +            'or in the upper triangle, of the',/,1X,'Hessian of ',
     +            'the Lagrangian computed by the user-supplied ',
     +            'subroutine EVALHL. It',/,1X,'will be ignored.')

 300  format(/,1X,'VEVALHL WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the',/,1X,'Hessian of the ',
     +            'Lagrangian computed by the user-supplied ',
     +            'subroutine EVALHL.')

 400  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ivevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

      include "dim.par"
      include "algparam.inc"

C     LOCAL SCALARS
      integer col,con,hnnztmp,i,ind,itmp,j,lin,nextj,rnnz
      double precision val

C     LOCAL ARRAYS
      integer hcon(hnnzmax),pos(nmax),rind(nmax),stlin(nmax)
      double precision rval(nmax)

C     EXTERNAL SUBROUTINES
      external vevalh,vevalhc

C     ==================================================================
C     COMPUTE HESSIANS
C     ==================================================================

C     COMPUTE HESSIAN OF THE OBJECTIVE FUNCTION

      if ( ignoref ) then
          hnnz = 0

      else
C         Compute Hessian of the objective function
          call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

C         For each element of the Hessian of the objective function,
C         set constraint index as zero
          do i = 1,hnnz
              hval(i) = hval(i) * sf
              hcon(i) = 0
          end do
      end if

      if ( m .eq. 0 ) return

C     COMPUTE HESSIANS OF THE CONSTRAINTS

      ind = 0

      do j = 1,m
C         Compute Hessian of constraint j
          call vevalhc(n,x,j,hlin(hnnz+ind+1),hcol(hnnz+ind+1),
     +    hval(hnnz+ind+1),hnnztmp,inform)
          if ( inform .lt. 0 ) return

C         For each element of the Hessian, set constraint as j
          do i = hnnz + ind + 1,hnnz + ind + hnnztmp
              hval(i) = hval(i) * sc(j)
              hcon(i) = j
          end do

          ind = ind + hnnztmp
      end do

      if ( ind .eq. 0 ) return

      hnnz = hnnz + ind

C     ==================================================================
C     SET ROW LINKED LISTS
C     ==================================================================

C     Initialize pointers to the first element of each row
      do i = 1,n
         stlin(i) = 0
      end do

C     Set row linked lists
      do i = 1,hnnz
         lin        = hlin(i)
         itmp       = stlin(lin)
         stlin(lin) = i
         hlin(i)    = itmp
      end do

C     ==================================================================
C     BUILD HESSIAN OF THE LAGRANGIAN ROW BY ROW
C     ==================================================================

C     Initialize array pos
      do i = 1,n
          pos(i) = 0
      end do

      do i = 1,n
C         Initialize the i-th row of the Hessian of the Lagrangian
          rnnz = 0

C         Process the i-th row of all the Hessians
          j = stlin(i)

 10       if ( j .ne. 0 ) then

C             Process element (i,hcol(j)) of the Hessian of constraint
C             hcon(j) (Hessian of the objective function if hcon(j)=0)

              col = hcol(j)
              con = hcon(j)
              if ( con .eq. 0 ) then
                  val = hval(j)
              else
                  val = hval(j) * lambda(con)
              end if

              if ( pos(col) .ne. 0 ) then
                  rval(pos(col)) = rval(pos(col)) + val

              else
                  rnnz           = rnnz + 1
                  pos(col)       = rnnz
                  rind(pos(col)) = col
                  rval(pos(col)) = val
              end if

C             Get next element in the i-th row linked list
              j = hlin(j)
              go to 10
          end if

C         Clean array pos
          do j = 1,rnnz
              pos(rind(j)) = 0
          end do

C         Set i-th row of hl (over the i-th rows of the Hessians)
C         and mark remaining elements to be deleted
          j = stlin(i)

 20       if ( j .ne. 0 ) then
              nextj = hlin(j)

              if ( rnnz .ne. 0 ) then
                  hlin(j) = i
                  hcol(j) = rind(rnnz)
                  hval(j) = rval(rnnz)
                  rnnz    = rnnz - 1
              else
                  hlin(j) = 0
              end if

              j = nextj
              go to 20
          end if

      end do

C     Eliminate remaining elements (marked with hlin(j)=0)
      j = 1

 30   if ( j .le. hnnz ) then
          if ( hlin(j) .eq. 0 ) then
              if ( j .ne. hnnz ) then
                  hlin(j) = hlin(hnnz)
                  hcol(j) = hcol(hnnz)
                  hval(j) = hval(hnnz)
              end if
              hnnz = hnnz - 1
          else
              j = j + 1
          end if

          go to 30
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

      include "dim.par"
      include "algparam.inc"
      include "counters.inc"
      include "outtyp.inc"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhlp,ivevalhlp,reperr

      if ( hlpcoded ) then
          call evalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,flag)

          ehlpcnt = ehlpcnt + 1

          if ( flag .ne. 0 ) then
              if ( iprintctl(3) ) then
C                  write(* ,100)
                  write(10,100)
              end if

              if ( safemode ) then
                  inform = - 90
                  call reperr(inform)
                  return
              end if
          end if

          do i = 1,n
              if ( .not. IsANumber(hp(i)) ) then
                  if ( iprintctl(3) ) then
C                      write(* ,200)
C                      write(* ,300) n,i,hp(i)
                      write(10,200)
                      write(10,300) n,i,hp(i)
                  end if

                  if ( safemode ) then
                      inform = - 90
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else if ( seconde ) then
          call ivevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHLP WARNING: A non-null flag was returned.',/)

 200  format(/,1X,'VEVALHLP WARNING: There is an element in the ',
     +            'product of the Hessian of the Lagrangian by a ',
     +            'computed by the user-supplied subroutine EVALHLP ',
     +            'whose value is +Inf, -Inf or NaN.')

 300  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ivevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

      include "dim.par"
      include "hessdat.inc"

C     LOCAL SCALARS
      integer col,i,lin
      double precision val

C     EXTERNAL SUBROUTINES
      external vevalhl

      if ( .not. gothl ) then
          gothl = .true.

          call vevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      do i = 1,n
          hp(i) = 0.0d0
      end do

      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          val = hval(i)

          hp(lin) = hp(lin) + p(col) * val

          if ( lin .ne. col ) then
              hp(col) = hp(col) + p(lin) * val
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine reperr(inform)

C     SCALAR ARGUMENTS
      integer inform

      include "outtyp.inc"

      if ( iprintctl(3) ) then
C          write(* ,100) inform
          write(10,100) inform
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'*** Safemode is enable. Terminating execution. ***',
     +       /,1X,'*** Error code = ',I3,'***',/)

      end

C     ******************************************************************
C     ******************************************************************

      logical function IsANumber(x)

      implicit none

C     SCALAR ARGUMENTS
      double precision x

      include "machconst.inc"

      IsANumber = .true.

      if ( .not. abs( x ) .le. bignum ) then
          IsANumber = .false.
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vsetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

      include "dim.par"
      include "graddat.inc"

C     EXTERNAL SUBROUTINES
      external setp

      gotc = .false.
      gotj = .false.

      call setp(n,x)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vunsetp()

      implicit none

      include "dim.par"
      include "graddat.inc"

C     EXTERNAL SUBROUTINES
      external unsetp

      gotc = .false.
      gotj = .false.

      call unsetp()

      end
