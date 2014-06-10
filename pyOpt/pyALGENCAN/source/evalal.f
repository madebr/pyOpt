C     *****************************************************************
C     *****************************************************************

      subroutine sevalal(n,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision al
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

      include "dim.par"
      include "graddat.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer j
      double precision f,p

      if ( innercall ) then
          call minsqf(n,x,al,inform)
          return
      end if

      if ( fccoded ) then

C         COMPUTE OBJECTIVE FUNTION AND CONSTRAINTS
          call sevalfc(n,x,f,m,c,inform)
          if ( inform .lt. 0 ) return

C         COMPUTES AL = f + sum_j P(c_j, rho_j, lambda_j)
          al = f

          do j = 1,m
C             ADD P(c_j, rho_j, lambda_j)
              call evalp(c(j),rho(j),lambda(j),equatn(j),p)
              al = al + p
          end do

      else if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then

C         COMPUTE OBJECTIVE FUNCTION
          call sevalf(n,x,f,inform)
          if ( inform .lt. 0 ) return

C         COMPUTES AL = f + sum_j P(c_j, rho_j, lambda_j)
          al = f

          do j = 1,m
C             COMPUTE j-TH CONSTRAINT
              call sevalc(n,x,j,c(j),inform)
              if ( inform .lt. 0 ) return

C             ADD P(c_j, rho_j, lambda_j)
              call evalp(c(j),rho(j),lambda(j),equatn(j),p)
              al = al + p
          end do

      end if

      gotc = .true.

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnl(n,x,m,lambda,equatn,linear,nl,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nl(n),x(n)

      include "dim.par"
      include "algparam.inc"
      include "graddat.inc"

C     LOCAL SCALARS
      logical alllin,allnol
      integer i,ind,j,jcnnz
      double precision dum

C     LOCAL ARRAYS
      integer jcfun(jcnnzmax)
      double precision p(mmax)

      if ( fccoded ) then

          if ( gjacpcoded ) then

C             COMPUTE p
              constrc = .false.
              do j = 1,m
                  if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then
                      p(j) = lambda(j)
                      constrc = .true.
                  else
                      p(j) = 0.0d0
                  end if
              end do

C             COMPUTE nl = Jacobian^T p
              call sevalgjacp(n,x,g,m,p,nl,'T',gotj,inform)
              if ( inform .lt. 0 ) return

C             COMPUTE gparc = Jacobian^t p IGNORING LINEAR CONSTRAINTS
              allnol = .true.
              alllin = .true.
              do j = 1,m
                  if ( linear(j) ) then
                      p(j) = 0.0d0
                      allnol = .false.
                  else
                      alllin = .false.
                  end if
              end do

              if ( allnol ) then
                  do i = 1,n
                      gparc(i) = nl(i)
                  end do
              else if ( alllin ) then
                  do i = 1,n
                      gparc(i) = 0.0d0
                  end do
              else
                  call sevalgjacp(n,x,dum,m,p,gparc,'t',gotj,inform)
                  if ( inform .lt. 0 ) return
              end if

C             ADD THE GRADIENT OF THE OBJECTIVE FUNCTION
              do i = 1,n
                  nl(i)    = nl(i)    + g(i)
                  gparc(i) = gparc(i) + g(i)
              end do

          else ! if ( gjaccoded ) then
C         In fact, this else is the choice for gjaccoded or 'not 
C         gjacpcoded and not gjaccoded', in which case the calling 
C         sequence starting with sevalgjac below will end up using 
C         finite differences to approximate first derivatives of the 
C         objective function and the constraints

C             COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AND 
C             JACOBIAN OF CONSTRAINTS
              call sevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
              if ( inform .lt. 0 ) return

C             CONVERT JACOBIAN OF CONSTRAINTS FROM COORDINATE FORMAT TO
C             COMPRESSED SPARSE ROW FORMAT
              call coo2csr(m,jcnnz,jcfun,jcvar,jcval,jclen,jcsta)

C             COMPUTE \nabla L = \nabla f + \sum_j lambda_j * \nabla c_j
              constrc = .false.

              do i = 1,n
                  nl(i)    = g(i)
                  gparc(i) = g(i)
              end do

              do j = 1,m
                  if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then
                      do i = jcsta(j),jcsta(j) + jclen(j) - 1
                          nl(jcvar(i)) = 
     +                    nl(jcvar(i)) + lambda(j) * jcval(i)
                      end do

                      if ( .not. linear(j) ) then
                          do i = jcsta(j),jcsta(j) + jclen(j) - 1
                              gparc(jcvar(i)) =
     +                        gparc(jcvar(i)) + lambda(j) * jcval(i)
                          end do
                      end if

                      constrc = .true.
                  end if
              end do

          end if

      else if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then

C         COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION
          call sevalg(n,x,g,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              nl(i)    = g(i)
              gparc(i) = g(i)
          end do

C         ADD nl = Jacobian^T lambda
          constrc = .false.

          ind = 0

          do j = 1,m
              if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then

                  jcsta(j) = ind + 1

C                 COMPUTE THE GRADIENT OF THE j-TH CONSTRAINT
                  call sevaljac(n,x,j,jcvar(ind+1),jcval(ind+1),
     +            jclen(j),inform)
                  if ( inform .lt. 0 ) return

                  ind = ind + jclen(j)

C                 ADD lambda_j * \nabla c_j
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      nl(jcvar(i)) = 
     +                nl(jcvar(i)) + lambda(j) * jcval(i)
                  end do

                  if ( .not. linear(j) ) then
                      do i = jcsta(j),jcsta(j) + jclen(j) - 1
                          gparc(jcvar(i)) =
     +                    gparc(jcvar(i)) + lambda(j) * jcval(i)
                      end do
                  end if

                  constrc = .true.

              end if
          end do

      end if

      gotj = .true.

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnal(n,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(n),rho(m),x(n)

      include "dim.par"
      include "graddat.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer j
      double precision dum

      if ( innercall ) then
          call minsqg(n,x,nal,inform)
          return
      end if

      if ( fccoded ) then
C         COMPUTE CONSTRAINTS
          if ( .not. gotc .and. m .gt. 0 ) then
              call sevalfc(n,x,dum,m,c,inform)
              if ( inform .lt. 0 ) return
          end if

      else if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then
C         COMPUTE CONSTRAINTS ONE BY ONE
          if ( .not. gotc ) then
              do j = 1,m
                  call sevalc(n,x,j,c(j),inform)
                  if ( inform .lt. 0 ) return
              end do
          end if
      end if

      gotc = .true.

C     COMPUTE dP/dc
      do j = 1,m
          call evaldpdy(c(j),rho(j),lambda(j),equatn(j),dpdc(j))
      end do

C     COMPUTE GRADIENT OF THE LAGRANGIAN WITH DPDC INSTEAD OF LAMBDA
      call sevalnl(n,x,m,dpdc,equatn,linear,nal,inform)
      if ( inform .lt. 0 ) return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalp(y,rho,lambda,equatn,p)

      implicit none

C     SCALAR ARGUMENTS
      logical equatn
      double precision lambda,p,rho,y

      if ( equatn ) then
          p  = y * ( lambda + 0.5d0 * rho * y )
      else
          if ( lambda + rho * y .ge. 0.0d0 ) then
              p = y * ( lambda + 0.5d0 * rho * y )
          else
              p = - 0.5d0 * lambda ** 2 / rho
          end if
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evaldpdy(y,rho,lambda,equatn,dpdy)

      implicit none

C     SCALAR ARGUMENTS
      logical equatn
      double precision y,rho,lambda,dpdy

      if ( equatn ) then
          dpdy = lambda + rho * y
      else
          dpdy = max( 0.0d0, lambda + rho * y )
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalnal(n,xp,m,lambda,rho,equatn,linear,iglin,nalp,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical iglin
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nalp(n),rho(m),xp(n)

C     This subroutine computes the gradient of the Augmented Lagrangian
C     function at a point xp, which is near to x, taking care of the
C     non-differentiability (i.e. considering the contribution of the 
C     same constraints that contributed to compute the gradient of the 
C     augmented Lagrangian at x). The augmented Lagrangian gradient must 
C     have been previously computed at x. If iglin = true, linear 
C     constraints are ignored.

      include "dim.par"
      include "algparam.inc"
      include "graddat.inc"

C     LOCAL SCALARS
      logical alllin,allnol
      integer i,j,jcpnnz
      double precision cpj,dpdcpj,dum

C     LOCAL ARRAYS
      integer jcpfun(jcnnzmax),jcplen(mmax),jcpsta(mmax),
     +        jcpvar(jcnnzmax)
      double precision cp(mmax),gp(nmax),jcpval(jcnnzmax),dpdcp(mmax)

      if ( innercall ) then
          call minsqg(n,xp,nalp,inform)
          return
      end if

      if ( fccoded ) then

C         COMPUTE nalp = gp + Jacobian^T dpdcp IGNORING THE LINEAR 
C         CONSTRAINTS IF iglin IS TRUE

          if ( gjacpcoded ) then

C             COMPUTE CONSTRAINTS AT xp
              if ( m .gt. 0 ) then
                  call sevalfc(n,xp,dum,m,cp,inform)
                  if ( inform .lt. 0 ) return
              end if

C             COMPUTE dpdcp
              do j = 1,m
                  if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and.
     +                 .not. ( iglin .and. linear(j) ) ) then
                      dpdcp(j) = lambda(j) + rho(j) * cp(j)
                  else
                      dpdcp(j) = 0.0d0
                  end if
              end do

C             COMPUTE nalp = Jacobian^T dpdcp
              call sevalgjacp(n,xp,gp,m,dpdcp,nalp,'T',gotj,inform)
              if ( inform .lt. 0 ) return

C             ADD THE GRADIENT OF THE OBJECTIVE FUNCTION
              do i = 1,n
                  nalp(i) = nalp(i) + gp(i)
              end do

          else ! if ( gjaccoded ) then
C         In fact, this else is the choice for gjaccoded or 'not 
C         gjacpcoded and not gjaccoded', in which case the calling 
C         sequence starting with sevalgjac below will end up using 
C         finite differences to approximate first derivatives of the 
C         objective function and the constraints

C             COMPUTE CONSTRAINTS AT xp
              if ( m .gt. 0 ) then
                  call sevalfc(n,xp,dum,m,cp,inform)
                  if ( inform .lt. 0 ) return
              end if

C             COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AND JACOBIAN
C             OF CONSTRAINTS AT xp
              call sevalgjac(n,xp,gp,m,jcpfun,jcpvar,jcpval,jcpnnz,
     +        inform)
              if ( inform .lt. 0 ) return

              do i = 1,n
                  nalp(i) = gp(i)
              end do

C             CONVERT JACOBIAN OF CONSTRAINTS FROM COODINATE FORMAT TO
C             COMPRESSED SPARSE ROW FORMAT
              call coo2csr(m,jcpnnz,jcpfun,jcpvar,jcpval,jcplen,jcpsta)

C             COMPUTE dpdcp AND ADD dPdc * dcdx
              do j = 1,m
                  if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and.
     +                 .not. ( iglin .and. linear(j) ) ) then
                      dpdcpj = lambda(j) + rho(j) * cp(j)

                      do i = jcpsta(j),jcpsta(j) + jcplen(j) - 1
                          nalp(jcpvar(i)) =
     +                    nalp(jcpvar(i)) + dpdcpj * jcpval(i)
                      end do
                  end if
              end do
          end if

      else if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then

C         COMPUTE GRADIENT OF THE OBJECTIVE FUNCTION
          call sevalg(n,xp,gp,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              nalp(i) = gp(i)
          end do

C         ADD Jacobian^T dpdcp IGNORING THE LINEAR CONSTRAINTS 
C         IF iglin IS TRUE
          do j = 1,m
              if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and.
     +             .not. ( iglin .and. linear(j) ) ) then

C                 COMPUTE j-TH CONSTRAINT
                  call sevalc(n,xp,j,cpj,inform)
                  if ( inform .lt. 0 ) return

                  dpdcpj = lambda(j) + rho(j) * cpj

C                 COMPUTE THE GRADIENT OF THE j-TH CONSTRAINT
                  call sevaljac(n,xp,j,jcpvar,jcpval,jcpnnz,inform)
                  if ( inform .lt. 0 ) return

C                 ADD dPdc * dcdx
                  do i = 1,jcpnnz
                      nalp(jcpvar(i)) = 
     +                nalp(jcpvar(i)) + dpdcpj * jcpval(i)
                  end do
              end if
          end do

      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhal(n,x,m,lambda,rho,equatn,linear,hallin,halcol,
     +halval,halnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer halnnz,inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      integer hallin(*),halcol(*)
      double precision lambda(m),halval(*),rho(m),x(n)

C     This subroutine computes the Hessian of the augmented Lagrangian.

      include "dim.par"
      include "graddat.inc"

C     LOCAL SCALARS
      integer i,lin,j,k,l,var

C     LOCAL ARRAYS
      integer stlin(nmax)
      double precision r(nmax)

      call sevalhl(n,x,m,dpdc,hallin,halcol,halval,halnnz,inform)
      if ( inform .lt. 0 ) return

      if ( m .eq. 0 ) return

C     PUT MATRIX INTO A ROW-LINKED LIST

      do i = 1,n
         r(i) = 0.0d0
      end do

      do i = 1,n
          stlin(i) = 0
      end do

      do i = 1,halnnz
          lin = hallin(i)
          k   = stlin(lin)
          stlin(lin) = i
          hallin(i)  = k
      end do

C     ADD \sum_j \rho_j * \nabla c_j \nabla c_j^t

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             ADD \rho_j * \nabla c_j \nabla c_j^t

              do k = jcsta(j),jcsta(j) + jclen(j) - 1

                  var = jcvar(k)

C                 PUT ROW jcvar(k) INTO A DENSE VECTOR

                  lin = stlin(var)
 10               if ( lin .ne. 0 ) then
                      r(halcol(lin)) = r(halcol(lin)) + halval(lin)
                      lin = hallin(lin)
                      go to 10
                  end if

C                 ADD VALUE

                  do l = jcsta(j),jcsta(j) + jclen(j) - 1
                      if ( jcvar(l) .le. var ) then
                          r(jcvar(l)) =
     +                    r(jcvar(l)) + rho(j) * jcval(l) * jcval(k)
                      end if
                  end do

C                 UPDATE VALUES IN HALVAL

                  lin = stlin(var)
 20               if ( lin .ne. 0 ) then
                      halval(lin) = r(halcol(lin))
                      r(halcol(lin)) = 0.0d0
                      lin = hallin(lin)
                      go to 20
                  end if

C                 INSERT NEW ELEMENTS IN HAL REPRESENTATION

                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      l = jcvar(i)
                      if ( r(l) .ne. 0.0d0 ) then
                          halnnz = halnnz + 1
                          halval(halnnz) = r(l)
                          halcol(halnnz) = l
                          hallin(halnnz) = stlin(var)
                          stlin(var) = halnnz
                          r(l) = 0.0d0
                      end if
                  end do

              end do

          end if

      end do

C     PUT MATRIX BACK INTO COORDINATE SQUEME

      do i = 1,n

          lin = stlin(i)
 30       if ( lin .ne. 0 ) then
              k = hallin(lin)
              hallin(lin) = i
              lin = k
              go to 30
          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,gothl,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hp(n),lambda(m),p(n),rho(m),x(n)

      include "dim.par"
      include "graddat.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      integer i,j
      double precision ajtp,dum

C     LOCAL ARRAYS
      double precision ap(mmax),atap(nmax)

      if ( innercall ) then
          call minsqhp(n,x,p,hp,gothl,inform)
          return
      end if

C     --------------------------------------------------------------
C     Hessian approximation
C     --------------------------------------------------------------

      if ( hptype .eq. 'HAPPRO' .and. constrc ) then

          call applyhapp(n,m,rho,equatn,gothl,p,hp)

C     --------------------------------------------------------------
C     Incremental quotients
C     --------------------------------------------------------------

      else if ( hptype .eq. 'INCQUO' .or. hptype .eq. 'HAPPRO' ) then

          call ievalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,inform)
          if ( inform .lt. 0 ) return

C     --------------------------------------------------------------
C     True Hessian
C     --------------------------------------------------------------

      else if ( hptype .eq. 'TRUEHP' ) then

C         Compute Hessian of Lagrangian times p using dpdc
C         instead of lambda

          call sevalhlp(n,x,m,dpdc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

C         Add rho A^T A

          if ( gjacpcoded ) then
              call sevalgjacp(n,x,dum,m,ap,p,'j',gotj,inform)
              if ( inform .lt. 0 ) return

              do j = 1,m
                  ap(j) = ap(j) * rho(j)
              end do

              call sevalgjacp(n,x,dum,m,ap,atap,'t',gotj,inform)
              if ( inform .lt. 0 ) return

              do i = 1,n
                  hp(i) = hp(i) + atap(i)
              end do

          else
              do j = 1,m
                  if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

                      ajtp = 0.0d0
                      do i = jcsta(j),jcsta(j) + jclen(j) - 1
                          ajtp = ajtp + jcval(i) * p(jcvar(i))
                      end do

                      ajtp = ajtp * rho(j)

                      do i = jcsta(j),jcsta(j) + jclen(j) - 1
                          hp(jcvar(i)) = hp(jcvar(i)) + ajtp * jcval(i)
                      end do

                  end if
              end do
          end if

      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hp(n),lambda(m),p(n),rho(m),x(n)

C     Computes an approximation of the product of the Hessian of the
C     Augmented Lagrangian times a vector using incremental quotients.

      include "machconst.inc"
      include "dim.par"
      include "graddat.inc"
      include "algparam.inc"

C     LOCAL SCALARS
      logical iglin
      integer i,j
      double precision ajtp,dum,psupn,step,xsupn

C     LOCAL  ARRAYS
      double precision ap(mmax),atap(nmax),gpparc(nmax),ptmp(nmax),
     +        xp(nmax)

C     ------------------------------------------------------------------
C     Set auxiliary point
C     ------------------------------------------------------------------

      xsupn = 0.0d0
      psupn = 0.0d0
      do i = 1,n
          xsupn = max( xsupn, abs( x(i) ) )
          psupn = max( psupn, abs( p(i) ) )
      end do

      step = macheps12 * max( xsupn / psupn, 1.0d0 )

      do i = 1,n
          xp(i) = x(i) + step * p(i)
      end do

C     ------------------------------------------------------------------
C     Compute gradient of the augmented Lagrangian at xp considering the
C     same constraints considered at x and ignoring linear constraints
C     ------------------------------------------------------------------

      call ssetp(n,xp)

      iglin = .true.
      call ievalnal(n,xp,m,lambda,rho,equatn,linear,iglin,gpparc,inform)
      if ( inform .lt. 0 ) return

C     ------------------------------------------------------------------
C     Compute gradients difference
C     ------------------------------------------------------------------

      do i = 1,n
          hp(i) = ( gpparc(i) - gparc(i) ) / step
      end do

C     ------------------------------------------------------------------
C     Add contribution of linear constraints
C     ------------------------------------------------------------------

      if ( gjacpcoded ) then
          do i = 1,n
              if ( linear(i) ) then
                  ptmp(i) = p(i)
              else
                  ptmp(i) = 0.0d0
              end if
          end do

          call sevalgjacp(n,x,dum,m,ap,ptmp,'j',gotj,inform)
          if ( inform .lt. 0 ) return

          do j = 1,m
              ap(j) = ap(j) * rho(j)
          end do

          call sevalgjacp(n,x,dum,m,ap,atap,'t',gotj,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              hp(i) = hp(i) + atap(i)
          end do

      else
          do j = 1,m
              if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then
                  if ( linear(j) ) then

C                     Compute inner product <a,p>
                      ajtp = 0.0d0
                      do i = jcsta(j),jcsta(j) + jclen(j) - 1
                          ajtp = ajtp + jcval(i) * p(jcvar(i))
                      end do

                      ajtp = ajtp * rho(j)

C                     Add rho * ajtp * a
                      do i = jcsta(j),jcsta(j) + jclen(j) - 1
                          hp(jcvar(i)) = hp(jcvar(i)) + ajtp * jcval(i)
                      end do

                  end if
              end if
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine coo2csr(m,nnz,alin,acol,aval,alen,asta)

      implicit none

C     SCALAR ARGUMENTS
      integer m,nnz

C     ARRAY ARGUMENTS
      integer acol(nnz),alen(m),alin(nnz),asta(m)
      double precision aval(nnz)

C     This subroutines converts a matrix from coordinate format to
C     compressed sparse row format.

C     LOCAL SCALARS
      integer i,j,col,coltmp,lin,lintmp
      double precision val,valtmp

      do i = 1,m
          alen(i) = 0
      end do

      do i = 1,nnz
          lin = alin(i)
          alen(lin) = alen(lin) + 1
      end do

      asta(1) = 1
      do i = 2,m
          asta(i) = asta(i-1) + alen(i-1)
      end do

      do i = 1,nnz

          val = aval(i)
          col = acol(i)
          lin = alin(i)

          alin(i) = - 1

 10       if ( lin .ge. 0 ) then

              j = asta(lin)
              asta(lin) = j + 1

              valtmp = aval(j)
              coltmp = acol(j)
              lintmp = alin(j)

              aval(j) = val
              acol(j) = col
              alin(j) = - 1

              val = valtmp
              col = coltmp
              lin = lintmp

              go to 10

          end if

      end do

      do i = 1,m
          asta(i) = asta(i) - alen(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      logical function sstop(n,x,m,lambda,rho,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

      include "algparam.inc"

C     EXTERNAL FUNCTIONS
      logical minsqstop

      if ( innercall ) then
          sstop = minsqstop(n,x,inform)
          return
      end if

      end
