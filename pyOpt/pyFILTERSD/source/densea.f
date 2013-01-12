christen this file    denseA.f

c  Copyright (C) 1996 Roger Fletcher

c  Current version dated 21 May 1998

c  ******************************************
c  Specification of A in dense matrix format
c  ******************************************

c  The matrix A contains gradients of the linear terms in the objective
c  function (column 0) and the general constraints (columns 1:m).
c  No explicit reference to simple bound constraints is required in A.
c  The information is set in the parameters a(*) and la.

c  In this dense case A is set in standard matrix format as a(la,0:m), where la
c  is the stride between columns. la is an integer which must be greater or
c  equal to n.

c  In the straightforward case that la=n, columns of A follow successively
c  in the space occupied by a(.).


      subroutine saipy(s,a,la,i,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),y(*)
c  saxpy with column i of A
      call mysaxpy(s,a(1,i),y,n)
      return
      end

      subroutine isaipy(s,a,la,i,y,n,lr,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),y(*),lr(*),li(*)
c  indirectly addressed saxpy with column i of A
      call isaxpy(s,a(1,i),lr,y,n)
      return
      end

      subroutine isaipy1(s,a,la,i,y,n,lr,li,m1)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),y(*),lr(*),li(*)
c  indirectly addressed saxpy with column i of A_1
      call isaxpy(s,a(1,i),lr,y,m1)
      return
      end

c     subroutine isaipy2(s,a,la,i,y,n,lr,li,m1)
c     implicit double precision (a-h,r-z), integer (i-q)
c     dimension a(la,0:*),y(*),lr(*),li(*)
c  indirectly addressed saxpy with column i of A_2
c     call isaxpy(s,a(1,i),lr(m1+1),y(m1+1),n-m1)
c     return
c     end

c     subroutine ssaipy(s,a,la,i,y,n)
c     implicit double precision (a-h,r-z), integer (i-q)
c     dimension a(la,0:*),y(*)
c  ssaxpy with column i of A
c     call ssaxpy(s,a(1,i),y,n)
c     return
c     end

c     subroutine ssaxpy(a,x,y,n)
c     implicit double precision (a-h,r-z), integer (i-q)
c     dimension x(*),y(*)
c  saxpy with squares of x
c     do i=1,n
c       y(i)=y(i)+a*x(i)**2
c     enddo
c     return
c     end

      function aiscpr(n,a,la,i,x,b)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*)
c  scalar product with column i of A
      aiscpr=scpr(b,a(1,i),x,n)
      return
      end

      function daiscpr(n,a,la,i,x,b)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*)
      DOUBLE PRECISION daiscpr,dscpr
      daiscpr=dscpr(b,a(1,i),x,n)
      return
      end

      function aiscpri(n,a,la,i,x,b,lr,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),lr(*),li(*)
c  indirectly addressed scalar product with column i of A
      aiscpri=scpri(b,a(1,i),lr,x,n)
      return
      end

      function daiscpri(n,a,la,i,x,b,lr,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),lr(*),li(*)
      DOUBLE PRECISION daiscpri,dscpri
      daiscpri=dscpri(b,a(1,i),lr,x,n)
      return
      end

      function aiscpri1(n,a,la,i,x,b,lr,li,m1)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),lr(*),li(*)
c  indirectly addressed scalar product with column i of A_1
      aiscpri1=scpri(b,a(1,i),lr,x,m1)
      return
      end

      function aiscpri2(n,a,la,i,x,b,lr,li,m1)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),lr(*),li(*)
c  indirectly addressed scalar product with column i of A_2
      aiscpri2=scpri(b,a(1,i),lr(m1+1),x(m1+1),n-m1)
      return
      end

      function ailen(n,a,la,i)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*)
c  L2 length of column i of A
      ailen=scpr(0.D0,a(1,i),a(1,i),n)
      ailen=sqrt(ailen)
      return
      end

      subroutine iscatter(a,la,i,li,an,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),li(*),an(*)
c  indirect scatter into vector an
      do j=1,n
        an(li(j))=a(j,i)
      enddo
      return
      end

      subroutine iunscatter(a,la,i,li,an,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),li(*),an(*)
c  included for compatibility with sparseA.f
      return
      end

      function aij(i,j,a,la)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*)
c  get element A(i,j)
      aij=a(i,j)
      return
      end

      subroutine setaij(aij,i,j,a,la)
      implicit double precision (a-h,o-z)
      dimension a(la,0:*)
c  set element A(i,j)
      aij=a(i,j)
      end

      subroutine isaxpy(a,x,lr,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension x(*),lr(*),y(*)
c  saxpy with x indirectly addressed
      if(a.eq.0.D0)return
      do i=1,n
        y(i)=y(i)+a*x(lr(i))
      enddo
      return
      end

      function dscpr(a,x,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension x(*),y(*)
      DOUBLE PRECISION dscpr
      dscpr=dble(a)
      do i=1,n
        dscpr=dscpr+dble(x(i))*dble(y(i))
      enddo
      return
      end

      function scpri(a,x,lr,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension x(*),lr(*),y(*)
c  scpr with x indirectly addressed
      scpri=a
      do i=1,n
        scpri=scpri+x(lr(i))*y(i)
      enddo
      return
      end

      function dscpri(a,x,lr,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension x(*),lr(*),y(*)
      DOUBLE PRECISION dscpri
      dscpri=dble(a)
      do i=1,n
        dscpri=dscpri+dble(x(lr(i)))*dble(y(i))
      enddo
      return
      end

      subroutine cscale(n,m,a,la,x,bl,bu,s,menu,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),bl(*),bu(*),s(*)

c     Constraint scaling procedure for use prior to calling bqpd when using
c     denseA.f

c     The user must set the parameter menu to control how the
c     x-variables are scaled (or equivalently how constraints i = 1:n
c     are scaled), as follows

c     menu = 1 indicates that a unit scaling applies to the x-variables

c     menu = 2 the user provides estimates s(i)>0 of the magnitude of
c              x(i) for i = 1:n. In this case the elements  x(i), bl(i), bu(i)
c              are divided by s(i) for i = 1:n.

c     In all cases, cscale goes on to scale the general constraints, in
c     such a way that the normal vector of each nontrivial constraint in
c     the scaled problem has an l_2 norm of unity. This scaling is also
c     applied to the right hand sides  bl(i), bu(i) for i = n+1:n+m.
c     The scaled data overwrites the original data.

c     cscale also scales the constant vector of the quadratic function,
c     which is found in a(1:n). However if a non-unit x-variable scaling
c     is used, it is necessary for the user to scale the Hessian matrix
c     G appropriately. This can be done by passing the x-variable scale
c     factors s(i) i = 1:n into the subroutine gdotx using the
c     parameter ws, and multiplying G(i,j) by s(i)*s(j) (possibly
c     implicitly).

c     cscale sets ifail = 1 to indicate that some s(i)< = 0,
c             and ifail = 2 to indicate an incorrect setting of menu.
c       Otherwise ifail = 0.

      ifail=2
      if(menu.lt.1.or.menu.gt.2)return
c     z=1.D0/log(2.D0)
      if(menu.eq.1)then
        do j=1,n
          s(j)=1.D0
        enddo
      else
        ifail=1
        do j=1,n
          if(s(j).le.0.D0)return
        enddo
c       if(menu.eq.2)then
c         do j=1,n
c           s(j)=2.D0**nint(log(s(j))*z)
c         enddo
c       endif
        do j=1,n
          if(s(j).ne.1.D0)then
            x(j)=x(j)/s(j)
            bl(j)=bl(j)/s(j)
            bu(j)=bu(j)/s(j)
            a(j,0)=a(j,0)*s(j)
          endif
        enddo
      endif
      do i=1,m
        t=0.D0
        do j=1,n
          a(j,i)=a(j,i)*s(j)
          t=t+a(j,i)**2
        enddo
        t=sqrt(t)
        if(t.eq.0.D0)then
          s(n+i)=1.D0
        else
c         t=2.D0**nint(log(t)*z)
          s(n+i)=t
          do j=1,n
            a(j,i)=a(j,i)/t 
          enddo
          bl(n+i)=bl(n+i)/t
          bu(n+i)=bu(n+i)/t
        endif
      enddo
      ifail=0
      return
      end
