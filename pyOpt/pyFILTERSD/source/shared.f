
c  this file contains subordinate routines shared by glcpd.f, qlcpd.f and
c  ql1lcpd.f

      block data defaults
      implicit double precision (a-h,o-z)
      common/epsc/eps,tol,emin
      common/repc/sgnf,nrep,npiv,nres
      common/refactorc/mc,mxmc
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      data  eps,    tol,   emin, sgnf, nrep, npiv, nres, mxmc, kk, ll
     * /1111.D-19, 1.D-12, 0.D0, 1.D-8,  2,    3,   2,   500,   0,  0/
      end

      subroutine optest(jmin,jmax,r,e,ls,rp,pj)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension r(*),e(*),ls(*)
      rp=0.D0
      do 1 j=jmin,jmax
        i=abs(ls(j))
        if(e(i).eq.0.D0)print *,'e(i).eq.0.D0: i =',i
        ri=r(i)/e(i)
        if(ri.ge.rp)goto1
c       if((i/2)*2.eq.i)then
c         ic=i-1
c  this change is needed when nm is odd
c         ic=i+1
c       else
c         ic=i+1
c         ic=i-1
c       endif
c       do jj=jmin,jmax
c         if(ls(jj).eq.ic)goto2
c       enddo
c       goto1
c   2   continue
        rp=ri
        pj=j
    1 continue
      return
      end

      subroutine form_Ats(jmin,jmax,n,plus,a,la,an,w,ls,tol)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),an(*),w(*),ls(*)
      logical plus
c  form At.s and denominators
      if(plus)then
        do j=jmin,jmax
          i=abs(ls(j))
          if(i.gt.n)then
            wi=aiscpr(n,a,la,i-n,an,0.D0)
          else
            wi=an(i)
          endif
          if(wi.ne.0.D0)then
            if(abs(wi).le.tol)then
              wi=0.D0
            elseif(ls(j).ge.0)then
              wi=-wi
            endif
          endif
          w(i)=wi
        enddo
      else
        do j=jmin,jmax
          i=abs(ls(j))
          if(i.gt.n)then
            wi=aiscpr(n,a,la,i-n,an,0.D0)
          else
            wi=an(i)
          endif
          if(wi.ne.0.D0)then
            if(abs(wi).le.tol)then
              wi=0.D0
            elseif(ls(j).lt.0)then
              wi=-wi
            endif
          endif
          w(i)=wi
        enddo
      endif
      return
      end

      subroutine zprod(k,n,a,la,an,r,w,ls,aa,ll)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),an(*),r(*),w(*),ls(*),aa(*),ll(*)
      common/noutc/nout
      do j=1,n-k
        w(abs(ls(j)))=0.D0
      enddo
      do j=n-k+1,n
        w(ls(j))=-r(ls(j))
      enddo
      call tfbsub(n,a,la,0,w,an,aa,ll,ep,.false.)
      return
      end

      subroutine signst(n,r,w,ls)
      implicit double precision (a-h,o-z)
      dimension r(*),w(*),ls(*)
c  transfer with sign change as necessary
        do j=1,n
          i=abs(ls(j))
          if(ls(j).ge.0)then
            r(i)=w(i)
          else
            r(i)=-w(i)
          endif
        enddo
      return
      end

      subroutine warm_start(n,nm,a,la,x,bl,bu,b,ls,aa,ll,an,vstep)
      implicit double precision (a-h,o-z)
      dimension a(*),la(*),x(*),bl(*),bu(*),b(*),ls(*),
     *  aa(*),ll(*),an(*)
      DOUBLE PRECISION daiscpr
      common/epsc/eps,tol,emin
      common/noutc/nout
      common/iprintc/iprint
      do j=1,n
        i=abs(ls(j))
        if(i.le.n)then
          b(i)=0.D0
        else
          if(ls(j).ge.0)then
            b(i)=daiscpr(n,a,la,i-n,x,-bl(i))
          else
            b(i)=daiscpr(n,a,la,i-n,x,-bu(i))
          endif
        endif
      enddo
c     print 3,'ls =',(ls(i),i=1,n)
    3 format(A/(15I5))
c     write(nout,1000)'x =',(x(i),i=1,n)
c     write(nout,1000)'b =',(b(i),i=1,nm)
c     write(nout,1000)'r =',(b(abs(ls(j))),j=1,n)
      call tfbsub(n,a,la,0,b,an,aa,ll,ep,.false.)
c     write(nout,1000)'d =',(an(i),i=1,n)
      call linf(n,an,vstep,i)
      if(iprint.ge.1)
     *  write(nout,*)'infinity norm of vertical step =',vstep
      do i=1,n
        x(i)=x(i)-an(i)
      enddo
      do j=1,n
        i=abs(ls(j))
        if(i.le.n)then
          if(x(i).ge.bu(i)-tol)then
            x(i)=bu(i)
            ls(j)=-i
          else
            ls(j)=i
            if(x(i).le.bl(i)+tol)x(i)=bl(i)
          endif
        endif
      enddo
c     write(nout,*)'x =',(x(i),i=1,n)
 1000 format(a/(e16.5,4e16.5))
c1001 format(a/(i4,1x,e11.5,4(i4,1x,e11.5)))
      return
      end

      subroutine residuals(n,jmin,jmax,a,la,x,bl,bu,r,ls,f,g,ninf)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),x(*),bl(*),bu(*),r(*),ls(*),g(*)
      common/epsc/eps,tol,emin
      DOUBLE PRECISION daiscpr,z
      f=0.D0
      ninf=0
      do i=1,n
        g(i)=0.D0
      enddo
      do j=jmin,jmax
        i=abs(ls(j))
        if(i.gt.n)then
          z=daiscpr(n,a,la,i-n,x,0.D0)
        else
          z=dble(x(i))
        endif
        ri=z-dble(bl(i))
        ro=dble(bu(i))-z
        if(ri.le.ro)then
          ls(j)=i
        else
          ri=ro
          ls(j)=-i
        endif
        if(abs(ri).le.tol)then
          ri=0.D0
        elseif(ri.lt.0.D0)then
          f=f-ri
          ninf=ninf+1
          if(i.gt.n)then
            call saipy(-sign(1.D0,dble(ls(j))),a,la,i-n,g,n)
          else
            g(i)=g(i)-sign(1.D0,dble(ls(j)))
          endif
        endif
        r(i)=ri
      enddo
      return
      end

      subroutine store_rg(k,ig,G,r,ls)
      implicit double precision (a-h,o-z)
      dimension G(k,*),r(*),ls(*)
      do j=1,k
        G(j,ig)=r(ls(j))
      enddo
      return
      end

      subroutine trid(d,e,n)
      implicit double precision (a-h,o-z)
      dimension d(*),e(*)
c  QL method for eigenvalues of a tridiagonal matrix.  d(i) i=1,..,n
c  is diagonal and e(i) i=1,..,n-1 is subdiagonal (storage for dummy e(n)
c  must be available). Eigenvalues returned in d
c  This routine is adapted from the EISPAK subroutine imtq11.
      if(n.eq.1)return
      e(n)=0.D0
      do 15 l=1,n
      iter=0
    1 continue
      do 12 m=l,n-1
      dd=abs(d(m))+abs(d(m+1))
      if(abs(e(m))+dd.eq.dd)goto2
   12 continue
      m=n
    2 continue
      if(m.eq.l.or.iter.eq.30)goto15
      iter=iter+1
      g=(d(l+1)-d(l))*5.D-1/e(l)
      r=sqrt(g**2+1.D0)
      g=d(m)-d(l)+e(l)/(g+sign(r,g))
      s=1.D0
      c=1.D0
      p=0.D0
      do 14 i=m-1,l,-1
      f=s*e(i)
      b=c*e(i)
      if(abs(f).ge.abs(g))then
        c=g/f
        r=sqrt(c**2+1.D0)
        e(i+1)=f*r
        s=1.D0/r
        c=c*s
      else
        s=f/g
        r=sqrt(s**2+1.D0)
        e(i+1)=g*r
        c=1.D0/r
        s=s*c
      endif
      g=d(i+1)-p
      r=(d(i)-g)*s+2.*c*b
      p=s*r
      d(i+1)=g+p
      g=c*r-b
   14 continue
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.D0
      goto1
   15 continue
      return
      end

      subroutine formR(nv,k,ig,maxg,a,b,c,d,e,G,R)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*),c(*),d(*),e(*),G(k,*),R(*)
c     print *,'G and ig',ig
c     do i=1,k
c       print 5,(G(i,j),j=1,6)
c     enddo
c     print 4,'a =',(a(i),i=1,6)
c     print 4,'b =',(b(i),i=1,6)
c     print 4,'c =',(c(i),i=1,6)
   10 ii=ig-nv
      if(ii.lt.0)ii=ii+maxg
      do i=1,nv
        ii=ii+1
        if(ii.gt.maxg)ii=1
        e(i)=a(ii)
        nvi=nv-i+1
        nvi1=nvi+1
        d(1)=b(ii)
        d(2)=c(ii)
        jj=ii+1
        if(jj.gt.maxg)jj=1
        do j=3,nvi1
          jj=jj+1
          if(jj.gt.maxg)jj=1
          d(j)=scpr(0.D0,G(1,jj),G(1,ii),k)
        enddo
        ij=i
        do j=1,i-1
          call mysaxpy(-R(ij),R(ij),d,nvi1)
          ij=ij+maxg-j
        enddo
        if(d(1).le.0.D0)then
c         print *,'R is singular'
          nv=i-1
          goto10
        endif
        t=sqrt(d(1))
        R(ij)=t
        do j=1,nvi
          R(ij+j)=d(j+1)/t
        enddo
c       print 4,'row of R =',(R(ij+j),j=0,nvi)
      enddo
c     print 4,'R matrix',(R(j),j=1,nv+1)
c     ii=maxg
c     do i=1,nv-1
c       print 5,(0.D0,j=1,i),(R(ii+j),j=1,nv-i+1)
c       ii=ii+maxg-i
c     enddo
      return
    2 format(A,5E13.5)
    4 format(A/(6E13.5))
    5 format((6E13.5))
      end

      subroutine formT(n,nmax,R,d,e)
      implicit double precision (a-h,o-z)
      dimension R(*),d(*),e(*)
c  forms the tridiagonal matrix  T = [R | Qt.gp].S.R^(-1) in d and e
      t=e(1)*R(2)/R(1)
      d(1)=e(1)-t
      ir=1
      irp=nmax+1
      do i=2,n
        im=i-1
        e(im)=-e(im)*R(irp)/R(ir)
        ir=irp
        irp=irp+nmax-im
        dii=t
        t=e(i)*R(ir+1)/R(ir)
        d(i)=dii+e(i)-t
      enddo
      return
      end

      subroutine insort(nv,v)
      implicit double precision (a-h,o-z)
      dimension v(*)
c  insertion sort into ascending order
      do i=2,nv
        t=v(i)
        do j=i-1,1,-1
          if(v(j).gt.t)then
            v(j+1)=v(j)
          else
            v(j+1)=t
            goto10
          endif
        enddo
        v(1)=t
   10   continue
      enddo
      return
      end

      subroutine checkT(n,nmax,R,a,d)
      implicit double precision (a-h,o-z)
      dimension R(*),a(*),d(*)
      dimension T(10,10)
      if(n.gt.9)print *,'Increase dimension of T'
      if(n.gt.9)stop
      do j=1,n+1
        T(1,j)=R(j)
      enddo
      ii=nmax
      do i=1,n-1
        do j=1,i
          T(i+1,j)=0.D0
        enddo
        do j=1,n-i+1
          T(i+1,j+i)=R(ii+j)
        enddo
        ii=ii+nmax-i
      enddo
c     print 4,'R matrix'
c     do i=1,n
c       print 5,(T(i,j),j=1,n+1)
c     enddo
c     print 2,'a =',(a(i),i=1,n)
      do i=1,n
        call mysaxpy(-1.D0,T(1,i+1),T(1,i),n)
        do j=1,n
          T(j,i)=T(j,i)*a(i)
        enddo
      enddo
c     print 4,'R*J matrix'
c     do i=1,n
c       print 5,(T(i,j),j=1,n)
c     enddo
      print 4,'T matrix'
      do i=1,n
        do j=1,n
          d(j)=T(i,j)
        enddo
        call rtsol(n,nn,nmax,R,d)
        print 5,(d(j),j=1,n)
      enddo
      return
    2 format(A,6E15.7)
    4 format(A/(5E15.7))
    5 format((5E15.7))
      end
