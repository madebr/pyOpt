christen this file denseL.f

c  Copyright (C) 1996 Roger Fletcher

c  Current version dated 4 October 2011

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

c***************** dense matrix routines for manipulating L ********************

c  ***************************************************************
c  Basis matrix routines for bqpd with dense matrices (block form)
c  ***************************************************************

c  These routines form and update L-Implicit-U factors LPB=U of a matrix B
c  whose columns are the normal vectors of the active constraints. In this
c  method only the unit lower triangular matrix L and the diagonal of U (in
c  addition to the row permutation P) is stored. B is represented in block form

c    | A_1  0 |    where the first m1 columns (A_1 and A_2) come from the
c    | A_2  I |    general constraint normals (columns of the matrix A in bqpd)

c  and the remaining unit columns come from simple bounds. The matrix A may be
c  specified in either dense or sparse format and the user is referred to the
c  files  denseA.f  or  sparseA.f. About m1*m1/2 locations are required to store
c  L-Implicit-U factors of B. The user MUST supply an upper bound on m1 by
c  setting mxm1 in the labelled common block

c     common/mxm1c/mxm1

c  Setting  mxm1=min(m+1,n)  is always sufficient.

c  Workspace
c  *********
c  denseL.f requires
c     mxm1*(mxm1+1)/2+3*n+mxm1   locations of real workspace, and
c     n+mxm1+n+m                 locations of integer workspace
c  These are stored at the end of the workspace arrays ws and lws in bqpd.
c  The user MUST set the lengths of these arrays in mxws and mxlws in
c     common/wsc/kk,ll,kkk,lll,mxws,mxlws
c  along with the values kk and ll of space to be used by gdotx.

c  Other information
c  *****************

c  L-Implicit-U factors are updated by a variant of the Fletcher-Matthews
c  method, which has proved very reliable in practice. The method is described
c  in the reference
c    Fletcher R., Dense Factors of Sparse Matrices, in "Approximation Theory
c    and Optimization. Tributes to M.J.D. Powell", (M.D. Buhmann and A. Iserles,
c    eds), Cambridge University Press (1997), pp. 145-166.

c  Steepest edge coefficients e(i) are also updated in these routines

c  The file contains routines for solving systems with B or its transpose
c  which might be of use in association with bqpd. These routines are
c  documented below.

      subroutine start_up(n,nm,nmi,a,la,nk,e,ls,aa,ll,mode,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),e(*),ls(*),aa(*),ll(*)
      common/noutc/nout
      common/wsc/kk,ll_,kkk,lll,mxws,mxlws
      common/epsc/eps,tol,emin
      common/densec/ns,ns1,nt,nt1,nu,nu1,mx1,lc,lc1,li,li1
      common/factorc/m0,m1,mm0,mm,mp,mq
      common/refactorc/nup,nfreq
      common/mxm1c/mxm1
      if(mxm1.le.0)then
        write(nout,*)'mxm1 =',mxm1,' is not set correctly'
        ifail=7
        return
      endif
      ns=kk+kkk+mxm1*(mxm1+1)/2+3*n+mxm1
      nt=ll_+lll+n+mxm1+nmi
      if(ns.gt.mxws.or.nt.gt.mxlws)then
        write(nout,*)'not enough real (ws) or integer (lws) workspace'
        write(nout,*)'you give values for mxws and mxlws as',mxws,mxlws
        write(nout,*)'minimum values for mxws and mxlws are',ns,nt
        ifail=7
        return
      endif
      nup=0
      small=max(1.D1*tol,sqrt(eps))
      smallish=max(eps/tol,1.D1*small)
c  set storage map for dense factors
      ns=mxm1*(mxm1+1)/2
      ns1=ns+1
      nt=ns+n
      nt1=nt+1
      nu=nt+n
      nu1=nu+1
      mx1=nu1+n
      lc=n
      lc1=lc+1
      li=lc+mxm1
      li1=li+1
c     write(nout,*)'ls',(ls(ij),ij=1,nk)
c     write(nout,*)'ls',(ls(ij),ij=nm+1,nmi)
      if(mode.ge.3)then
        call re_factor(n,nm,a,la,aa,aa(ns1),aa(nt1),ll,ll(lc1),ll(li1))
        call check_L(n,aa,ifail)
        if(ifail.eq.1)then
          mode=2
          goto1
        endif
        if(nk.eq.n)return
c  reset ls from e
        do j=1,nk
          i=-ls(j)
          if(i.gt.0)e(i)=-e(i)
        enddo
        j=0
        nk=nmi
        do i=1,nmi
          if(e(i).ne.0.D0)then
            j=j+1
            if(e(i).gt.0.D0)then
              ls(j)=i
            else
              ls(j)=-i
              e(i)=-e(i)
            endif
          else
            ls(nk)=i
            nk=nk-1
          endif
        enddo
        if(j.ne.n)then
          write(nout,*)'malfunction in reset sequence in start_up'
          stop
        endif
        ifail=0
        return
      endif
    1 continue
      if(emin.eq.0.D0)then
c  set a lower bound on e(i)
        emin=1.D0
        do i=1,nmi-n
          emin=max(emin,ailen(n,a,la,i))
        enddo
        emin=1.D0/emin
      endif
      do i=1,n
        e(i)=1.D0
        ll(i)=i
      enddo
      do i=n+1,nmi
        e(i)=0.D0
        ll(li+i)=0
      enddo
c  shift designated bounds to end
      nn=n
      do j=nk,1,-1
        i=abs(ls(j))
        if(i.eq.0.or.i.gt.nmi)then
          write(nout,*)
     *      'ls(j) is zero, or greater in modulus than n+m, for j =',j
          ifail=4
          return
        endif
        if(i.le.n)then
          ls(j)=ls(nk)
          nk=nk-1
          call iexch(ll(nn),ll(i))
          nn=nn-1
        endif
      enddo
      do i=1,n
        ll(li+ll(i))=i
      enddo
      m0=(max(mxm1-nk,0))/2
      mm0=m0*(m0+1)/2
      m1=0
      mm=mm0
      j=1
    2 continue
        if(j.gt.nk)goto3
        q=abs(ls(j))
c  extend factors
        call aqsol(n,a,la,q,aa,aa(nt1),aa(mx1),aa,ll,ll(lc1),ll(li1))
        m1p=m1+1
        call linf(nn-m1,aa(nt+m1p),z,iz)
        iz=iz+m1
        if(z.le.tol)then
c         write(nout,*)'reject c/s',q
          nk=nk-1
          do ij=j,nk
            ls(ij)=ls(ij+1)
          enddo
          goto2
        endif
        if(m1p.gt.mxm1)then
          write(nout,*)'mxm1 =',mxm1,'  is insufficient'
          ifail=7
          return
        endif
        if(iz.gt.m1p)then
c  pivot interchange
          ll(li+ll(m1p))=iz
          call iexch(ll(m1p),ll(iz))
          call rexch(aa(nt+m1p),aa(nt+iz))
          ll(li+ll(m1p))=m1p
        endif
        p=ll(m1p)
        tp=aa(nt+m1p)
        call eptsol(n,a,la,p,a,aa,aa(ns1),aa(nt1),ll,ll(lc1),ll(li1))
        aa(ns+m1p)=1.D0
c  update steepest edge coefficients
        ep=e(p)
c       eq=ep/tp
        eq=abs(ep/tp)
        tp_=tp/ep
        tpsq=tp_**2
        call aqsol(n,a,la,-1,a,aa(nu1),aa(mx1),aa,ll,ll(lc1),ll(li1))
        do i=1,m1p
          aa(nu+i)=aa(ns+i)/ep
        enddo
        do i=m1p+1,n
          aa(nu+i)=0.D0
        enddo
        e(p)=0.D0
        do i=1,nmi
          if(e(i).gt.0.D0)then
            ij=ll(li+i)
            ei=e(i)
c           ti=aa(nt+ij)*eq/ei
c           e(i)=max(emin,ei*sqrt(max(1.D0-ti*(2.D0*aa(nu+ij)/ei-ti),0.D0)))
            ti=aa(nt+j)/ei
            e(i)=max(emin,
     *        ei*sqrt(max(tpsq-ti*(2.D0*tp*aa(nu+j)/ei-ti),0.D0))*eq)
          endif
        enddo
c       e(q)=max(emin,abs(eq))
        e(q)=max(emin,eq)
        m1=m1p
        mm=mm+m0
        do ij=1,m1
          aa(mm+ij)=aa(ns+ij)
        enddo
        ll(lc+m1)=q
        ll(li+q)=m1
        mm=mm+m1
        aa(mm)=tp
        j=j+1
        goto2
    3 continue
c  complete the vector ls
      do i=nn+1,n
        nk=nk+1
        ls(nk)=ll(i)
      enddo
      j=nk
      do i=m1+1,nn
        j=j+1
        ls(j)=ll(i)
      enddo
      do j=nm+1,nmi
        e(abs(ls(j)))=1.D0
      enddo
      j=n
      do i=1,nmi
        if(e(i).eq.0.D0)then
          j=j+1
          ls(j)=i
        endif
      enddo
      do j=nm+1,nmi
        e(abs(ls(j)))=0.D0
      enddo
      if(mode.gt.2)then
        z=sqrt(eps)
        do j=1,n
          i=abs(ls(j))
          e(i)=max(z,e(i))
        enddo
        do j=n+1,nmi
          i=abs(ls(j))
          e(i)=0.D0
        enddo
      endif
c     write(nout,*)'e =',(e(ij),ij=1,nmi)
c     write(nout,*)'PAQ factors'
c     ij=mm0+m0
c     do ii=1,m1
c       write(nout,*)(aa(ij+j),j=1,ii)
c       ij=ij+m0+ii
c     enddo
c     write(nout,*)'m0,mm0,m1,mm',m0,mm0,m1,mm
c     write(nout,*)'ls',(ls(ij),ij=1,nmi)
c     write(nout,*)'row perm',(ll(ij),ij=1,n)
c     write(nout,*)'column perm',(ll(lc+ij),ij=1,m1)
c     write(nout,*)'inverse perm',(ll(li+ij),ij=1,nmi)
c     call checkout(n,a,la,aa,ll,ll(lc1),ll(li1))
      mp=-1
      mq=-1
      ifail=0
      return
      end

      subroutine refactor(n,nm,a,la,aa,ll,ifail)
      implicit double precision (a-h,o-z)
      dimension a(*),la(*),aa(*),ll(*)
      common/densec/ns,ns1,nt,nt1,nu,nu1,mx1,lc,lc1,li,li1
      common/factorc/m0,m1,mm0,mm,mp,mq
c     write(nout,*)'refactor'
      call re_factor(n,nm,a,la,aa,aa(ns1),aa(nt1),ll,ll(lc1),ll(li1))
      call check_L(n,aa,ifail)
      return
      end

      subroutine pivot(p,q,n,nm,a,la,e,aa,ll,ifail,info)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),e(*),aa(*),ll(*),info(*)
      common/noutc/nout
      common/iprintc/iprint
      common/densec/ns,ns1,nt,nt1,nu,nu1,mx1,lc,lc1,li,li1
      common/factorc/m0,m1,mm0,mm,mp,mq
      common/mxm1c/mxm1
      common/refactorc/nup,nfreq
      common/epsc/eps,tol,emin
c     write(nout,*)'pivot: p,q =',p,q
      ifail=0
      if(p.ne.mp)then
        call eptsol(n,a,la,p,a,aa,aa(ns1),aa(nt1),ll,ll(lc1),ll(li1))
        e(p)=sqrt(scpr(0.D0,aa(ns1),aa(ns1),m1+1))
        mp=p
      endif
      if(q.ne.mq)then
        call aqsol(n,a,la,q,a,aa(nt1),aa(mx1),aa,ll,ll(lc1),ll(li1))
        mq=q
      endif
c  update steepest edge coefficients
      tp=aa(nt+ll(li+p))
      if(tp.eq.0.D0)tp=eps
      ep=e(p)
c     eq=ep/tp
      eq=abs(ep/tp)
      tp=tp/ep
      tpsq=tp**2
      do i=1,m1+1
        aa(nu+i)=aa(ns+i)/ep
      enddo
      do i=m1+2,n
        aa(nu+i)=0.D0
      enddo
      call aqsol(n,a,la,-1,a,aa(nu1),aa(mx1),aa,ll,ll(lc1),ll(li1))
c     write(nout,*)'row perm',(ll(ij),ij=1,n)
c     write(nout,*)'column perm',(ll(lc+ij),ij=1,m1)
c     write(nout,*)'s =',(aa(ns+ij),ij=1,n)
c     write(nout,*)'t =',(aa(nt+ij),ij=1,n)
c     write(nout,*)'u =',(aa(nu+ij),ij=1,n)
      e(p)=0.D0
      do i=1,nm
        if(e(i).gt.0.D0)then
          j=ll(li+i)
          ei=e(i)
c         ti=aa(nt+j)*eq/ei
c         e(i)=max(emin,ei*sqrt(max(1.D0-ti*(2.D0*aa(nu+j)/ei-ti),0.D0)))
          ti=aa(nt+j)/ei
          e(i)=max(emin,
     *      ei*sqrt(max(tpsq-ti*(2.D0*tp*aa(nu+j)/ei-ti),0.D0))*eq)
        endif
      enddo
c     e(q)=max(emin,abs(eq))
      e(q)=max(emin,eq)
      info(1)=info(1)+1
      if(nup.ge.nfreq)then
c  refactorize L
        ip=ll(li+p)
        if(p.gt.n)then
          qq=ll(lc+m1)
          ll(lc+ip)=qq
          ll(li+qq)=ip
          m1=m1-1
          ll(li+p)=0
        else
          m1p=m1+1
          ll(ip)=ll(m1p)
          ll(li+ll(ip))=ip
          ll(m1p)=p
          ll(li+p)=m1p
        endif
        if(q.gt.n)then
          if(m1.eq.mxm1)then
            ifail=7
            return
          endif
          m1=m1+1
          ll(lc+m1)=q
          ll(li+q)=m1
        else
          iq=ll(li+q)
          m1p=m1+1
          ll(iq)=ll(m1p)
          ll(li+ll(iq))=iq
          ll(m1p)=q
          ll(li+q)=m1p
        endif
        call re_factor(n,nm,a,la,aa,aa(ns1),aa(nt1),ll,ll(lc1),ll(li1))
      else
c  update L
        nup=nup+1
        if(p.le.n)then
          if(m1.eq.mxm1)then
            ifail=7
            return
          endif
          call linf(m1,aa(ns1),z,iz)
          if(z.le.4.D0)then
            if(m0+m1.eq.mxm1)then
c             write(nout,*)'m0 + m1 = mxm1:  re-centre triangle'
              ii=mm0
              mo=m0
              m0=m0/2
              mm0=m0*(m0+1)/2
              mm=mm0
              do i=1,m1
                ii=ii+mo+i
                mm=mm+m0+i
                do j=1-i,0
                  aa(mm+j)=aa(ii+j)
                enddo
              enddo
            endif
            do i=1,m1
              aa(mm+m0+i)=aa(ns+i)
            enddo
            goto1
          endif
        endif
        call c_flma(n,a,la,p,aa,ll,ll(lc1),ll(li1))
    1   continue
        if(q.le.n)then
          call r_flma(n,a,la,q,aa,ll,ll(lc1),ll(li1))
        else
          m1=m1+1
          mm=mm+m0+m1
          aa(mm)=1.D0
          aa(mm)=aiscpri1(n,a,la,q-n,aa(mm-m1+1),0.D0,ll,ll(li1),m1)
          if(abs(aa(mm)).le.eps)aa(mm)=eps
          ll(lc+m1)=q
          ll(li+q)=m1
        endif
        mp=-1
        mq=-1
      endif
      call check_L(n,aa,ifail)
c     write(nout,*)'PAQ factors'
c     ij=m0+mm0
c     do ii=1,m1
c       write(nout,*)(aa(ij+j),j=1,ii)
c       ij=ij+m0+ii
c     enddo
c     write(nout,*)'m0,mm0,m1,mm',m0,mm0,m1,mm
c     write(nout,*)'row perm',(ll(ij),ij=1,n)
c     write(nout,*)'column perm',(ll(lc+ij),ij=1,m1)
c     write(nout,*)'inverse perm',(ll(li+ij),ij=1,nm)
c     call checkout(n,a,la,aa,ll,ll(lc1),ll(li1))
c     write(nout,*)'steepest edge coefficients',(e(ij),ij=1,nm)
c     emax=0.D0
c     do i=1,nm
c       if(e(i).gt.0.D0)then
c         call eptsol(n,a,la,i,a,aa,aa(ns1),aa(nt1),ll,ll(lc1),ll(li1))
c         ei=sqrt(scpr(0.D0,aa(ns1),aa(ns1),n))
c         emax=max(emax,abs(ei-e(i)))
c       endif
c     enddo
c     if(emax.ge.tol)
c    *  write(nout,*)'error in steepest edge coefficients =',emax
      return
      end

      subroutine fbsub(n,jmin,jmax,a,la,q,b,x,ls,aa,ll,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
      dimension a(*),la(*),b(*),x(*),ls(*),aa(*),ll(*)

c  solves a system  B.x=b

c  Parameter list
c  **************
c   n   number of variables (as for bqpd)
c   a,la   specification of QP problem data (as for bqpd)
c   jmin,jmax  (see description of ls below)
c   q   an integer which, if in the range 1:n+m, specifies that the rhs vector
c       b is to be column q of the matrix A of general constraint normals.
c       In this case the parameter b is not referenced by fbsub.
c       If q=0 then b is taken as the vector given in the parameter b.
c   b(n)  must be set to the r.h.s. vector b (but only if q=0)
c   x(n+m)  contains the required part of the solution x, set according to the
c       index number of that component (in the range 1:n for a simple bound and
c       n+1:n+m for a general constraint)
c   ls(*)  an index vector, listing the components of x that are required.
c       Only the absolute value of the elements of ls are used (this allows
c       the possibility of using of the contents of the ls parameter of bqpd).
c       Elements of x in the range abs(ls(j)), j=jmin:jmax are set by fbsub.
c       These contortions allow bqpd to be independent of the basis matrix code.
c   aa(*)  real storage used by the basis matrix code (supply the vector
c       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
c   ll(*)  integer storage used by the basis matrix code (supply the vector
c       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
c   save   indicates if fbsub is to save its copy of the solution for possible
c       future use. We suggest that the user only sets save = .false.

      common/noutc/nout
      common/densec/ns,ns1,nt,nt1,nu,nu1,mx1,lc,lc1,li,li1
      common/factorc/m0,m1,mm0,mm,mp,mq
c     write(nout,*)'fbsub  q =',q
      if(save)then
        if(q.ne.mq)then
          call aqsol(n,a,la,q,b,aa(nt1),aa(mx1),aa,ll,ll(lc1),ll(li1))
          mq=q
        endif
        do j=jmin,jmax
          i=abs(ls(j))
          x(i)=aa(nt+ll(li+i))
        enddo
      else
        call aqsol(n,a,la,q,b,aa(nu1),aa(mx1),aa,ll,ll(lc1),ll(li1))
        do j=jmin,jmax
          i=abs(ls(j))
          x(i)=aa(nu+ll(li+i))
        enddo
      endif
      return
      end

      subroutine tfbsub(n,a,la,p,b,x,aa,ll,ep,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
      dimension a(*),la(*),b(*),x(*),aa(*),ll(*)

c  solves a system  Bt.x=b

c  Parameter list
c  **************
c   n   number of variables (as for bqpd)
c   a,la   specification of QP problem data (as for bqpd)
c   p    an integer which, if in the range 1:n+m, specifies that the rhs vector
c        b is a unit vector appropriate to the position of p in the current
c        ordering. In this case b is not referenced by tfbsub.
c   b(n+m)  If p=0, this must be set to the r.h.s. vector b. Only the components
c        of b need be set, according to the index number of each component (in
c        the range 1:n for a simple bound and n+1:n+m for a general constraint)
c   x(n)  contains the solution x (in natural ordering)
c   aa(*)  real storage used by the basis matrix code (supply the vector
c       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
c   ll(*)  integer storage used by the basis matrix code (supply the vector
c       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
c   ep  if p.ne.0 and save is true, ep contains the l_2 length of x on exit
c   save  indicates if tfbsub is to save its copy of the solution for possible
c       future use. We suggest that the user only sets save = .false.

      common/noutc/nout
      common/densec/ns,ns1,nt,nt1,nu,nu1,mx1,lc,lc1,li,li1
      common/factorc/m0,m1,mm0,mm,mp,mq
c     write(nout,*)'tfbsub  p =',p
      if(save)then
        if(p.ne.mp)then
          call eptsol(n,a,la,p,b,aa,aa(ns1),aa(nt1),ll,ll(lc1),ll(li1))
          mp=p
        endif
        do i=1,n
          x(ll(i))=aa(ns+i)
        enddo
        if(p.gt.0)ep=sqrt(scpr(0.D0,aa(ns1),aa(ns1),m1+1))
      else
        call eptsol(n,a,la,p,b,aa,aa(nu1),aa(nt1),ll,ll(lc1),ll(li1))
        do i=1,n
          x(ll(i))=aa(nu+i)
        enddo
      endif
c     write(nout,*)'x =',(x(i),i=1,n)
      return
      end

      subroutine newg
      common/factorc/m0,m1,mm0,mm,mp,mq
      mq=-1
      return
      end

c******** The following routines are internal to denseL.f **************

      subroutine re_factor(n,nm,a,la,T,sn,tn,lr,lc,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),T(*),sn(*),tn(*),lr(*),lc(*),li(*)
      common/noutc/nout
      common/iprintc/iprint
      common/refactorc/nup,nfreq
      common/factorc/m0,m1,mm0,mm,mp,mq
      common/mxm1c/mxm1
      common/epsc/eps,tol,emin
c     write(nout,*)'re_factor'
      nup=0
      if(m1.eq.0)return
      m0=(mxm1-m1)/2
      mm0=m0*(m0+1)/2
c     write(nout,*)'row perm',(lr(ij),ij=1,n)
c     write(nout,*)'column perm',(lc(ij),ij=1,m1)
      do i=1,m1
        sn(i)=0.D0
      enddo
      mm=mm0
      do i=1,m1-1
        mm=mm+m0+i
        im=i-1
        i1=mm-im
        q=lc(i)-n
        if(q.le.0)goto1
c  form L.a_q
        call iscatter(a,la,q,li,sn,n)
c       write(nout,*)'aq =',(sn(ij),ij=1,m1)
        jj=mm
        j1=i1
        do j=i,m1
          tn(j)=scpr(sn(j),T(j1),sn,im)
          j1=jj+m0+1
          jj=j1+j
        enddo
        call iunscatter(a,la,q,li,sn,n)
c       write(nout,*)'L.aq =',(tn(ij),ij=i,m1)
        call linf(m1-im,tn(i),z,iz)
        if(iz.gt.1)then
c  pivot interchange
          iz=iz-1
          call vexch(T(i1),T(i1+iz*(m0+i)+iz*(iz-1)/2),im)
          iz=iz+i
          call rexch(tn(i),tn(iz))
          li(lr(i))=iz
          call iexch(lr(i),lr(iz))
          li(lr(i))=i
        endif
        if(tn(i).eq.0.D0)tn(i)=eps
c  update L
        j1=i1+m0+i
        zz=-tn(i)
        do j=i+1,m1
          z=tn(j)/zz
          call mysaxpy(z,T(i1),T(j1),i-1)
          T(j1+im)=z
c         write(nout,*)'L(j) =',(T(ij),ij=j1,j1+im)
          j1=j1+m0+j
        enddo
        T(mm)=-zz
      enddo
      mm=mm+m0+m1
      q=lc(i)-n
      if(q.le.0)goto1
      call iscatter(a,la,q,li,sn,n)
      T(mm)=scpr(sn(m1),T(mm-m1+1),sn,m1-1)
      if(T(mm).eq.0.D0)T(mm)=eps
c     write(nout,*)'PAQ factors'
c     ij=mm0+m0
c     do ii=1,m1
c       write(nout,*)(T(ij+j),j=1,ii)
c       ij=ij+m0+ii
c     enddo
c     write(nout,*)'m0,mm0,m1,mm',m0,mm0,m1,mm
c     write(nout,*)'row perm',(lr(ij),ij=1,n)
c     write(nout,*)'column perm',(lc(ij),ij=1,m1)
c     write(nout,*)'inverse perm',(li(ij),ij=1,nm)
c     call checkout(n,a,la,T,lr,lc1,li)
      mp=-1
      mq=-1
      return
    1 continue
      write(nout,*)'malfunction in re_factor:  i,lc(i) =',i,q+n
      stop
      end

      subroutine check_L(n,T,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension T(*)
      common/noutc/nout
      common/factorc/m0,m1,mm0,mm,mp,mq
      common/epsc/eps,tol,emin
c     write(nout,*)'check_L'
      ifail=1
      kk=mm0
c     dmin=1.D37
      do k=1,m1
        kk=kk+m0+k
c       dmin=min(dmin,abs(T(kk)))
        if(abs(T(kk)).le.tol)return
      enddo
c     write(nout,*)'dmin =',dmin
      ifail=0
      return
      end

      subroutine aqsol(n,a,la,q,b,tn,xm,T,lr,lc,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),tn(*),xm(*),T(*),lr(*),lc(*),li(*)
      common/noutc/nout
      common/factorc/m0,m1,mm0,mm,mp,mq
c     write(nout,*)'aqsol  q =',q
      if(q.gt.0)then
        do i=1,n
          tn(i)=0.D0
        enddo
        if(q.le.n)then
          tn(li(q))=1.D0
        else
c         call isaipy(1.D0,a,la,q-n,tn,n,lr,li)
          call iscatter(a,la,q-n,li,tn,n)
        endif
      elseif(q.eq.0)then
        do i=1,n
          tn(li(i))=b(i)
        enddo
      endif
c     write(nout,*)'tn =',(tn(i),i=1,n)
      ii=mm
      do i=m1,1,-1
        xm(i)=(scpr(tn(i),T(ii-i+1),tn,i-1))/T(ii)
        call isaipy(-xm(i),a,la,lc(i)-n,tn,n,lr,li)
        ii=ii-m0-i
      enddo
      do i=1,m1
        tn(i)=xm(i)
      enddo
c     write(nout,*)'tn =',(tn(i),i=1,n)
      return
      end

      subroutine eptsol(n,a,la,p,b,T,sn,tn,lr,lc,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),T(*),sn(*),tn(*),lr(*),lc(*),li(*)
      common/noutc/nout
      common/iprintc/iprint
      common/epsc/eps,tol,emin
      common/factorc/m0,m1,mm0,mm,mp,mq
c     write(nout,*)'eptsol  p =',p
c     if(p.eq.9)then
c         write(nout,9)'row perm',(lr(ij),ij=1,n)
c         write(nout,9)'column perm',(lc(ij),ij=1,m1)
c         write(nout,9)'inverse perm',(li(ij),ij=1,p)
c   9     format(A/(15I5))
c     endif
      if(p.gt.n)then
        pr=li(p)
        if(pr.le.0)print *,'here1'
        if(pr.le.0)goto1
        if(pr.ne.m1)then
          z=tn(pr)
          call r_shift(tn(pr),m1-pr,1)
          tn(m1)=z
          call c_flma(n,a,la,p,T,lr,lc,li)
          m1=m1+1
          mm=mm+m0+m1
          li(p)=m1
          lc(m1)=p
          T(mm)=1.D0
          T(mm)=aiscpri1(n,a,la,p-n,T(mm-m1+1),0.D0,lr,li,m1)
          if(T(mm).eq.0.D0)T(mm)=eps
c         write(nout,*)'PAQ factors'
c         ij=m0+mm0
c         do ii=1,m1
c           write(nout,*)(T(ij+j),j=1,ii)
c           ij=ij+m0+ii
c         enddo
c         write(nout,*)'m0,mm0,m1,mm',m0,mm0,m1,mm
c         write(nout,*)'row perm',(lr(ij),ij=1,n)
c         write(nout,*)'column perm',(lc(ij),ij=1,m1)
c         write(nout,*)'inverse perm',(li(ij),ij=1,p)
c         call checkout(n,a,la,T,lr,lc,li)
        endif
        ii=mm-m1
        z=1.D0/T(mm)
        do i=1,m1-1
          sn(i)=T(ii+i)*z
        enddo
        sn(m1)=z
        do i=m1+1,n
          sn(i)=0.D0
        enddo
      else
        ii=m0+mm0
        if(p.eq.0)then
          do i=1,m1
            sn(i)=0.D0
          enddo
          do i=m1+1,n
            sn(i)=b(lr(i))
          enddo
          do i=1,m1
            ii=ii+i
            j=lc(i)
            sn(i)=-aiscpri(n,a,la,j-n,sn,-b(j),lr,li)/T(ii)
            call mysaxpy(sn(i),T(ij),sn,i-1)
            ii=ii+m0
            ij=ii+1
          enddo
        else
          pr=li(p)
          if(pr.le.m1)print *,'here2'
          if(pr.le.m1)goto1
          m1p=m1+1
          call iexch(lr(pr),lr(m1p))
          call iexch(li(lr(pr)),li(lr(m1p)))
          call rexch(tn(pr),tn(m1p))
          do i=1,n
            sn(i)=0.D0
          enddo
          sn(m1p)=1.D0
          do i=1,m1
            ii=ii+i
            sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,0.D0,lr,li)/T(ii)
            call mysaxpy(sn(i),T(ij),sn,i-1)
            ii=ii+m0
            ij=ii+1
          enddo
        endif
      endif
c     write(nout,*)'sn =',(sn(i),i=1,n)
      return
    1 continue
      write(nout,*)'malfunction detected in eptsol: p =',p
      stop
      end

      subroutine c_flma(n,a,la,q,T,lr,lc,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),T(*),lr(*),lc(*),li(*)
      common/noutc/nout
      common/mxm1c/mxm1
      common/epsc/eps,tol,emin
      common/factorc/m0,m1,mm0,mm,mp,mq
      double precision l21
c     write(nout,*)'c_flma: q =',q
      qc=li(q)
      if(q.gt.n)then
        if(qc.le.0)goto1
        call ishift(lc(qc),m1-qc,1)
        do j=qc,m1-1
          li(lc(j))=j
        enddo
        li(q)=0
        mm=mm-m1-m0
        m1=m1-1
      else
        if(qc.le.m1)goto1
        call iexch(lr(qc),lr(m1+1))
        call iexch(li(lr(qc)),li(lr(m1+1)))
        call ishift(lr(2),m1,-1)
        lr(1)=q
        do i=1,m1+1
          li(lr(i))=i
        enddo
        if(m0.eq.0)then
c         write(nout,*)'m0 = 0:  re-centre triangle'
          m0=(mxm1+1-m1)/2
          mm0=m0*(m0+1)/2
          ii=mm
          mm=(m0+m1)*(m0+m1+1)/2
          ii=ii-mm
          ij=mm+m0+1
          do i=m1,1,-1
            ij=ij-m0-i
            call r_shift(T(ij),i,ii)
            ii=ii+m0
          enddo
        endif
        mm=mm-m0-m1
        m0=m0-1
        do i=1,m1
          mm0=mm0+m0+i
          T(mm0)=0.D0
        enddo
        mm0=m0*(m0+1)/2
        qc=1
      endif
      iswap=0
      ii=(qc+m0)*(qc+m0+1)/2
      do i=qc,m1
        im=i+m0
        ii1=ii+m0+1
        iip=ii1+i
        T(ii)=1.D0
        u21=T(iip)
        u11=aiscpri1(n,a,la,lc(i)-n,T(ii1-im),0.D0,lr,li,i)
        ij=ii+im-iswap
c       write(nout,*)'i,im,ii,iip,iswap,ij',i,im,ii,iip,iswap,ij
        l21=T(ij)
        if(abs(l21).le.eps)l21=0.D0
        if(iswap.gt.0)call r_shift(T(ij),iswap,1)
        del=u21-l21*u11
c       write(nout,*)'l21,u11,u21,del =',l21,u11,u21,del
c       write(nout,*)'old row =',(T(j),j=ii1-im,ii)
c       write(nout,*)'new row =',(T(j),j=ii1,ii+im)
        if(abs(del).le.abs(u11)*max(1.D0,abs(l21)))then
c         if(u11.eq.0.D0)then
c           r=0.D0
c         else
            if(u11.eq.0.D0)u11=eps
            r=-u21/u11
            if(abs(r).le.eps)r=0.D0
            call mysaxpy(r,T(ii1-im),T(ii1),i-1)
c         endif
          T(ii)=u11
          T(ii+im)=l21+r
          if(iswap.gt.0)then
            do j=im+1,m0+m1
              ij=ij+j
              r=T(ij)
              call r_shift(T(ij),iswap,1)
              T(ij+iswap)=r
            enddo
          endif
          iswap=0
        else
          r=-u11/del
          if(abs(r).le.eps)r=0.D0
          call permop(T(ii1-im),T(ii1),r,-l21,i-1)
          T(ii)=del
          T(ii+im)=r
          call iexch(lr(i),lr(i+1))
          call iexch(li(lr(i)),li(lr(i+1)))
          iswap=iswap+1
        endif
        ii=iip
      enddo
      return
    1 continue
      write(nout,*)'malfunction detected in c_flma: q =',q
      stop
      end

      subroutine r_flma(n,a,la,p,T,lr,lc,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),T(*),lr(*),lc(*),li(*)
      common/noutc/nout
      common/epsc/eps,tol,emin
      common/factorc/m0,m1,mm0,mm,mp,mq
      double precision l11
c     write(nout,*)'r_flma: p =',p
      pr=li(p)
      if(pr.gt.m1)then
        if(pr.eq.m1+1)return
        write(nout,*)'malfunction detected in r_flma: p =',p
        stop
      endif
      ii=(pr+m0)*(pr+m0+1)/2
      u11=T(ii)
      T(ii)=1.D0
      ip=ii
      do i=pr,m1-1
        im=i+m0
        ii1=ii+m0+1
        iip=ii1+i
        u22=T(iip)
        l11=-T(ip+im)/T(ip)
        if(abs(l11).le.eps)l11=0.D0
        u12=aiscpri1(n,a,la,lc(i+1)-n,T(ii1-im),0.D0,lr,li,i)
        del=l11*u12+u22
c       write(nout,*)'l11,u11,u12,u22,del',l11,u11,u12,u22,del
c       write(nout,*)'old row =',(T(j),j=ii1-im,ii)
c       write(nout,*)'new row =',(T(j),j=ii1,ii+im)
        if(abs(del).le.abs(l11)*max(abs(u11),abs(u12)))then
          call saxpyx(l11,T(ii1-im),T(ii1),i)
          u11=l11*u11
          if(u11.eq.0.D0)u11=eps
          T(iip)=1.D0
        else
          r=-u12/del
          if(abs(r).le.eps)r=0.D0
          call permop(T(ii1-im),T(ii1),r,l11,i)
          call iexch(lc(i),lc(i+1))
          call iexch(li(lc(i)),li(lc(i+1)))
          T(iip)=r
          u22=u11*u22/del
          u11=del
        endif
        call r_shift(T(ip),i-pr,1)
        T(ii)=u11
        u11=u22
        ip=ip+im
        ii=iip
      enddo
      call ishift(lr(pr),m1-pr+1,1)
      lr(m1+1)=p
      do j=pr,m1+1
        li(lr(j))=j
      enddo
c     if(T(ip).eq.0.D0)T(ip)=eps
      l11=-T(ip+m0+m1)/T(ip)
      call saxpyx(l11,T(mm-m1+1),T(mm+m0+1),m1)
      call r_shift(T(ip),m1-pr,1)
      T(mm)=l11*u11
      if(T(mm).eq.0.D0)T(mm)=eps
      return
      end

      subroutine permop(v1,v2,r,s,n)
      implicit double precision (a-h,o-z)
      dimension v1(*),v2(*)
      common/noutc/nout
      if(s.eq.0)then
        if(r.eq.0)then
          call vexch(v1,v2,n)
        else
          do i=1,n
            z=v2(i)
            v2(i)=v1(i)+r*z
            v1(i)=z
          enddo
        endif
      else
        if(r.eq.0)then
          do i=1,n
            z=v1(i)
            v1(i)=v2(i)+s*z
            v2(i)=z
          enddo
        else
          do i=1,n
            z=v1(i)
            v1(i)=v2(i)+s*z
            v2(i)=z+r*v1(i)
          enddo
        endif
      endif
      return
      end

      subroutine checkout(n,a,la,T,lr,lc,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(*),T(*),lr(*),lc(*),li(*)
      common/noutc/nout
      common/mxm1c/mxm1
      common/epsc/eps,tol,emin
      common/factorc/m0,m1,mm0,mm,mp,mq
      emax=0.D0
      gmax=0.D0
      ii=mm0
      do i=1,m1
        ii1=ii+m0+1
        ii=ii+m0+i
        d=T(ii)
        T(ii)=1.D0
        do j=1,i-1
          e=aiscpri1(n,a,la,lc(j)-n,T(ii1),0.D0,lr,li,i)
          emax=max(emax,abs(e))
          gmax=max(gmax,abs(T(ii+m0+j)))
        enddo
        e=aiscpri1(n,a,la,lc(i)-n,T(ii1),-d,lr,li,i)
        emax=max(emax,abs(e))
        T(ii)=d
      enddo
c     if(emax.gt.tol.or.gmax.gt.1.D1)
c    *  write(nout,*)'error in LA=U is ',emax,'  growth in L =',gmax
      write(nout,*)'error in LA=U is ',emax,'  growth in L =',gmax
      return
      end
