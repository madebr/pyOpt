
christen this file l1sold.f
cut here >>>>>>>>>>>>>>>>>>

c  Copyright (C) 2010 Roger Fletcher

      subroutine l1sold(n,m,k,kmax,maxg,a,la,x,bl,bu,f,g,r,w,e,ls,
     * alp,lp,mlp,peq,ws,lws,cws,v,nv,rgtol,ifail,iprint,nout)
      implicit double precision (a-h,r-z), integer (i-q)

c  This routine is a post-processor for qlcpd and glcpd. In the case that the
c  constraint set is infeasible (ifail=3), l1sold finds a best l1 solution of
c  the general constraints, subject to the simple bounds being satisfied.

c  Parameters are a subset of those for glcpd and must be passed through
c  unchanged. For qlcpd a dummy parameter cws must be included.

c  A modified form of Wolfe's method is used to resolve degeneracy.
c  If the solution is degenerate, there may be inactive constraints with
c  zero residual and multiplier 1. Such constraints are marked on exit by
c  setting their residual value (in r(*)) to -eps, (see common/epsc for eps)

      parameter (ainfty=1.D100)
      dimension a(*),la(*),x(*),bl(*),bu(*),g(*),r(*),w(*),e(*),ls(*),
     *  alp(*),lp(*),ws(*),lws(*),v(*)
      character cws(*)
      character*32 spaces
      common/lcpdc/na,na1,nb,nb1,krg,krg1,kr,kr1,
     *  ka,ka1,kb,kb1,kc,kc1,kd,kd1,ke,ke1,lu1,ll1
      common/epsc/eps,tol,emin
      common/infoc/vstep,iter,npv,nfn,ngr
      common/repc/sgnf,nrep,npiv,nres
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/refactorc/nup,nfreq
      common/alphac/alpha,rp,pj,qqj,qqj1
      logical plus

    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(15I5))
    4 format(A/(5E15.7))
    5 format((6E15.7))

c     if(iprint.eq.3)print 4,'a =',(a(i),i=1,110)
      spaces='         '
      n1=n+1
      nm=n+m
      lp(1)=nm
      lev=1
      npv=0
c  collect simple bound equations
      peq=0
      do j=peq+1,n-k
        i=abs(ls(j))
        if(i.le.n.and.bl(i).eq.bu(i))then
          peq=peq+1
          call iexch(ls(j),ls(peq))
        endif
      enddo
      gnorm=sqrt(scpr(0.D0,g,g,n))
      gtol=sgnf*gnorm
      rgtol=max(rgt0l,gtol)
      if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *    ''  level = 1    f ='',E16.8)')npv,f
c     print 4,'gradient =',(g(i),i=1,n)
      goto20
c  start of major iteration
   10 continue
      if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *    ''  level = 1    f ='',E16.8)')npv,f
c  calculate multipliers
c     print 4,'gradient =',(g(i),i=1,n)
      do i=1,nm
        w(i)=0.D0
      enddo
      call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(n,r,w,ls)

   20 continue
      if(iprint.ge.3)then
        write(nout,1001)'costs vector and indices',
     *    (ls(j),r(abs(ls(j))),j=1,n)
c       write(nout,1000)'steepest edge coefficients',
c    *    (e(abs(ls(j))),j=1,n)
        write(nout,1)'# of bound equations and free variables = ',peq,k
      endif
c     if(iprint.eq.3)print 4,'gradient =',(g(i),i=1,n)
c     if(iprint.eq.3)write(nout,1001)'residual vector and indices',
c    *    (ls(j),r(abs(ls(j))),j=n1,lp(1))
c     call check1(n,lp(1),nm,k,kmax,g,a,la,x,bl,bu,r,ls,lp,ws(nb1),f,
c    *  ws,lws,cws,1,p,rp)

   21 continue
c  l1 optimality test
      call optest1(peq,k,n,bl,bu,r,e,ls,rp,pj)

   22 continue
      if(rp.le.gtol)then
c  allow for changes to norm(g)
        gnorm=sqrt(scpr(0.D0,g,g,n))
        gtol=sgnf*gnorm
      endif

      if(rp.le.gtol)then
        do j=peq+1,n-k
          i=abs(ls(j))
          if(i.le.n)then
            if(ls(j).ge.0)then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            endif
          endif
        enddo
        do i=1,n
          x(i)=max(min(x(i),bu(i)),bl(i))
        enddo
        do j=n1,nm
          i=abs(ls(j))
          if(r(i).eq.0.D0.and.i.le.n)then
            if(ls(j).ge.0)then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            endif
          endif
        enddo
c       write(nout,1)'# of simple bound equations = ',peq
        do j=peq+1,n-k
          i=abs(ls(j))
          if(bl(i).eq.bu(i))then
            peq=peq+1
            call iexch(ls(j),ls(peq))
          endif
        enddo
c       write(nout,1001)'costs vector and indices',
c    *    (ls(j),r(abs(ls(j))),j=1,n)
c       write(nout,1)'# of active equations and free variables = ',peq,k
        if(iprint.ge.2)then
          write(nout,*)'OPTIMAL l1 solution'
          if(iprint.ge.3)then
c           write(nout,1000)'x variables',(x(i),i=1,n)
            write(nout,1001)'residual vector and indices',
     *        (ls(j),r(abs(ls(j))),j=n1,nm)
          endif
        endif
        ifail=0
        return
      endif

      p=abs(ls(pj))
      if(iprint.ge.2)write(nout,*)'CHOOSE p,pj =',ls(pj),pj,r(p)
c  compute +/- Steepest Edge search direction s in an(.)
      call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1),
     *  e(p),.true.)

        rp=scpr(0.D0,ws(na1),g,n)
        if(ls(pj).lt.0)rp=-rp
        if(rp*r(p).le.0.D0)then
          print 2,'3rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
          r(p)=0.D0
          goto98
        endif

      if(pj.gt.n-k)then
        rp=-abs(r(p))
        plus=ls(pj).ge.0.eqv.r(p).lt.0.D0
      elseif(r(p).gt.0.D0)then
        rp=1.D0-r(p)
        plus=ls(pj).lt.0
      else
        rp=r(p)
        plus=ls(pj).ge.0
      endif
      snorm=e(p)
c     print 4,'s (or -s if .not.plus) =',(ws(i),i=na1,na+n)

c  form At.s and denominators
      call form_Ats(n1,lp(1),n,plus,a,la,ws(na1),w,ls,snorm*tol)

c  return from degeneracy code
   30 continue
      if(iprint.ge.3)then
c       write(nout,1000)'x variables',(x(i),i=1,n)
        write(nout,1001)'residual vector and indices',
     *    (ls(j),r(abs(ls(j))),j=n1,lp(1))
        write(nout,1000)'denominators',(w(abs(ls(j))),j=n1,lp(1))
      endif
c     print 2,'slope for r(169) =',aiscpr(n,a,la,169-n,ws(na1),0.D0)
c     print 2,'slope for r(217) =',aiscpr(n,a,la,217-n,ws(na1),0.D0)
c     print *,'plus =',plus

   40 continue
c  level 1 ratio tests
      amax=ainfty
      qj=0
      do 41 j=n-k+1,n
        i=ls(j)
        si=ws(na+i)
        t=abs(si)
        if(t.le.tol)goto41
        if(si.gt.0.D0.eqv.plus)then
          z=bu(i)-x(i)
          if(abs(z).lt.tol)then
            z=0.D0
            x(i)=bu(i)
          else
            z=z/t
          endif
        else
          z=x(i)-bl(i)
          if(abs(z).lt.tol)then
            z=0.D0
            x(i)=bl(i)
          else
            z=z/t
          endif
        endif
        if(z.gt.amax)goto41
        amax=z
        qj=j
   41 continue
      if(pj.le.n-k.and.r(p).lt.0.D0.and.bu(p)-bl(p).lt.amax)then
        amax=bu(p)-bl(p)
        qj=pj
      endif
      alpha=amax
      do 44 j=n1,lp(1)
        i=abs(ls(j))
        wi=w(i)
        if(wi.eq.0.D0)goto44
        ri=r(i)
        if(ri.eq.-eps)then
          if(bl(i).eq.bu(i).and.wi.lt.0.D0)then
            alpha=0.D0
            q=i
            qj=j
            goto45
          endif
          goto44
        elseif(ri.lt.0.D0)then
          if(wi.gt.0.D0)goto44
          z=(ri-tol)/wi
        elseif(wi.gt.0.D0)then
          z=(ri+tol)/wi
        else
          z=(bl(i)-bu(i)+ri-tol)/wi
        endif
        if(z.ge.alpha)goto44
        alpha=z
        qj=j
   44 continue
      q=abs(ls(qj))
      if(qj.gt.n)then
        if(r(q).lt.0.D0.eqv.w(q).lt.0.D0)then
          alpha=r(q)/w(q)
        else
          alpha=(bl(q)-bu(q)+r(q))/w(q)
        endif
      endif
   45 continue
      if(iprint.ge.2)then
        write(nout,2)'r(q),w(q) =',r(q),w(q)
        write(nout,*)'alpha =',alpha,'   q =',q
      endif

      if(alpha.eq.0.D0.and.pj.le.n-k)then
        if(iprint.ge.2)
     *    write(nout,*)'degeneracy block at level 1'
        if(w(q).lt.0.D0.and.r(q).ne.-eps)then
          w(q)=-w(q)
          ls(qj)=-ls(qj)
        endif
        alp(1)=f
        plev=n
        do j=n1,lp(1)
          i=abs(ls(j))
          if(r(i).eq.-eps)then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            r(i)=-1.D0
          elseif(r(i).eq.0.D0)then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            r(i)=1.D0
          endif
        enddo
        lp(2)=plev
        lev=2
c       write(nout,1001)'costs vector and indices',
c    *    (ls(j),r(abs(ls(j))),j=1,n)
c       write(nout,1)'# of bound equations and free variables = ',peq,k
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,''     level = 2'',
     *    ''    f ='',E16.8)')npv,f
        goto86
      endif

      if(alpha.gt.0.D0)then
        ff=f
        f=f+alpha*rp
        if(f.ge.ff)then
          if(pj.gt.n-k.or.r(p).lt.0.D0)then
            r(p)=0.D0
          else
            r(p)=1.D0
          endif
          goto20
        endif
        if(plus)then
          call mysaxpy(alpha,ws(na1),x,n)
        else
          call mysaxpy(-alpha,ws(na1),x,n)
        endif
c  update r for inactive c/s
        do 61 j=n1,lp(1)
          i=abs(ls(j))
          if(w(i).eq.0.D0)goto61
          ri=r(i)-alpha*w(i)
          if(abs(ri).le.tol)ri=0.D0
          if(r(i).lt.0.D0.and.ri.ge.0.D0)then
c  remove contribution to gradient
            call saipy(sign(1.D0,dble(ls(j))),a,la,i-n,g,n)
            call newg
          endif
          if(w(i).lt.0.D0)then
            ro=(bu(i)-bl(i))-ri
            if(abs(ro).le.tol)ro=0.D0
            if(ro.lt.ri)then
              ri=ro
c             w(i)=-w(i)
              ls(j)=-ls(j)
            endif
          endif
          if(ri.eq.0.D0.and.i.le.n)then
            if(ls(j).ge.0)then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            endif
          endif
          r(i)=ri
   61   continue
      endif

   70 continue
      if(qj.ne.pj)then
c  pivot interchange
        if(iprint.ge.2)write(nout,*)'replace',p,' by',q
        call pivot(p,q,n,nm,a,la,e,ws(lu1),lws(ll1),ifail,npv)
        if(ifail.ge.1)then
          if(iprint.ge.1)write(nout,*)'near singularity in pivot (1)'
          goto98
        endif
        if(pj.gt.n-k)then
          r(p)=x(p)-bl(p)
          rpu=bu(p)-x(p)
          if(rpu.lt.r(p))then
            r(p)=rpu
            ls(pj)=-p
          else
            ls(pj)=p
          endif
          if(r(p).le.tol)r(p)=0.D0
        elseif(r(p).gt.0.D0)then
          call saipy(-sign(1.D0,dble(ls(pj))),a,la,p-n,g,n)
          call newg
          r(p)=-alpha
        else
          rpu=bu(p)-bl(p)-alpha
          if(abs(rpu).le.tol)rpu=0.D0
          if(alpha.le.rpu)then
            r(p)=alpha
          else
            r(p)=rpu
            ls(pj)=-ls(pj)
          endif
        endif
        if(pj.gt.n-k)then
          k=k-1
          call iexch(ls(pj),ls(n-k))
          pj=n-k
        endif
        call iexch(ls(pj),ls(qj))
        if(q.le.n.and.bl(q).eq.bu(q))then
          peq=peq+1
          call iexch(ls(pj),ls(peq))
        endif
        goto10
      endif
c  opposite bound comes active
c     r(p)=-rp
      if(pj.le.n-k)then
        ls(pj)=-ls(pj)
        goto10
      endif
c  free variable reaches its bound
      if(r(p).lt.0.D0)ls(pj)=-p
      k=k-1
      call iexch(ls(pj),ls(n-k))
      goto10

c  recursive code for resolving degeneracy (Wolfe's method)
   80 continue
c  calculate multipliers
      call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(n,r,w,ls)
c  reset multiplier loop
   82 continue
      if(iprint.ge.3)then
         write(nout,1001)'costs vector and indices',
     *    (ls(j),r(abs(ls(j))),j=1,n)
c       write(nout,1000)'steepest edge coefficients',
c    *    (e(abs(ls(j))),j=1,n)
        write(nout,1)'# of bound equations and free variables = ',peq,k
      endif

   84 continue
c     call check1(n,lp(1),nm,k,kmax,g,a,la,x,bl,bu,r,ls,lp,ws(nb1),f,
c    *  ws,lws,cws,lev,p,rp)

   85 continue
      call optest1(peq,k,n,bl,bu,r,e,ls,rp,pj)

      if(rp.le.gtol.or.pj.gt.n-k)then
        if(iprint.ge.2)write(nout,*)'return to level 1'
        do j=n1,lp(2)
          i=abs(ls(j))
          if(r(i).lt.0.D0)then
            r(i)=-eps
          else
            r(i)=0.D0
          endif
        enddo
        lev=1
        f=alp(1)
        goto22
      endif
      p=abs(ls(pj))
      if(iprint.ge.2)write(nout,*)'CHOOSE p,pj =',p,pj
c  compute +/- Steepest Edge (SE) search direction s in an(.)
      call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1),
     *  e(p),.true.)

        rp=scpr(0.D0,ws(na1),g,n)
        if(ls(pj).lt.0)rp=-rp
        if(rp*r(p).le.0.D0)then
          print 2,'4rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
          do j=n1,lp(2)
            i=abs(ls(j))
            if(r(i).lt.0.D0)then
              r(i)=-eps
            else
              r(i)=0.D0
            endif
          enddo
          f=alp(1)
          goto98
        endif

      if(r(p).gt.0.D0)then
        rp=1.D0-r(p)
        plus=ls(pj).lt.0
      else
        rp=r(p)
        plus=ls(pj).ge.0
      endif
      snorm=e(p)
c     print 4,'s (or -s if .not.plus) =',(ws(i),i=na1,na+n)

c  form At.s and denominators
      call form_Ats(n1,lp(lev),n,plus,a,la,ws(na1),w,ls,snorm*tol)
c     print 2,'slope for r(169) =',aiscpr(n,a,la,169-n,ws(na1),0.D0)
c     print 2,'slope for r(217) =',aiscpr(n,a,la,217-n,ws(na1),0.D0)
c     print *,'plus =',plus
   86 continue
      if(iprint.ge.3)then
        write(nout,1001)'residual vector and indices',
     *    (ls(j),r(abs(ls(j))),j=n1,lp(lev))
        write(nout,1000)'denominators',(w(abs(ls(j))),j=n1,lp(lev))
      endif
   88 continue
c  ratio test at higher levels
      alpha=ainfty
      qj=0
      do 90 j=n1,lp(lev)
        i=abs(ls(j))
        wi=w(i)
        if(wi.eq.0.D0)goto90
        ri=r(i)
        if(ri.eq.-eps)then
          if(bl(i).eq.bu(i).and.wi.lt.0.D0)then
            alpha=0.D0
            q=i
            qj=j
            goto91
          endif
          goto90
        elseif(ri.lt.0.D0)then
          if(wi.gt.0.D0)goto90
          z=(ri-tol)/wi
        elseif(wi.gt.0.D0)then
          z=(ri+tol)/wi
        else
          goto90
c         if(bl(i).lt.bu(i))goto90
c         z=(ri-2.D0-tol)/wi
        endif
        if(z.ge.alpha)goto90
        alpha=z
        qj=j
   90 continue
      if(qj.eq.0)then
        do j=n1,lp(lev)
          i=abs(ls(j))
          if(r(i).lt.0.D0)then
            r(i)=-eps
          else
            r(i)=0.D0
          endif
        enddo
        call form_Ats(lp(lev)+1,lp(lev-1),n,plus,a,la,ws(na1),
     *    w,ls,snorm*tol)
        lev=lev-1
        f=alp(lev)
        if(iprint.ge.2)write(nout,*)'UNBOUNDED:   p =',p,
     *    '   return to level',lev
        if(lev.gt.1)goto86
        if(iprint.ge.3)then
          write(nout,1001)'costs vector and indices',
     *      (ls(j),r(abs(ls(j))),j=1,n)
        write(nout,1)'# of bound equations and free variables = ',peq,k
        endif
        goto30
      endif
      q=abs(ls(qj))
      alpha=r(q)/w(q)
c     if(alpha.lt.0.D0)then
c       r(q)=r(q)-2.D0
c       w(q)=-w(q)
c       ls(qj)=-ls(qj)
c     endif
c     print *,'alpha =',alpha
   91 continue
      if(iprint.ge.2)then
        write(nout,*)'alpha =',alpha,'   p =',p,'   q =',q
        write(nout,2)'r(p),r(q),w(q) =',r(p),r(q),w(q)
      endif

      if(alpha.eq.0.D0)then
        if(iprint.ge.2)write(nout,1)
     *    'degeneracy block at level',lev
        if(lev+2.gt.mlp)then
          ifail=5
          return
        endif
        if(w(q).lt.0.D0.and.r(q).ne.-eps)then
          w(q)=-w(q)
          ls(qj)=-ls(qj)
        endif
c       r(q)=0.D0
c       alp(lev)=f
        plev=n
        do j=n1,lp(lev)
          i=abs(ls(j))
          if(r(i).eq.-eps)then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            r(i)=-1.D0
          elseif(r(i).eq.0.D0)then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            r(i)=1.D0
          endif
        enddo
        lev=lev+1
        lp(lev)=plev
        if(iprint.ge.2)write(nout,*)
     *    'degeneracy: increase level to ',lev       
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,A,''level ='',I2,
     *    ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
          goto86
      endif
c  update r and f
      if(alpha.gt.0.D0)then
c       ff=f
c       f=f+alpha*rp
c       if(f.ge.ff)then
c         if(r(p).gt.0.D0)then
c           r(p)=1.D0
c         else
c           r(p)=0.D0
c         endif
c         goto85
c       endif
        do 92 j=n1,lp(lev)
          i=abs(ls(j))
          if(w(i).eq.0.D0)goto92
          ri=r(i)-alpha*w(i)
          if(abs(ri).le.tol)ri=0.D0
          if(r(i).lt.0.D0.and.ri.ge.0.D0)then
c  remove contribution to gradient
            call saipy(sign(1.D0,dble(ls(j))),a,la,i-n,g,n)
            call newg
          endif
c         if(w(i).lt.0.D0.and.bl(i).eq.bu(i))then
c           ro=2.D0-ri
c           if(abs(ro).le.tol)ro=0.D0
c           if(ro.lt.ri)then
c             ri=ro
c             w(i)=-w(i)
c             ls(j)=-ls(j)
c           endif
c         endif
          r(i)=ri
   92   continue
      endif
      if(iprint.ge.2)write(nout,*)'replace',p,' by',q
      call pivot(p,q,n,nm,a,la,e,ws(lu1),lws(ll1),ifail,npv)
      if(ifail.ge.1)then
        if(ifail.ge.2)return
c       call iexch(ls(pj),ls(qj))
        if(iprint.ge.1)write(nout,*)'near singularity in pivot (4)'
        goto98
      endif
      if(r(p).gt.0.D0)then
c  add contribution to g from c/s p
        call saipy(-sign(1.D0,dble(ls(pj))),a,la,p-n,g,n)
        call newg
        r(p)=-alpha
      else
        rpu=2.D0-alpha
        if(alpha.le.rpu.or.bl(p).lt.bu(p))then
          r(p)=alpha
        else
          r(p)=rpu
          ls(pj)=-ls(pj)
        endif
      endif
      if(abs(r(p)).le.tol)r(p)=0.D0
c  exchange a constraint
      call iexch(ls(pj),ls(qj))
      if(q.le.n.and.bl(q).eq.bu(q))then
        peq=peq+1
        call iexch(ls(pj),ls(peq))
      endif
      if(iprint.ge.1)write(nout,'(''pivots ='',I5,A,''level ='',I2,
     *  ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
      goto80
c  restart sequence
   98 continue
      return
 1000 format(a/(e16.5,4e16.5))
 1001 format(a/(i4,1x,e11.5,4(i4,1x,e11.5)))
c1000 format(a/(e18.8,3e19.8))
c1001 format(a/(i3,1x,e14.8,3(i4,1x,e14.8)))
      end

      subroutine optest1(peq,k,n,bl,bu,r,e,ls,rp,pj)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension bl(*),bu(*),r(*),e(*),ls(*)
      rp=0.D0
      do 1 j=peq+1,n-k
        i=abs(ls(j))
        ri=-r(i)/e(i)
        if(i.le.n)then
          if(ri.le.rp)goto1
        elseif(ri.lt.0.D0)then
          ri=(r(i)-1.D0)/e(i)
          if(ri.le.rp)goto1
        elseif(ri.gt.rp)then
          if(bl(i).eq.bu(i))then
            ri=(-r(i)-1.D0)/e(i)
            if(ri.le.rp)goto1
            r(i)=-r(i)
            ls(j)=-ls(j)
          endif
        else
          goto1
        endif
        rp=ri
        pj=j
    1 continue
c  additional test for free variables
      do j=n-k+1,n
        i=abs(ls(j))
        ri=abs(r(i))/e(i)
        if(ri.ge.rp)then
          rp=ri
          pj=j
        endif
      enddo
      return
      end

      subroutine check1(n,nm,nmi,k,kmax,g,a,la,x,bl,bu,r,ls,lp,an,f,
     *  ws,lws,cws,lev,p,alp2)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension g(*),a(*),la(*),x(*),bl(*),bu(*),r(*),ls(*),lp(*),
     *  an(*),ws(*),lws(*)
      character cws(*)
      common/noutc/nout
      common/epsc/eps,tol,emin
c     print *,'ENTER check1'
      e=0.D0
      if(lev.eq.1)then
        do j=n+1,nm
          i=abs(ls(j))
          if(i.le.n)then
            s=x(i)
          else
            s=aiscpr(n,a,la,i-n,x,0.D0)
          endif
          if(ls(j).gt.0)then
c           print *,'i,s,r(i),bl(i)',i,s,r(i),bl(i)
            s=r(i)-s+bl(i)
          else
            s=r(i)+s-bu(i)
          endif
          if(abs(s).le.tol*max(1.D0,abs(r(i))))s=0.D0
          if(abs(s).gt.e)then
            e=abs(s)
            ie=i
          endif
        enddo
      else
        do j=n+1,lp(2)
          i=abs(ls(j))
          if(i.le.n)then
            s=x(i)
          else
            s=aiscpr(n,a,la,i-n,x,0.D0)
          endif
c         print 2,'s,bl(i),bu(i) =',s,bl(i),bu(i)
          if(ls(j).gt.0)then
            s=-s+bl(i)
          else
            s=s-bu(i)
          endif
          if(abs(s).le.tol)s=0.D0
          if(abs(s).gt.e)then
            e=abs(s)
            ie=i
          endif
        enddo
      endif
      if(e.gt.tol)write(nout,*)'inactive c/s residual error = ',e,ie
c     if(e.gt.tol)stop
      if(e.gt.1.D-6)print 2,'r(ie)',r(ie)
      if(e.gt.1.D-6)stop
      if(lev.eq.1)then
        ff=0.D0
        do j=n+1,nm
          i=abs(ls(j))
          if(r(i).lt.0.D0)ff=ff-r(i)
        enddo
        e=abs(ff-f)
        if(e.gt.tol*max(1.D0,abs(f)))write(nout,*)'function error = ',e,
     *    '   f(x) =',ff
        if(e.gt.tol*max(1.D0,abs(f)))stop
      endif  
c       print 4,'g =',(g(i),i=1,n)
c       print 4,'an =',(an(i),i=1,n)
c       err=0.D0
c       do i=1,n
c         err=err+abs(g(i)-an(i))
c       enddo
c       print 2,'check err =',err
      do i=1,n
        an(i)=0.D0
      enddo
      do j=n+1,nm
        i=abs(ls(j))
        if(r(i).lt.0.D0)then
          if(i.gt.n)then
            call saipy(-sign(1.D0,dble(ls(j))),a,la,i-n,an,n)
          else
            an(i)=an(i)-sign(1.D0,dble(ls(j)))
          endif
        endif
      enddo
      gnm=sqrt(scpr(0.D0,an,an,n))
      e=0.D0
      do j=1,n
c       write(nout,*)'an =',(an(i),i=1,n)
        i=abs(ls(j))
        s=sign(1.D0,dble(ls(j)))
        if(i.le.n)then
          an(i)=an(i)-s*r(i)
          if(j.gt.n-k)then
            s=max(0.D0,bl(i)-x(i),x(i)-bu(i))
          elseif(ls(j).gt.0)then
            s=x(i)-bl(i)
          else
            s=bu(i)-x(i)
          endif
        else
          call saipy(-s*r(i),a,la,i-n,an,n)
          if(ls(j).gt.0)then
            s=aiscpr(n,a,la,i-n,x,-bl(i))
          else
            s=-aiscpr(n,a,la,i-n,x,-bu(i))
          endif
        endif
        if(abs(s).gt.e)then
          e=abs(s)
          ie=i
        endif
      enddo
      if(e.gt.tol)write(nout,*)'active c/s residual error = ',e,ie
c     if(e.gt.tol)stop
c     if(e.gt.1.D-4)print 4,'x =',(x(i),i=1,n)
      if(e.gt.1.D-4)stop
      e=0.D0
      do j=1,n
        if(abs(an(j)).gt.e)then
          e=abs(an(j))
          ie=ls(j)
          je=j
        endif
      enddo
      if(e.gt.gnm*tol)write(nout,*)'KT condition error = ',e,je,ie,gnm
c     if(e.gt.gnm*tol)write(nout,4)'KT cond_n errors = ',(an(i),i=1,n)
c     if(e.gt.gnm*tol)stop
      if(e.gt.1.D-4)stop
    1 format(A,10I5)
    2 format(A,5E15.7)
    4 format(A/(5E15.6))
      return
      end
