christen this file glcpd.f
cut here >>>>>>>>>>>>>>>>>

c  Copyright (C) 2010 Roger Fletcher

c  Current version dated 27 March 2013

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

      subroutine glcpd(n,m,k,kmax,maxg,a,la,x,bl,bu,f,fmin,g,r,w,e,ls,
     * alp,lp,mlp,peq,ws,lws,cws,v,nv,rgtol,m0de,ifail,mxgr,iprint,nout)
      implicit double precision (a-h,r-z), integer (i-q)

c  This routine finds a KT point for the General LCP (Linearly Constrained
c  Problem)

c       minimize    f(x)

c       subject to  l <= [I : A]t.x <= u                  (t = transpose)

c  where f(x) is a given function of n variables x, to be determined.
c  Lower and upper bound constraints on the variables x and the linear
c  functions At.x may be supplied, where A is an n*m matrix.
c  A recursive form of an active set method is used, using Wolfe's method to
c  resolve degeneracy. A limited memory reduced gradient sweep method is used
c  for minimization in the null space, so usually the KT point is a local
c  minimizer. Matrix information is made available and processed by calls to
c  external subroutines. Details of these are given in an auxiliary file
c  named either 'denseL.f' or 'schurQR.f'. (schurQR.f is a more recent
c  replacement for the file sparseL.f)

c  parameter list  (variables in a line starting with C must be set on entry)
c  **************

C  n     number of variables
C  m     number of general constraints (columns of A)
c  k     dimension of the null space obtained by eliminating the active
c        constraints (only to be set if mode>=2). The number of constraints in
c        the active set is n-k
C  kmax  maximum value of k (kmax <= n)
C  maxg  max number of reduced gradient vectors stored in sweep method:
C        (1 < maxg <= kmax+1 when kmax>0), typically maxg = min(6,kmax+1)
C  a(*)  storage of reals associated with A. This storage may be provided
C        in either dense or sparse format. Refer to either denseA.f or sparseA.f
C        for information on how to set a(*) and la(*). The vector c referred to
C        in these files should be set to the zero vector.
C  la(*) storage of integers associated with c and A
C  x(n)  contains the vector of variables. Initially an estimate of the solution
C        must be set, replaced by the solution (if it exists) on exit.
C  bl(n+m)  vector of lower bounds for variables and general constraints
C  bu(n+m)  vector of upper bounds (use numbers less than about 1.e30, and
C        where possible supply realistic bounds on the x variables)
c  f     returns the value of f(x) when x is a feasible solution
c        Otherwise f stores the sum of constraint infeasibilities
C  fmin  set a strict lower bound on f(x) (used to identify an unbounded LCP)
c  g(n)  returns the gradient vector of f(x) when x is feasible
c  r(n+m) workspace: stores constraint residuals (or multipliers if the
c        constraint is active). The sign convention is such that these are
c        nonnegative at a solution (except multipliers of equality constraints)
c  w(n+m) workspace: stores denominators for ratio tests
c  e(n+m) stores steepest-edge normalization coefficients: if mode>2 then
c        information in this vector from a previous call should not be changed.
c        (In mode 3 these values provide approximate coefficients)
c  ls(n+m) stores indices of the active constraints in locations 1:n and of
c        the inactive constraints in locations n+1:n+m. The simple bounds
c        on the variables are indexed by 1:n and the general constraints by
c        n+1:n+m. The sign of ls(j) indicates whether the lower bound (+) or
c        the upper bound (-) of constraint ls(j) is currently significant.
c        Within the set of active constraints, locations 1:peq store the indices
c        of any equality constraints, locations peq+1:n-k store the indices of
c        any inequality constraints, and locations n-k+1:n store the indices of 
c        any free variables (variables not on a bound, which are used to
c        parametrise the null space: ls(j) is always positive in this range)
c          If mode>=2, the first n-k elements of ls must be set on entry
c  alp(mlp) workspace associated with recursion
c  lp(mlp)  list of pointers to recursion information in ls
C  mlp   maximum number of levels of recursion allowed (mlp>2: typically
C        mlp=50 would usually be adequate but mlp=m is an upper bound)
c  peq   pointer to the end of equality constraint indices in ls
c  ws(*) real workspace for gdotx (see below), qlcpd and denseL.f (or schurQR.f)
c          Set the total number in mxws (see "Common" below).
c  lws(*) integer workspace for gdotx, qlcpd and denseL.f (or schurQR.f). 
c          Set the total number in mxlws (see "Common" below).
c        The storage maps for ws and lws are set by the routine stmap below
c  cws(*) character workspace (if any) needed by funct
C  v(maxg) set nv estimates of the eigenvalues of the reduced Hessian of f(x)
C          (for example from a previous run of glcpd). Set nv=1 and v(1)=1.D0
C          in absence of other information. New values of v are left on exit
C  nv    Number of estimates in v
C  rgtol required accuracy in the reduced gradient l2 norm: it is advisable not 
C        to seek too high accuracy - rgtol may be increased by the code if it
c        is deemed to be too small, see the definition of sgnf below
C  m0de  mode of operation (larger numbers imply extra information):
C          0 = cold start (no other information available, takes simple
C                bounds for the initial active set)
C          1 = as 0 but includes all equality constraints in initial active set
C          2 = user sets n-k active constraint indices in ls(j), j=1,..,n-k.
c                For a general constraint the sign of ls(j) indicates which
c                bound to use. For a simple bound the current value of x is used
C          3 = takes active set and other information from a previous call.
C                Steepest edge weights are approximated using previous values.
C          4 = as 3 but it is also assumed that columns of A are unchanged
c                so that factors of the basis matrix stored in ws and lws are
c                valid (changes in f(x) and the vectors l and u are allowed)
C        A local copy (mode) of m0de is made and may be changed by glcpd
c  ifail   outcome of the process
c              0 = solution obtained
c              1 = unbounded problem (f(x)<fmin has occurred: note grad is not
c                    evaluated in this case)
c              2 = bl(i) > bu(i) for some i
c              3 = infeasible problem detected in Phase 1
c              4 = line search cannot improve f (possibly increase rgtol)
c              5 = mxgr gradient calls exceeded (this test is only carried
c                    out at the start of each iteration)
c              6 = incorrect setting of m, n, kmax, maxg, mlp, m0de or tol
c              7 = not enough space in ws or lws
c              8 = not enough space in lp (increase mlp)
c              9 = dimension of reduced space too large (increase kmax)
c             10 = maximum number of unsuccessful restarts taken
c            >10= possible use by later sparse matrix codes
C  mxgr  maximum number of gradient calls
C  iprint  switch for diagnostic printing (0 = off, 1 = summary,
C                 2 = scalar information, 3 = verbose)
C  nout  channel number for output

c  Storage Allocation
c  ******************
c  User information about the lengths of ws and lws is supplied to glcpd in
c    common/wsc/kk,ll,kkk,lll,mxws,mxlws
c  kk and ll refer to the lengths of ws and lws needed by the user subroutines.
c  kkk and lll are the numbers of locations used by glcpd and are set by glcpd.
c  The rest of ws and lws is used by the files denseL.f or schurQR.f
c  mxws and mxlws must be set to the total lengths of ws and lws available: a
c  message will be given if more storage is needed.

c  User subroutines
c  ****************

c  The user must provide two subroutines as follows

c      subroutine funct(n,x,f,ws,lws,cws)
c      implicit double precision (a-h,o-z)
c      dimension x(*),ws(*),lws(*)
c      character cws(*)
c      ...
c      statements to compute f(x) from x
c      ...
c      return
c      end

c      subroutine grad(n,x,g,ws,lws,cws)
c      implicit double precision (a-h,o-z)
c      dimension x(*),ws(*),lws(*)
c      character cws(*)
c      ...
c      statements to compute grad.f(x) in g from x (the user
c      may assume that a call of grad immediately follows one
c      of funct with the same vector x.)
c      ...
c      return
c      end

c  The parameters ws, lws and cws in the above subroutines enables data to be
c  passed from the user's calling program to these subroutines

c  Tolerances, accuracy and diagnostics
c  ************************************
c  glcpd uses tolerance and accuracy information stored in
c     common/epsc/eps,tol,emin
c     common/repc/sgnf,nrep,npiv,nres
c     common/refactorc/mc,mxmc
c     common/infoc/rgnorm,vstep,iter,npv,nfn,ngr
c  eps must be set to the machine precision (unit round-off) and tol is a
c  tolerance such that numbers whose absolute value is less than tol are
c  truncated to zero. This tolerance strategy in the code assumes that the
c  problem is well-scaled. The parameter sgnf is used to measure the maximum
c  allowable relative error in gradient values. If at any stage the accuracy
c  requirement rgtol < sgnf*rgnorm then rgtol is increased to sgnf*rgnorm
c    The code allows one or more refinement steps after the
c  calculation has terminated, to improve the accuracy of the solution,
c  and a fixed number nrep of such repeats is allowed. However the code
c  terminates without further repeats if no more than npiv pivots are taken.
c    In case of any breakdown, the code is restarted in mode 0.
c  The maximum number of unsuccessful restarts allowed is set in nres.
c    The basis matrix may be refactorised on occasions, for example to prevent
c  build-up of round-off in the factors or (when using schurQR.f) to limit
c  the growth in the Schur complement. The maximum interval between
c  refactorizations (or size of Schur complement) is set in mxmc.
c    Default values are set in block data but can be reset by the user.
c    infoc returns information about the progress of the method: rgnorm is the
c  norm of the reduced gradient on exit, and vstep is the length of the vertical
c  step in the warm start process. iter is the total number of iterations taken,
c  npv is the number of pivots, nfn is the number of function evaluations, and
c  ngr is the number of gradient evaluations.

      parameter (ainfty=1.D100)
      dimension a(*),la(*),x(*),bl(*),bu(*),g(*),r(*),w(*),e(*),ls(*),
     *  alp(*),lp(*),ws(*),lws(*),v(*)
      character cws(*)
      character*32 spaces
      common/lcpdc/na,na1,nb,nb1,krg,krg1,kr,kr1,
     *  ka,ka1,kb,kb1,kc,kc1,kd,kd1,ke,ke1,lu1,ll1
      common/epsc/eps,t0l,emin
c     common/epsc/eps,tol,emin
      common/infoc/rgnorm,vstep,iter,npv,nfn,ngr
      common/repc/sgnf,nrep,npiv,nres
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/refactorc/mc,mxmc
      common/alphac/alpha,rp,pj,qqj,qqj1
      logical plus

    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(15I5))
    4 format(A/(5E15.7))
    5 format((6E15.7))
    6 format(A,I5,2E15.7)

      spaces='         '
      mode=m0de
      tol=t0l
      iter=0
      npv=0
      if(m.lt.0.or.n.le.0.or.mlp.lt.2.or.mode.lt.0.or.mode.gt.4.or.
     *  kmax.lt.0.or.(kmax.gt.0.and.maxg.le.1).or.tol.le.0.D0)then
        ifail=6
        return
      endif
      rgt0l=rgtol
      n1=n+1
      nm=n+m
      nmi=nm
      nfn=0
      ngr=0
      nv0=nv
      if(iprint.ge.3)then
        write(nout,1000)'lower bounds',(bl(i),i=1,nm)
        write(nout,1000)'upper bounds',(bu(i),i=1,nm)
      endif
      irep=0
      ires=0
      do i=1,nm
        t=bu(i)-bl(i)
        if(t.lt.-tol)then
          print *,'i,bl(i),bu(i)',i,bl(i),bu(i)
          ifail=2
          return
        elseif(t.le.tol)then
          bl(i)=5.D-1*(bl(i)+bu(i))
          bu(i)=bl(i)
        endif
      enddo
      vmax=0.D0
      do i=1,n
        x(i)=min(bu(i),max(bl(i),x(i)))
        vmax=max(vmax,bu(i)-bl(i))
      enddo
      if(mode.le.2)then
        call stmap(n,nm,kmax,maxg)
        if(mode.eq.0)then
          nk=0
        elseif(mode.eq.1)then
c  collect equality c/s
          nk=0
          do i=1,nm
            if(bu(i).eq.bl(i))then
              nk=nk+1
              ls(nk)=i
            endif
          enddo
c         write(nout,*)'number of eqty c/s =',nk
        else
          nk=n-k
        endif
      endif
c  restarts loop
    7 continue
      lp(1)=nm
      lev=1
      if(mode.le.3)then
c  set up factors of basis matrix and permutation vectors
        ifail=mode
        call start_up(n,nm,nmi,a,la,nk,e,ls,ws(lu1),lws(ll1),mode,ifail)
        if(ifail.gt.0)return
      endif
    8 continue
      peq=0
      ig=0
c  refinement step loop
      mpiv=iter+npiv
      ninf=0
      do i=1,n
        g(i)=0.D0
      enddo
      if(mode.gt.0)then
        call warm_start(n,nm,a,la,x,bl,bu,r,ls,ws(lu1),
     *    lws(ll1),ws(na1),vstep)
c       print *,'vstep,vmax',vstep,vmax
        if(vstep.gt.2.D0*vmax)then
          mpiv=0
          mode=0
          nk=0
          do i=1,n
            x(i)=min(bu(i),max(bl(i),x(i)))
          enddo
          goto7
        endif
        if(vstep.gt.tol)mpiv=0
      endif
      k=0
c  collect free variables
      do j=n,1,-1
        i=abs(ls(j))
        if(i.le.n.and.x(i).gt.bl(i).and.x(i).lt.bu(i))then
          call iexch(ls(j),ls(n-k))
          k=k+1
        endif
      enddo
      if(mode.eq.0)then
        do j=1,n-k
          i=ls(j)
          if(x(i).eq.bu(i))ls(j)=-i
        enddo
        lp(1)=n
        goto9
      endif
      phase=0
c  move inactive general c/s to the end
      do j=nm,n1,-1
        i=abs(ls(j))
        if(i.gt.n)then
          call iexch(ls(j),ls(lp(1)))
          lp(1)=lp(1)-1
        endif
      enddo
      call residuals(n,n1,lp(1),a,la,x,bl,bu,r,ls,f,g,ninf)
      if(ninf.gt.0)then
        gnorm=sqrt(dble(ninf))
        gtol=sgnf*gnorm
        rgtol=max(rgt0l,gtol)
        goto15
      endif
    9 continue
c  enter phase 1
      phase=1
c  collect active equality c/s
      do j=1,n-k
        i=abs(ls(j))
        if(bu(i).eq.bl(i))then
          peq=peq+1
          call iexch(ls(j),ls(peq))
        endif
      enddo
      call residuals(n,lp(1)+1,nm,a,la,x,bl,bu,r,ls,f,g,ninf)
      lp(1)=nm
      if(ninf.gt.0)then
        gnorm=sqrt(scpr(0.D0,g,g,n))
        gtol=sgnf*gnorm
        rgtol=max(rgt0l,gtol)
        goto15
      endif
   10 continue
      phase=2
      if(iprint.ge.1)write(nout,*)'FEASIBILITY OBTAINED at level 1'
      n_inf=0
      call funct(n,x,f,ws,lws,cws)
      nfn=nfn+1
      if(f.lt.fmin)goto75
      call grad(n,x,g,ws,lws,cws)
      ngr=ngr+1
c     write(nout,4)'x =',(x(i),i=1,n)
c     write(nout,4)'g =',(g(i),i=1,n)
      call newg
      gnorm=sqrt(scpr(0.D0,g,g,n))
      gtol=sgnf*gnorm
      rgtol=max(rgt0l,gtol)
      alpha=1.D0
      ig=0
      if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *  ''  level = 1    f ='',E16.8)')npv,f
      goto16
c  start of major iteration
   15 continue
      if(iprint.ge.1)then
        if(ninf.eq.0)then
          if(k.gt.0)then
c           write(nout,'(''pivots ='',I5,
c    *        ''  level = 1    f ='',E16.8,''   k ='',I4)')npv,f,k
            write(nout,'(''pivots ='',I5,
     *        ''  level = 1    f ='',E16.8,''   rg ='',E12.4,
     *        ''  k ='',I4)')npv,f,rgnorm,k
          else
            write(nout,'(''pivots ='',I5,
     *        ''  level = 1    f ='',E16.8)')npv,f
          endif
        elseif(phase.eq.0)then
          write(nout,'(''pivots ='',I5,''  level = 1    f ='',
     *      E16.8,''   ninfb ='',I4)')npv,f,ninf
        else
          write(nout,'(''pivots ='',I5,''  level = 1    f ='',
     *      E16.8,''   ninf ='',I4)')npv,f,ninf
        endif
      endif
   16 continue
c  calculate multipliers
c     print 4,'gradient =',(g(i),i=1,n)
      do i=1,nm
        w(i)=0.D0
      enddo
      call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(n,r,w,ls)
c  opposite bound or reset multiplier loop
   20 continue
      if(iprint.ge.3)then
        write(nout,1001)'costs vector and indices',
     *    (ls(j),r(abs(ls(j))),j=1,n)
c       write(nout,1000)'steepest edge coefficients',
c    *    (e(abs(ls(j))),j=1,n)
        if(peq.gt.0.or.k.gt.0)write(nout,1)
     *    '# active equality c/s and free variables = ',peq,k
      endif
c     call check(n,lp(1),nmi,kmax,g,a,la,x,bl,bu,r,ls,ws(nb1),f,
c    *  ws,lws,cws,ninf,peq,k,1,p,rp)

   21 continue
      call optest(peq+1,n-k,r,e,ls,rp,pj)
      if(phase.eq.0)then
c  possibly choose an active general c/s to relax (marked by rp>0)
        t=-1.D1*rp
        do 13 j=1,n
          i=abs(ls(j))
          if(i.le.n)goto13
          if(bu(i).eq.bl(i).and.r(i).lt.0.D0)then
            r(i)=-r(i)
            ls(j)=-ls(j)
          endif
          if(r(i)/e(i).le.t)goto13
          rp=r(i)
          t=rp/e(i)
          pj=j
   13   continue
      endif

      if(ig.eq.0)then
        gg=0.D0
        do j=n-k+1,n
          i=ls(j)
          gg=gg+r(i)**2
        enddo
        rgnorm=sqrt(gg)
      endif
c     print 2,'rgtol,rgnorm,rp',rgtol,rgnorm,rp

   25 continue
      if(rgnorm.le.rgtol.and.abs(rp).le.gtol)then
c  allow for changes to norm(g)
        gnorm=sqrt(scpr(0.D0,g,g,n))
        gtol=sgnf*gnorm
        rgtol=max(rgt0l,gtol)
      endif
      if(iprint.eq.3)print 2,'gtol,rgtol,rgnorm,rp',gtol,rgtol,rgnorm,rp

      if((rgnorm.le.rgtol.and.abs(rp).le.gtol).or.ngr.gt.mxgr)then
c  optimal at current level: first tidy up x
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
          if(r(i).le.tol.and.i.le.n)then
            r(i)=0.D0
            if(ls(j).ge.0)then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            endif
          endif
        enddo
        if(ngr.gt.mxgr)then
          ifail=5
          return
        endif
        if(iprint.ge.2)then
          write(nout,*)'OPTIMAL at level 1'
          if(iprint.ge.3)then
c           write(nout,1000)'x variables',(x(i),i=1,n)
            write(nout,1001)'residual vector and indices',
     *        (ls(j),r(abs(ls(j))),j=n1,nm)
          endif
        endif
        irep=irep+1
        if(irep.le.nrep.and.iter.gt.mpiv)then
          if(iprint.ge.1)write(nout,*)'refinement step #',irep
          mode=4
          goto8
        endif
        if(iprint.ge.2.and.nrep.gt.0)
     *     write(nout,*)'total number of restarts =',ires
        if(ninf.gt.0)then
          ifail=3
          return
        endif
        nv=nv0
        ifail=0
        return
      endif

      if(rgnorm.ge.abs(rp))then
c  ignore the multiplier of c/s p and set up or continue SD steps
        p=0
      else
        p=abs(ls(pj))
        if(iprint.ge.2)print 1,'CHOOSE p =',p
        rp=r(p)
        call iexch(ls(pj),ls(n-k))
        pj=n-k
        ig=0
      endif

c     if(k.eq.0.or.p.gt.n)then
      if(p.gt.0)then
c  compute +/- Steepest Edge (SE) search direction s in an(.)
        call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1),
     *    e(p),.true.)
        rp=scpr(0.D0,ws(na1),g,n)
        if(ls(pj).lt.0)rp=-rp
        if(rp*r(p).le.0.D0)then
          r(p)=0.D0
          goto21
        endif
        if(abs(rp-r(p)).gt.5.D-1*max(abs(rp),abs(r(p))))then
c       if(abs(rp-r(p)).gt.1.D-1*gnorm)then
          print 2,'1rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
          goto98
        endif
        snorm=e(p)
        plus=ls(pj).ge.0.eqv.rp.lt.0.D0
        f0=f
        ig=0
      else
        if(ig.eq.0)then
c  start up the limited memory sweep method
c         if(p.gt.0)then
c  transfer c/s p into Z
c           if(ls(pj).lt.0)then
c             r(p)=-r(p)
c             ls(pj)=-ls(pj)
c           endif
c           k=k+1
c           gg=gg+r(p)**2
c         endif
          ig=1
          ngv=1
          f0=f
          ws(kb1)=gg
          rgnorm=sqrt(gg)
c         print 2,'initial rg =',(r(ls(j)),j=n-k+1,n)
          if(k*ngv.gt.kmax*maxg)then
            ifail=9
            return
          endif
          call store_rg(k,ig,ws(krg1),r,ls(n-k+1))
        endif
c  compute Steepest Descent (SD) search direction s = -Z.rg in an(.)
        call zprod(k,n,a,la,ws(na1),r,w,ls,ws(lu1),lws(ll1))
        rp=scpr(0.D0,ws(na1),g,n)
        if(abs(gg+rp).gt.5.D-1*max(gg,abs(rp)))then
c       if(abs(gg+rp).gt.1.D-2*max(gg,abs(rp)))then
          print 2,'gg,rp,gg+rp',gg,rp,gg+rp
          goto98
        endif
        snorm=sqrt(scpr(0.D0,ws(na1),ws(na1),n))
        plus=.true.
      endif
c     print 4,'s (or -s if .not.plus) =',(ws(i),i=na1,na+n)

c  form At.s and denominators
      call form_Ats(n1,lp(1),n,plus,a,la,ws(na1),w,ls,snorm*tol)

c  return from degeneracy code
   30 continue

      if(iprint.ge.3)then
        write(nout,1000)'x variables',(x(i),i=1,n)
        write(nout,1001)'residual vector and indices',
     *    (ls(j),r(abs(ls(j))),j=n1,lp(1))
        write(nout,1000)'denominators',(w(abs(ls(j))),j=n1,lp(1))
      endif

   40 continue
c  level 1 ratio tests
      amax=ainfty
      qj=0
      qj1=0
      do 41 j=n-k+1,n
        i=ls(j)
        if(i.le.0)print *,'i.le.0'
        if(i.le.0)goto98
        si=ws(na+i)
        if(si.eq.0.D0)goto41
        t=abs(si)
c       if(t.le.tol)goto41
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
      if(ig.eq.0.and.rp.lt.0.D0.and.bu(p)-bl(p).lt.amax)then
        amax=bu(p)-bl(p)
        qj=pj
      endif
      if(ninf.gt.0)then
        alpha1=ainfty
        do 42 j=n1,lp(1)
          i=abs(ls(j))
          wi=w(i)
          if(wi.eq.0.D0)goto42
          ri=r(i)
          if(wi.gt.0.D0)then
            if(ri.lt.0.D0)goto42
            z=(ri+tol)/wi
          else
            if(ri.lt.0.D0)then
              z=ri/wi
              if(z.lt.alpha1)then
                alpha1=z
                qj1=j
              endif
            endif
            z=(bl(i)-bu(i)+ri-tol)/wi
          endif
          if(z.ge.amax)goto42
          amax=z
          qj=j
   42   continue
        if(qj1.gt.0.and.alpha1.le.amax)then
c  find feasible step that zeros most infeasible c/s
          do 43 j=n1,lp(1)
            i=abs(ls(j))
            wi=w(i)
            if(wi.ge.0.D0)goto43
            ri=r(i)
            if(ri.lt.0.D0)then
              z=ri/wi
              if(z.gt.alpha1.and.z.le.amax)then
                alpha1=z
                qj1=j
              endif
            endif
   43     continue
          amax=alpha1
          qj=qj1
        else
          qj1=0
        endif
      else
        do 44 j=n1,lp(1)
          i=abs(ls(j))
          wi=w(i)
          if(wi.eq.0.D0)goto44
          ri=r(i)
          if(wi.gt.0.D0)then
            z=(ri+tol)/wi
          else
            z=(bl(i)-bu(i)+ri-tol)/wi
          endif
          if(z.ge.amax)goto44
          amax=z
          qj=j
   44   continue
      endif
      q=abs(ls(qj))
      if(iprint.ge.2.and.q.ne.p.and.qj.gt.n)
     *  write(nout,*)'q,r(q),w(q) =',q,r(q),w(q)
      if(qj.gt.n.and.qj1.eq.0)then
        if(w(q).gt.0.D0)then
          amax=r(q)/w(q)
        else
          amax=(bl(q)-bu(q)+r(q))/w(q)
        endif
      endif

      if(amax.eq.0.D0.and.rp.le.0.D0)then
        alpha=0.D0
c  potential degeneracy block at level 1
        if(p.eq.0)goto65
        if(bu(q).eq.bl(q))goto70
        plev=n
        do j=n1,lp(1)
          i=abs(ls(j))
          if(r(i).eq.0.D0)then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            if(bu(i).gt.bl(i))r(i)=1.D0
          endif
        enddo
        if(plev.gt.n1)then
          lp(2)=plev
          lev=2
          alp(1)=f
          f=0.D0
          qj=pj
          q=p
          if(iprint.ge.1)write(nout,'(''pivots ='',I5,''     level = 2'',
     *      ''    f ='',E16.8)')npv,f
          goto86
        endif
        qj=n1
        r(q)=0.D0
c       print *,'only one degenerate c/s'
        goto70
      endif

      if(ninf.gt.0)then
        alpha=amax
        if(plus)then
          call mysaxpy(alpha,ws(na1),x,n)
        else
          call mysaxpy(-alpha,ws(na1),x,n)
        endif
      else
c  take a Ritz value off the stack
c       print 4,'Ritz values =',(v(i),i=1,nv)
        if(nv.gt.0.and.v(nv).gt.0.D0)then
          alpha=min(1.D0/v(nv),amax)
          nv=nv-1
        else
          alpha=amax
          nv=0
        endif
c  line search
        alphar=amax
        alphal=0.D0
        dalpha=alpha
        fi=f
        fr=ainfty
        ggo=gg
        gs=rp
        gsi=gs
c       print 2,'f0,fi,gsi,amax =',f0,fi,gsi,amax
   51   continue
c  calculate new x
        if(plus)then
          call saxpyz(alpha,ws(na1),x,ws(nb1),n)
        else
          call saxpyz(-alpha,ws(na1),x,ws(nb1),n)
        endif
        call funct(n,ws(nb1),fp,ws,lws,cws)
        if(fp.lt.fmin)goto75
        nfn=nfn+1
        df=f-fp
c  check for lack of improvement
        if(fp.ge.f0)then
c         print 2,'alphal,alpha,fp =',alphal,alpha,fp
          if(dalpha.lt.1.D-10.and.df.lt.-dalpha*gs)then
c           print *,'alpha too small'
            if(alphal.gt.0.D0)goto52
            ifail=4
            return
          endif
          fr=fp
          alphar=alpha
          z=5.D-1/(1.D0+df/(gs*dalpha))
c         print 2,'df,z =',df,z
          dalpha=dalpha*max(1.D-1,z)
          alpha=alphal+dalpha
          nv=0
          goto51
        endif
        f=fp
        call grad(n,ws(nb1),g,ws,lws,cws)
        ngr=ngr+1
c       print 4,'new g =',(g(i),i=1,n)
        call newg
        gps=scpr(0.D0,g,ws(na1),n)
        if(.not.plus)gps=-gps
c       print 2,'fp,gps',fp,gps
c       print 2,'alphal,alpha,alphar',alphal,alpha,alphar
c  check for non-positive curvature
        if(alpha.lt.amax.and.(gps.le.gsi.or.(gps.lt.25.D-2*gsi.and.
     *    (alphal.gt.0.D0.or.fr.lt.ainfty))))then
c       if(alpha.lt.amax.and.gps.le.gsi)then
          alphal=alpha
          if(fr.eq.ainfty)then
            alpha=min(alpha*5.D0,amax)
            dalpha=alpha-alphal
          else
            dalpha=alphar-alpha
            z=max(2.D-1,5.D-1/(1.D0+(f-fr)/(gps*dalpha)))
            dalpha=dalpha*z
            alpha=min(alpha+dalpha,amax)
          endif
          gs=gps
          nv=0
          goto51
        endif
c  end of line search
   52   continue
        do i=1,n
          x(i)=ws(nb+i)
        enddo
        if(ig.eq.0)goto60
        ig1=ig+1
        if(ig1.gt.maxg)ig1=1
        call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
c       print 4,'new rg =',(w(ls(j)),j=n-k+1,n)
        if(ngv.lt.maxg)ngv=ngv+1
        if(k*ngv.gt.kmax*maxg)then
          ifail=9
          return
        endif
        call store_rg(k,ig1,ws(krg1),w,ls(n-k+1))
        gpg=0.D0
        gg=0.D0
        do j=n-k+1,n
          i=ls(j)
          gpg=gpg+r(i)*w(i)
          gg=gg+w(i)**2
        enddo
        rgnorm=sqrt(gg)
c       print 2,'gpg,gg',gpg,gg
c       print 2,'f =',f
        call signst(n,r,w,ls)
        ws(ka+ig)=1.D0/alpha
        ws(kb+ig1)=gg
        ws(kc+ig)=gpg
        if(nv.eq.0.or.gg.gt.ggo)then
c  compute new Ritz values
          if(ngv.eq.0)then
            nv=1
            v(1)=1.D0/alpha
          else
            nv=min(ngv-1,k)
            if(nv.le.0)print 1,'ngv,k,ig,nv =',ngv,k,ig,nv
            if(nv.le.0)goto98
c           print 1,'ngv,k,ig,nv =',ngv,k,ig,nv
c           print 4,'G =',(ws(krg+i),i=1,k*ngv)
c           print 4,'a =',(ws(ka+i),i=1,ngv)
c           print 4,'b =',(ws(kb+i),i=1,ngv+1)
c           print 4,'c =',(ws(kc+i),i=1,ngv)
            call formR(nv,k,ig,maxg,ws(ka1),ws(kb1),ws(kc1),ws(kd1),
     *        ws(ke1),ws(krg1),ws(kr1))
c           call checkT(nv,maxg,ws(kr1),ws(ke1),ws(kd1))
            call formT(nv,maxg,ws(kr1),v,ws(ke1))
c           print 4,'T matrix',(v(i),i=1,nv)
c             if(nv.gt.1)print 5,(ws(ke+i),i=1,nv-1)
            call trid(v,ws(ke1),nv)
c           print 4,'eigenvalues of T',(v(i),i=1,nv)
            call insort(nv,v)
c           print 4,'sorted eigenvalues of T',(v(i),i=1,nv)
          endif
          nv0=nv
          f0=f
        endif
        ig=ig1
      endif

   60 continue
      if(alpha.gt.0.D0)then
c  update r for inactive c/s
        iter=iter+1
        if(ninf.gt.0)then
          n_inf=0
          ff=f
          f=0.D0
          do 61 j=n1,lp(1)
            i=abs(ls(j))
            if(w(i).eq.0.D0)then
              if(r(i).ge.0.D0)goto61
              n_inf=n_inf+1
              f=f-r(i)
              goto61
            endif
            ri=r(i)-alpha*w(i)
            if(abs(ri).le.tol)ri=0.D0
            if(r(i).lt.0.D0)then
              if(ri.ge.0.D0)then
c  remove contribution to gradient
                if(i.gt.n)then
                  call saipy(sign(1.D0,dble(ls(j))),a,la,i-n,g,n)
                else
                  g(i)=0.D0
                endif
              else
                n_inf=n_inf+1
                f=f-ri
              endif
            endif
            if(w(i).lt.0.D0)then
              ro=(bu(i)-bl(i))-ri
              if(abs(ro).le.tol)ro=0.D0
              if(ro.lt.ri)then
                ri=ro
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
   61     continue
          if(n_inf.ne.ninf)then
            call iexch(ninf,n_inf)
            call newg
c         elseif(f.ge.ff)then
          elseif(f.ge.eps*ff+ff)then
            goto98
          endif
        else
          n_inf=0
          do 62 j=n1,lp(1)
            i=abs(ls(j))
            if(w(i).eq.0.D0)goto62
            ri=r(i)-alpha*w(i)
            if(w(i).lt.0.D0)then
              ro=(bu(i)-bl(i))-ri
              if(ro.lt.ri)then
                ri=ro
                w(i)=-w(i)
                ls(j)=-ls(j)
              endif
            endif
            if(ri.le.tol)then
              ri=0.D0
              if(i.le.n)then
                if(ls(j).ge.0)then
                  x(i)=bl(i)
                else
                  x(i)=bu(i)
                endif
              endif
            endif
            r(i)=ri
   62     continue
        endif
      endif

      if(alpha.lt.amax)then
        if(ig.gt.0)then
c  continue limited memory SD iterations
          if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *      ''  level = 1    f ='',E16.8,''   rg ='',E12.4,
     *      ''  k ='',I4)')npv,f,rgnorm,k
          if(alpha.gt.0.D0)goto20
          print *,'alpha.le.0'
          goto98
        endif
c  Cauchy step with SE iteration
        k=k+1
        if(p.le.n)then
          ls(pj)=p
          goto15
        endif
c  case p>n: find best inactive simple bound to replace p in ls(pj)
        t=0.D0
        do j=n1,lp(1)
          i=abs(ls(j))
          if(i.le.n)then
            ti=abs(ws(na+i))
            if(ti.gt.t)then
              t=ti
              qj=j
            endif
          endif
        enddo
        if(t.le.snorm*tol)then
          print *,'no suitable simple bound available'
          goto98
        endif
        q=abs(ls(qj))
        ls(qj)=q
        if(iprint.ge.2)write(nout,1)'New free variable',q
        goto70
      endif

   65 continue
      if(iprint.ge.2)
     *  write(nout,*)'New active c/s:  alpha =',alpha,'   q =',q
      if(ig.gt.0)then
c  case alpha=amax and SD step: find best free variable to relax
        k=k-1
        if(qj.le.n)then
c  case: q is a free variable
          if(ws(na+q).gt.0.D0)ls(qj)=-q
          call iexch(ls(qj),ls(n-k))
          ig=0
          if(n_inf.gt.0.and.ninf.eq.0)goto10
          goto15
        endif
        call fbsub(n,n-k,n,a,la,q,w,w,ls,ws(lu1),lws(ll1),.false.)
c       print 4,'w(n-k:n) =',(w(ls(j)),j=n-k,n)
        t=0.D0
        do j=n-k,n
          i=ls(j)
          ti=abs(w(i))/e(i)
          if(ti.gt.t)then
            t=ti
            pj=j
          endif
        enddo
        if(t.le.tol)then
          print *,'no suitable free variable to relax'
          goto98
        endif
        p=ls(pj)
        call iexch(ls(pj),ls(n-k))
        pj=n-k
        if(iprint.ge.2)write(nout,*)'relax free variable',p
      endif


c  return from degeneracy with an equality c/s
   70 continue
      if(qj.ne.pj)then
c  pivot interchange
        if(iprint.ge.2)write(nout,*)'replace',p,' by',q
        if(p.eq.0)print *,'p.eq.0'
        if(p.eq.0)goto98
        call pivot(p,q,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,npv)
        if(ifail.ge.1)then
c         if(ifail.ge.2)return
          if(ifail.eq.7)return
          if(iprint.ge.1)write(nout,*)'failure detected in pivot (1)'
c         print *,'r(q),w(q),q',r(q),w(q),q
          goto98
        endif
        if(rp.gt.0.D0)then
          if(phase.gt.0)print *,'phase =',phase
          call iexch(ls(pj),ls(qj))
          call iexch(ls(lp(1)),ls(qj))
          lp(1)=lp(1)-1
          if(ninf.gt.0)goto15
          goto9
        endif
        if(ig.gt.0)then
          ri=x(p)-bl(p)
          ro=bu(p)-x(p)
          if(ro.lt.ri)then
            ri=ro
            ls(pj)=-p
          endif
          if(ri.le.tol)ri=0.D0
          r(p)=ri
          ig=0
        else
          rpu=max(bu(p)-bl(p)-alpha,0.D0)
          if(alpha.le.rpu)then
            rpu=alpha
          else
            ls(pj)=-ls(pj)
          endif
          if(abs(rpu).le.tol)rpu=0.D0
          r(p)=rpu
        endif
c       print 2,'r(p)',r(p)
        call iexch(ls(pj),ls(qj))
        if(phase.gt.0.and.bu(q).eq.bl(q))then
          peq=peq+1
          call iexch(ls(pj),ls(peq))
        endif
        if(ninf.eq.0)then
          if(phase.eq.0)goto9
          if(phase.eq.1)goto10
        endif
        goto15
      endif
c  opposite bound comes active
      if(ninf.eq.0)then
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *    ''  level = 1    f ='',E16.8)')npv,f
      elseif(phase.eq.0)then
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *    ''  level = 1    f ='',E16.8,''   ninfb ='',I4)')
     *    npv,f,ninf
      else
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *    ''  level = 1    f ='',E16.8,''   ninf ='',I4)')
     *    npv,f,ninf
      endif
      ls(pj)=-ls(pj)
      if(ninf.eq.0.or.ninf.ne.n_inf)goto16
      r(p)=-rp
      goto20

c  unbounded solution case
   75 continue
      irep=irep+1
      if(irep.le.nrep.and.iter.gt.mpiv)then
        mode=4
        if(iprint.ge.1)write(nout,*)
     *    'unbounded solution identified: refinement step #',irep
        goto8
      endif
      ifail=1
c  tidy up x
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
      nv=nv0
      return

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
        if(peq.gt.0.or.k.gt.0)write(nout,1)
     *    '# active equality c/s and free variables = ',peq,k
      endif

   84 continue
      call optest(peq+1,n-k,r,e,ls,rp,pj)

      if(-rp.le.gtol)then
        if(iprint.ge.2)write(nout,*)'return to level 1'
        lev=1
        f=alp(1)
        do j=n1,lp(2)
          r(abs(ls(j)))=0.D0
        enddo
        lev=1
        if(rp.eq.0.D0.and.phase.gt.0)goto25
        goto20
      endif
      call iexch(ls(pj),ls(n-k))
      pj=n-k
      plus=ls(pj).ge.0
      p=abs(ls(pj))
      rp=r(p)
c  compute search direction s in an(.)
      call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1),
     *  e(p),.true.)

        rp=scpr(0.D0,ws(na1),g,n)
        if(ls(pj).lt.0)rp=-rp
        if(rp*r(p).le.0.D0)then
          r(p)=0.D0
          goto84
        endif
        if(abs(rp-r(p)).gt.5.D-1*max(abs(rp),abs(r(p))))then
c       if(abs(rp-r(p)).gt.1.D-1*gnorm)then
          print 2,'2rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
          goto98
        endif

      snorm=e(p)
c  form At.s and denominators
      call form_Ats(n1,lp(lev),n,plus,a,la,ws(na1),w,ls,snorm*tol)
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
        if(wi.le.0.D0)goto90
        if(r(i).lt.0.D0)goto90
        z=(r(i)+tol)/wi
        if(z.ge.alpha)goto90
        alpha=z
        qj=j
   90 continue
      if(qj.eq.0)then
        do j=n1,lp(lev)
          i=abs(ls(j))
          w(i)=min(w(i),0.D0)
          r(i)=0.D0
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
          if(peq.gt.0.or.k.gt.0)print 1,
     *      '# active equality c/s and free variables = ',peq,k
        endif
c       call check(n,lp(1),nmi,kmax,g,a,la,x,bl,bu,r,ls,ws(nb1),f,
c    *    ws,lws,cws,ninf,peq,k,1,p,rp)
        goto30
      endif
      q=abs(ls(qj))
      alpha=r(q)/w(q)
      ff=f+alpha*rp
      if(iprint.ge.2)then
        write(nout,*)'alpha =',alpha,'   p =',p,'   q =',q
        write(nout,2)'r(p),r(q),w(q) =',r(p),r(q),w(q)
      endif
c  test for equality c/s
      if(bu(q).eq.bl(q))then
        do j=n1,lp(2)
          r(abs(ls(j)))=0.D0
        enddo
        lev=1
        f=alp(1)
        alpha=0.D0
        if(iprint.ge.2)write(nout,*)'EQTY:   p =',p,'   q =',q,
     *    '   return to level 1'
        goto70
      endif
      if(alpha.eq.0.D0)then
c  potential degeneracy block at level lev
        if(lev+2.gt.mlp)then
          ifail=8
          return
        endif
        r(q)=0.D0
        plev=n
        do j=n1,lp(lev)
          i=abs(ls(j))
          if(r(i).eq.0.D0)then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            if(bu(i).gt.bl(i))r(i)=1.D0
          endif
        enddo
        if(plev.gt.n1)then
          lev=lev+1
          lp(lev)=plev
          alp(lev)=f
          f=0.D0
          if(iprint.ge.2)write(nout,*)
     *      'degeneracy: increase level to ',lev       
          if(iprint.ge.1)write(nout,'(''pivots ='',I5,A,''level ='',I2,
     *      ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
          goto86
        endif
        qj=n1
      endif
      iter=iter+1
      if(iprint.ge.2)write(nout,*)'replace',p,' by',q
      call pivot(p,q,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,npv)
      if(ifail.ge.1)then
c       if(ifail.ge.2)return
        if(ifail.eq.7)return
c       call iexch(ls(pj),ls(qj))
        if(iprint.ge.1)write(nout,*)'failure detected in pivot (2)'
c       print *,'r(q),w(q),q',r(q),w(q),q
        goto98
      endif
c  update r and f
      do j=n1,lp(lev)
        i=abs(ls(j))
        ri=r(i)-alpha*w(i)
        if(abs(ri).le.tol)ri=0.D0
        r(i)=ri
      enddo
      f=ff
c  exchange a constraint
      r(p)=alpha
      if(r(p).le.tol)r(p)=0.D0
      call iexch(ls(pj),ls(qj))
      if(iprint.ge.1)write(nout,'(''pivots ='',I5,A,''level ='',I2,
     *  ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
      goto80
c  restart sequence
   98 continue
      do i=1,n
        x(i)=min(bu(i),max(bl(i),x(i)))
      enddo
      nk=peq
      do j=peq+1,n-k
        i=abs(ls(j))
        if(i.gt.n)then
          nk=nk+1
          ls(nk)=ls(j)
        endif
      enddo
      k=n-nk
      mode=2
      ires=ires+1
      if(iprint.ge.1)write(nout,*)'major restart #',ires
      tol=1.D1*tol
      if(ires.le.nres)goto7
      ifail=10
      return
 1000 format(a/(e16.5,4e16.5))
 1001 format(a/(i4,1x,e11.5,4(i4,1x,e11.5)))
c1000 format(a/(e18.8,3e19.8))
c1001 format(a/(i3,1x,e14.8,3(i4,1x,e14.8)))
      end

c     block data defaults
c     implicit double precision (a-h,o-z)
c     common/epsc/eps,tol,emin
c     common/repc/sgnf,nrep,npiv,nres
c     common/refactorc/mc,mxmc
c     common/wsc/kk,ll,kkk,lll,mxws,mxlws
c     data  eps,    tol,   emin, sgnf, nrep, npiv, nres, mxmc, kk, ll
c    * /1111.D-19, 1.D-12, 0.D0, 1.D-8,  2,    3,   2,   500,   0,  0/
c     end

      subroutine stmap(n,nm,kmax,maxg)
c  set storage map for workspace in glcpd and auxiliary routines
      implicit double precision (a-h,r-z), integer (i-q)
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/lcpdc/na,na1,nb,nb1,krg,krg1,kr,kr1,
     *  ka,ka1,kb,kb1,kc,kc1,kd,kd1,ke,ke1,lu1,ll1
c  double precision storage (ws)
c  locations 1:kk are user workspace for funct and grad
c  scratch slots of length n+m and n
      na=kk
      na1=kk+1
      nb=na+nm
      nb1=nb+1
c  workspace of length kmax*maxg for reduced gradient vectors
      krg=nb+n
      krg1=krg+1
c  a slot of length maxg*(maxg+1)/2 and 5 slots of length maxg for sweep method
      kr=krg+kmax*maxg
      kr1=kr+1
      ka=kr+maxg*(maxg+1)/2
      ka1=ka+1
      kb=ka+maxg
      kb1=kb+1
      kc=kb+maxg
      kc1=kc+1
      kd=kc+maxg
      kd1=kd+1
      ke=kd+maxg
      ke1=ke+1
c  remaining space for use by denseL.f or schurQR.f
      lu1=ke1+maxg
c  total number of double precision locations required by glcpd
      kkk=nm+n+maxg*(maxg+1)/2+maxg*(kmax+5)
c  integer storage (lws)
c  locations 1:ll are user workspace for funct and grad
c  number of integer locations required by glcpd
      lll=0
c  remaining space for use by denseL.f or schurQR.f
      ll1=ll+1
      return
      end

      subroutine check(n,nm,nmi,kmax,g,a,la,x,bl,bu,r,ls,an,f,
     *  ws,lws,cws,ninf,peq,k,lev,p,alp2)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension g(*),a(*),la(*),x(*),bl(*),bu(*),r(*),ls(*),
     *  an(*),ws(*),lws(*)
      character cws(*)
      common/noutc/nout
      common/epsc/eps,tol,emin
c     if(lev.eq.2)then
c       do i=1,n
c         an(i)=g(i)
c       enddo
c       e=alp2*sign(1.D0,dble(p))
c       i=abs(p)
c       if(i.le.n)then
c         an(i)=an(i)-e
c       else
c         call saipy(-e,a,la,i-n,an,n)
c       endif
c       goto1
c     endif
      j=nmi*(nmi+1)/2
      do i=1,nmi
        j=j-abs(ls(i))
      enddo
      if(j.ne.0)write(nout,*)'indexing error'
      if(j.ne.0)stop
      do j=1,peq
        i=abs(ls(j))
        if(bu(i).gt.bl(i))then
          write(nout,*)'non-equality constraint i =',i
          stop
        endif
      enddo
      do j=n-k+1,n
        i=ls(j)
        if(i.le.0.or.i.gt.n)then
          write(nout,*)'faulty free variable: i, j =',i,j
          stop
        endif
      enddo
      e=0.D0
      do j=n+1,nm
        i=abs(ls(j))
        if(i.le.n)then
          s=x(i)
        else
          s=aiscpr(n,a,la,i-n,x,0.D0)
        endif
        if(ls(j).gt.0)then
c         print *,'i,s,r(i),bl(i)',i,s,r(i),bl(i)
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
      if(e.gt.tol)write(nout,*)'residual error at level 1 = ',e,ie
c     if(e.gt.tol)stop
      if(e.gt.1.D-6)stop
      if(ninf.eq.0)then
        call funct(n,x,ff,ws,lws,cws)
        call grad(n,x,an,ws,lws,cws)
      else
        do i=1,n
          an(i)=0.D0
        enddo
        ff=0.D0
        do j=n+1,nm
          i=abs(ls(j))
          if(r(i).lt.0.D0)then
            ff=ff-r(i)
            if(i.gt.n)then
              call saipy(-sign(1.D0,dble(ls(j))),a,la,i-n,an,n)
            else
              an(i)=an(i)-sign(1.D0,dble(ls(j)))
            endif
          endif
        enddo
      endif
      gnm=sqrt(scpr(0.D0,an,an,n))
      if(lev.eq.1.and.max(abs(f),abs(ff)).lt.1.D20)then
        e=abs(ff-f)
        if(e.gt.tol*max(1.D0,abs(f)))write(nout,*)'function error = ',e,
     *    '   f(x) =',ff
c     if(e.gt.tol)stop
c       if(e.gt.tol*max(1.D0,abs(f)))print 4,'x =',(x(j),j=1,n)
        if(e.gt.tol*max(1.D0,abs(f)))stop
      endif
    1 continue
      e=0.D0
      do j=1,n
c       write(nout,*)'an =',(an(i),i=1,n)
        i=abs(ls(j))
        s=sign(1.D0,dble(ls(j)))
        if(i.le.n)then
c         print *,'i,s,r(i)',i,s,r(i)
          an(i)=an(i)-s*r(i)
          if(j.gt.n-k)then
            s=max(0.D0,bl(i)-x(i),x(i)-bu(i))
          elseif(ls(j).gt.0)then
            s=x(i)-bl(i)
          else
            s=bu(i)-x(i)
          endif
        else
c         print *,'i,s,r(i)',i,s,r(i)
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
      if(e.gt.tol)write(nout,*)'residual error at level 2 = ',e,ie
c     if(e.gt.tol)stop
c     if(e.gt.1.D-6)print 4,'x =',(x(i),i=1,n)
      if(e.gt.1.D-6)stop
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
    2 format(A,5E15.7)
    4 format(A/(5E15.6))
      return
      end
