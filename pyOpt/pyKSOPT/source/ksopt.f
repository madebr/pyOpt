      subroutine ksopt (isel,x,obj,g,df,dg,nomax,ngmax,work)
      implicit double precision (a-h,o-z)
      common /kscomm/ rdf   ,adf   ,fdl   ,fdm   ,rho   ,drho  ,rhomax,
     1                fun0  ,slope ,delx  ,alpha ,alpmax,a1    ,a2    ,
     2                a3    ,a4    ,f1    ,f2    ,f3    ,f4    ,alim  ,
     3                atest ,ftest ,ifscl ,ifoff ,isx   ,isx0  ,isxlb ,
     4                isxub ,iscl  ,ig0   ,idf   ,islp  ,iobj0 ,iy    ,
     5                ip    ,ih    ,ihess ,iside ,isact ,idobj ,idg   ,
     6                itmp1 ,itmp2 ,inext ,jnext ,jsel  ,itcnt ,icntr ,
     7                icnta ,isdflg,isdrst,ifncl ,nunit ,ndv   ,ncon  ,
     8                nobj  ,nside ,nscale,iprnt ,itmax ,igrad ,limit
      dimension x(*),obj(*),g(*),df(nomax,*),dg(ngmax,*),work(*)
c
c          KSOPT is a general purpose multiobjective optimization
c          program which performs a constrained to unconstrained
c          transformation on the given optimization problem.
c          The problem is then solved using a Davidon-Fletcher-Powell
c          unconstrained optimization method.
c
c          input parameters are
c               obj   -- current set of objective function values
c               g     -- current set of constraints
c               df    -- derivatives of objective funcs. w.r.t. x-vector
c               dg    -- derivatives of constraints w.r.t. x-vector
c               nomax -- first dimension of df matrix
c               ngmax -- first dimension of dg matrix
c               work  -- working array from routine ksinit
c                        which must not be altered during optimization.
c
c          output parameters are
c               isel  -- flag requesting user supplied information
c                        =0 -- minimization terminated
c                        =1 -- user supplies obj and g
c                        =2 -- user supplies df and dg
c               x     -- new set of design variables
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 19 July 1996
c
c
c          choose a path based on the value of jsel
c
      call kscomg (work(1))
      iprnt3 = iprnt / 100
      ileft  = mod(iprnt,100)
      iprnt2 = ileft / 10
      iprnt1 = mod(ileft,10)
      nodim  = nobj
      if (nodim .le. 0) nodim = 1
      ncdim  = ncon
      if (ncdim .le. 0) ncdim = 1
c
      isel = 0
      if (jsel .le. 0) return
c
      go to (100,200,300,400,500,600) jsel
c
c          begin an iteration
c
  100 continue
      isel = 1
      jsel = 2
      inext = 0
      if (itcnt .eq. 1) go to 900
c
c          compute each objective function's scale and offset
c
  200 continue
      call ksando (obj,g,work(ifscl),work(ifoff),ncon,nobj)
c
c          calculate initial unconstrained function
c
      call ksfun (fun0,obj,g,rho,work(ifscl),work(ifoff),ncon,nobj,
     1            work(itmp1))
      if (igrad .eq. 1) go to 350
c
c          get finite difference gradients of obj and constraints
c
  300 continue
      jsel = 3
      isel = 1
      call ksgrad (inext,work(isx),work(isx0),work(isxlb),work(isxub),
     1             g,work(ig0),obj,work(iobj0),work(idobj),work(idg),
     2             work(iscl),delx,ndv,ncon,nobj,nside,fdl,fdm,
     3             nodim,ncdim)
      if (inext .ne. 0) go to 900
c
      do 310 i = 1,nobj
        obj(i) = work(iobj0 + i - 1)
  310 continue
c
      if (ncon .le. 0) go to 330
      do 320 i = 1,ncon
        g(i) = work(ig0 + i - 1)
  320 continue
  330 continue
      go to 440
c
c          get analytical gradients of obj and constraints
c
  350 continue
      do 360 i = 1,nobj
        work(iobj0 + i - 1) = obj(i)
  360 continue
c
      if (ncon .eq. 0) go to 380
      do 370 i = 1,ncon
        work(ig0 + i - 1) = g(i)
  370 continue
  380 continue
      isel = 2
      jsel = 4
      go to 900
c
c          scale analytical gradients of obj and constraints
c
  400 continue
      do 430 i = 1,ndv
        ss = work(iscl + i - 1)
c
        do 410 l = 1,nobj
          j = (i - 1) * nobj + l
          work(idobj + j - 1) = df(l,i) * ss
  410   continue
c
        if (ncon .eq. 0) go to 430
        do 420 k = 1,ncon
          j = (i - 1) * ncon + k
          work(idg + j - 1) = dg(k,i) * ss
  420   continue
  430 continue
c
c          compute gradient of unconstrained function
c          to be minimizied
c
  440 continue
      call ksdfun (work(idf),work(iobj0),work(ifscl),work(ifoff),
     1             work(idobj),work(ig0),work(idg),rho,ndv,ncon,nobj,
     2             work(itmp1),work(itmp2),nodim,ncdim)
c
c          compute an initial or restart hessian matrix
c
      call kshess (work(ihess),work(iobj0),work(ifscl),work(ifoff),
     1             work(idobj),work(ig0),work(idg),rho,
     2             ndv,ncon,nobj,work(itmp1),work(itmp2),nodim,ncdim)
c
c          check for side constraint violation
c
      call ksside (work(isx),work(isxlb),work(isxub),work(iside),
     1             work(idf),ndv,nside)
c
c          compute search direction
c
      call ksdfp (work(isx),work(iside),work(isact),work(idf),ndv,
     1            work(islp),slope,work(iy),work(ip),work(ih),
     2            work(ihess),isdflg)
c
c          print beginning of iteration information
c
      ipflag = 1
      call ksprnt (ipflag,iprnt1,iprnt2,x,work(iobj0),work(ig0),
     1             work(idobj),work(idg),work(iside),work(iscl),
     2             nodim,ncdim,work(itmp1),work(1))
c
c          compute unconstrained function to be minimized
c
  500 continue
      if (inext .eq. 0) go to 510
      call ksfun (fun,obj,g,rho,work(ifscl),work(ifoff),ncon,nobj,
     1            work(itmp1))
  510 continue
c
c          perform one-dimensional search
c
      jsel = 5
      isel = 1
      call ksoned (inext,jnext,work(isx),work(isx0),work(isxlb),
     1             work(isxub),fun,fun0,work(islp),slope,alpha,alpmax,
     2             ndv,a1,a2,a3,a4,f1,f2,f3,f4,alim,atest,ftest,
     3             nside,limit,nunit,iprnt3,work(iscl),work(itmp1),
     4             isdflg)
      if (inext .ne. 0) go to 900
      jsel = 6
      isel = 1
      go to 900
c
c          test for termination criteria
c
  600 continue
      call ksfun (fun,obj,g,rho,work(ifscl),work(ifoff),ncon,nobj,
     1            work(itmp1))
      if (itcnt .ge. itmax) go to 700
      af    = abs(fun)
      adelf = abs(fun - fun0)
      rdelf = adelf
      if (af .gt. 1.0e-12) rdelf = adelf / af
      icntr = icntr + 1
      icnta = icnta + 1
      if (rdelf .ge. rdf) icntr = 0
      if (adelf .ge. adf) icnta = 0
      if (icntr .ge. 3 .or. icnta .ge. 3) go to 700
c
c          go to next iteration
c
      itcnt = itcnt + 1
c
c          increment rho when necessary
c
      if (icntr .eq. 0 .and. icnta .eq. 0) go to 650
      rho = rho + drho
      if (rho .gt. rhomax) rho = rhomax
  650 continue
c
c          set flag for hessian matrix restart
c
      isdflg = isdflg + 1
      if (isdflg .gt. isdrst) isdflg = 0
c
c          re-scale when necessary
c
      iscale = 0
      if (nscale .gt. 0) iscale = nscale
      if (iscale .ne. 0) itflag = mod(itcnt - 1,iscale)
      if (iscale .eq. 0) itflag = 1
c
c          currently re-scale only when hessian matrix is reset
c          28 March 1991  g.a.w.
c
      itflag = 1
      if (isdflg .eq. 0) itflag = 0
      if (itflag .eq. 0) call ksscal (work(isx),work(isx0),work(isxlb),
     1                                work(isxub),work(iscl),
     2                                ndv,nside,nscale)
      go to 100
c
c          terminate
c
  700 continue
      isel = 0
c
      ipflag = 2
      call ksprnt (ipflag,iprnt1,iprnt2,x,obj,g,work(idobj),work(idg),
     1             work(iside),work(iscl),nodim,ncdim,work(itmp1),
     2             work(1))
c
c          un-scale design variables and return to calling routine
c
  900 continue
      call ksxlim (work(isx),work(isxlb),work(isxub),ndv,nside)
      call ksunsc (x,work(isx),work(iscl),ndv)
      if (isel .eq. 1) ifncl = ifncl + 1
      call kscomp (work(1))
      return
c
      end
