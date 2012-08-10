      subroutine ksinit (x,xlb,xub,scale,work,mdv,mcon,mobj,mside,
     1                   mscale,jprnt,jtmax,jgrad,jsdrst,rdfun,adfun,
     2                   fdelt,fdmin,rhomn,rhomx,rhodl,munit,ireq)
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
      dimension x(*),xlb(*),xub(*),scale(*),work(*)
c
c          routine to initialize optimization parameters
c          for routine ksopt
c
c          input parameters are
c
c               x     -- initial design variable values
c
c               xlb   -- lower bounds on design variables
c
c               xub   -- upper bounds on design variables
c
c               scale -- scale factors for design variables
c
c               work  -- scratch array
c
c               mdv   -- number of design variables
c
c               mcon  -- number of constraints
c
c               mobj  -- number of objective functions
c
c               mside -- flag  =0 if no side constraints
c
c               mscale-- flag selecting design variable scaling
c                        =0 -- no scaling
c                        <0 -- user supplied scaling in vector scale
c                        >0 -- automatic scaling by ksopt every mscale
c                              iterations
c
c               jprnt -- print level flag (3 digit)
c                        each print level includes all lower print levels
c
c                        hundreds digit - one-dimensional search print
c                                         0 = no print
c                                         1 = alpha and k-s function
c                                         2 = proposed d.v. vector
c
c                        tens     digit - gradient print
c                                         0 = no print
c                                         1 = df and dg
c
c                        ones     digit - d.v. and constraint print
c                                         0 = no print
c                                         1 = initial and final iterations
c                                         2 = all iterations
c                                         3 = slope and search direction
c                                         4 = hessian matrix
c
c               jtmax -- maximum number of iterations
c                        (default is 20)
c
c               jgrad -- flag selecting user supplied information
c                        =0 -- user supplies function and constraints
c                              with gradients computed by finite
c                              differences from within ksopt
c                        =1 -- user supplies all function, constraint,
c                              and gradient information
c                        (default is 0)
c
c               jsdrst-- number of iterations before restarting search
c                        direction with steepest descent
c                        (default is ndv+1)
c
c               rdfun -- relative function change for termination
c                        (default is 0.01)
c
c               adfun -- absolute function change for termination
c                        (default is 0.0)
c
c               fdelt -- step size for computing finite differences
c                        derivatives
c                        (default is 0.01)
c
c               fdmin -- minimum difference for computing
c                        finite differences
c                        (default is 0.0001)
c
c               rhomn -- minimum multiplier for ks function
c                        (default is 5.0)
c
c               rhomx -- maximum multiplier for ks function
c                        (default is 100.0)
c
c               rhodl -- increment for rho
c                        (default is computed internally)
c
c               munit -- fortran unit number for output file
c                        (default is 0)
c
c          output parameters are
c
c               work  -- array containing information required
c                        by routine ksopt.  this array must not
c                        be altered during execution of ksopt.
c
c               ireq  -- required length of work array
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 19 July 1996
c
c          these subroutines are compatible with fortran 77
c
      jsel   = 1
      ndv    = mdv
      ncon   = mcon
      nobj   = mobj
      nside  = 0
      nscale = mscale
      iprnt  = 1
      itmax  = 20
      igrad  = 0
      isdrst = jsdrst
      rdf    = 0.01
      adf    = 0.0
      fdl    = 0.01
      fdm    = 0.0001
      rho    = 5.0
      rhomax = 100.0
      drho   = 0.0
      nunit  = 0
      itcnt  = 1
      icntr  = 0
      icnta  = 0
      isdflg = 0
      ifncl  = 0
c
      if (mdv   .le. 0)   ndv   = 0
      if (mcon  .le. 0)   ncon  = 0
      if (mobj  .le. 0)   nobj  = 0
      if (mside .ne. 0)   nside = 1
      if (jprnt .ge. 0)   iprnt = jprnt
      if (jtmax .gt. 0)   itmax = jtmax
      if (jgrad .ne. 0)   igrad = 1
      if (jsdrst.le. 0)   isdrst= ndv + 1
      if (rdfun .gt. 0.0) rdf   = rdfun
      if (adfun .gt. 0.0) adf   = adfun
      if (fdelt .gt. 0.0) fdl   = fdelt
      if (fdmin .gt. 0.0) fdm   = fdmin
      if (rhomn .gt. 0.0) rho   = rhomn
      if (rhomx .ge. rho) rhomax= rhomx
      if (rhodl .gt. 0.0) drho  = rhodl
      if (munit .gt. 0)   nunit = munit
c
      ntmp  = 2 * ndv
      if (ntmp .lt. (nobj + ncon)) ntmp = nobj + ncon
      ifscl = 64
      ifoff = ifscl + nobj
      isx   = ifoff + nobj
      isx0  = isx   + ndv
      isxlb = isx0  + ndv
      isxub = isxlb + ndv
      iscl  = isxub + ndv
      ig0   = iscl  + 2 * ndv
      idf   = ig0   + ncon
      islp  = idf   + ndv
      iobj0 = islp  + ndv
      iy    = iobj0 + nobj
      ip    = iy    + ndv
      ih    = ip    + ndv
      ihess = ih    + ndv * (ndv + 1) / 2
      iside = ihess + ndv * (ndv + 1) / 2
      isact = iside + ndv
      idobj = isact + ndv
      idg   = idobj + nobj * ndv
      itmp1 = idg   + ncon * ndv
      itmp2 = itmp1 + ntmp
      ireq  = itmp2 + ntmp - 1
c
      if (ndv .eq. 0) go to 40
      do 30 i = 1,ndv
        j = i - 1
        xx  = x(i)
        if (nside .eq. 0) go to 20
        xbl = xlb(i)
        xbu = xub(i)
        if (xbl .le. xbu) go to 10
        xb = (xbl + xbu) / 2.0
        xbl = xb
        xbu = xb
   10 continue
        if (xx .lt. xbl) xx = xbl
        if (xx .gt. xbu) xx = xbu
        x(i) = xx
   20 continue
        ss = 1.0
        if (nscale .gt. 0) ss = abs(xx)
        if (nscale .lt. 0) ss = abs(scale(i))
        if (ss .lt. 1.0e-04) ss = 1.0
        work(iscl  + j) = ss
        work(isx   + j) = xx / ss
        work(isx0  + j) = xx / ss
        if (nside .eq. 0) go to 30
        work(isxlb + j) = xbl / ss
        work(isxub + j) = xbu / ss
   30 continue
c
      if (drho .gt. 0.0) go to 40
      drho = (rhomax - rho) / 5.0
      if (drho .lt. 10.0) drho = 10.0
      if (drho .gt. 40.0) drho = 40.0
c
   40 continue
c
      if (iprnt .eq. 0 .and. ndv .gt. 0 .and. nobj .gt. 0) go to 60
      write (nunit,70)
      write (nunit,80) ireq
      if (ndv .eq. 0) write (nunit,90)
      if (nobj .eq. 0) write (nunit,100)
      write (nunit,110) ndv,ncon,nobj,nside,nscale,iprnt,itmax,igrad
      write (nunit,120) rdf,adf,fdl,fdm,rho,rhomax,drho
      if (nside .ne. 0) go to 50
      write (nunit,130)
      write (nunit,140) (i,x(i),work(iscl+i-1),i=1,ndv)
      go to 60
   50 continue
      write (nunit,150)
      write (nunit,160) (i,x(i),work(iscl+i-1),
     1                   work(isxlb+i-1)*work(iscl+i-1),
     2                   work(isxub+i-1)*work(iscl+i-1),i=1,ndv)
   60 continue
      if (ndv .eq. 0 .or. nobj .eq. 0) jsel = 0
      call kscomp (work(1))
      return
c
   70 format (1h1,
     1        20x,40h========================================/
     1        21x,40h=                                      =/
     2        21x,40h=    KSOPT Multiobjective Optimizer    =/
     3        21x,40h=                                      =/
     4        21x,40h=    Version 2.8         19 July 96    =/
     5        21x,40h=                                      =/
     6        21x,40h=    Written by    Gregory A. Wrenn    =/
     7        21x,40h=                                      =/
     8        21x,40h========================================///)
   80 format (10x,43hThe work array must be dimensioned at least,
     1        i10,12h words long.//)
   90 format (/53h   Number of design variables is zero -- optimization,
     1        16h is not possible//)
  100 format (/44h   Number of objective functions is zero -- ,
     1        28hoptimization is not possible//)
  110 format (/10h      ndv=,i5,10h     ncon=,i5,10h     nobj=,i5,
     1         10h    nside=,i5,//10h   nscale=,i5,10h   iprint=,i5,
     2         10h    itmax=,i5,10h    igrad=,i5)
  120 format (/10h    rdfun=,e14.7,2x,10h    adfun=,e14.7/
     1        /10h    fdelt=,e14.7,2x,10h    fdmin=,e14.7/
     2        /10h   rhomin=,e14.7,2x,10h   rhomax=,e14.7,2x,
     3         10h     drho=,e14.7)
  130 format (/42(1h-)/
     1        5x,6h d.v. ,6x,7hinitial,9x,6h scale/
     2        5x,6hnumber,6x,7h value ,9x,6hfactor/
     3        1x,42(1h-))
  140 format (4x,i5,4x,e12.5,4x,e12.5)
  150 format (/73(1h-)/
     1        5x,6h d.v. ,6x,7hinitial,9x,6h scale,10x,5hlower,
     2        11x,5hupper/
     3        5x,6hnumber,6x,7h value ,9x,6hfactor,10x,5hbound,
     4        11x,5hbound/
     5        1x,73(1h-))
  160 format (4x,i5,4x,e12.5,4x,e12.5,4x,e12.5,4x,e12.5)
      end
