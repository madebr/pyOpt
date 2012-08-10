      subroutine ksoned (inext,jtry,x,x0,xlb,xub,fun,fun0,s,slope,
     1                   alpha,alpmax,ndv,a1,a2,a3,a4,f1,f2,f3,f4,
     2                   alim,atest,ftest,nside,limit,nunit,iprnt3,
     3                   scale,temp,isdflg)
      implicit double precision (a-h,o-z)
      character*4 ip1
      dimension x(*),x0(*),xlb(*),xub(*),s(*),scale(ndv,2),temp(*)
      data ip1 /'(   '/
      data rtoler /0.00001/
      data atoler /0.001/
      data jmax /20/
c
c          routine to perform one-dimensional search
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 19 July 1996
c
      if (inext .gt. 0) go to 50
c
c          initialization
c
      do 10 i = 1,ndv
        x0(i) = x(i)
   10 continue
      alpha = 0.0
      alim = 0.0
      a4 = -1.0
      atest = 0.0
      if (slope .ge. 0.0) go to 500
c
c          estimate initial alpha to limit obj to a small change
c
      alpha  = abs((0.10 * fun0) / slope)
      otest1 = abs(0.10 / slope)
      if (alpha .gt. otest1) alpha = otest1
c
c          estimate initial alpha to limit x(i) to a small change
c
      do 20 i = 1,ndv
        si = abs(s(i))
        if (si .lt. 1.0e-08) go to 20
        xtest  = abs((0.10 * x(i)) / si)
        xtest1 = abs(0.10 / si)
        if (xtest .lt. xtest1) xtest = xtest1
        if (alpha .gt. xtest) alpha = xtest
   20 continue
c
c          estimate initial alpha to prevent a bounds violation
c
      alpmax = 1.0e+10
      if (nside .eq. 0) go to 40
      do 30 i = 1,ndv
        if (abs(s(i)) .lt. 1.0e-08) go to 30
        btest = (xub(i) - x(i)) / s(i)
        if (btest .gt. 0.0 .and. btest .lt. alpmax) alpmax = btest
        btest = (xlb(i) - x(i)) / s(i)
        if (btest .gt. 0.0 .and. btest .lt. alpmax) alpmax = btest
   30 continue
c
      if (alpha .gt. alpmax) alpha = alpmax
      if (alpha .lt. 1.0e-06) alpha = 1.0e-06
c
   40 continue
c
      if (iprnt3 .gt. 0) write (nunit,1010) alpha,alpmax,fun0
c
      a1 = 0.0
      f1 = fun0
      jtry = 0
      inext = 1
      go to 600
c
c          store fun in appropriate location
c
   50 continue
      if (iprnt3 .gt. 0) write (nunit,1050) alpha,fun
      go to (100,200,300,400) inext
c
c          store trial point a2
c
  100 continue
      a2 = alpha
      f2 = fun
      if (f2 .gt. f1) go to 150
      if (limit .eq. 1) go to 500
      alpha = 2.0 * a2
      inext = 3
      go to 600
c
  150 continue
      a3 = a2
      f3 = f2
      alpha = 0.1 * (9.0 * a1 + a3)
      if (alpha .le. 1.0e-08) go to 500
      inext = 2
      go to 600
c
c          point between a1 and a3
c
  200 continue
      jtry = jtry + 1
      if (jtry .gt. jmax) go to 500
c
      a2 = alpha
      f2 = fun
      if (f2 .gt. f1) go to 150
      go to 350
c
c          point after a2
c
  300 continue
c
      jtry = jtry + 1
      if (jtry .gt. jmax) go to 500
c
      a3 = alpha
      f3 = fun
      if (f3 .lt. f2) go to 310
      f12 = f1 - f2
      f32 = f3 - f2
      if (f12 .lt. 1.e-6) f12 = 1.e-6
      if (f32 .lt. 10.0 * f12) go to 350
      alpha = 0.5 * (a2 + a3)
      alim = a3
      inext = 3
      go to 600
  310 continue
      if (limit .eq. 1) go to 500
      call ksales (a1,a2,a3,f1,f2,f3,alim,alpha)
      a1 = a2
      f1 = f2
      a2 = a3
      f2 = f3
      inext = 3
      go to 600
c
c          perform quadratic or cubic interpolation
c
  350 continue
      if (a4 .le. 0.0) call ksquad (a1,a2,a3,f1,f2,f3,astar,fstar)
      if (a4 .gt. 0.0) call kscube (a1,a2,a3,a4,f1,f2,f3,f4,astar,fstar)
      alpha = astar
      inext = 4
      go to 600
c
c          check for convergence
c
  400 continue
c
      jtry = jtry + 1
      if (jtry .gt. jmax) go to 500
c
      if (atest .le. 0.0) go to 450
      adiv = alpha
      fdiv = abs(fun)
      if (adiv .lt. 1.e-06) adiv = 1.0
      if (fdiv .lt. 1.e-06) fdiv = 1.0
      aa = abs(atest - alpha)
      ff = abs(ftest - fun)
      atol = atoler * abs(alpha)
      ftol = atoler * abs(fun)
      if (atol .lt. 1.e-06) atol = 1.0e-06
      if (ftol .lt. 1.e-06) ftol = 1.0e-06
      if (aa .lt. atol .or. ff .lt. ftol) go to 500
      aa = aa / adiv
      ff = ff / fdiv
      if (aa .lt. rtoler .or. ff .lt. rtoler) go to 500
  450 continue
      atest = alpha
      ftest = fun
      call ksqmin (a1,a2,a3,a4,alpha,f1,f2,f3,f4,fun)
      go to 350
c
c          end of one-dimensional search
c
  500 continue
      inext = 0
      if (jtry .gt. jmax .or. limit .eq. 1) isdflg = -1
c
c          compute x-vector and get new function
c
  600 continue
      limit = 0
      if ( alpha .lt. alpmax) go to 610
      alpha = alpmax
      limit = 1
  610 continue
      do 620 i = 1,ndv
        x(i) = x0(i) + alpha * s(i)
  620 continue
      if (iprnt3 .lt. 1) return
      if (inext .gt. 0) write (nunit,1020) inext,jtry,alpha
      if (inext .eq. 0) write (nunit,1080) inext,jtry,alpha
      if (limit .eq. 1) write (nunit,1060)
      if (iprnt3 .lt. 2) return
      call ksunsc (temp,x,scale,ndv)
      if (inext .ne. 0) write (nunit,1030)
      if (inext .eq. 0) write (nunit,1070)
      write (nunit,1000) (ip1,i,temp(i),i=1,ndv)
      write (nunit,1040)
      return
c
 1000 format (3(3x,a1,i5,2h)=,e14.7))
 1010 format (/43h       one-dimensional minimization started/
     1        /32h       initial alpha estimate = ,e12.6/
     2         32h       maximum alpha allowed  = ,e12.6/
     3         32h       initial k-s function   = ,e12.6/)
 1020 format (15h       inext = ,i2,2x,8h jtry = ,i2,2x,
     1        18h proposed alpha = ,e12.6)
 1030 format (/35h       proposed design variable set/)
 1040 format (/)
 1050 format (/23h       current alpha = ,e12.6,14h   function = ,e12.6)
 1060 format (/51h       upper or lower bounds limit has been reached)
 1070 format (/32h       final design variable set/)
 1080 format (15h       inext = ,i2,2x,8h jtry = ,i2,2x,
     1        15h final alpha = ,e12.6)
      end
