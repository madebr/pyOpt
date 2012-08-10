      subroutine ksscal (x,x0,xlb,xub,scale,ndv,nside,nscale)
      implicit double precision (a-h,o-z)
      dimension x(*),x0(*),xlb(*),xub(*),scale(ndv,2)
c
c          routine to compute new scaling vector and re-scale
c          design variables and lower and upper bounds
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 17 July 1996
c
      do 10 i = 1,ndv
        xx    = x(i)
        xx0   = x0(i)
        sold  = scale(i,1)
        snew  = sold
        if (nscale .gt. 0) snew = abs(xx * sold)
        if (nscale .lt. 0) snew = abs(scale(i,2))
        if (snew .lt. 1.0e-08) snew = 1.0
        scale(i,1) = snew
        x(i)    = xx  * sold / snew
        x0(i)   = xx0 * sold / snew
        if (nside .le. 0) go to 10
        xlb(i)  = xlb(i) * sold / snew
        xub(i)  = xub(i) * sold / snew
   10 continue
c
      return
      end
