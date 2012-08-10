      subroutine ksxlim (x,xlb,xub,ndv,nside)
      implicit double precision (a-h,o-z)
      dimension x(*),xlb(*),xub(*)
c
c          routine to insure x-vector does not violate
c          upper or lower bounds
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 17 July 1996
c
      if (nside .eq. 0) go to 20
      do 10 i = 1,ndv
        xx  = x(i)
        xl  = xlb(i)
        xu  = xub(i)
        if (xx .lt. xl) xx = xl
        if (xx .gt. xu) xx = xu
        x(i) = xx
   10 continue
   20 continue
c
      return
      end
