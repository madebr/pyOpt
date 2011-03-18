      subroutine ks (f,g,ng,rho)
      implicit double precision (a-h,o-z)
      dimension g(ng)
      data toler /-40.0/
c
c          routine to compute k-s function (type 2)
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification -  9 July 1988
c
      sum = 0.0
      gmax = g(1)
      if (ng .lt. 2) go to 30
      do 10 i = 2,ng
        if (g(i) .gt. gmax) gmax = g(i)
   10 continue
      do 20 i = 1,ng
        val = rho * (g(i) - gmax)
        if (val .lt. toler) go to 20
        sum = sum + exp(val)
   20 continue
   30 continue
      f = gmax
      if (ng .gt. 1) f = gmax + log(sum) / rho
c
      return
      end
