      subroutine ksd (df,g,dgdx,ng,rho)
      implicit double precision (a-h,o-z)
      dimension g(ng),dgdx(ng)
      data toler /-40.0/
c
c          routine to compute partial k-s w.r.t. x (type 2).
c          partial of g(i) w.r.t. x (e.g.  by finite difference)
c          must be supplied in array dgdx(i).
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification -  9 July 1988
c
      sum1 = 0.0
      sum2 = 0.0
      gmax = g(1)
      if (ng .lt. 2) go to 30
      do 10 i = 2,ng
        if (g(i) .gt. gmax) gmax = g(i)
   10 continue
      do 20 i = 1,ng
        val = rho * (g(i) - gmax)
        if (val .lt. toler) go to 20
        sum1 = sum1 + exp(val) * dgdx(i)
        sum2 = sum2 + exp(val)
   20 continue
   30 continue
      df = dgdx(1)
      if (ng .gt. 1) df = sum1 / sum2
c
      return
      end
