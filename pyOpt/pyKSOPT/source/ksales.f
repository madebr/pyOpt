      subroutine ksales (a1,a2,a3,f1,f2,f3,alim,alpha)
      implicit double precision (a-h,o-z)
c
c          routine to estimate the next one-dimensional
c          step size alpha, based on the slopes of the
c          previous two trial points.
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 21 April 1989
c
      alpha = 0.0
      alpmin = 2.0 * a3 - a2
      alpmax = 3.0 * a3
c
      if (abs(a2 - a1) .lt. 1.0e-08) go to 10
      x1 = a2
      y1 = (f2 - f1) / (a2 - a1)
c
      if (abs(a3 - a2) .lt. 1.0e-08) go to 10
      x2 = a3
      y2 = (f3 - f2) / (a3 - a2)
c
      a = (y2 - y1) / (x2 - x1)
      b = y1 - a * x1
c
      if (abs(a) .ge. 1.0e-08) alpha = -b / a
   10 continue
      if (alpha .le. 0.0   ) alpha = alpmax
      if (alpha .lt. alpmin) alpha = alpmin
      if (alpha .gt. alpmax) alpha = alpmax
      if (alim .gt. 0.0 .and. alpha .ge. alim) alpha = 0.5 * (a3 + alim)
c
      return
      end

