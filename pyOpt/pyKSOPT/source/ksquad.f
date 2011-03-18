      subroutine ksquad (a1,a2,a3,f1,f2,f3,astar,fstar)
      implicit double precision (a-h,o-z)
c
c          routine to compute the minimum of a quadratic
c          curve through points a1,f1 a2,f2 and a3,f3.
c          f2 must be less than f1 and f3.
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 27 April 1989
c
      f21 = f2 - f1
      a21 = a2 - a1
      if (a21 .lt. 1.0e-12) go to 10
      f31 = f3 - f1
      a31 = a3 - a1
      a3121 = a31 / a21
      a21s = a2 * a2 - a1 * a1
      a31s = a3 * a3 - a1 * a1
c
      d = a31s - a21s * a3121
      if (d .lt. 1.0e-12) go to 10
c
      a = (f31 - f21 * a3121) / d
      b = (f21 - a * a21s) / a21
      c = f1 - b * a1 - a * a1 * a1
c
      if (a .lt. 1.0e-12) go to 10
      astar = -b / (2.0 * a)
      if (astar .lt. a1 .or. astar .gt. a3) go to 10
      fstar = a * astar * astar + b * astar + c
      go to 30
c
   10 continue
      astar = a1
      fstar = f1
      if (f2 .ge. fstar) go to 20
      astar = a2
      fstar = f2
   20 continue
      if (f3 .ge. fstar) go to 30
      astar = a3
      fstar = f3
c
   30 continue
      return
      end
