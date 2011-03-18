      subroutine kscube (a1,a2,a3,a4,f1,f2,f3,f4,astar,fstar)
      implicit double precision (a-h,o-z)
c
c          routine to compute the minimum of a cubic
c          curve through points a1,f1 a2,f2 a3,f3 and a4,f4.
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 27 April 1989
c
      a1c = a1 * a1 * a1
      a2c = a2 * a2 * a2
      a3c = a3 * a3 * a3
      a4c = a4 * a4 * a4
c
      q1 = a3c * (a2 - a1) - a2c * (a3 - a1) + a1c * (a3 - a2)
      q2 = a4c * (a2 - a1) - a2c * (a4 - a1) + a1c * (a4 - a2)
      q3 = (a3 - a2) * (a2 - a1) * (a3 - a1)
      q4 = (a4 - a2) * (a2 - a1) * (a4 - a1)
      q5 = f3 * (a2 - a1) - f2 * (a3 - a1) + f1 * (a3 - a2)
      q6 = f4 * (a2 - a1) - f2 * (a4 - a1) + f1 * (a4 - a2)
c
      d1 = q2 * q3 - q1 * q4
      d2 = a2 - a1
      if (d1 .lt. 1.0e-12 .or. d2 .lt. 1.0e-12) go to 10
c
      b3 = (q3 * q6 - q4 * q5) / d1
      b2 = (q5 - b3 * q1) / q3
      b1 = (f2 - f1) / d2 - b3 * (a2c - a1c) / d2 - b2 * (a1 + a2)
      b0 = f1 - b1 * a1 - b2 * a1 * a1 - b3 * a1c
c
      bb = b2 * b2 - 3.0 * b1 * b3
      if (bb .lt. 0.0 .or. abs(b3) .lt. 1.0e-12) go to 10
      astar = (-b2 + sqrt(bb)) / (3.0 * b3)
      if (astar .lt. a1 .or. astar .gt. a4) go to 10
      as2 = astar * astar
      as3 = astar * as2
      fstar = b0 + b1 * astar + b2 * as2 + b3 * as3
      go to 40
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
   30 continue
      if (f4 .ge. fstar) go to 40
      astar = a4
      fstar = f4
c
   40 continue
      return
      end
 