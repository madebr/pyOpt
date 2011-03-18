      subroutine ksqmin (a1,a2,a3,a4,alpha,f1,f2,f3,f4,fun)
      implicit double precision (a-h,o-z)
c
c          routine to place alpha and fun in the appropriate
c          position relative to the three or four previously
c          defined points a1,f1 through a4,f4, dropping the
c          most distant point from the current minimum.
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 26 April 1989
c
      if (a4 .gt. 0.0) go to 20
c
      if (alpha .lt. a2) go to 10
      a4 = a3
      f4 = f3
      a3 = alpha
      f3 = fun
      go to 50
   10 continue
      a4 = a3
      f4 = f3
      a3 = a2
      f3 = f2
      a2 = alpha
      f2 = fun
      go to 50
   20 continue
      if (alpha .le. a2) go to 10
      if (alpha .lt. a3) go to 30
      a1 = a2
      f1 = f2
      a2 = a3
      f2 = f3
      a3 = alpha
      f3 = fun
      go to 50
   30 continue
      if ( (alpha - a1) .gt. (a4 - alpha) ) go to 40
      a4 = a3
      f4 = f3
      a3 = alpha
      f3 = fun
      go to 50
   40 continue
      a1 = a2
      f1 = f2
      a2 = alpha
      f2 = fun
c
   50 continue
      return
      end
