      subroutine ksside (x,xlb,xub,iside,dfun,ndv,nside)
      implicit double precision (a-h,o-z)
      dimension x(1),xlb(1),xub(1),iside(1),dfun(1)
c
c          routine to compute a vector of flags iside(i),i=1,ndv
c          iside(i) =  0  -- d.v. i is not at a side constraint
c                   = -1  -- d.v. i is at a lower bound
c                   = +1  -- d.v. i is at an upper bound
c                   = 999 -- d.v. i is at both upper and lower bounds
c
c          also zero components of dfun that would violate
c          side constraints
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification -  9 July 1988
c
      do 10 i = 1,ndv
        iside(i) = 0
        if (nside .eq. 0) go to 10
c
        xl = abs(xlb(i))
        if (xl .lt. 1.0) xl = 1.0
        xx = (xlb(i) - x(i)) / xl
        if (xx .ge. -1.0e-6) iside(i) = -1
        if (iside(i) .lt. 0 .and. dfun(i) .gt. 0.0) dfun(i) = 0.0
c
        xu = abs(xub(i))
        if (xu .lt. 1.0) xu = 1.0
        xx = (x(i) - xub(i)) / xu
        if (xx .ge. -1.0e-6 .and. iside(i) .ne. 0) iside(i) = 999
        if (xx .ge. -1.0e-6 .and. iside(i) .eq. 0) iside(i) = 1
        if (iside(i) .gt. 0 .and. dfun(i) .lt. 0.0) dfun(i) = 0.0
   10 continue
c
      return
      end
