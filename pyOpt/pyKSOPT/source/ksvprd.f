      subroutine ksvprd (x,y,prod,nrow)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
c
c          routine to multiply vectors x and y
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 17 July 1996
c
      prod = 0.0
      do 10 i = 1,nrow
        prod = prod + x(i) * y(i)
   10 continue
c
      return
      end
