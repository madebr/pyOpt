      subroutine kshmul (a,b,x,nrow)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*),x(*)
c
c          routine to perform matrix multiplication a * b = x
c          for a triangular matrix a
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 17 July 1996
c
      k = 0
      do 10 i = 1,nrow
        x(i) = 0.0
   10 continue
      do 20 i = 1,nrow
      do 20 j = 1,i
        k = k + 1
        x(i) = x(i) + a(k) * b(j)
        if (j .lt. i) x(j) = x(j) + a(k) * b(i)
   20 continue
c
      return
      end
