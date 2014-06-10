C     COMMON SCALARS
      logical rmfixv

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,rmfixv
      save   /fixvar/
