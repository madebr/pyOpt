      subroutine ksunsc (x,sx,scale,ndv)
      implicit double precision (a-h,o-z)
      dimension x(*),sx(*),scale(ndv,2)
c
c          routine to un-scale design variables before returning
c          to the real world
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 17 July 1996
c
      do 10 i = 1,ndv
        ss = scale(i,1)
        x(i) = sx(i) * ss
   10 continue
c
      return
      end
