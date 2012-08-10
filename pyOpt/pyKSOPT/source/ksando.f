      subroutine ksando (obj,g,fscale,offset,ncon,nobj)
      implicit double precision (a-h,o-z)
      dimension obj(*),g(*),fscale(*),offset(*)
      data eps/1.0e-6/
c
c          routine to compute scale and offset for obj
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 17 July 1996
c
      if (ncon .le. 0) go to 20
      gmax = g(1)
      if (ncon .eq. 1) go to 20
      do 10 i = 2,ncon
        if (g(i) .gt. gmax) gmax = g(i)
   10 continue
c
   20 continue
      do 30 j = 1,nobj
        fscl = abs(obj(j))
        foff = abs(obj(j))
        if (fscl .lt. 1.0) fscl = 1.0
        if (foff .gt. 1.0) foff = 1.0
        fscale(j) = fscl
        offset(j) = 0.0
c
c          compute offset if ncon > 0
c
        if (ncon .le. 0) go to 30
        if (obj(j) .ge. 0.0) offset(j) = -foff - gmax - eps
        if (obj(j) .lt. 0.0) offset(j) =  foff - gmax + eps
   30 continue
c
      return
      end
