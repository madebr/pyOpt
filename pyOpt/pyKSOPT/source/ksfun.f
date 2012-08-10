      subroutine ksfun (fun,obj,g,rho,fscale,offset,ncon,nobj,temp)
      implicit double precision (a-h,o-z)
      dimension obj(*),g(*),fscale(*),offset(*),temp(*)
c
c          routine to compute unconstrained function to be minimized
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 17 July 1996
c
      do 10 i = 1,nobj
        temp(i) = obj(i) / fscale(i) + offset(i)
   10 continue
      j = nobj
      if (ncon .le. 0) go to 30
      do 20 i = 1,ncon
        j = j + 1
        temp(j) = g(i)
   20 continue
   30 continue
c
      call ks (fun,temp,j,rho)
c
      return
      end
