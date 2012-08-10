      subroutine ksgrad (inext,x,x0,xlb,xub,g,g0,obj,obj0,df,dg,scale,
     1                   delx,ndv,ncon,nobj,nside,fdelt,fdmin,nodim,
     2                   ncdim)
      implicit double precision (a-h,o-z)
      dimension x(*),x0(*),xlb(*),xub(*),g(*),g0(*)
      dimension obj(*),obj0(*),df(nodim,*),dg(ncdim,*),scale(ndv,2)
c
c          routine to determine df and dg by finite differences
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 19 July 1996
c
      if (inext .gt. 0) go to 40
c
c          initialize obj0, x0, and g0
c
      do 10 i = 1,nobj
        obj0(i) = obj(i)
   10 continue
c
      do 20 i = 1,ndv
        x0(i) = x(i)
   20 continue

      if (ncon .le. 0) go to 40
      do 30 i = 1,ncon
        g0(i) = g(i)
   30 continue
c
c          calculate gradients and restore design variables
c
c          note - these gradients are already scaled by the
c                 design variable scale factors
c
   40 continue
      if (inext .eq. 0) go to 80
c
      do 50 i = 1,nobj
        df(i,inext) = (obj(i) - obj0(i)) / delx
   50 continue
c
      if (ncon .le. 0) go to 70
      do 60 i = 1,ncon
        dg(i,inext) = (g(i) - g0(i)) / delx
   60 continue
   70 continue
      x(inext) = x0(inext)
c
c          increment design variable inext and get obj and g
c
   80 continue
      inext = inext + 1
      if (inext .gt. ndv) go to 100
      xhere = x0(inext)
      delmin = fdmin / scale(inext,1)
      delx = abs(xhere * fdelt)
      if (delx .lt. delmin) delx = delmin
c
      if (nside .le. 0) go to 90
c
      xtest = xhere + delx
      if (xtest .le. xub(inext)) go to 90
c
c          adjust delta-x for upper bound
 
      delx1 = delx
      delx = abs(xub(inext) - xhere)
      if (delx .ge. delmin) go to 90
c
c          adjust delta-x for lower bound
c
      delx = -delx1
      xtest = xhere + delx
      if (xtest .lt. xlb(inext)) delx = xlb(inext) - xhere
c
   90 continue
      x(inext) = x0(inext) + delx
      go to 110
c
c          reset inext to zero when finished
c
  100 continue
      inext = 0
  110 continue
c
      return
      end
