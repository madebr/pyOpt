      subroutine ksdfun (dfun,obj,fscale,offset,df,g,dg,rho,ndv,ncon,
     1                   nobj,temp1,temp2,nodim,ncdim)
      implicit double precision (a-h,o-z)
      dimension dfun(*),obj(*),fscale(*),offset(*),df(nodim,*)
      dimension g(*),dg(ncdim,*),temp1(*),temp2(*)
c
c          routine to compute gradients of function to be minimized
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 19 July 1996
c
      do 10 i = 1,nobj
        temp1(i) = obj(i) / fscale(i) + offset(i)
   10 continue
      j = nobj
      if (ncon .le. 0) go to 30
      do 20 i = 1,ncon
        j = j + 1
        temp1(j) = g(i)
   20 continue
   30 continue
c
      do 70 n = 1,ndv
        do 40 i = 1,nobj
          temp2(i) = df(i,n) / fscale(i)
   40   continue
        j = nobj
        if (ncon .le. 0) go to 60
        do 50 i = 1,ncon
          j = j + 1
          temp2(j) = dg(i,n)
   50   continue
   60   continue
c
      call ksd (dfun(n),temp1,temp2,j,rho)
c
   70 continue
c
      return
      end
