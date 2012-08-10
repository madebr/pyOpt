      subroutine kshess (hess,obj,fscale,offset,df,g,dg,rho,
     1                   ndv,ncon,nobj,temp1,temp2,nodim,ncdim)
      implicit double precision (a-h,o-z)
      dimension hess(*),obj(*),fscale(*),offset(*),df(nodim,*)
      dimension g(*),dg(ncdim,*),temp1(*),temp2(2,*)
      data toler /-40.0/
c
c          routine to compute approximate hessian matrix from
c          first order gradient information
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
      kk = 0
      do 140 m = 1,ndv
c
        do 40 i = 1,nobj
          temp2(1,i) = df(i,m) / fscale(i)
   40   continue
        j = nobj
        if (ncon .le. 0) go to 60
        do 50 i = 1,ncon
          j = j + 1
          temp2(1,j) = dg(i,m)
   50   continue
   60   continue
c
        do 130 n = 1,m
c
          do 70 i = 1,nobj
            temp2(2,i) = df(i,n) / fscale(i)
   70     continue
          j = nobj
          if (ncon .le. 0) go to 90
          do 80 i = 1,ncon
            j = j + 1
            temp2(2,j) = dg(i,n)
   80     continue
   90     continue
c
          sum1 = 0.0
          sum2 = 0.0
          sum3 = 0.0
          sum4 = 0.0
          ng = nobj + ncon
          gmax = temp1(1)
          if (ng .lt. 2) go to 110
          do 100 i = 2,ng
            if (temp1(i) .gt. gmax) gmax = temp1(i)
  100     continue
  110     continue
          do 120 i = 1,ng
            d2grad = 0.0
            if (n .eq. m) d2grad = 1.0
            val = rho * (temp1(i) - gmax)
            if (val .lt. toler) go to 120
            sum1 = sum1 + exp(val)
            sum2 = sum2 + exp(val) * temp2(1,i)
            sum3 = sum3 + exp(val) * rho * temp2(2,i)
            sum4 = sum4 + exp(val) * (rho * temp2(1,i) * temp2(2,i)
     1                                    + d2grad)
  120     continue
c
          kk = kk + 1
          if (n .eq. m) hess(kk) = (sum4 / sum1) -
     1                             (sum2 * sum3) / (sum1 * sum1)
          if (n .ne. m) hess(kk) = 0.0
c
  130   continue
  140 continue
c
      return
      end
