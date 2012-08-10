      subroutine ksdfp (x,side,act,dfun,ndv,s,slope,y,p,h,hess,
     1                  isdflg)
      implicit double precision (a-h,o-z)
      dimension x(*),side(*),act(*),dfun(*),s(*),y(*),p(*),h(*)
      dimension hess(*)
c
c          routine to compute search direction using the
c          davidon-fletcher-powell method
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 19 July 1996
c
      if (isdflg .gt. 0) go to 40
c
c          reset approximate hessian matrix to initial condition
c
   10 continue
c
c          set h matrix to hess matrix
c
      isdflg = 0
      nh = ndv * (ndv + 1) / 2
      do 20 i = 1,nh
        h(i) = hess(i)
   20 continue
      do 30 i = 1,ndv
        act(i) = side(i)
   30 continue
      go to 80
c
c          dfp algorithm
c
   40 continue
c
c          check for validity of hessian matrix
c          in the case of side constraints
c
      do 50 i = 1,ndv
        if (side(i) .ne. act(i)) go to 10
   50 continue
c
c          compute vectors p and y
c
      do 60 i = 1,ndv
        y(i) = dfun(i) - y(i)
        p(i) = x(i) - p(i)
   60 continue
c
c          compute update to hessian matrix
c
      call kshmul (h,y,s,ndv)
      call ksvprd (p,y,sigma,ndv)
      call ksvprd (y,s,tau,ndv)
      if (abs(sigma) .lt. 1.0e-08) sigma = 1.0e-08
      if (abs(tau) .lt. 1.0e-08) tau = 1.0e-08
      sigma = 1.0 / sigma
      tau = 1.0 / tau
c
      k = 0
      do 70 i = 1,ndv
      do 70 j = 1,i
        k = k + 1
        h(k) = h(k) + sigma * p(i) * p(j) - tau * s(i) * s(j)
   70 continue
c
   80 continue
      call kshmul (h,dfun,s,ndv)
c
c          negate s-vector
c
      smax = 1.0e-08
      do 90 i = 1,ndv
        s(i) = -s(i)
        si = abs(s(i))
        if (si .gt. smax) smax = si
   90 continue
c
c          normalize search direction to unit maximum
c          and save previous gradients in y-vector
c          and previous design variables in p-vector
c
      do 100 i = 1,ndv
        s(i) = s(i) / smax
        y(i) = dfun(i)
        p(i) = x(i)
  100 continue
c
c          compute slope
c
      call ksvprd (s,dfun,slope,ndv)
      if (slope .gt. 0.0 .and. isdflg .gt. 0) go to 10
      if (slope .gt. 0.0) slope = 0.0
c
      return
      end
