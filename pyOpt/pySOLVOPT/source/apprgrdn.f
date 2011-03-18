      subroutine apprgrdn(n,g,x,f,fun,deltax,obj)
c Subroutine APPRGRDN performs the finite difference approximation 
c of the gradient <g> at a point <x>.
c f      is the calculated function value at a point <x>,
c <fun>  is the name of a subroutine that calculates function values,
c deltax is an array of the relative stepsizes.
c obj    is the flag indicating whether the gradient of the objective
c        function (1) or the constraint function (0) is to be calculated. 
c
      implicit none
      double precision x(*), g(*), f, deltax(*)
      double precision lowbndobj,lowbndcnt,d,y,fi,one,ten,half
      integer n, i, j
      logical center,obj
      external fun
      data lowbndobj /2.d-10/, lowbndcnt /5.d-15/
      data one /1.d0/, ten /1.d1/, half /.5d0/
      do i=1,n
         y=x(i)
         d=dmax1(lowbndcnt,dabs(y))
         d=deltax(i)*d
         if (obj) then
           if (dabs(d).lt.lowbndobj) then
              d=lowbndobj*dsign(one,deltax(i))
              center=.true.
           else
              center=.false.
           endif
         else
           if (dabs(d).lt.lowbndcnt) then
              d=lowbndcnt*dsign(one,deltax(i))
           endif
         endif     
         x(i)=y+d
         call fun(x,fi)
         if (obj) then
          if (fi.eq.f) then
           do j=1,3
              d=d*ten
              x(i)=y+d
              call fun(x,fi)
              if (fi.ne.f) exit
           enddo
          endif 
         endif     
         g(i)=(fi-f)/d
         if (obj) then
           if (center) then
              x(i)=y-d
              call fun(x,fi)
              g(i)=half*(g(i)+(f-fi)/d)
           endif
         endif     
         x(i)=y
      enddo    
      end
