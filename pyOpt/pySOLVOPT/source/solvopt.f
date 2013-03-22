      subroutine solvopt(n,x,f,fun,flg,grad,options,flfc,func,flgc,
     1     gradc,B,g,g0,g1,gt,gc,z,x1,xopt,xrec,grec,xx,deltax,idx,
     2     iout,ifile)
c-----------------------------------------------------------------------------
c The subroutine SOLVOPT performs a modified version of Shor's r-algorithm in
c order to find a local minimum resp. maximum of a nonlinear function
c defined on the n-dimensional Euclidean space 
c or 
c a local minimum for a nonlinear constrained problem: 
c min { f(x): g(x) (<)= 0, g(x) in R(m), x in R(n) }.
c Arguments:
c n       is the space dimension (integer*4),
c x       is the n-vector, the coordinates of the starting point
c         at a call to the subroutine and the optimizer at regular return
c         (double precision),
c f       returns the optimum function value
c         (double precision),
c fun     is the entry name of a subroutine which computes the value 
c         of the function <fun> at a point x, should be declared as external
c         in a calling routine,
c         synopsis: fun(x,f)
c grad    is the entry name of a subroutine which computes the gradient 
c         vector of the function <fun> at a point x, should be declared as 
c         external in a calling routine,
c         synopsis: grad(x,g)
c func    is the entry name of a subroutine which computes the MAXIMAL 
c         RESIDIAL!!! (a scalar) for a set of constraints at a point x,
c         should be declared as external in a calling routine,
c         synopsis: func(x,fc)
c gradc   is the entry name of a subroutine which computes the gradient
c         vector for a constraint with the MAXIMAL RESIDUAL at a point x,
c         should be declared as external in a calling routine,
c         synopsis: gradc(x,gc)
c flg,    (logical) is a flag for the use of a subroutine <grad>:
c         .true. means gradients are calculated by the user-supplied routine.
c flfc,   (logical) is a flag for a constrained problem:
c         .true. means the maximal residual for a set of constraints
c         is calculated by <func>.
c flgc,   (logical) is a flag for the use of a subroutine <gradc>:
c         .true. means gradients of the constraints are calculated 
c         by the user-supplied routine.
c options is a vector of optional parameters (double precision):
c     options(1)= H, where sign(H)=-1 resp. sign(H)=+1 means minimize resp. 
c         maximize <fun> (valid only for an unconstrained problem) and 
c         H itself is a factor for the initial trial step size 
c         (options(1)=-1.d0 by default),
c     options(2)= relative error for the argument in terms of the infinity-norm
c         (1.d-4 by default),
c     options(3)= relative error for the function value (1.d-6 by default),
c     options(4)= limit for the number of iterations (1.5d4 by default),
c     options(5)= control of the display of intermediate results and error 
c         resp. warning messages (default value is 0.d0, i.e., no intermediate 
c         output but error and warning messages, see the manual for more),
c     options(6)= maximal admissible residual for a set of constraints
c         (options(6)=1.d-8 by default, see the manual for more),
c    *options(7)= the coefficient of space dilation (2.5d0 by default),
c    *options(8)= lower bound for the stepsize used for the difference
c        approximation of gradients (1.d-11 by default,see the manual for more).
c   (* ... changes should be done with care)
c Returned optional values:
c     options(9),  the number of iterations, if positive,
c         or an abnormal stop code, if negative (see manual for more),
c                -1: allocation error,
c                -2: improper space dimension,
c                -3: <fun> returns an improper value,
c                -4: <grad> returns a zero vector or improper value at the
c                    starting point,
c                -5: <func> returns an improper value,
c                -6: <gradc> returns an improper value,
c                -7: function is unbounded,
c                -8: gradient is zero at the point, 
c                    but stopping criteria are not fulfilled,
c                -9: iterations limit exceeded,
c               -11: Premature stop is possible,
c               -12: Result may not provide the true optimum,
c               -13: Function is flat: result may be inaccurate 
c                    in view of a point.
c               -14: Function is steep: result may be inaccurate 
c                    in view of a function value,
c       options(10), the number of objective function evaluations, and
c       options(11), the number of gradient evaluations.
c       options(12), the number of constraint function evaluations, and
c       options(13), the number of constraint gradient evaluations.
c ____________________________________________________________________________
c
      implicit none
      include 'messages.inc'
      
      logical flg,flgc,flfc, constr, app, appconstr
      logical FsbPnt, FsbPnt1, termflag, stopf
      logical stopping, dispwarn, Reset, ksm,knan,obj
      integer n, kstore, ajp,ajpp,knorms, k, kcheck, numelem
      integer dispdata, ld, mxtc, termx, limxterm, nzero, krerun
      integer warnno, kflat, stepvanish, i,j,ni,ii, kd,kj,kc,ip
      integer iterlimit, kg,k1,k2, kless
      integer m1
      double precision options(13),doptions(13) 
      double precision x(n),f
      double precision nsteps(3), gnorms(10), kk, nx
      double precision ajb,ajs, des, dq,du20,du10,du03
      double precision n_float, cnteps
      double precision low_bound, ZeroGrad, ddx, y
      double precision lowxbound, lowfbound, detfr, detxr, grbnd
      double precision fp,fp1,fc,f1,f2,fm,fopt,frec,fst, fp_rate
      double precision PenCoef, PenCoefNew
      double precision gamma,w,wdef,h1,h,hp
      double precision dx,ng,ngc,nng,ngt,nrmz,ng1,d,dd, laststep
      double precision zero,one,two,three,four,five,six,seven
      double precision eight,nine,ten,hundr
      double precision infty, epsnorm,epsnorm2,powerm12
      double precision B,g,g0,g1,gt,gc,z
      double precision x1,xopt,xrec,grec,xx,deltax
      integer idx
      dimension B(n,n),g(n),g0(n),g1(n),gt(n),gc(n),z(n)
      dimension x1(n),xopt(n),xrec(n),grec(n),xx(n),deltax(n)
      dimension idx(n)
      
      character endwarn*100
      external fun,grad,func,gradc
      
      integer iout
      character*(*) ifile
      
      data 
     1    zero/0.d0/, one/1.d0/, two/2.d0/, three/3.d0/, four/4.d0/,
     2    five/5.d0/, six/6.d0/, seven/7.d0/, eight/8.d0/, nine/9.d0/, 
     3    ten/1.d1/,  hundr/1.d2/, powerm12/1.d-12/,   
     4    infty /1.d100/, epsnorm /1.d-15/,  epsnorm2 /1.d-30/
     
      if (options(5).gt.-one) then
        open(unit=iout,file=ifile(1:len_trim(ifile)),status='unknown')
      endif
     
c Check the dimension:
      if (n.lt.2) then
		if (options(5).gt.-one) then
          write(iout,*) errmes
          write(iout,*) error2
        endif
        options(9)=-one
        goto 999
      endif  
      n_float=dble(n)

c store flags:
      app=.not.flg
      constr=flfc
      appconstr=.not.flgc
c Default values for options:
      call soptions(doptions)
      do i=1,8
            if (options(i).eq.zero) then
               options(i)=doptions(i)
            elseif (i.eq.2.or.i.eq.3.or.i.eq.6) then
               options(i)=dmax1(options(i),powerm12)
               options(i)=dmin1(options(i),one)
               if (i.eq.2)options(i)=dmax1(options(i),options(8)*hundr)
            elseif (i.eq.7) then
               options(7)=dmax1(options(i),1.5d0)
            endif
      enddo
               
c WORKING CONSTANTS AND COUNTERS ----{
                        
      options(10)=zero    !! counter for function calculations 
      options(11)=zero    !! counter for gradient calculations
      options(12)=zero    !! counter for constraint function calculations 
      options(13)=zero    !! counter for constraint gradient calculations
      iterlimit=idint(options(4))
      if (constr) then
        h1=-one           !! NLP: restricted to minimization 
        cnteps=options(6)
      else 
        h1=dsign(one,options(1))  !! Minimize resp. maximize a function
      endif
      k=0                         !! Iteration counter
      wdef=one/options(7)-one     !! Default space transf. coeff.

c Gamma control ---{
      ajb=one+1.d-1/n_float**2    !! Base I
      ajp=20  
      ajpp=ajp                    !! Start value for the power 
      ajs=1.15d0                  !! Base II
      knorms=0
      do i=1,10
       gnorms(i)=zero
      enddo  
c---}
c Display control ---{
      if (options(5).le.zero) then
         dispdata=0
         if (options(5).eq.-one) then
            dispwarn=.false.
         else 
            dispwarn=.true.
         endif
      else 
         dispdata=idnint(options(5))
         dispwarn=.true.
      endif
      ld=dispdata
c---}

c Stepsize control ---{
      dq=5.1d0           !! Step divider (at f_{i+1}>gamma*f_{i})
      du20=two
      du10=1.5d0
      du03=1.05d0        !! Step multipliers (at certain steps made)
      kstore=3
      do i=1,kstore
       nsteps(i)=zero    !! Steps made at the last 'kstore' iterations 
      enddo 
      if (app) then
        des=6.3d0        !! Desired number of steps per 1-D search
      else
        des=3.3d0
      endif
      mxtc=3             !! Number of trial cycles (steep wall detect)
c---}
      termx=0
      limxterm=50        !! Counter and limit for x-criterion
c stepsize for gradient approximation 
      ddx=dmax1(1.d-11,options(8))     

      low_bound=-one+1.d-4     !! Lower bound cosine used to detect a ravine
      ZeroGrad=n_float*1.d-16  !! Lower bound for a gradient norm
      nzero=0                  !! Zero-gradient events counter
c Low bound for the values of variables to take into account
      lowxbound=dmax1(options(2),1.d-3)  
c Lower bound for function values to be considered as making difference
      lowfbound=options(3)**2
      krerun=0                 !! Re-run events counter
      detfr=options(3)*hundr   !! Relative error for f/f_{record}
      detxr=options(2)*ten     !! Relative error for norm(x)/norm(x_{record})
      warnno=0                 !! the number of warn.mess. to end with
      kflat=0                  !! counter for points of flatness
      stepvanish=0             !! counter for vanished steps
      stopf=.false.
c ----}  End of setting constants
c ----}  End of the preamble
c--------------------------------------------------------------------
c COMPUTE THE FUNCTION  ( FIRST TIME ) ----{
      call fun(n,x,f)
c	  write(iout,*) 'f'
c	  write(iout,*) f
      options(10)=options(10)+one
      if (dabs(f).ge.infty) then
         if (dispwarn) then
            write(iout,*) errmes
            write(iout,*) error32
            write(iout,*) error6
         endif   
         options(9)=-three
         goto 999
      endif
      do i=1,n
        xrec(i)=x(i)
      enddo  
      frec=f     !! record point and function value
c Constrained problem   
      if (constr)  then
          kless=0
          fp=f
          call func(n,x,fc)
c		  write(iout,*) 'g'
c		  write(iout,*) fc
          options(12)=options(12)+one 
          if (dabs(fc).ge.infty) then
             if (dispwarn) then
                write(iout,*) errmes
                write(iout,*) error52
                write(iout,*) error6
             endif
             options(9)=-five
             goto 999
          endif   
        PenCoef=one          !! first rough approximation
        if (fc.le.cnteps) then  
         FsbPnt=.true.       !! feasible point  
         fc=zero             
        else
         FsbPnt=.false. 
        endif
        f=f+PenCoef*fc
      endif   
c ----}
c COMPUTE THE GRADIENT ( FIRST TIME ) ----{
      if (app) then
        do i=1,n
         deltax(i)=h1*ddx
        enddo
        obj=.true.
        if (constr) then
           call apprgrdn(n,g,x,fp,fun,deltax,obj)
        else
           call apprgrdn(n,g,x,f,fun,deltax,obj)  
        endif
        options(10)=options(10)+n_float
      else
        call grad(n,x,g)
c	    write(iout,*) 'gf'
c	    write(iout,*) g
        options(11)=options(11)+one
      endif
      ng=zero      
      do i=1,n
         ng=ng+g(i)*g(i)
      enddo
      ng=dsqrt(ng)
      if (ng.ge.infty) then
         if (dispwarn) then
            write(iout,*) errmes
            write(iout,*) error42
            write(iout,*) error6
         endif
         options(9)=-four
         goto 999
      elseif (ng.lt.ZeroGrad) then
         if (dispwarn) then
            write(iout,*) errmes
            write(iout,*) error43
            write(iout,*) error6
         endif
         options(9)=-four
         goto 999
      endif
      if (constr) then
       if (.not.FsbPnt) then
         if (appconstr) then
            do j=1,n
              if (x(j).ge.zero) then
                 deltax(j)=ddx
              else
                 deltax(j)=-ddx
              endif
            enddo
            obj=.false.     
            call apprgrdn(n,gc,x,fc,func,deltax,obj)
         else
            call gradc(n,x,gc)
c	        write(iout,*) 'dg'
c	        write(iout,*) gc
         endif
         ngc=zero      
         do i=1,n
           ngc=ngc+gc(i)*gc(i)
         enddo
         ngc=dsqrt(ngc)
         if (ng.ge.infty) then
            if (dispwarn) then
               write(iout,*) errmes
               write(iout,*) error62
               write(iout,*) error6
            endif
            options(9)=-six
            goto 999
         elseif (ng.lt.ZeroGrad) then
            if (dispwarn) then
               write(iout,*) errmes
               write(iout,*) error63
            endif
            options(9)=-six
            goto 999
         endif
         do i=1,n
           g(i)=g(i)+PenCoef*gc(i)
         enddo
         ng=zero      
         do i=1,n
           ng=ng+g(i)*g(i)
           grec(i)=g(i) 
         enddo
         ng=dsqrt(ng)
       endif
      endif
      do i=1,n
        grec(i)=g(i) 
      enddo
      nng=ng
c ----}
c INITIAL STEPSIZE
      d=zero
      do i=1,n
        if (d.lt.dabs(x(i))) d=dabs(x(i))
      enddo  
      h=h1*dsqrt(options(2))*d                  !! smallest possible stepsize
      if (dabs(options(1)).ne.one) then 
        h=h1*dmax1(dabs(options(1)),dabs(h))    !! user-supplied stepsize
      else  
          h=h1*dmax1(one/dlog(ng+1.1d0),dabs(h)) !! calculated stepsize
      endif

C RESETTING LOOP ----{
      do while (.true.)
        kcheck=0                       !! Set checkpoint counter.
        kg=0                           !! stepsizes stored
        kj=0                           !! ravine jump counter
        do i=1,n
          do j=1,n
            B(i,j)=zero
          enddo
          B(i,i)=one                   !! re-set transf. matrix to identity 
          g1(i)=g(i)
        enddo     
        fst=f
        dx=0 
c ----}    

c MAIN ITERATIONS ----{
   
        if (options(5).gt.zero) then
          write(iout,*) 
     1        'Iteration # ..... Function Value ..... ',
     2        'Step Value ..... Gradient Norm'
        endif
        
        do while (.true.)
          k=k+1
          kcheck=kcheck+1
          laststep=dx
c ADJUST GAMMA --{
           gamma=one+dmax1(ajb**((ajp-kcheck)*n),two*options(3))
           gamma=dmin1 ( gamma,ajs**dmax1(one,dlog10(nng+one)) )
c --}      
       ngt=zero
       ng1=zero
       dd=zero
       do i=1,n
         d=zero
         do j=1,n
            d=d+B(j,i)*g(j)
         enddo
         gt(i)=d
         dd=dd+d*g1(i)
         ngt=ngt+d*d
         ng1=ng1+g1(i)*g1(i)
       enddo
       ngt=dsqrt(ngt)      
       ng1=dsqrt(ng1)
       dd=dd/ngt/ng1      
       
       w=wdef       
c JUMPING OVER A RAVINE ----{      
       if (dd.lt.low_bound) then
        if (kj.eq.2) then
          do i=1,n
           xx(i)=x(i)
          enddo
        endif   
        if (kj.eq.0) kd=4
        kj=kj+1
        w=-.9d0              !! use large coef. of space dilation 
        h=h*two  
        if (kj.gt.2*kd) then
          kd=kd+1
          warnno=1
          endwarn=endwarn1  
          do i=1,n
            if (dabs(x(i)-xx(i)).lt.epsnorm*dabs(x(i))) then
             if (dispwarn)  then
                write(iout,*) wrnmes
                write(iout,*) warn08
             endif
            endif
          enddo
        endif  
       else
        kj=0 
       endif
c ----}
c DILATION ----{      
       nrmz=zero
       do i=1,n
         z(i)=gt(i)-g1(i)
         nrmz=nrmz+z(i)*z(i)
       enddo  
       nrmz=dsqrt(nrmz)
       if (nrmz.gt.epsnorm*ngt) then 
        do i=1,n
         z(i)=z(i)/nrmz               
        enddo 
c New direction in the transformed space: g1=gt+w*(z*gt')*z and
c new inverse matrix: B = B ( I + (1/alpha -1)zz' )
        d = zero
        do i=1,n
          d=d+z(i)*gt(i)
        enddo
        ng1=zero
        d = d*w
        do i=1,n
          dd=zero
          g1(i)=gt(i)+d*z(i)
          ng1=ng1+g1(i)*g1(i)
          do j=1,n
             dd=dd+B(i,j)*z(j)
          enddo 
          dd=w*dd
          do j=1,n
            B(i,j)=B(i,j)+dd*z(j)
          enddo
        enddo      
        ng1=dsqrt(ng1)
       else
        do i=1,n
         z(i)=zero
         g1(i)=gt(i)  
        enddo
        nrmz=zero
       endif
       do i=1,n
           gt(i)=g1(i)/ng1
       enddo
        do i=1,n
          d=zero
            do j=1,n
               d=d+B(i,j)*gt(j)
            enddo  
          g0(i)=d
        enddo
c ----}
c RESETTING ----{
        if (kcheck.gt.1) then
           numelem=0
           do i=1,n
              if (dabs(g(i)).gt.ZeroGrad) then
                 numelem=numelem+1
                 idx(numelem)=i
              endif   
           enddo   
           if (numelem.gt.0) then
              grbnd=epsnorm*dble(numelem**2)
              ii=0
              do i=1,numelem
                 j=idx(i)
                 if (dabs(g1(j)).le.dabs(g(j))*grbnd) ii=ii+1
              enddo   
              if (ii.eq.n .or. nrmz.eq.zero) then
                if (dispwarn) then
                  write(iout,*) wrnmes
                  write(iout,*) warn20
                endif
                if (dabs(fst-f).lt.dabs(f)*1.d-2) then
                   ajp=ajp-10*n
                else
                   ajp=ajpp
                endif
                h=h1*dx/three
                k=k-1 
                exit
              endif
           endif 
        endif
c ----}
c STORE THE CURRENT VALUES AND SET THE COUNTERS FOR 1-D SEARCH 
        do i=1,n
         xopt(i)=x(i)
        enddo 
        fopt=f   
        k1=0
        k2=0
        ksm=.false.
        kc=0
        knan=.false.
        hp=h
        if (constr) Reset=.false.
c 1-D SEARCH ----{ 
        do while (.true.)
         do i=1,n
          x1(i)=x(i)
         enddo 
         f1=f   
         if (constr) then
           FsbPnt1=FsbPnt
           fp1=fp
         endif
c NEW POINT         
         do i=1,n
            x(i)=x(i)+hp*g0(i)
         enddo
           ii=0
           do i=1,n   
            if (dabs(x(i)-x1(i)).lt.dabs(x(i))*epsnorm) ii=ii+1
           enddo    
c FUNCTION VALUE         
         call fun(n,x,f)
c	     write(iout,*) 'f'
c	     write(iout,*) f
         options(10)=options(10)+one 
         if (h1*f.ge.infty) then
            if (dispwarn) then
              write(iout,*) errmes 
              write(iout,*) error5
            endif
            options(9)=-seven
            goto 999
         endif
         if (constr) then
           fp=f
           call func(n,x,fc)
c	       write(iout,*) 'g'
c	       write(iout,*) fc
           options(12)=options(12)+one
           if (dabs(fc).ge.infty) then 
               if (dispwarn) then
                  write(iout,*) errmes
                  write(iout,*) error52
                  write(iout,*) error6
               endif
               options(9)=-five
               goto 999
           endif
           if (fc.le.cnteps) then
              FsbPnt=.true.
              fc=zero
           else
              FsbPnt=.false.
              fp_rate=fp-fp1 
              if (fp_rate.lt.-epsnorm) then
               if (.not.FsbPnt1) then
                d=zero
                do i=1,n
                  d=d+(x(i)-x1(i))**2
                enddo
                d=dsqrt(d)  
                PenCoefNew=-1.5d1*fp_rate/d
                if (PenCoefNew.gt.1.2d0*PenCoef) then
                  PenCoef=PenCoefNew
                  Reset=.true.
                  kless=0 
                  f=f+PenCoef*fc
                  exit
                endif
               endif 
              endif 
           endif
           f=f+PenCoef*fc
         endif
         if (dabs(f).ge.infty) then
             if (dispwarn) then
               write(iout,*) wrnmes
               write(iout,*) error32
             endif
             if (ksm.or.kc.ge.mxtc) then
                options(9)=-three
                goto 999
             else 
                k2=k2+1
                k1=0 
                hp=hp/dq 
                do i=1,n
                 x(i)=x1(i)
                enddo 
                f=f1 
                knan=.true. 
                if (constr) then
                  FsbPnt=FsbPnt1 
                  fp=fp1 
                endif
             endif
c STEP SIZE IS ZERO TO THE EXTENT OF EPSNORM
         elseif (ii.eq.n) then
                stepvanish=stepvanish+1
                if (stepvanish.ge.5) then
                    options(9)=-ten-four
                    if (dispwarn) then 
                       write(iout,*) termwarn1
                       write(iout,*) endwarn4 
                    endif
                    goto 999
                else 
                    do i=1,n
                     x(i)=x1(i)
                    enddo 
                    f=f1 
                    hp=hp*ten 
                    ksm=.true.
                    if (constr) then
                       FsbPnt=FsbPnt1 
                       fp=fp1 
                    endif
                endif
c USE SMALLER STEP
         elseif (h1*f.lt.h1*gamma**idint(dsign(one,f1))*f1) then
             if (ksm) exit  
             k2=k2+1
             k1=0 
             hp=hp/dq
             do i=1,n
              x(i)=x1(i)
             enddo 
             f=f1 
             if (constr) then
                FsbPnt=FsbPnt1 
                fp=fp1 
             endif
             if (kc.ge.mxtc) exit
c 1-D OPTIMIZER IS LEFT BEHIND
         else   
             if (h1*f.le.h1*f1) exit
c USE LARGER STEP
             k1=k1+1 
             if (k2.gt.0) kc=kc+1 
             k2=0
             if (k1.ge.20) then 
                 hp=du20*hp 
             elseif (k1.ge.10) then 
                 hp=du10*hp
             elseif (k1.ge.3) then
                 hp=du03*hp
             endif
         endif
        enddo
c ----}  End of 1-D search
c ADJUST THE TRIAL STEP SIZE ----{
        dx=zero
        do i=1,n
           dx=dx+(xopt(i)-x(i))**2
        enddo
        dx=dsqrt(dx)   
        if (kg.lt.kstore)  kg=kg+1
        if (kg.ge.2)  then
           do i=kg,2,-1
             nsteps(i)=nsteps(i-1)
           enddo
        endif
        d=zero   
        do i=1,n
           d=d+g0(i)*g0(i)
        enddo
        d=dsqrt(d)
        nsteps(1)=dx/(dabs(h)*d)
        kk=zero
        d=zero
        do i=1,kg
           dd=dble(kg-i+1) 
           d=d+dd
           kk=kk+nsteps(i)*dd
        enddo
        kk=kk/d   
        if     (kk.gt.des) then
             if (kg.eq.1) then
                h=h*(kk-des+one)
             else   
                h=h*dsqrt(kk-des+one) 
             endif
        elseif (kk.lt.des) then
             h=h*dsqrt(kk/des)  
        endif

        if (ksm) stepvanish=stepvanish+1
c ----}
c COMPUTE THE GRADIENT ----{
        if (app) then
          do j=1,n
            if (g0(j).ge.zero) then
               deltax(j)=h1*ddx
            else
               deltax(j)=-h1*ddx
            endif
          enddo
          obj=.true.     
          if (constr)  then
             call apprgrdn(n,g,x,fp,fun,deltax,obj)
          else 
             call apprgrdn(n,g,x,f,fun,deltax,obj)
          endif
          options(10)=options(10)+n_float
        else
          call grad(n,x,g)
c	      write(iout,*) 'df'
c	      write(iout,*) g
          options(11)=options(11)+one
        endif
        ng=zero
        do i=1,n
          ng=ng+g(i)*g(i)
        enddo
        ng=dsqrt(ng)
        if (ng.ge.infty) then
         if (dispwarn) then
           write(iout,*) errmes
           write(iout,*) error42
         endif
         options(9)=-four
         goto 999
        elseif (ng.lt.ZeroGrad) then
         if (dispwarn) then
           write(iout,*) wrnmes
           write(iout,*) warn1
         endif
         ng=ZeroGrad
        endif
c Constraints:      
        if (constr) then
         if (.not.FsbPnt) then
           if (ng.lt.1.d-2*PenCoef) then
              kless=kless+1
              if (kless.ge.20) then
                 PenCoef=PenCoef/ten
                 Reset=.true.
                 kless=0
              endif   
           else
              kless=0
           endif      
           if (appconstr) then
                 do j=1,n
                   if (x(j).ge.zero) then
                      deltax(j)=ddx
                   else
                      deltax(j)=-ddx
                   endif
                 enddo
                 obj=.false.     
                 call apprgrdn(n,gc,x,fc,func,deltax,obj)
                 options(12)=options(12)+n_float 
           else
                 call gradc(n,x,gc)
c			     write(iout,*) 'dg'
c	             write(iout,*) gc
                 options(13)=options(13)+one 
           endif
           ngc=zero
           do i=1,n
              ngc=ngc+gc(i)*gc(i)
           enddo
           ngc=dsqrt(ngc)
           if (ngc.ge.infty) then
                  if (dispwarn) then
                     write(iout,*) errmes
                     write(iout,*) error62
                  endif
                  options(9)=-six
                  goto 999
           elseif (ngc.lt.ZeroGrad .and. .not.appconstr) then 
                  if (dispwarn) then
                     write(iout,*) errmes
                     write(iout,*) error63
                  endif
                  options(9)=-six 
                  goto 999
           endif
           do i=1,n
             g(i)=g(i)+PenCoef*gc(i) 
           enddo
           ng=zero
           do i=1,n
              ng=ng+g(i)*g(i)
           enddo
           ng=dsqrt(ng)
           if (Reset) then
              if (dispwarn) then
                 write(iout,*) wrnmes
                 write(iout,*) warn21 
              endif
              h=h1*dx/three
              k=k-1
              nng=ng
              exit
           endif
         endif 
        endif
        if (h1*f.gt.h1*frec) then
          frec=f 
          do i=1,n
            xrec(i)=x(i)
            grec(i)=g(i)
          enddo 
        endif
c ----}
       if (ng.gt.ZeroGrad) then
        if (knorms.lt.10)  knorms=knorms+1
        if (knorms.ge.2)  then
          do i=knorms,2,-1
           gnorms(i)=gnorms(i-1) 
          enddo
        endif  
        gnorms(1)=ng
        nng=one  
          do i=1,knorms
            nng=nng*gnorms(i)
          enddo
        nng=nng**(one/dble(knorms))
       endif
c Norm X:
       nx=zero
       do i=1,n
        nx=nx+x(i)*x(i)
       enddo
       nx=dsqrt(nx)         
  
c DISPLAY THE CURRENT VALUES ----{
       if (k.eq.ld) then
c         write(iout,*) 
c     1        'Iteration # ..... Function Value ..... ',
c     2        'Step Value ..... Gradient Norm'
         write(iout,'(5x,i5,7x,f13.5,6x,f13.5,7x,f13.5)') k,f,dx,ng
         do i=1,n,6
          m1=min0(n,i+5)
          write(iout,'(3x,a,2x,6e13.5)') '  X =',(x(j),j=1,m1)
         enddo
         ld=k+dispdata
       endif
c----}
c CHECK THE STOPPING CRITERIA ----{
      termflag=.true.
      if (constr) then
        if (.not.FsbPnt) termflag=.false.
      endif
      if(kcheck.le.5.or.kcheck.le.12.and.ng.gt.one)termflag=.false.
      if(kc.ge.mxtc .or. knan)termflag=.false.
c ARGUMENT
       if (termflag) then
           ii=0
           stopping=.true.
           do i=1,n
             if (dabs(x(i)).ge.lowxbound) then
                ii=ii+1
                idx(ii)=i
                if (dabs(xopt(i)-x(i)).gt.options(2)*dabs(x(i))) then
                  stopping=.false.
                endif  
             endif
           enddo
           if (ii.eq.0 .or. stopping)  then
                stopping=.true.
                termx=termx+1
                d=zero
                do i=1,n
                  d=d+(x(i)-xrec(i))**2
                enddo
                d=dsqrt(d)  
c FUNCTION
                if(dabs(f-frec).gt.detfr*dabs(f) .and.
     1             dabs(f-fopt).le.options(3)*dabs(f) .and.
     2             krerun.le.3 .and. .not. constr) then
                   stopping=.false.
                   if (ii.gt.0) then
                    do i=1,ii
                     j=idx(i)
                     if (dabs(xrec(j)-x(j)).gt.detxr*dabs(x(j))) then
                       stopping=.true.
                       exit
                     endif
                    enddo
                   endif     
                   if (stopping) then
                      if (dispwarn) then
                        write(iout,*) wrnmes
                        write(iout,*) warn09
                      endif
                      ng=zero
                      do i=1,n
                       x(i)=xrec(i)
                       g(i)=grec(i)
                       ng=ng+g(i)*g(i)
                      enddo
                      ng=dsqrt(ng)  
                      f=frec
                      krerun=krerun+1
                      h=h1*dmax1(dx,detxr*nx)/dble(krerun)
                      warnno=2
                      endwarn=endwarn2
                      exit
                   else 
                      h=h*ten
                   endif
                elseif(dabs(f-frec).gt.options(3)*dabs(f) .and.
     1             d.lt.options(2)*nx. and. constr) then
                   continue
                elseif  (dabs(f-fopt).le.options(3)*dabs(f) .or.
     1              dabs(f).le.lowfbound .or.
     2              (dabs(f-fopt).le.options(3).and.
     3               termx.ge.limxterm )) then
                  if (stopf) then
                   if (dx.le.laststep) then
                    if (warnno.eq.1 .and. ng.lt.dsqrt(options(3))) then
                       warnno=0 
                    endif
                    if (.not.app) then
                      do i=1,n
                       if (dabs(g(i)).le.epsnorm2) then
                         warnno=3
                         endwarn=endwarn3
                         exit
                       endif
                      enddo
                    endif  
                    if (warnno.ne.0) then
                       options(9)=-dble(warnno)-ten
                       if (dispwarn) then 
                         write(iout,*) termwarn1
                         write(iout,*) endwarn
                         if (app) write(iout,*) appwarn
                       endif    
                    else 
                       options(9)=dble(k)
                       if (dispwarn) write(iout,*) termwarn0
                    endif
                    goto 999
                   endif
                  else
                   stopf=.true.
                  endif 
                elseif (dx.lt.powerm12*dmax1(nx,one) .and. 
     1                  termx.ge.limxterm ) then
                     options(9)=-four-ten
                     if (dispwarn) then 
                       write(iout,*) termwarn1 
                       write(iout,*) endwarn4
                       if (app) write(iout,*) appwarn
                       f=frec
                       do i=1,n
                        x(i)=xrec(i)
                       enddo
                     endif
                     goto 999
                endif
           endif
       endif 
c ITERATIONS LIMIT
            if(k.eq.iterlimit) then
                options(9)=-nine
                if (dispwarn) then
                  write(iout,*) wrnmes
                  write(iout,*) warn4 
                endif
                goto 999
            endif
c ----}
c ZERO GRADIENT ----{
          if (constr) then 
            if (ng.le.ZeroGrad) then
                if (dispwarn) then  
                  write(iout,*) termwarn1
                  write(iout,*) warn1
                endif
                options(9)=-eight
                goto 999
            endif
          else  
            if (ng.le.ZeroGrad) then
             nzero=nzero+1
             if (dispwarn) then
               write(iout,*) wrnmes
               write(iout,*) warn1
             endif
             if (nzero.ge.3) then
               options(9)=-eight
               goto 999
             endif  
             do i=1,n
               g0(i)=-h*g0(i)/two
             enddo  
             do i=1,10
               do j=1,n
                x(j)=x(j)+g0(j)
               enddo                
               call fun(n,x,f)
c			   write(iout,*) 'f'
c			   write(iout,*) f
               options(10)=options(10)+one
               if (dabs(f).ge.infty) then 
                 if (dispwarn) then
                   write(iout,*) errmes
                   write(iout,*) error32
                 endif
                 options(9)=-three
                 goto 999
               endif
               if (app) then
                   do j=1,n
                     if (g0(j).ge.zero) then
                        deltax(j)=h1*ddx
                     else
                        deltax(j)=-h1*ddx
                     endif
                   enddo
                   obj=.true.     
                   call apprgrdn(n,g,x,f,fun,deltax,obj)
                   options(10)=options(10)+n_float
               else
                   call grad(n,x,g)
c	               write(iout,*) 'df'
c	               write(iout,*) g
                   options(11)=options(11)+one
               endif    
               ng=zero
               do j=1,n
                  ng=ng+g(j)*g(j)
               enddo
               ng=dsqrt(ng)
               if (ng.ge.infty) then
                    if (dispwarn) then
                      write(iout,*) errmes
                      write(iout,*) error42
                    endif
                    options(9)=-four
                    goto 999
               endif
               if (ng.gt.ZeroGrad) exit
             enddo
             if (ng.le.ZeroGrad) then
                if (dispwarn) then  
                  write(iout,*) termwarn1
                  write(iout,*) warn1
                endif
                options(9)=-eight
                goto 999
             endif
             h=h1*dx
             exit
            endif
          endif  
c ----}
c FUNCTION IS FLAT AT THE POINT ----{
          if (.not.constr .and. 
     1        dabs(f-fopt).lt.dabs(fopt)*options(3) .and.
     2        kcheck.gt.5  .and. ng.lt.one ) then
          
           ni=0
           do i=1,n
             if (dabs(g(i)).le.epsnorm2) then
               ni=ni+1
               idx(ni)=i
             endif
           enddo    
           if (ni.ge.1 .and. ni.le.n/2 .and. kflat.le.3) then
             kflat=kflat+1
             if (dispwarn) then
                write(iout,*) wrnmes
                write(iout,*) warn31
             endif
             warnno=1
             endwarn=endwarn1
             do i=1,n
               x1(i)=x(i)
             enddo  
             fm=f 
             do i=1,ni
              j=idx(i)
              f2=fm
              y=x(j)
              if (y.eq.zero) then
                x1(j)=one
              elseif (dabs(y).lt.one) then
                x1(j)=dsign(one,y)
              else
                x1(j)=y
              endif
              do ip=1,20
               x1(j)=x1(j)/1.15d0
               call fun(n,x1,f1)
c			   write(iout,*) 'f1'
c			   write(iout,*) f1
               options(10)=options(10)+one
               if (dabs(f1).lt.infty) then
                 if (h1*f1.gt.h1*fm) then
                   y=x1(j)
                   fm=f1
                 elseif (h1*f2.gt.h1*f1) then
                   exit
                 elseif (f2.eq.f1) then
                   x1(j)=x1(j)/1.5d0   
                 endif  
                 f2=f1   
               endif
              enddo
              x1(j)=y
             enddo
             if (h1*fm.gt.h1*f) then
              if (app) then
                do j=1,n
                  deltax(j)=h1*ddx
                enddo
                obj=.true.     
                call apprgrdn(n,gt,x1,fm,fun,deltax,obj)
                options(10)=options(10)+n_float
              else
                call grad(n,x1,gt)
c			    write(iout,*) 'df'
c			    write(iout,*) gt
                options(11)=options(11)+one
              endif
              ngt=zero
              do i=1,n
                ngt=ngt+gt(i)*gt(i)
              enddo  
              if (ngt.gt.epsnorm2 .and. ngt.lt.infty) then  
                if (dispwarn) write(iout,*) warn32
                do i=1,n
                 x(i)=x1(i)
                 g(i)=gt(i)
                enddo  
                ng=ngt
                f=fm
                h=h1*dx/three
                options(3)=options(3)/five
                exit
              endif   !! regular gradient
             endif   !! a better value has been found
           endif   !! function is flat
          endif   !! pre-conditions are fulfilled
c ----}
       enddo   !! iterations
      enddo   !! restart

999   continue
      
      end
