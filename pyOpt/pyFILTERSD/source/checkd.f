christen this file checkd.f

c  Copyright (C) 2010 Roger Fletcher

c  Current version 20 January 2011

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

      subroutine checkd(n,m,x,al,ws,lws,maxa,maxla,maxu,maxiu,
     *  mxws,mxlws,tol,iprint,nout,ifail)
      implicit double precision (a-h,o-z)

c  first derivative checking subroutine for use with filterSD code

c  Parameters
c  **********
c  n    set number of variables
c  m    set number of constraints
c  x    set any vector of variables in x(i), i=1,n
c  al   set a vector of differencing increments in al(i), i=1,n
c        (al(i) is the difference interval for x(i) and must be nonzero,
c        say 0.01 times the typical magnitude of x(i))
c  ws   double precision workspace (the amount required for filterSD.f should
c       be plenty)
c  lws  integer workspace (use the amount required for filterSD.f)
c  maxa set as maxa parameter for filterSD
c  maxla set as maxla parameter for filterSD
c  maxu  set as maxu parameter for filterSD
c  maxiu  set as maxiu parameter for filterSD
c  mxws  set the length of ws as provided in the driver
c  mxlws  set the length of lws as provided in the driver
c  tol  tolerace (>0) for reporting inconsistencies in derivatives (eg 1.D-12)
C  iprint  verbosity of printing (0=none, 1=one line per iteration,
C            2=additional text information given)
C  nout    output channel for write(nout,*)ing
c  ifail   returns failure indication as follows
c              -1 = ws not large enough
c              -2 = lws not large enough
c              -3 = derivative inconsistency in constraint/variable

c  Usage
c  *****
c  The user must write subroutines 'functions' and 'gradients' as for filterSD
c  Write a driver program for your problem,  but replace the call of filterSD
c  by a call of checkd (having set differencing increments in al).
c    The program will report any inconsistencies in the derivatives.
c  If the difference quotient estimate lies between the derivatives
c  at x and x+h (h is the perturbation stored in in al) then the
c  derivative is assumed to be correct. Small errors in this
c  comparison may be ignored. If no errors are reported then the
c  call of filter.. may be restored.

      dimension x(*),al(*),ws(*),lws(*)
      if (iprint.gt.0) write(nout,*) 'entering checkd'
      m1=m+1
c  set real storage map for ws
c  first maxu locations are user storage for functions and gradients
c  vectors required by checkd: two slots of length maxa for a(*)
      last1=maxu+1
      next1=last1+maxa
c  slot of length m+1 for f,c at x
      ncx0=next1+maxa
      ncx1=ncx0+1
c  slot of length m+1 for f,c at x + h.e_i
      ncxd0=ncx0+m1
      ncxd1=ncxd0+1
c  total length of ws used is
      kk=ncxd0+m
      if(kk.gt.mxws)then
        if (iprint.gt.0)then
          write(nout,1) 'ws not large enough: kk, mxws =',kk,mxws
        endif
        ifail = -1
        return
      endif

c  set integer storage map for lws
c  first maxiu locations are user storage for functions and gradients
c  storage of length maxla for la(0:*)
      nla1=maxiu+1
c  total storage needed is
      ll=nla1+maxla-1
      if(ll.gt.mxlws)then
        if (iprint.gt.0)then
          write(nout,1) 'lws not large enough: ll, mxlws =',ll,mxlws
        endif
        ifail = -2
        return
      endif


      call functions(n,m,x,ws(ncx0),ws(ncx1),ws,lws)
      call gradients(n,m,x,ws(last1),ws,lws)
c     if (iprint.gt.0)
c       write(nout,4) 'ws_0',(ws(j),j=last1,last1+7)
c     endif 
      do i=1,n
        xi=x(i)
        x(i)=x(i)+al(i)
        call functions(n,m,x,ws(ncxd0),ws(ncxd1),ws,lws,iflag)
        call gradients(n,m,x,ws(next1),ws,lws,iflag)
        do 10 j=0,m
          dfi=(ws(ncxd0+j)-ws(ncx0+j))/al(i)
          a_ij=aij(i,j,ws(last1),lws(nla1))
          ah_ij=aij(i,j,ws(next1),lws(nla1))
          if((dfi.ge.a_ij-tol.and.dfi.le.ah_ij+tol).or.
     *      (dfi.ge.ah_ij-tol.and.dfi.le.a_ij+tol))goto10
          if (iprint.gt.0)then
            write(nout,1) 
     *        'derivative inconsistency in constraint/variable',j,i
            write(nout,*) 'deriv at x, diff quotient, deriv at x+h =',
     *          a_ij,dfi,ah_ij
            write(nout,*) 'c at x+h, c at x',ws(ncxd0+j),ws(ncx0+j)
c           write(nout,1) 'la(0) =',la(0)
c           write(nout,1) 'c/s j pointer =',la(la(0)+j)
c           write(nout,3) 'c/s j indices',
c     *        (la(k),k=la(la(0)+j),la(la(0)+j+1)-1)
c           write(nout,4) 'c/s j entries',
c     *        (a(k),k=la(la(0)+j),la(la(0)+j+1)-1)
          endif
          ifail = -3
          return
   10   continue
        x(i)=xi
      enddo
      if (iprint.gt.0) write(nout,*) 'derivative check successful'
      if (iprint.gt.0) write(nout,*) 'exiting checkd'
    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(20I4))
    4 format(A/(6E15.7))
      return
      end
