      subroutine filtersd_wrap(n,m,xx,xl,xu,lm,f,g,
     *  ainf,ubgv,rho,htol,rgtol,maxit,maxgr,dchk,
     *  dtol,iprint,iout,ifile,ifail,nfs,ngs)
c
      implicit double precision (a-h, o-z)
c
      parameter (mbar=5,maxla=1,maxu=0,maxiu=0)
c
      integer n,m,maxit,maxgr,dchk,iprint,iout,ifail,nfs,ngs
      double precision xx(n),xl(n),xu(n),lm(m),g(m)
      double precision x(n+m),bl(n+m),bu(n+m),al(n+m),v(mbar)
      double precision ws(maxu+10*n+7*m+2*n*(m+1)+3*max(m,50)+
     *  (min(n,mbar)+1)*((min(n,mbar)+1)+1)/2+(min(n,mbar)+1)*
     *  (n+5)+min(m+1,n)*(min(m+1,n)+1)/2+3*n+min(m+1,n))
      integer lws(maxiu+1+maxla+n+m+max(m,50)-1+n+min(m+1,n)+n+m)
      character cstype(m)
      character*(*) ifile
c
c  commons
c
      common/defaultc/ainfty,ubd,mlp,mxf
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/statsc/dnorm,h,hJt,hJ,ipeq,k,itn,nft,ngt
      common/ngrc/mxgr
      common/mxm1c/mxm1
c
c  set values
c
      mxm1=min(m+1,n)
      maxa=n*(m+1)
      maxg=min(n,mbar)+1
      kmax=n
      mlp=max(m,50)
      mxf=max(m,50)
c
      mxws_filterSD = maxu+1+maxa+maxa-1+n+m+n+n+m+1+m+n+
     *  n+m+n+m+n+n+m+mlp+mxf+mxf-1
      mxws_glcpd = n+m+n+maxg*(maxg+1)/2+maxg*(kmax+5)
      mxws_denseL = mxm1*(mxm1+1)/2+3*n+mxm1
      mxws=mxws_filterSD+mxws_glcpd+mxws_denseL
      mxlws=maxiu+1+maxla+n+m+mlp-1+n+min(m+1,n)+n+m
c
      do i=maxu+1,maxu+maxa
        ws(maxa+i)=ws(i)
      enddo
      lws(1)=n
c
      do i=1,n
        x(i)=xx(i)
        bl(i)=xl(i)
        bu(i)=xu(i)
      enddo
      do i=n+1,n+m
        bl(i)=-ainf
        bu(i)=0.D0
        al(i)=lm(i)
      enddo
c
      fmin=-ainf
      ainfty=ainf
      ubd=ubgv
      mxgr=maxgr
      nft=nfs
      ngt=ngs
c
c  open output file
c
      if(iprint.gt.0)then
        open(unit=iout,file=ifile(1:len_trim(ifile)),
     *    status='unknown')
      endif
c   
c  check gradients
c
      if(dchk.gt.0)then
        do i=1,n
          al(i)=1.D-6
        enddo
        call checkd(n,m,x,al,ws,lws,maxa,maxla,maxu,maxiu,
     *    mxws,mxlws,dtol,iprint,iout,ifail)
        do i=1,n
          al(i)=0.0
        enddo
      endif
      if(ifail.lt.0) return
c  
c  call filterSD
c
c   10 continue
c
      call filterSD(n,m,x,al,f,fmin,cstype,bl,bu,ws,lws,v,nv,
     *  maxa,maxla,maxu,maxiu,kmax,maxg,rho,htol,rgtol,maxit,
     *  iprint,iout,ifail)
c
      if (iprint.gt.0)then
        write(iout,1) 'ifail =',ifail
      endif
c      if(ifail.eq.4.and.h.gt.ubd)then
c        ubd=11.D-1*h
c        if (iprint.gt.0)then
c          write(iout,2) 'initial x has h > ubd, ubd updated to ',ubd
c        endif
c        goto10
c      endif
c
c  assign output values
c
      nfs=nft
      ngs=ngt
      do i=1,n
        xx(i)=x(i)
      enddo
      do i=n+1,n+m
        lm(i)=al(i)
      enddo
c
c  print final results (h,rgnorm,k,itn,)
c
      call functions(n,m,x,f,g,ws,lws)
      if (ifail.ge.0.and.ifail.le.6)then
        if (iprint.gt.0)then
          write(iout,1) 'number of free variables =',k
          write(iout,1) 'number of function and gradient calls =',
     *      nft,ngt
          write(iout,4) 'x =',(x(i),i=1,n)
          write(iout,4) 'f =',f
          write(iout,4) 'g =',(g(i),i=1,m)
          write(iout,4) 'al =',(al(i),i=1,n+m)
        endif
      endif
c
      return
c  ------------------------------------------------------------------
c  FORMATS
c  ------------------------------------------------------------------
    1 format(A,15I5)
    2 format(A,5E15.7)
    3 format(A/(20I4))
    4 format(A/(5E15.7))
c
      end
