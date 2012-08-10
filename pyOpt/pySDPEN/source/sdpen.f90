  subroutine sdpen(n,m,x,f,bl,bu,alfa_stop,nf_max,ni,nf,&
      iprint,iout,istop,qq,eps)

  implicit none

  logical :: cambio_eps
  integer :: n,m,i,j,i_corr,nf,ni,nf_max,qq
  integer :: num_fal,istop
  integer :: iprint,iout,i_corr_fall

  real*8 :: x(n),z(n),d(n)
  real*8 :: alfa_d(n),alfa,alfa_max
  real*8 :: f,fz 
  real*8 :: bl(n),bu(n),eps(qq-1),alfa_stop,maxeps 
  real*8 :: fstop(n+1)

  num_fal=0

  istop = 0

  fstop=0.d0

  ! ---- scelta iniziale dei passi lungo le direzioni --------
  do i=1,n
    alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i))))
    if(iprint.ge.2) then
      write(iout,*) ' alfainiz(',i,')=',alfa_d(i)
    endif
  end do
  ! -----------------------------------------------------------

  !  ---- scelta iniziale delle direzioni ----------------------
  do i=1,n      
    d(i)=1.d0 
  end do
  ! -----------------------------------------------------------  
     
  call funct(n,m,x,eps,f)
  nf=nf+1

  i_corr=1

  fstop(i_corr)=f

  do i=1,n
    z(i)=x(i)
  end do

  if(iprint.ge.2) then
    write(iout,*) ' ----------------------------------'
    write(iout,*) ' finiz =',f
    do i=1,n
      write(iout,*) ' xiniz(',i,')=',x(i)
    enddo
  endif

  !---------------------------   
  !     ciclo principale
  !---------------------------
  do 

    if(iprint.ge.1) then
      write(iout,*) '----------------------------------------------'
      write(iout,100) ni,nf,f,alfa_max
      100 format(' ni=',i4,'  nf=',i5,'   f=',d12.5,'   alfamax=',d12.5)
    ! pause
    endif
    if(iprint.ge.2) then
      do i=1,n
        write(iout,*) ' x(',i,')=',x(i)
      enddo
    endif

    !-------------------------------------
    ! campionamento lungo asse i_corr
    !-------------------------------------
    call linesearchbox(n,m,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                alfa_max,i_corr_fall,iprint,iout,bl,bu,ni,nf,eps)
              
    if(dabs(alfa).ge.1.d-12) then
      x(i_corr) = x(i_corr)+alfa*d(i_corr)
      f=fz
      fstop(i_corr)=f
      
      num_fal=0
      ni=ni+1
    else
      if(i_corr_fall.lt.2) then 
        fstop(i_corr)=fz         
        
        num_fal=num_fal+1
        ni=ni+1
      endif
    end if

    z(i_corr) = x(i_corr)

    if(i_corr.lt.n) then
      i_corr=i_corr+1
    else
      i_corr=1
    end if 

    call stop(n,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max)

    if (istop.ge.1) exit

    !------------------------------------------------
    ! Aggiornamento parametro di smoothing eps
    !------------------------------------------------
    cambio_eps=.false.
    maxeps = maxval(eps)
    do i = 1,qq-1
      if(eps(i) == maxeps) then
        if(eps(i) > 1.0d-2*sqrt(alfa_max)) then
          if(iprint.ge.0) then
            write(iout,*) '**************************************'
            write(iout,*) '*********** updating eps *************'
          endif
          eps(i) =min(1.d-2*eps(i), 1.0d-1*sqrt(alfa_max))
          cambio_eps=.true.
          call funct(n,m,x,eps,f)
        endif
      endif
    enddo
    if(cambio_eps) then 
      do i=1,n 
        alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i)))) 
      enddo
    endif
  enddo

  return
  end
