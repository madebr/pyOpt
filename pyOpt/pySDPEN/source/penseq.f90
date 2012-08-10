  subroutine penseq(n,m,x,lb,ub,fob,constr,alfa_stop,nf_max,&
      iprint,iout,ifile,istop,num_funct)

  implicit none

  integer :: n, m
  integer :: num_funct, num_iter, i, qq
  integer :: nf_max, iprint, iout, istop

  real*8 :: f, alfamax, fob
  real*8 :: violiniz, finiz, alfa_stop
  real*8 :: x(n), lb(n), ub(n)
  real*8 :: constr(m)

  real*8, allocatable :: epsiniz(:),eps(:)
  
  character*(*) ifile
  
  common /num/f
  common /calfamax/alfamax

  !---------------------------------------
  ! allocate eps
  !---------------------------------------
  qq = m+1
  
  allocate(eps(qq-1), epsiniz(qq-1))
  
  !---------------------------------------
  !  punto iniziale, valore ottimo
  !---------------------------------------
  num_funct   = 0 
  num_iter    = 0
  
  call fobcon(n,m,x,fob,constr)
  num_funct   = 1
  
  do i = 1,qq-1
    if(max(0.d0,constr(i)) < 1.d-0) then
      eps(i) = 1.d-3
    else
      eps(i) = 1.d-1
    endif
  enddo
  
  epsiniz     = eps
  finiz       = fob
  violiniz    = max(0.d0,maxval(constr))
  
  !---------------------------------------
  ! print initial info 
  !---------------------------------------
  if (iprint.ge.0) then
    
    open(unit=iout,file=ifile,status='unknown')
    
    write(iout,*) '------------------------'
    write(iout,*) '---- Initial values ----'
    write(iout,*) '------------------------'
    write(iout,*) 'fob = ',fob
    write(iout,*)
    do i=1,n
      write(iout,*) 'x(',i,')=',x(i)
    enddo
    write(iout,*)
    do i=1,m
      write(iout,*) 'con(',i,')=',constr(i),' eps(',i,')=',eps(i)
    enddo
    write(iout,*) '------------------------'
    write(iout,*)
    write(iout,*) 'Start the optimizer:'
  
  endif

  !---------------------------------------
  ! 
  !---------------------------------------
  call sdpen(n,m,x,f,lb,ub,alfa_stop,nf_max,num_iter,&
        num_funct,iprint,iout,istop,qq,eps)
  
  !---------------------------------------
  ! 
  !---------------------------------------
  call fobcon(n,m,x,fob,constr)
  num_funct = num_funct + 1
  
  !---------------------------------------
  ! print final info 
  !---------------------------------------
  if (iprint.ge.0) then
  
    write(iout,*) 'Done'
    write(iout,*)
    write(iout,*) '------------------------'
    write(iout,*) '----   Final values ----'
    write(iout,*) '------------------------'
    write(iout,*)
    write(iout,*) 'fob  = ',fob
    write(iout,*)
    do i=1,n
      write(iout,*) 'x(',i,')=',x(i)
    enddo
    write(iout,*)
    do i=1,m
      write(iout,*) 'con(',i,')=',constr(i),' eps(',i,')=',eps(i)
    enddo
    write(iout,*) '------------------------'
    write(iout,*)
    write(iout,1010) n,m,num_funct, &
          finiz,violiniz,fob,max(0.d0,maxval(constr))
    1010 format(1x,2(' & ', i2),1(' & ',i7), ' & ',es10.3,' & ',es8.1,' & ',es10.3,' & ',es8.1,'\\\\') 
    write(iout,*) '------------------------'
  
  endif
  
  return
  end
