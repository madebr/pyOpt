  subroutine stop(n,alfa_d,istop,alfa_max,nf,&
      ni,fstop,f,alfa_stop,nf_max)

  implicit none

  integer :: n,istop,i,nf,ni,nf_max
  real*8 :: alfa_d(n),alfa_max,fstop(n+1),ffstop,ffm,f,alfa_stop

  istop=0

  alfa_max=alfa_d(1)
  do i=1,n
    if(alfa_d(i).gt.alfa_max) then
      alfa_max=alfa_d(i)
    end if
  end do

  if(ni.ge.(n+1)) then
    ffm=f
    do i=1,n
      ffm=ffm+fstop(i)
    enddo
    ffm=ffm/dfloat((n+1))

    ffstop=(f-ffm)*(f-ffm)
    do i=1,n
      ffstop=ffstop+(fstop(i)-ffm)*(fstop(i)-ffm)
    enddo

    ffstop=dsqrt(ffstop/dfloat(n+1))

  ! if(ffstop.le.alfa_stop) then
  !   istop = 1
  ! end if

  endif

  if(alfa_max.le.alfa_stop) then
    istop = 1
  end if

  if(nf.gt.nf_max) then
    istop = 2
  end if

  return
  end
