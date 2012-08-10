  subroutine funct(n,m,x,eps,f)

  implicit none

  integer n,m,i
  real*8 x(n),f,fob,fmax
  real*8 constr(m), eps(m)

  call fobcon(n,m,x,fob,constr)

  fmax = 0.d0

  do i = 1,m
    fmax = fmax + (max(0.d0,constr(i))**1.1d0)/eps(i)
  enddo

  f = fob + fmax

  return
  end
