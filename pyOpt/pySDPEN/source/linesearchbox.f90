  subroutine linesearchbox(n,m,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                    alfa_max,i_corr_fall,iprint,iout,bl,bu,ni,nf,eps)

  implicit none

  integer :: n,m,i_corr,nf
  integer :: i,j
  integer :: ni,num_fal
  integer :: iprint,iout,i_corr_fall
  integer :: ifront,ielle
  real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n)
  real*8 :: f,alfa,alfa_max,alfaex, fz,gamma
  real*8 :: delta,delta1,fpar,fzdelta
  real*8 :: eps(m)

  gamma=1.d-6

  delta =0.5d0
  delta1 =0.5d0

  i_corr_fall=0

  ifront=0

  ! indice della direzione corrente

  j=i_corr

  if(iprint.ge.1) then
    write(iout,*) ' j =',j,'    d(j) =',d(j)
  endif

  if(dabs(alfa_d(j)).le.1.d-3*dmin1(1.d0,alfa_max)) then
    alfa=0.d0
    if(iprint.ge.1) then
      write(iout,*) '  alfa piccolo'
      write(iout,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
    endif
    return
  endif

  do ielle=1,2
    if(d(j).gt.0.d0) then
      if((alfa_d(j)-(bu(j)-x(j))).lt.(-1.d-6)) then                 
        alfa=dmax1(1.d-24,alfa_d(j))
      else
        alfa=bu(j)-x(j)
        ifront=1
      endif
    else
      if((alfa_d(j)-(x(j)-bl(j))).lt.(-1.d-6)) then
        alfa=dmax1(1.d-24,alfa_d(j))
      else
        alfa=x(j)-bl(j)
        ifront=1
      endif
    endif
    
    if(dabs(alfa).le.1.d-3*dmin1(1.d0,alfa_max)) then
      d(j)=-d(j)
      i_corr_fall=i_corr_fall+1
      alfa=0.d0
      ifront=0
    
      if(iprint.ge.1) then
        write(iout,*) ' direzione opposta per alfa piccolo'
        write(iout,*) ' j =',j,'    d(j) =',d(j)
        write(iout,*) ' alfa=',alfa,'    alfamax=',alfa_max
      endif
    
      cycle
    
    endif
    
    alfaex=alfa
    
    z(j) = x(j)+alfa*d(j)
   
    call funct(n,m,z,eps,fz)
    nf=nf+1

    if(iprint.ge.1) then
      write(iout,*) ' fz =',fz,'   alfa =',alfa
    endif
    if(iprint.ge.2) then
      do i=1,n
        write(iout,*) ' z(',i,')=',z(i)
      enddo
    endif

    fpar= f-gamma*alfa*alfa

    if(fz.lt.fpar) then

      ! espansione

      do
      
        if((ifront.eq.1).or.(num_fal.gt.n-1)) then
  !     if((ifront.eq.1)) then
          alfa_d(j)=delta*alfa
          return
        end if

        if(d(j).gt.0.d0) then 
          if((alfa/delta1-(bu(j)-x(j))).lt.(-1.d-6)) then
            alfaex=alfa/delta1
          else
            alfaex=bu(j)-x(j)
            ifront=1
            if(iprint.ge.1) then
              write(iout,*) ' punto espan. sulla front.'
            endif
          end if
        else
          if((alfa/delta1-(x(j)-bl(j))).lt.(-1.d-6)) then
            alfaex=alfa/delta1
          else
            alfaex=x(j)-bl(j)
            ifront=1
            if(iprint.ge.1) then
              write(iout,*) ' punto espan. sulla front.'
            endif
          end if
        endif
             
        z(j) = x(j)+alfaex*d(j)    
    
        call funct(n,m,z,eps,fzdelta)
        nf=nf+1

        if(iprint.ge.1) then
          write(iout,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
        endif
        if(iprint.ge.2) then
          do i=1,n
            write(iout,*) ' z(',i,')=',z(i)
          enddo
        endif

        fpar= f-gamma*alfaex*alfaex

        if(fzdelta.lt.fpar) then
          fz=fzdelta
          alfa=alfaex
        else               
          alfa_d(j)=delta*alfa
          return
        end if
      enddo
    else 
      d(j)=-d(j)
      ifront=0

      if(iprint.ge.1) then
        write(iout,*) ' direzione opposta'
        write(iout,*) ' j =',j,'    d(j) =',d(j)
      endif
    endif
  enddo

  if(i_corr_fall.eq.2) then
     alfa_d(j)=alfa_d(j)
  else
     alfa_d(j)=delta*alfa_d(j)
  end if

  alfa=0.d0

  return         
  end
