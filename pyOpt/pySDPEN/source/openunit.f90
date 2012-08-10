  subroutine openunit(unitnum,filename,filestatus,fileaction,ierror)

  implicit none

  integer unitnum
  character*(*) filename
  character*(*) filestatus
  character*(*) fileaction
  integer ierror

  open(unit=unitnum,file=filename,status=filestatus,&
      access=fileaction,iostat=ierror)

  return
  end


  subroutine pyflush(unitnum)

  implicit none

  integer unitnum

  call flush(unitnum)

  return
  end
