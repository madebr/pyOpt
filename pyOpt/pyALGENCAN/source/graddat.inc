C     COMMON SCALARS
      logical constrc,gotc,gotj

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),gparc(nmax),
     +        jcval(jcnnzmax),xset(nmax)

C     COMMON BLOCKS
      common /gdata/ xset,g,gparc,c,dpdc,jcval,jcvar,jcsta,jclen,
     +               constrc,gotc,gotj
      save   /gdata/
