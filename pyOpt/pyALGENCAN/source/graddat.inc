C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),gparc(nmax),
     +        jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,gparc,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
      save   /gdata/
