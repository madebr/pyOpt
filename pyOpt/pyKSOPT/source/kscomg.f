      subroutine kscomg (work)
      implicit double precision (a-h,o-z)
      common /kscomm/ rdf   ,adf   ,fdl   ,fdm   ,rho   ,drho  ,rhomax,
     1                fun0  ,slope ,delx  ,alpha ,alpmax,a1    ,a2    ,
     2                a3    ,a4    ,f1    ,f2    ,f3    ,f4    ,alim  ,
     3                atest ,ftest ,ifscl ,ifoff ,isx   ,isx0  ,isxlb ,
     4                isxub ,iscl  ,ig0   ,idf   ,islp  ,iobj0 ,iy    ,
     5                ip    ,ih    ,ihess ,iside ,isact ,idobj ,idg   ,
     6                itmp1 ,itmp2 ,inext ,jnext ,jsel  ,itcnt ,icntr ,
     7                icnta ,isdflg,isdrst,ifncl ,nunit ,ndv   ,ncon  ,
     8                nobj  ,nside ,nscale,iprnt ,itmax ,igrad ,limit
      dimension work(*)
c
c          copy variables from work array to common block kscomm
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 19 July 1996
c
      rdf    = work(1)
      adf    = work(2)
      fdl    = work(3)
      fdm    = work(4)
      rho    = work(5)
      drho   = work(6)
      rhomax = work(7)
      fun0   = work(8)
      slope  = work(9)
      delx   = work(10)
      alpha  = work(11)
      alpmax = work(12)
      a1     = work(13)
      a2     = work(14)
      a3     = work(15)
      a4     = work(16)
      f1     = work(17)
      f2     = work(18)
      f3     = work(19)
      f4     = work(20)
      alim   = work(21)
      atest  = work(22)
      ftest  = work(23)
      ifscl  = int(work(24))
      ifoff  = int(work(25))
      isx    = int(work(26))
      isx0   = int(work(27))
      isxlb  = int(work(28))
      isxub  = int(work(29))
      iscl   = int(work(30))
      ig0    = int(work(31))
      idf    = int(work(32))
      islp   = int(work(33))
      iobj0  = int(work(34))
      iy     = int(work(35))
      ip     = int(work(36))
      ih     = int(work(37))
      ihess  = int(work(38))
      iside  = int(work(39))
      isact  = int(work(40))
      idobj  = int(work(41))
      idg    = int(work(42))
      itmp1  = int(work(43))
      itmp2  = int(work(44))
      inext  = int(work(45))
      jnext  = int(work(46))
      jsel   = int(work(47))
      itcnt  = int(work(48))
      icntr  = int(work(49))
      icnta  = int(work(50))
      isdflg = int(work(51))
      isdrst = int(work(52))
      ifncl  = int(work(53))
      nunit  = int(work(54))
      ndv    = int(work(55))
      ncon   = int(work(56))
      nobj   = int(work(57))
      nside  = int(work(58))
      nscale = int(work(59))
      iprnt  = int(work(60))
      itmax  = int(work(61))
      igrad  = int(work(62))
      limit  = int(work(63))
c
      return
      end
