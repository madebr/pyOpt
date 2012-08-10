      subroutine kscomp (work)
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
c          copy variables from common block kscomm to work array
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 19 July 1996
c
      work(1)  = rdf
      work(2)  = adf
      work(3)  = fdl
      work(4)  = fdm
      work(5)  = rho
      work(6)  = drho
      work(7)  = rhomax
      work(8)  = fun0
      work(9)  = slope
      work(10) = delx
      work(11) = alpha
      work(12) = alpmax
      work(13) = a1
      work(14) = a2
      work(15) = a3
      work(16) = a4
      work(17) = f1
      work(18) = f2
      work(19) = f3
      work(20) = f4
      work(21) = alim
      work(22) = atest
      work(23) = ftest
      work(24) = float(ifscl)
      work(25) = float(ifoff)
      work(26) = float(isx)
      work(27) = float(isx0)
      work(28) = float(isxlb)
      work(29) = float(isxub)
      work(30) = float(iscl)
      work(31) = float(ig0)
      work(32) = float(idf)
      work(33) = float(islp)
      work(34) = float(iobj0)
      work(35) = float(iy)
      work(36) = float(ip)
      work(37) = float(ih)
      work(38) = float(ihess)
      work(39) = float(iside)
      work(40) = float(isact)
      work(41) = float(idobj)
      work(42) = float(idg)
      work(43) = float(itmp1)
      work(44) = float(itmp2)
      work(45) = float(inext)
      work(46) = float(jnext)
      work(47) = float(jsel)
      work(48) = float(itcnt)
      work(49) = float(icntr)
      work(50) = float(icnta)
      work(51) = float(isdflg)
      work(52) = float(isdrst)
      work(53) = float(ifncl)
      work(54) = float(nunit)
      work(55) = float(ndv)
      work(56) = float(ncon)
      work(57) = float(nobj)
      work(58) = float(nside)
      work(59) = float(nscale)
      work(60) = float(iprnt)
      work(61) = float(itmax)
      work(62) = float(igrad)
      work(63) = float(limit)
c
      return
      end
