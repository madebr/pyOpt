      subroutine ksinfo (iout,rout)
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
      dimension iout(15),rout(25)
c
c          This user called routine returns the current values of
c          several common block variables to satisfy the curiosity
c          of the user.  This routine may be called at any time after
c          routine ksinit is executed.
c
c          author   - Gregory A. Wrenn
c          location - Lockheed Engineering and Sciences Co.
c                     144 Research Drive
c                     Hampton, Va. 23666
c
c          last modification - 30 August 1991
c
      iout(1)  = inext
      iout(2)  = jnext
      iout(3)  = jsel
      iout(4)  = itcnt
      iout(5)  = icntr
      iout(6)  = icnta
      iout(7)  = isdflg
      iout(8)  = isdrst
      iout(9)  = ifncl
      iout(10) = limit
      iout(11) = -1
      iout(12) = -1
      iout(13) = -1
      iout(14) = -1
      iout(15) = -1
c
      rout(1)  = rho
      rout(2)  = drho
      rout(3)  = fun0
      rout(4)  = slope
      rout(5)  = delx
      rout(6)  = alpha
      rout(7)  = alpmax
      rout(8)  = a1
      rout(9)  = a2
      rout(10) = a3
      rout(11) = a4
      rout(12) = f1
      rout(13) = f2
      rout(14) = f3
      rout(15) = f4
      rout(16) = alim
      rout(17) = atest
      rout(18) = ftest
      rout(19) = -1.0
      rout(20) = -1.0
      rout(21) = -1.0
      rout(22) = -1.0
      rout(23) = -1.0
      rout(24) = -1.0
      rout(25) = -1.0
c
      return
      end
