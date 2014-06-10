C     COMMON SCALARS
      character * 6 hptype
      character * 2 innslvr
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,gjacpcoded,
     +        innercall,jaccoded,jacpcoded,hcoded,hccoded,hlcoded,
     +        hlpcoded,firstde,seconde,truehpr,ignoref,skipacc,sclsys,
     +        useustp,safemode

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,jacpcoded,
     +                  hccoded,hlcoded,hlpcoded,fccoded,gjaccoded,
     +                  gjacpcoded,firstde,seconde,truehpr,ignoref,
     +                  skipacc,sclsys,innercall,useustp,innslvr,hptype,
     +                  safemode
      save   /algparam/
