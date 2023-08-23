!COMTR8-internal variables for chord in TRACK
!  x0-x position of step along chord (m)
!  y0-y position of step along chord (m)
!  z0-z position of step along chord (m)
!  gx-dx/ds along chord (-)
!  gy-dy/ds along chord (-)
!  gz-dz/ds along chord (-)
!  d0-distance along chord from beginning (m)
!  tolinv-tolerence for R,Z inversion (m)
      REAL           x0,                      y0,
     #               z0,                      gx,
     #               gy,                      gz,
     #               d0,                      tolinv
      COMMON/COMTR8/ x0,                      y0,
     #               z0,                      gx,
     #               gy,                      gz,
     #               d0,                      tolinv
     
