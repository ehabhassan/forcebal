!COMTR7-internal variables for TRACK flux surface intersections
!  rhohit-rho of last flux surface hit (-)
!  rhflx-rho at intersection or tangency in step (-)
!  thflx-theta at intersection or tangency in step (radians)
!  rhflx0-rho at beginning of step (-)
!  thflx0-theta at beginning of step (radians)
!  rhflxs-rho at end of attempted step (-)
!  thflxs-theta at end of attempted step (radians)
!  nrho-number of flux surfaces in plasma (-)
!  iflx-internal option switch for setting which point to use (-) 
      REAL           rhohit,                  rhflx,
     #               thflx,                   rhflx0,
     #               thflx0,                  rhflxs,
     #               thflxs
      INTEGER        nrho,                    iflx
      COMMON/COMTR7/ rhohit,                  rhflx,
     #               thflx,                   rhflx0,
     #               thflx0,                  rhflxs,
     #               thflxs,
     #               nrho,                    iflx
     
