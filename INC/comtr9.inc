!COMTR9-internal variables for TRACK
!  mxit-max Jacobian iterations for R,Z, inversion (-)
      INTEGER        mxit
      PARAMETER      (mxit=100)
!  dthmax-max theta step (radians)
!  rdbg-record of iterated values of R (m)
!  zdbg-record of iterated values of Z (m)
!  toltra-tolerance for convergence (-)
!  maxstp-max steps per chord segment (-)
!  maxsrf-max surfaces intersected by chord segment (-)
!  kwtra-option diagnostic output (-)
!  iter-number of Jacobian iterations (-)
      REAL           dthmax,                  rdbg,
     #               zdbg,                    toltra
      INTEGER        maxstp,                  maxsrf,
     #               kwtra,                   iter
      COMMON/COMTR9/ dthmax,                  rdbg(mxit),
     #               zdbg(mxit),              toltra,
     #               maxstp,                  maxsrf,
     #               kwtra,                   iter
     
