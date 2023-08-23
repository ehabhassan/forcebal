!PAMXNC-parameters for radial cells in plasma transport 
!  mxnc-max radial cells (-)
!  mxnc1-max radial cells plus boundary node (-)
!  mx4nc2-max nodes for double pass chords (-)
      INTEGER        mxnc,                    mxnc1,
     #               mx4nc2
      PARAMETER      (mxnc=60,                mxnc1=mxnc+1,
     #               mx4nc2=4*mxnc+8)
      
