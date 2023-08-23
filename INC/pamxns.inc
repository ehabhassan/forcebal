!PAMXNS-parameters for plasma species
!  mxns-max ion species  (-)
!  mxmh-max impurity species (-)
!  mxmi-max e+i+im species (-)
!  mxsp-max e+i+im species counting all charge states (-)
      INTEGER        mxns,                    mxmh,
     #               mxmi,                    mxsp 
      PARAMETER      (mxns=6,                 mxmh=2,
     #               mxmi=1+mxns+mxmh,        mxsp=40)
     
