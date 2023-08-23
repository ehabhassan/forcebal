!COMM31-poloidal and toroidal arrays for surface integrals - MHD (VMEC)
!  phfcn-toroidal node (phi) weighting factors (-)
!  thfcn-poloidal node (theta) weighting factors (-)
!  phnorm-toroidal (phi) normalization (-)
!  thnorm-poloidal (theta) normalization (-)
!  cth-cosine function for each mode and node (-)
!  sth-sine function for each mode and node (-)
!  mp3d-number of toroidal (phi) nodes (-)
!  mt3d-number of poloidal (theta) nodes (-)
!  nper-toroidal field periods (0 for axisymmetry) (-)
      REAL           phfcn,                   thfcn,
     #               phnorm,                  thnorm,
     #               cth,                     sth
      INTEGER        mp3d,                    mt3d,
     #               nper
      COMMON/COMM31/ phfcn(mxmp),             thfcn(mxmt),
     #               phnorm,                  thnorm,
     #               cth(mxmt,mxmp,mxkm),     sth(mxmt,mxmp,mxkm),
     #               mp3d,                    mt3d,
     #               nper
     
