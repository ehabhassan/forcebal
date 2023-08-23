!COMM30-toroidal geometry - MHD (VMEC)
!  fiota-iota-bar, 1/q (-)
!  fiotab,fiotac,fiotad-cubic spline coefficients for fiota (-)
!  rhm3d-radial grid normalized to a0*sqrt(tor flux) (m)
!  xm3d-radial grid normalized to sqrt(tor flux) (-)
!  emn-lambda (magnetic stream function) coefficients (-)
!  emna-lambda coefficients normalized to xm3d**m3d (-)
!  emnb,emnc,emnd-cubic spline coefficients for emna (-)
!  rmn-R coefficients (m)
!  rmna-R coefficients normalized to xm3d**m3d (m)
!  rmnb,rmnc,rmnd-cubic spline coefficients for rmna (m)
!  zmn-Z coefficients (m)
!  zmna-Z coefficients normalized to xm3d**m3d (m)
!  zmnb,zmnc,zmnd-cubic spline coefficients for zmna (m)
!  phitot-total toroidal flux (Wb)
!  km1n0-mode number for m=1,n=0 component (-)
!  nk3d-number of modes for lambda expansion (-)
!  nk3drz-number of modes for R and Z expansions (-)
!  nx3d-number of radial nodes (-)
!  m3d-poloidal mode numbers (-)
!  n3d-toroidal mode numbers (-)
      REAL           fiota,                   fiotab,
     #               fiotac,                  fiotad,
     #               rhm3d,                   xm3d,
     #               emn,                     emna,
     #               emnb,                    emnc,
     #               emnd,                    rmn,
     #               rmna,                    rmnb,
     #               rmnc,                    rmnd,
     #               zmn,                     zmna,
     #               zmnb,                    zmnc,
     #               zmnd,
     #               phitot
      INTEGER        km1n0,                   nk3d,
     #               nk3drz,                  nx3d,
     #               m3d,                     n3d
      COMMON/COMM30/ fiota(mxnx),             fiotab(mxnx),
     #               fiotac(mxnx),            fiotad(mxnx),
     #               rhm3d(mxnx),             xm3d(mxnx),
     #               emn(mxnx,mxkm),          emna(mxnx,mxkm),
     #               emnb(mxnx,mxkm),         emnc(mxnx,mxkm),
     #               emnd(mxnx,mxkm),         rmn(mxnx,mxkm),
     #               rmna(mxnx,mxkm),         rmnb(mxnx,mxkm),
     #               rmnc(mxnx,mxkm),         rmnd(mxnx,mxkm),
     #               zmn(mxnx,mxkm),          zmna(mxnx,mxkm),
     #               zmnb(mxnx,mxkm),         zmnc(mxnx,mxkm),
     #               zmnd(mxnx,mxkm),
     #               phitot,
     #               km1n0,                   nk3d,
     #               nk3drz,                  nx3d,
     #               m3d(mxkm),               n3d(mxkm)
     
