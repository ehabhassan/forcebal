MODULE FLUXAV_MOD
!-------------------------------------------------------------------------------
!FLUXAV_MOD is an F90 module of routines that determine flux surface averages
!  and other quantities from MHD equilibrium solutions given by Psi(R,Z).
!
!References:
!
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routines:
!
!  FLUXAV_LOAD         -loads MHD equilibria in private data
!  FLUXAV              -calculates flux surface averages
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE LINEAR1_MOD
USE SPLINE1_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  FLUXAV_EC,           & !generates a set of points lying on a flux surface
                         !  called from FLUXAV_SURF
  FLUXAV_ELN,          & !equation of line through axis
                         !  called from FLUXAV_EC
  FLUXAV_ENP,          & !approximates next point on flux surface
                         !  called from FLUXAV_EC
  FLUXAV_FIT,          & !calculates bicubic spline coefficients
                         !  called from FLUXAV_SURF
  FLUXAV_NUC,          & !provides the nucleus of FLUXAV_FIT
                         !  called from FLUXAV_FIT
  FLUXAV_EVAL,         & !evaluates bicubic splines
                         !  called from FLUXAV_EC
                         !  called from FLUXAV_SURF
  FLUXAV_SEARCH,       & !searches to find indices of an array that bound a value
                         !  called from FLUXAV_EVAL
  FLUXAV_SPLN,         & !
  FLUXAV_SURF            !generates a set of flux surfaces
                         !  called from FLUXAV_LOAD

!-------------------------------------------------------------------------------
! Private data
!-------------------------------------------------------------------------------
!(x,y)=(R,Z) used interchangeably
!Logical switches
LOGICAL, PRIVATE, SAVE :: &
  l_flip_f               !option to flip sign of pol magnetic flux [logical]

!Constants
INTEGER, PRIVATE, SAVE :: &
  nlim_f,              & !number of points on limiter [-]
  nr_f,                & !number of radial points [-]
  np_f,                & !number of poloidal points [-]
  nx_f,                & !number of x points on psi(x,y) grid [-]
  ny_f                   !number of y points on psi(x,y) grid [-]

!0-D variables
REAL(KIND=rspec), PRIVATE, SAVE :: &
  cur0_f,              & !toroidal plasma current [A]
  r0_f,                & !reference major radius [m]
  bmag_f,              & !B at magnetic axis [T]
  rmag_f,              & !major radius of magnetic axis [m]
  xmaxlim_f,           & !maximum x of limiter points [m]
  xminlim_f,           & !minimum x of limiter points [m]
  ymaxlim_f,           & !maximum y of limiter opints [m]
  yminlim_f,           & !minimum y of limiter points
  zmag_f                 !vertical position of magnetic axis [m]

!Arrays
REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  f_r_f(:),            & !R*B_t on radial grid [m*T]
  ffp_r_f(:),          & !F*dF/dpsi on radial grid [rad*T]
  psi_r_f(:),          & !pol magnetic flux on radial grid [Wb/rad]
  rhop_r_f(:),         & !radial grid uniform in poloidal flux [-]
  q_r_f(:)               !safety factor on radial grid [-]

REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  x_f(:),              & !x grid for psi(x,y) [m]
  y_f(:),              & !y grid for psi(x,y) [m]
  psi_xy_f(:,:),       & !pol magnetic flux on (x,y) grid [Wb/rad]
  cspln_xy_f(:,:,:)      !spline coefficients for psi(x,y)

REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  xlim_f(:),           & !x nodes for limiter position [m]
  ylim_f(:)              !y nodes for limiter position [m]

REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  arclen_pr_f(:,:),    & !distance along surface from start [m]
  bp_pr_f(:,:),        & !pol magnetic field on surfaces [T]
  x_pr_f(:,:),         & !x nodes on surfaces, uniform in dl/Bp [m]
  y_pr_f(:,:),         & !y nodes on surfaces, uniform in dl/Bp [m]
  x_sr_f(:,:),         & !x nodes on surfaces, uniform in dl [m]
  y_sr_f(:,:),         & !y nodes on surfaces, uniform in dl [m]
  btot_pr_f(:,:),      & !tot magnetic field on surfaces [T]
  ctheta_pr_f(:,:),    & !Theta on surfaces [-]
  xndgrb_pr_f(:,:),    & !n.grad(B) on surfaces [T/m]
  xnthi_pr_f(:,:)        !|Btot|/Bp on surfaces [-]

REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  dlobp_f(:),          & !integral dl/Bp on radial grid [m/T]
  grth_f(:),           & !<n.grad(theta)> on radial grid [/m]
  phit_f(:),           & !tor magnetic flux on radial grid [Wb]
  rhot_f(:),           & !radial grid prop to sqrt(tor flux) [rho]
  rin_f(:),            & !major radius of inside intersections on radial grid [m]
  rm2_f(:),            & !<1/R**2> on radial grid [/m**2]
  rzmax_f(:),          & !major radius at zmax on surfaces [m]
  rzmin_f(:),          & !major radius at zmin on surfaces [m]
  vol_f(:),            & !plasma volume inside surfaces [m**3]
  zmax_f(:),           & !max Z on surfaces [m]
  zmin_f(:)              !min Z on surfaces [m]

!Constants
REAL(KIND=rspec), PRIVATE, PARAMETER :: &
  z_mu0=1.2566e-06_rspec,        & !permeability of free space [H/m]
  z_pi=3.141592654_rspec           !pi [-]

REAL(KIND=rspec), PRIVATE, PARAMETER :: &
  one=1.0_rspec,       & !REAL 1
  zero=0.0_rspec         !REAL 0

!-------------------------------------------------------------------------------
! Procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE FLUXAV_LOAD(cur0,r0,rmag,zmag,psimag,psilim, &
                       nx_xy,ny_xy,x_xy,y_xy,psi_xy,f_x,ffp_x,rhop_x,q_x, &
                       n_lim,x_lim,y_lim, &
                       iflag,message)
!-------------------------------------------------------------------------------
!FLUXAV_LOAD loads an MHD equilibrium of Psi(R,Z) form into private data and
!  constructs other dependent information
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  nx_xy,               & !number of x points on psi(x,y) grid [-]
  ny_xy,               & !number of y points on psi(x,y) grid [-]
  n_lim                  !number of points on limiter [-]

REAL(KIND=rspec), INTENT(IN) :: &
  cur0,                & !toroidal plasma current [A]
  r0,                  & !reference major radius [m]
  rmag,                & !major radius of magnetic axis [m]
  zmag,                & !vertical position of magnetic axis [m]
  psimag,              & !poloidal flux at magnetic axis/(2*pi) [Wb/rad]
  psilim,              & !poloidal flux at limiter/(2*pi) [Wb/rad]
  x_xy(:),             & !vertical grid for 2-D poloidal flux [m]
  y_xy(:),             & !horizontal grid for 2-D poloidal flux [m]
  psi_xy(:,:),         & !poloidal flux/(2*pi) on x-y grid [Wb/rad]
  x_lim(:),            & !radial positions of limiter points [m]
  y_lim(:),            & !vertical positions of limiter points [m]
  f_x(:),              & !R*B_t on equilibrium grid [m*T]
  ffp_x(:),            & !F*dF/dpsi on equilibrium grid [rad*T]
  rhop_x(:),           & !equilibrium grid uniform in poloidal flux [-]
  q_x(:)                 !safety factor on equilibrium grid [-]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,nset1,nset2

!-------------------------------------------------------------------------------
!Store 0-D information
!-------------------------------------------------------------------------------
cur0_f=cur0
r0_f=r0
rmag_f=rmag
zmag_f=zmag
bmag_f=f_x(1)/rmag_f

!Number of flux surfaces for internal calculations
nr_f=41
np_f=41
nx_f=nx_xy
ny_f=ny_xy
nlim_f=n_lim

!-------------------------------------------------------------------------------
!Store 1-D information
!-------------------------------------------------------------------------------
!Allocate arrays
IF(.NOT. ALLOCATED(psi_r_f)) THEN

  ALLOCATE(f_r_f(nr_f), &
           ffp_r_f(nr_f), &
           psi_r_f(nr_f), &
           rhop_r_f(nr_f), &
           q_r_f(nr_f))

ENDIF

!Null arrays
f_r_f(:)=zero
ffp_r_f(:)=zero
psi_r_f(:)=zero
rhop_r_f(:)=zero
q_r_f(:)=zero

!1D grid uniform in Psi^0.5 (_r)
psi_r_f(1:nr_f)=psimag+(psilim-psimag)*(/ (i-1,i=1,nr_f) /)**2/(nr_f-1)**2
rhop_r_f(1:nr_f)=(psi_r_f(1:nr_f)-psi_r_f(1))/(psi_r_f(nr_f)-psi_r_f(1))

!F=R*Bt
CALL LINEAR1_INTERP(nx_f,rhop_x,f_x,nr_f,rhop_r_f, &
                    f_r_f,iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='FLUXAV_LOAD(1)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

!F*F'
CALL LINEAR1_INTERP(nx_f,rhop_x,ffp_x,nr_f,rhop_r_f, &
                    ffp_r_f,iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='FLUXAV_LOAD(2)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

!q
CALL LINEAR1_INTERP(nx_f,rhop_x,q_x,nr_f,rhop_r_f, &
                    q_r_f,iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='FLUXAV_LOAD(3)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Store 2-D information
!-------------------------------------------------------------------------------
IF(ALLOCATED(psi_xy_f)) THEN

  nset1=SIZE(x_f,1)
  nset2=SIZE(y_f,1)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nx_f .OR. &
     nset2 /= ny_f) THEN

    !Reallocate
    DEALLOCATE(x_f, &
               y_f, &
               psi_xy_f, &
               cspln_xy_f)

    ALLOCATE(x_f(nx_f), &
             y_f(ny_f), &
             psi_xy_f(nx_f,ny_f), &
             cspln_xy_f(2,nx_f,2*ny_f))

  ENDIF

ELSE

  !Initial allocation
  ALLOCATE(x_f(nx_f), &
           y_f(ny_f), &
           psi_xy_f(nx_f,ny_f), &
           cspln_xy_f(2,nx_f,2*ny_f))

ENDIF

x_f(1:nx_f)=x_xy(1:nx_f)
y_f(1:ny_f)=y_xy(1:ny_f)
psi_xy_f(1:nx_f,1:ny_f)=psi_xy(1:nx_f,1:ny_f)
cspln_xy_f(:,:,:)=zero

!-------------------------------------------------------------------------------
!Store limiter information
!-------------------------------------------------------------------------------
IF(ALLOCATED(xlim_f)) THEN

  nset1=SIZE(xlim_f,1)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nlim_f) THEN

    !Reallocate
    DEALLOCATE(xlim_f, &
               ylim_f)

    ALLOCATE(xlim_f(nlim_f), &
             ylim_f(nlim_f))

  ENDIF

ELSE

  !Initial allocation
  ALLOCATE(xlim_f(nlim_f), &
           ylim_f(nlim_f))

ENDIF

xlim_f(1:nlim_f)=x_lim(1:nlim_f)
ylim_f(1:nlim_f)=y_lim(1:nlim_f)

!-------------------------------------------------------------------------------
!Allocate 1-D flux surface arrays
!-------------------------------------------------------------------------------
IF(ALLOCATED(vol_f)) THEN

  nset1=SIZE(vol_f,1)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nr_f) THEN

    !Reallocate
    DEALLOCATE(dlobp_f, &
               grth_f, &
               phit_f, &
               rhot_f, &
               rin_f, &
               rm2_f, &
               rzmax_f, &
               rzmin_f, &
               vol_f, &
               zmax_f, &
               zmin_f)

    ALLOCATE(dlobp_f(nr_f), &
             grth_f(nr_f), &
             phit_f(nr_f), &
             rhot_f(nr_f), &
             rin_f(nr_f), &
             rm2_f(nr_f), &
             rzmax_f(nr_f), &
             rzmin_f(nr_f), &
             vol_f(nr_f), &
             zmax_f(nr_f), &
             zmin_f(nr_f))

  ENDIF

ELSE

  !Initial allocation
  ALLOCATE(dlobp_f(nr_f), &
           grth_f(nr_f), &
           phit_f(nr_f), &
           rhot_f(nr_f), &
           rin_f(nr_f), &
           rm2_f(nr_f), &
           rzmax_f(nr_f), &
           rzmin_f(nr_f), &
           vol_f(nr_f), &
           zmax_f(nr_f), &
           zmin_f(nr_f))

ENDIF

dlobp_f(:)=zero
grth_f(:)=zero
phit_f(:)=zero
rhot_f(:)=zero
rin_f(:)=zero
rm2_f(:)=zero
rzmax_f(:)=zero
rzmin_f(:)=zero
vol_f(:)=zero
zmax_f(:)=zero
zmin_f(:)=zero

!-------------------------------------------------------------------------------
!Allocate 2-D flux surface arrays
!-------------------------------------------------------------------------------
IF(ALLOCATED(bp_pr_f)) THEN

  nset1=SIZE(bp_pr_f,2)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nr_f) THEN

    !Reallocate
    DEALLOCATE(arclen_pr_f, &
               bp_pr_f, &
               x_pr_f, &
               y_pr_f, &
               x_sr_f, &
               y_sr_f, &
               btot_pr_f, &
               ctheta_pr_f, &
               xndgrb_pr_f, &
               xnthi_pr_f)

    ALLOCATE(arclen_pr_f(np_f,nr_f), &
             bp_pr_f(np_f,nr_f), &
             x_pr_f(np_f,nr_f), &
             y_pr_f(np_f,nr_f), &
             x_sr_f(np_f,nr_f), &
             y_sr_f(np_f,nr_f), &
             btot_pr_f(np_f,nr_f), &
             ctheta_pr_f(np_f,nr_f), &
             xndgrb_pr_f(np_f,nr_f), &
             xnthi_pr_f(np_f,nr_f))

  ENDIF

ELSE

  !Initial allocation
  ALLOCATE(arclen_pr_f(np_f,nr_f), &
           bp_pr_f(np_f,nr_f), &
           x_pr_f(np_f,nr_f), &
           y_pr_f(np_f,nr_f), &
           x_sr_f(np_f,nr_f), &
           y_sr_f(np_f,nr_f), &
           btot_pr_f(np_f,nr_f), &
           ctheta_pr_f(np_f,nr_f), &
           xndgrb_pr_f(np_f,nr_f), &
           xnthi_pr_f(np_f,nr_f))

ENDIF

arclen_pr_f(:,:)=zero
bp_pr_f(:,:)=zero
x_pr_f(:,:)=zero
y_pr_f(:,:)=zero
x_sr_f(:,:)=zero
y_sr_f(:,:)=zero
btot_pr_f(:,:)=zero
ctheta_pr_f(:,:)=zero
xndgrb_pr_f(:,:)=zero
xnthi_pr_f(:,:)=zero

!-------------------------------------------------------------------------------
!Generate a set of flux surfaces and fill the 1-D and 2-D arrays
!-------------------------------------------------------------------------------
CALL FLUXAV_SURF(iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='FLUXAV_LOAD(4)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FLUXAV_LOAD

SUBROUTINE FLUXAV(k_grid,nr_r,rho_r, &
                  a0,b2_r,bm2_r,bpout_r,btout_r,elong_r,triang_r,f_r,fhat_r, &
                  fm_r,ftrap_r,gph_r,gr2bm2_r,grho1_r,grho2_r,grth_r,gth_r, &
                  phit_r,psi_r,q_r,r2_r,rin_r,rm2_r,rout_r,vol_r,vp_r,iflag, &
                  message)
!-------------------------------------------------------------------------------
!FLUXAV calculates flux surface quantities from EFIT equilibria
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  Modified from routines in ONETWO created by H. StJohn
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_grid,              & !rho_r grid definition option [-]
                         !=1 rho_r ~ poloidal flux
                         !=else rho_r ~ sqrt(toroidal flux)
  nr_r                   !number of radial nodes [-]

REAL(KIND=rspec), INTENT(IN) :: &
  rho_r(:)               !radial grid for output [-]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 no warnings or errors
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  a0                     !minor radius, half diameter of boundary flux surface [m]

REAL(KIND=rspec), INTENT(OUT) :: &
  b2_r(:),             & !<B**2> [T**2]
  bm2_r(:),            & !<1/B**2> [/T**2]
  bpout_r(:),          & !poloidal field at rout_r(i) [T]
  btout_r(:),          & !toroidal field at rout_r(i) [T]
  elong_r(:),          & !plasma elongation [-]
  f_r(:),              & !poloidal current=2pi*RB_t/mu0 [A]
  fhat_r(:),           & !RB_t/(dpsi_r/drho) [rho/m]
  fm_r(:,:),           & !geometric factor [-]
  ftrap_r(:),          & !trapped particle fraction [-]
  gph_r(:),            & !<1/R**2>V'/(2*pi)**2 [/rho]
  gr2bm2_r(:),         & !<|grad(rho)|**2/B**2> [rho**2/m**2/T**2]
  grho1_r(:),          & !<|grad(rho)|> [rho/m]
  grho2_r(:),          & !<|grad(rho)|**2> [rho**2/m**2]
  grth_r(:),           & !<n.grad(theta)> [/m]
  gth_r(:),            & !<gtt/sqrt(g)>-theta average [-]
  phit_r(:),           & !toroidal magnetic flux [Wb]
  psi_r(:),            & !poloidal magnetic flux/2pi [Wb/rad]
  q_r(:),              & !safety factor [-]
  r2_r(:),             & !<R**2> [m**2]
  rin_r(:),            & !major radius of intersection of surface on inside [m]
  rm2_r(:),            & !<1/R**2> [/m**2]
  rout_r(:),           & !major radius of intersection of surface on outside [m]
  triang_r(:),         & !average upper/lower triangularity [-]
  vol_r(:),            & !volume enclosed [m**3]
  vp_r(:)                !dV/drho [m**3/rho]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j,k

!Temporary profiles
INTEGER :: &
  k_vopt(3)

REAL(KIND=rspec) :: &
  avbp(nr_f),b2(nr_f),dpsidrho(nr_f),rinterp(nr_f),f(4,nr_f)

REAL(KIND=rspec) :: &
  value(3,nr_r)

!Data on a contour
REAL(KIND=rspec), DIMENSION(np_f) :: &
  hlin,ydum

REAL(KIND=rspec) :: &
  bmax,fmc,fms,ftl,ftu,h2a,ha,hca

REAL(KIND=rspec) :: &
  c1(2),c2(2)

!-------------------------------------------------------------------------------
!Scalar quantities
!-------------------------------------------------------------------------------
a0=(x_pr_f(1,nr_f)-rin_f(nr_f))/2

!-------------------------------------------------------------------------------
!Grid mapping
!-------------------------------------------------------------------------------
k_vopt(1)=1
k_vopt(2:3)=0
f(:,:)=zero
value(:,:)=zero

IF(k_grid == 1) THEN

  !Grid proportional to poloidal flux
  rinterp(:)=rhop_r_f(:)
  dpsidrho(1:nr_f)=ABS(psi_r_f(nr_f))

ELSE

  !Grid proportional to sqrt(toroidal flux)
  rinterp(:)=rhot_f(:)
  dpsidrho(1:nr_f)=ABS(phit_f(nr_f)*rhot_f(1:nr_f)/z_pi/q_r_f(1:nr_f))

ENDIF

!-------------------------------------------------------------------------------
!b2=<B**2>
!-------------------------------------------------------------------------------
f(:,:)=zero

DO j=2,nr_f

  ydum(:)=btot_pr_f(:,j)**2/bp_pr_f(:,j)
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(1)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  b2(j)=c2(2)/dlobp_f(j)
  f(1,j)=b2(j)

ENDDO

f(1,1)=bmag_f**2
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(2)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

b2_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!bm2=<1/B**2>
!-------------------------------------------------------------------------------
f(:,:)=zero

DO j=2,nr_f

  ydum(:)=one/btot_pr_f(:,j)**2/bp_pr_f(:,j)
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(3)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  f(1,j)=c2(2)/dlobp_f(j)

ENDDO

f(1,1)=one/bmag_f**2
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(4)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

bm2_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!bpout=outside poloidal field
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=bp_pr_f(1,:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(5)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

bpout_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!btout=outside toroidal field
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,1)=bmag_f
f(1,2:nr_f)=f_r_f(2:nr_f)/x_pr_f(1,2:nr_f)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(6)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

btout_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!elong=plasma elongation
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,2:nr_f)=(zmax_f(2:nr_f)-zmin_f(2:nr_f))/(x_pr_f(1,2:nr_f)-rin_f(2:nr_f))
f(1,1)=f(1,2)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(7)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

elong_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!f=poloidal current
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=f_r_f(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(8)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

f_r(1:nr_r)=2*z_pi*value(1,1:nr_r)/z_mu0

!-------------------------------------------------------------------------------
!gr2bm2=<|grad(rho)|**2/B**2>
!-------------------------------------------------------------------------------
f(:,:)=zero

DO j=2,nr_f

  ydum(:)=bp_pr_f(:,j)*(x_pr_f(:,j)/btot_pr_f(:,j))**2
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(9)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  f(1,j)=c2(2)/dlobp_f(j)/dpsidrho(j)**2

ENDDO

f(1,nr_f)=f(1,nr_f-1)

CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(10)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

gr2bm2_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!grho1=<|grad(rho)|>
!-------------------------------------------------------------------------------
f(:,:)=zero

DO j=2,nr_f-1

  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),x_pr_f(:,j),2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(11)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  f(1,j)=ABS(c2(2)/dlobp_f(j)/dpsidrho(j))

ENDDO

f(1,1)=f(1,2)
f(1,nr_f)=f(1,nr_f-1)

CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(12)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

grho1_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!grho2=<|grad(rho)|**2>
!-------------------------------------------------------------------------------
f(:,:)=zero

DO j=2,nr_f-1

  ydum(:)=x_pr_f(:,j)*x_pr_f(:,j)*bp_pr_f(:,j)
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(13)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  f(1,j)=c2(2)/dlobp_f(j)/dpsidrho(j)**2

ENDDO

f(1,1)=f(1,2)
f(1,nr_f)=f(1,nr_f-1)

CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(14)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

grho2_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!grth=<n.grad(theta)>
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=grth_f(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(15)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

grth_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!gth=<|grad(rho)|**2/R**2>
!-------------------------------------------------------------------------------
f(:,:)=zero

DO j=2,nr_f

  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),bp_pr_f(:,j),2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(16)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  f(1,j)=c2(2)/(2*z_pi*dpsidrho(j))

ENDDO

CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(17)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

gth_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!phit=toroidal flux
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=phit_f(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(18)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

phit_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!psi=poloidal flux
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=psi_r_f(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(19)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

psi_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!q=safety factor
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=q_r_f(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(20)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

q_r(1:nr_r)=value(1,1:nr_r)

!-----------------------------------------------------------------------------
!r2=<R**2>:
!-----------------------------------------------------------------------------
f(:,:)=zero

DO j=2,nr_f-1

  ydum(:)=x_pr_f(:,j)**2/bp_pr_f(:,j)
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

    IF(iflag /= 0) THEN

      iflag=1
      message='FLUXAV(21)/'//message
      IF(iflag == 1) GOTO 9999

    ENDIF

  f(1,j)=c2(2)/dlobp_f(j)

ENDDO

f(1,1)=rmag_f**2
f(1,nr_f)=f(1,nr_f-1)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(22)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

r2_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!rin=inside major radii 
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=rin_f(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(23)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

rin_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!rm2=<1/R**2> 
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=rm2_f(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(24)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

rm2_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!rout=outside major radii 
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,1)=rmag_f
f(1,2:nr_f)=x_pr_f(1,2:nr_f)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(25)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

rout_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!triang=average upper/lower triangularity
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,2:nr_f)=(x_pr_f(1,2:nr_f)+rin_f(2:nr_f)-rzmax_f(2:nr_f)-rzmin_f(2:nr_f))/2/a0
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(26)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

triang_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!vol=volume
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=vol_f(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(27)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

vol_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!vp=dV/drho
!-------------------------------------------------------------------------------
f(:,:)=zero
f(1,:)=ABS(2*z_pi*dlobp_f(:))*dpsidrho(:)
CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(28)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

vp_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!fm=poloidal expansion factors Fm(1...3), <n dot grad(small theta)>
!-------------------------------------------------------------------------------
DO j=2,nr_f !Over flux surfaces

  ydum(:)=btot_pr_f(:,j)/xnthi_pr_f(:,j)/bp_pr_f(:,j)
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(29)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  avbp(j)=2*z_pi*c2(2)/c1(2)/dlobp_f(j)

ENDDO !Over flux surfaces

DO k=1,3 !Over modes

  f(:,:)=zero

  DO j=2,nr_f !Over flux surfaces

    ydum(:)=SIN(k*ctheta_pr_f(:,j))*xndgrb_pr_f(:,j)/bp_pr_f(:,j)
    c1(1)=arclen_pr_f(1,j)
    c1(2)=arclen_pr_f(np_f,j)
    CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                       c2,iflag,message)

    IF(iflag /= 0) THEN

      iflag=1
      message='FLUXAV(30)/'//message
      IF(iflag == 1) GOTO 9999

    ENDIF

    fms=c2(2)/dlobp_f(j)
    ydum(:)=ydum(:)*grth_f(j)*ABS(btot_pr_f(:,j))
    CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                       c2,iflag,message)

    IF(iflag /= 0) THEN

      iflag=1
      message='FLUXAV(31)/'//message
      IF(iflag == 1) GOTO 9999

    ENDIF

    fms=fms*c2(2)/dlobp_f(j)
    ydum(:)=COS(k*ctheta_pr_f(:,j))*xndgrb_pr_f(:,j)/bp_pr_f(:,j)
    CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                       c2,iflag,message)

    IF(iflag /= 0) THEN

      iflag=1
      message='FLUXAV(32)/'//message
      IF(iflag == 1) GOTO 9999

    ENDIF

    fmc=c2(2)/dlobp_f(j)
    ydum(:)=ydum(:)*grth_f(j)*ABS(btot_pr_f(:,j))
    CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                       c2,iflag,message)

    IF(iflag /= 0) THEN

      iflag=1
      message='FLUXAV(33)/'//message
      IF(iflag == 1) GOTO 9999

    ENDIF

    fmc=fmc*c2(2)/dlobp_f(j)
    f(1,j)=2/(b2(j)*avbp(j))*(fms+fmc)

  ENDDO !Over flux surfaces

  CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                      value,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(34)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  fm_r(k,1:nr_r)=value(1,1:nr_r)

ENDDO !Over modes

!-------------------------------------------------------------------------------
!ftrap=trapped fraction, use approximation of Lin-Liu and Miller
!-------------------------------------------------------------------------------
f(:,:)=zero

DO j=2,nr_f !Over flux surfaces

  !Max total B field on contour
  bmax=MAXVAL(ABS(btot_pr_f(:,j)))

  !h=|B|/Bmax
  DO i=1,np_f

    hlin(i)=ABS(btot_pr_f(i,j))/bmax
    hlin(i)=MIN(hlin(i),one)   

  ENDDO

  !<h**2>
  ydum(:)=(hlin(:)**2)/bp_pr_f(:,j)
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(35)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  h2a=c2(2)/dlobp_f(j)

  !<h>   
  ydum(:)=hlin(:)/bp_pr_f(:,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(36)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  ha=c2(2)/dlobp_f(j)

  !Messy integral of Lin-liu & Miller  
  ydum(:)=(one-(SQRT(one-hlin(:))*(one+hlin(:)/2)))/bp_pr_f(:,j)/(hlin(:)**2)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),ydum,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV(37)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  hca=c2(2)/dlobp_f(j)
  ftl=one-h2a*hca
  ftu=one-(one-(one+ha/2)*SQRT(one-ha))*h2a/(ha*ha)
  f(1,j)=0.75*ftu+0.25*ftl

ENDDO !Over flux surfaces

CALL SPLINE1_INTERP(k_vopt,nr_f,rinterp,f,nr_r,rho_r, &
                    value,iflag,message)

IF(iflag /= 0) THEN

  iflag=1
  message='FLUXAV(38)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

ftrap_r(1:nr_r)=value(1,1:nr_r)

!-------------------------------------------------------------------------------
!gph=<1/R**2>V'/(2*pi)**2
!-------------------------------------------------------------------------------
gph_r(1:nr_r)=rm2_r(1:nr_r)*vp_r(1:nr_r)/(2*z_pi)**2

!-------------------------------------------------------------------------------
!fhat=RB_t/(dpsi_r/drho)
!-------------------------------------------------------------------------------
DO j=2,nr_r-1

  fhat_r(j)=z_mu0*f_r(j)/2/z_pi*(rho_r(j+1)-rho_r(j-1))/(psi_r(j+1)-psi_r(j-1))

ENDDO
fhat_r(1)=zero
fhat_r(nr_r)=fhat_r(nr_r-1)

!Force sign to be correct on fhat
fhat_r(1:nr_r)=SIGN(fhat_r(1:nr_r),bmag_f/cur0_f)

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FLUXAV

SUBROUTINE FLUXAV_EC(k_conv,mxncon,arcl,dpsia, &
                     l_auto,psivl,bperr, &
                     xp,yp,mp,bp,iflag,message)
!-------------------------------------------------------------------------------
!FLUXAV_EC generates a set of points lying on a flux surface (contour)
!
!References:
!  W.A.Houlberg, H.St.John 3/2000
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  Not intended to be a general purpose contouring routine
!  Designed specifically for finding plasma flux surfaces assuming that
!    - The plasma surface is convex (bean shapes may be problematic)
!    - The function to be contoured is monotonically decreasing in the
!      region of interest
!  Uses a bicubic spline representation of psi to get psi values at
!    arbitrary points
!  The spline coefficients must exist (in cspln) before entry
!  The contour must fully encircle (rmag,zmag)
!  Given (rmag,zmag) and a psi value, psivl, generate a contour of
!    ordered points,(xp(i),yp(i)), i = 1,mp, which has (rmag,zmag) as
!    as an interior point
!  The search is limited to a rectangle bounded by (xmin,xmax,ymin,ymax)
!  Define the search box (xmin,xmax,ymin,ymax) by a rectangle, within
!    which the plasma must reside under all circumstances:
!    - If the search box is set to the entire MHD grid, this will
!      normally not work satisfactorily because there are local minima
!      and maxima around the f coils which confuse the contour tracing
!    - Use the rectangle defined by the extremes of the limiters
!  If the number of elements (xp,yp) exceeds the limit (mxncon) bperr is
!    relaxed (up to twice) by the increment dbperr before an error exit
!    - bperr=0.03 is suggested
!  There are three ways to run:
!    a) l_auto=.FALSE. k_conv not used
!       Should be ok for all interior plasma surfaces
!       The value of psi to be traced, psivl, is not altered
!       Either the routine returns with the desired contour or a failure
!         is indicated with iflag when the surface passes outside box
!    b) l_auto=.TRUE. k_conv=0
!       Can be used on the plasma boundary near a separatrix
!       The value of psi to be traced, psivl, may be altered
!       If the contour with the value of psivl could not be traced or
!         passes outside the box the value of psivl will be repeatedly
!         brought closer to the value on the magnetic axis, either by
!         some fraction of dpsia if it is set or some fraction of
!         psiaxis-psilim if it is not
!       The first contour succesfully traced will be returned
!    c) l_auto=.TRUE. k_conv=1
!       Same as b) except that a binary search is done to find the
!         plasma surface
!       For this purpose the plasma surface is defined as the largest
!         closed flux surface that can be found inside the search box
!       The actual search box used is not critical for diverted plasmas
!         but is very important for limited plasmas
!  If iflag<0:
!    - The error condition on bperr has been relaxed before exit
!    - To continue this routine could be called again with a larger arcl
!    - This has not to be a problem if mxncon is large enough to begin
!      with to give a reasonable representation of the contour
!  If the value of psi to be traced leads to an error and l_auto=.TRUE. then
!    dpsia is used to control the amount by which the value of psi is
!    adjusted, otherwise dpsia is not used
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_conv,              & !convergence option when l_auto=.TRUE. [-]
                         !=0 not used
                         !=1 contour is converged to the outermost contour that
                         !   never leaves the search box given by
                         !   xmin,xmax,ymin,ymax

  mxncon                 !max number of points to be returned in vectors xp,yp [-]

REAL(KIND=rspec), INTENT(IN) :: &
  arcl,                & !arc length of step [m]
  dpsia                  !psi grid spacing at edge [Wb/rad]
                         !=psi(penultimate value)-psi(edge value)


!Declaration of input/output variables
LOGICAL, INTENT(INOUT) :: &
  l_auto                 !option for error recovery [-]
                         !=.TRUE.  automatic error recovery
                         !         may be set to 0 in if a failure occurs
                         !=.FALSE. no recovery


REAL(KIND=rspec), INTENT(INOUT) :: &
  bperr,               & !relative change in poloidal b field between (xp(i),yp(i))
                         ! and (xp(i+1),yp(i+1))
  psivl                  !psi value for which contour is desired [Wb/rad]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag,               & !iflag-error and warning flag
                         !=-1 warning
                         !=0 no warnings or errors
                         !=1 error
  mp                     !number of points on surface [-]

REAL(KIND=rspec), INTENT(OUT) :: &
  bp(:),               & !poloidal bfield at xp,yp [T]
  xp(:),               & !array of x values on contour [m]
  yp(:)                  !array of y values on contour [m]

!-------------------------------------------------------------------------------
!Declaration of local variables
LOGICAL :: &
  l_reset_psi

INTEGER :: &
  iflage,iflg,iqpt,isgn,isgnsave,itry,newti,nsteps

INTEGER :: &
  nitr8=10

REAL(KIND=rspec) :: &
  a,arclmin,b,bp1,bp2,bperrsave,bpmintol,cmult,cost,dbperr,derrt, &
  dth,dtharcl,dthmin,dx,dxmin,dxx,dy,dymin,dyy,psi1,psi2,psiin, &
  psimag,psiout,psisave,serr,serrt,sint,step,theta,thend,thnew, &
  thstart,x1,x2,xmult,xn,xns,xsave,y1,y2,ymult,yn,yns,ysave

REAL(KIND=rspec) :: &
  pds(6)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Number of tries to pull in flux surface to avoid limiter
itry=0

psiin=zero
bperrsave=bperr
dbperr=bperr/2
l_reset_psi=.FALSE.

!If min bp drops below bpmintol and l_auto=.TRUE., the contour is searched for
!  for an x point
bpmintol=0.05_rspec

!Set minimum arc length for determining minimum angle
arclmin=0.10*arcl

!Get psi at (rmag,zmag)
CALL FLUXAV_EVAL(x_f,nx_f,y_f,ny_f,cspln_xy_f,rmag_f,zmag_f,1, &
                 pds,iflag,message)

IF(iflag == -1) THEN

  message='FLUXAV_EC(1)/'//message

ELSEIF(iflag == 1) THEN

  !Magnetic axis is not in the (x,y) domain
  message='FLUXAV_EC(2)/ERROR:mag axis not in domain'
  GOTO 9999

ENDIF

psimag=pds(1)

IF(psimag < psivl) THEN

  iflag=1
  message='FLUXAV_EC(3)/ERROR:psimag less that psivl'
  GOTO 9999

ENDIF

10  mp=0
    thstart=1.25*z_pi
    thend=thstart+2*z_pi
    theta=thstart

    !Minimum increment in theta between rays
    !May be increased later, after the first point is found
    dthmin=1.0e-5_rspec
    dth=zero

    !Set step to twice grid spacing
    dxx=x_f(3)-x_f(1)
    dyy=y_f(3)-y_f(1)

    !Set minimum acceptable step
    dxmin=0.05*(x_f(2)-x_f(1))
    dymin=0.05*(y_f(2)-y_f(1))
    dx=dxx
    dy=dyy

    !Set absolute error convergence criteria for Newton's method
    serrt=3.5e-06_rspec
    derrt=0.5e-07_rspec

!Loop over theta from thstart to thend (=thstart+2*pi)
20  theta=theta+dth

    IF(theta > 2*z_pi .AND. &
       thstart /= zero) THEN

      theta=theta-2*z_pi
      thend=thend-2*z_pi

    ENDIF

    theta=MIN(theta,thend)

    IF(theta < thend) THEN

      dxx=dx
      dyy=dy
      iqpt=0

      !Get equation of ray emanating from (rmag_f,zmag)
      CALL FLUXAV_ELN(rmag_f,zmag_f,theta,a,b,iflg,isgn)

      !Now have y=a*x+b (iflg=0) or x=a*y+b (iflg=1)
      !Start search from axis point
      x1=rmag_f
      y1=zmag_f
      psi1=psimag
      cost=COS(theta)
      sint=SIN(theta)

      !Sliding interval search, max width of interval ~1.41*(dx or dy)
      nsteps=0
      cmult=one
      xsave=-1.0e20_rspec

30    IF(iflg == 0) THEN

        !Search in x
        xmult=one
        y2=yminlim_f-dymin

        DO WHILE(y2 < yminlim_f .OR. &
                 y2 > ymaxlim_f)

          !Increment x
          x2=x1+isgn*dxx*xmult*cmult

          !Limit the search
          x2=MAX(x2,xminlim_f)
          x2=MIN(x2,xmaxlim_f)

          IF(x2 == x1) THEN

            !If x values are the same, psivl is not within reach
            l_reset_psi=.TRUE.
            GOTO 9999

          ENDIF

          !Get corresponding value of y
          y2=a*x2+b
          xmult=xmult/2

        ENDDO

      ELSE

        !Search in y
        ymult=1
        x2=xminlim_f-dxmin

        DO WHILE(x2 < xminlim_f .OR. &
                 x2 > xmaxlim_f)

          !Increment y
          y2=y1+isgn*dyy*ymult*cmult

          !Limit the search
          y2=MAX(y2,yminlim_f)
          y2=MIN(y2,ymaxlim_f)

          !If y values are the same, psivl is not within reach
          IF(y2 == y1) THEN

            l_reset_psi=.TRUE.
            GOTO 9999

          ENDIF

          !Get corresponding value of x
          x2=a*y2+b
          ymult=ymult/2

        ENDDO

      ENDIF

      CALL FLUXAV_EVAL(x_f,nx_f,y_f,ny_f,cspln_xy_f,x2,y2,1, &
                       pds,iflag,message)

      IF(iflag == -1) THEN

        message='FLUXAV_EC(4)/'//message

      ELSEIF(iflag == 1) THEN

        !(x2,y2) not in domain
        message='FLUXAV_EC(5)/ERROR:(x2,y2) is not in domain'
        GOTO 9999

      ENDIF

      psi2=pds(1)
      IF((psivl-psi1)*(psivl-psi2) > zero) THEN

        !Since psi is increasing, we are in the vicinity of an x-point
        !To be sure that we stay on the contour which envelopes the
        !  plasma, we first search for the minimum in psi along the ray
        !We guarantee that the point we find is on the contour by not
        !  allowing the ray to extend past this minimum
        !Use Newton's method to search for the minimum -- the point
        !  along the ray where the directional derivative of psi along
        !  the ray is zero

        IF(cmult < 0.99_rspec) THEN

          cmult=cmult/2
          GOTO 30

        ENDIF

        IF(psi2-psi1 > zero) THEN    !psi is increasing

          xsave=x1   ! save current guess for possible later restore
          ysave=y1
          psisave=psi1
          isgnsave=isgn
          xn=x1
          yn=y1
          newti=0
          step=serrt+2

          DO WHILE(ABS(step) > serrt .AND. newti <= nitr8)

            CALL FLUXAV_EVAL(x_f,nx_f,y_f,ny_f,cspln_xy_f,xn,yn,6, &
                             pds,iflag,message)

            IF(iflag == -1) THEN

              message='FLUXAV_EC(6)/'//message

            ELSEIF(iflag == 1) THEN

              !(xn,yn) not in domain        
              message='FLUXAV_EC(7)/ERROR:(xn,yn) not in domain'
              GOTO 9999

            ENDIF

            step=-(pds(2)*cost+pds(3)*sint)/(cost*cost*pds(5) &
                 +2*cost*sint*pds(4)+sint*sint*pds(6))
            xn=xn+step*cost
            yn=yn+step*sint
            newti=newti+1

          ENDDO

          IF(newti <= nitr8) THEN

            !Newton's method converged
            isgn=-isgn
            IF(iflg == 1) dyy=dyy/2
            IF(iflg == 0) dxx=dxx/2
            x2=xn
            y2=yn
            psi2=pds(1)

            !The minimum in psi along the ray must be less than or
            !  equal to the value for which we are searching, psivl
            !If this is not the case then we can't find psivl on the
            !  ray so abandon the search
            IF(psi2 > psivl) THEN

              l_reset_psi=.TRUE.
              GOTO 9999

            ENDIF

          ELSE

            !Set cmult<1
            cmult=0.5_rspec
            x1=xsave
            y1=ysave
            psi1=psisave

            !Restore search direction
            isgn=isgnsave
            GOTO 30

          ENDIF

        ENDIF

        nsteps=nsteps+1
        x1=x2
        y1=y2
        psi1=psi2
        GOTO 30

      ENDIF

      !Now have psivl between psi1 and psi2, converge using Newton-Raphson
      newti=0

      IF(iflg == 0) THEN

        xn=x1+isgn*dxx/2
        yn=a*xn+b

      ELSE

        yn=y1+isgn*dyy/2
        xn=a*yn+b

      ENDIF

40    iflag=0
      CALL FLUXAV_EVAL(x_f,nx_f,y_f,ny_f,cspln_xy_f,xn,yn,3, &
                       pds,iflag,message)

      IF(iflag > 0) THEN

        !(xn,yn) is not in the (x,y) domain
        !Set number of Newton iterations to maximum to restart
        newti=nitr8
        iflag=0

      ELSE

        serr=-(pds(1)-psivl)/(pds(2)*cost+pds(3)*sint)
        IF(ABS(serr) < serrt) GOTO 50

        IF(psivl == zero) THEN

          IF(ABS(pds(1)-psivl) < derrt) GOTO 50

        ELSE

          IF(ABS((pds(1)-psivl)/psivl) < derrt) GOTO 50

        ENDIF

        xn=xn+serr*cost
        yn=yn+serr*sint
        newti=newti+1

      ENDIF

      IF(newti >= nitr8) THEN

        IF(dxx <= dxmin) GOTO 50
        IF(dyy <= dymin) GOTO 50
        IF(iflg == 0) dxx=dxx/2
        IF(iflg == 1) dyy=dyy/2
        IF(iqpt == 0) GOTO 30
        theta=theta-dth
        GOTO 20

      ENDIF

      GOTO 40

      !End of Newton iteration
      !Check for sufficient accuracy in point spacing as determined by theta
      !Accuracy test is based on a relative error in poloidal b field of bperr
50    bp2=SQRT(pds(2)**2+pds(3)**2)/xn

      IF(theta.ne.thstart) THEN

        !Not the first point of the contour
        IF(ABS(bp2-bp1)/MAX(bp2,bp1) >= bperr) THEN

          !Relative change in bp too large
          IF(dth.ne.dthmin) THEN

            !Spacing too large for grad psi, decrease theta and retry
            dth=dth/2

            !dth lower limit must be observed
            dth=MAX(dth,dthmin)

          ENDIF

        ENDIF

      ENDIF

      !Acceptable point found, collect it and set up for next point
      mp=mp+1
      IF(mp == mxncon) THEN

        !Ran out of storage for the contour points
        !If mxncon is reasonable then bperr may be too small
        IF(bperr-bperrsave >= 2*dbperr) THEN

          !Increasing bperr has been tried and failed
          iflag=-2
          message='FLUXAV_EC(8)/WARNING:bperr increased and failed'
          GOTO 9999

        ELSE

          !Increase bperr by dbperr and retry to generate the contour
          bperr=bperr+dbperr
          GOTO 10

        ENDIF

      ENDIF

      !Found new point - save it
      xp(mp)=xn
      yp(mp)=yn
      bp(mp)=bp2

      !Set up for next ray
      IF(mp == 1) THEN

        !Estimate the radius and relate the theta increment to arc length
        dthmin=arclmin/SQRT((xn-rmag_f)**2+(yn-zmag_f)**2)
        dtharcl=arcl/arclmin*dthmin

      ENDIF

      IF(ABS(bp2-bp1)/MAX(bp2,bp1) <= 0.25*bperr) dth=2*dth

      !Set bp1 for error test on next ray
      bp1=bp2

      !Don't let dth get out of range
      dth=MIN(dth,dtharcl)
      dth=MAX(dth,dthmin)

      !For contours sufficiently far removed from rmag_f,zmag, use the
      !  following approach:
      !  - An approximation to the new point,(xns,yns) is found by moving
      !    along the tangent line (after 2nd step), for a distance arcl
      !  - If near an x point, as signaled by bp < bpmintol, then switch
      !    to the ray method with small dx,dy increments to avoid crossing
      !    the x point
      IF(nsteps > 2 .AND. &
         bp(mp) > bpmintol) THEN

        !Use tangent line method
        iflage=0
        CALL FLUXAV_ENP(pds,arcl,xn,yn,theta, &
                        sint,cost,xns,yns,thnew,iflage)

        !If this failed, go back and use the ray method
        IF(iflage > 0) GOTO 20

        dth=thnew-theta
        !Do not use thnew here
        theta=theta+dth

        IF(theta > 2*z_pi .AND. &
           thstart /= zero) THEN

          theta=theta-2*z_pi
          thend=thend-2*z_pi

        ENDIF

        theta=MIN(theta,thend)

        IF(theta < thend) THEN

          sint=SIN(theta)
          cost=COS(theta)
          xn=xns
          yn=yns
          newti=0
          iqpt=1

          !Skip directly to Newton's method
          GOTO 40

        ENDIF

      ELSE

        !Use ray method
        GOTO 20

      ENDIF

    ENDIF

    !Normal exit
    mp=mp+1
    xp(mp)=xp(1)
    yp(mp)=yp(1)
    bp(mp)=bp(1)
    bperr=bperrsave

    IF(itry.ne.0 .AND. &
       k_conv == 1) THEN

        psiin=psivl
        IF(ABS((psiout+psiin)/2-psivl) < 1.0e-6_rspec) GOTO 9999
        psivl=(psiout+psiin)/2
        GOTO 10

    ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

!Check for open contour with possible recovery
IF(l_reset_psi) THEN

  SELECT CASE (l_auto)

  CASE (.TRUE.)

    psiout=psivl

    !Try to correct improper value of psivl
    IF(dpsia > zero) THEN

      psivl=psivl+dpsia/FLOAT(99)

    ELSE

      psivl=psivl+0.0005*(psimag-psivl)

    ENDIF

    IF((itry > 0) .AND. &
       (psiin /= zero) .AND. &
       (k_conv == 1)) THEN

      psivl=(psiin+psiout)/2

    ENDIF

    itry=itry+1

    IF(itry <= 100) THEN

      !Go back and try again
      l_reset_psi=.FALSE.
      GOTO 10

    ELSE

      !Allow error exit to be taken
      iflag=1
      message='FLUXAV_EC(10)/ERROR:no recovery'

    ENDIF

  CASE (.FALSE.)

    !No recovery
    iflag=1
    message='FLUXAV_EC(9)/ERROR:no recovery'

  END SELECT

ENDIF

END SUBROUTINE FLUXAV_EC

SUBROUTINE FLUXAV_ELN(xaxis,yaxis,theta,a,b,keqxy,keqdir)
!-------------------------------------------------------------------------------
!FLUXAV_ELN finds the equation of a straight line passing through
!  (xaxis,yaxis) and with slope whose angle is theta
!
!References:
!  W.A.Houlberg, H.St.John 3/1998
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  theta,               & !poloidal angle defining line direction [rad]
  xaxis,               & !starting x coordinate, nominally the magnetic axis [m]
  yaxis                  !starting y coordinate, nominally the magnetic axis [m]

!Declaration of output variables
INTEGER, INTENT(OUT) :: &
  keqdir,              & !switch indicating direction to increment abscissa to
                         !move along the line from (xaxis,yaxis) to the boundary
                         !=+1 x or y increases from axis to surface
                         !=-1 x or y decreases from axis to surface
  keqxy                  !switch indicating form of equation for line
                         !=0 y=a*x+b
                         !=1 x=a*y+b

REAL(KIND=rspec), INTENT(OUT) :: &
  a,                   & !slope of line [-]
  b                      !intercept of line [m]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec) :: &
  theta1

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Get equation of ray emanating from (xaxis,yaxis)
IF(((z_pi/4 <= theta) .AND. (theta <= 3*z_pi/4)) .OR. &
   ((5*z_pi/4 <= theta) .AND. (theta <= 7*z_pi/4))) THEN

  !Express x as a function of y, x=a*y+b
  keqxy=1

  IF(theta > z_pi) THEN

    !Bottom of plasma
    keqdir=-1
    theta1=3*z_pi/2-theta
    IF(theta > 3*z_pi/2) theta1=z_pi-ABS(theta1)

  ELSE

    !Top of plasma
    keqdir=1
    theta1=z_pi/2-theta
    IF(theta > z_pi/2) theta1=2*z_pi-ABS(theta1)

  ENDIF

  a=TAN(theta1)
  b=xaxis-a*yaxis

ELSE

  !Express y as a function of x, y=a*x+b
  keqxy=0

  IF((theta < z_pi/4) .OR. &
     (theta > 7*z_pi/4)) THEN

    !Outside of plasma
    keqdir=1

  ELSE

    !Inside of plasma
    keqdir=-1

  ENDIF

  a=TAN(theta)
  b=yaxis-a*xaxis

ENDIF

END SUBROUTINE FLUXAV_ELN

SUBROUTINE FLUXAV_ENP(pds,arcl,xn,yn,theta,sint,cost, &
                      xns,yns,thnew,iflag)
!-------------------------------------------------------------------------------
!FLUXAV_ENP makes an approximation to get the next point on a flux
!  surface contour given the current point and the gradient in psi
!
!References:
!  W.A.Houlberg, H.St.John 3/1998
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The tangent line at the point (xn,yn) is df = 0 = (df/dx)*dx+(df/dy)*dy
!  A circle of radius arcl, centered at (xn,yn) is arcl**2 = dr**2+dz**2
!  This gives us two equations in two unkowns: dr and dz .
!  The new point is xn +/- dx and yn +/- dy where the signs have to
!    be picked so that theta increases.
!  Note the special treatment required as thnew crosses zero (outside
!    this routine)
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  arcl,                & !radius of circle at xn,yn [m]
  cost,                & !cos(theta) [-]
  sint,                & !sin(theta) [-]
  theta,               & !reference poloidal angle [rad]
  xn,                  & !x value of reference point [m]
  yn                     !y value of reference point [m]

REAL(KIND=rspec), INTENT(IN) :: &
  pds(:)                 !bicubic spline coefficients for polidal flux [arb]

!Declaration of output variables
INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 no warnings or errors
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  thnew,               & !value of theta at new point [rad]
  xns,                 & !x value of new point [m]
  yns                    !y value of new point [m]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec) :: &
  alpha,dpx,dpy,dx,dy

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Error flag
iflag=0

dpx=pds(2)
dpy=pds(3)

IF(ABS(dpx) > ABS(dpy)) THEN

  IF(dpx == zero) THEN

    iflag=1
    GOTO 9999

  ENDIF

  alpha=dpy/dpx
  dy=arcl/SQRT(one+alpha*alpha)
  dx=-alpha*dy

ELSE

  IF(dpy == zero) THEN

    iflag=1
    GOTO 9999

  ENDIF

  alpha=dpx/dpy
  dx=arcl/SQRT(one+alpha*alpha)
  dy=-alpha*dx

ENDIF

!The sign on dx,dy must be taken so that theta increases
!A unit vector in the direction of increasing theta (i.e., theta
!  counterclockwise) is (-SIN (theta),COS (theta))
!The displacement vector is (dx,dy)
!Its projection on the above vector must be positive and equals
!  -dx*SIN (theta)+dy*COS (theta)
IF(-dx*sint+dy*cost < zero) THEN

  dx=-dx
  dy=-dy

ENDIF

xns=xn+dx
yns=yn+dy

!Make sure the point is inside the boundaries
IF(xns < xminlim_f .OR. &
   xns > xmaxlim_f .OR. &
   yns < yminlim_f .OR. &
   yns > ymaxlim_f) THEN

  iflag=1
  GOTO 9999

ENDIF

thnew=ATAN2(yns-zmag_f,xns-rmag_f)

IF(thnew < zero) THEN

  thnew=thnew+2*z_pi

ELSEIF(theta > 7*z_pi/4) THEN

  thnew=thnew+2*z_pi

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FLUXAV_ENP

SUBROUTINE FLUXAV_FIT(f,x,nx,y,ny,c,iflag,message)
!-------------------------------------------------------------------------------
!FLUXAV_FIT calculates bicubic spline coefficients
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  nx,                  & !number of x values  >=  4 [-]
  ny                     !number of y values  >=  4 [-]

REAL(KIND=rspec), INTENT(IN) :: &
  x(:),                & !array of x values
  y(:),                & !array of y values
  f(:,:)                 !array of values at (x(i),y(j))

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 no warnings or errors
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  c(:,:,:)               !array of spline coefficients
                         !c(1,i,j)=s
                         !c(2,i,j)=ds/dx
                         !c(1,i,j+ny)=ds/dy
                         !c(2,i,j+ny)=d(ds/dx)/dy
                         !where s(x,y) is the spline approximation


!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j

REAL(KIND=rspec) :: &
  wk(ny,2,nx),f2(ny,2*nx),wk2(2*nx,2,ny)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Check for minumum number of x points
IF(nx < 4) THEN

  iflag=1
  message='FLUXAV_FIT(1)/ERROR:nx < 4'
  GOTO 9999

ENDIF

!Check for minumum number of y points
IF(ny < 4) THEN

  iflag=1
  message='FLUXAV_FIT(2)/ERROR:ny < 4'
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Do spline fits
!-------------------------------------------------------------------------------
!x direction
CALL FLUXAV_NUC(x,f,nx,ny,wk,iflag,message)

!Check messages
IF(iflag > 0) THEN

  message='FLUXAV_FIT(3)/'//message
  GOTO 9999

ENDIF

!y direction
DO i=1,nx

  f2(1:ny,2*i-1)=wk(1:ny,1,i)
  f2(1:ny,2*i)=wk(1:ny,2,i)

ENDDO

CALL FLUXAV_NUC(y,f2,ny,2*nx,wk2,iflag,message)

!Check messages
IF(iflag > 0) THEN

  message='FLUXAV_FIT(4)/'//message
  GOTO 9999

ENDIF

DO i=1,nx

  DO j=1,ny

    c(1,i,2*j-1)=wk2(2*i-1,1,j)
    c(2,i,2*j-1)=wk2(2*i,1,j)
    c(1,i,2*j)  =wk2(2*i-1,2,j)
    c(2,i,2*j)  =wk2(2*i,2,j)

  ENDDO

ENDDO

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FLUXAV_FIT

SUBROUTINE FLUXAV_NUC(tau,gtau,n,m,vs,iflag,message)
!-------------------------------------------------------------------------------
!FLUXAV_NUC provides the nucleus of the bicubic spline routine FLUXAV_FIT
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n,                   & !
  m                      !

REAL(KIND=rspec), INTENT(IN) :: &
  tau(:),              & !
  gtau(:,:)              !

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 no warnings or errors
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  vs(:,:,:)              !

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j,jj,jm1,jp1,k

REAL(KIND=rspec) :: &
  aa,bb,c1,c2,cc,dd,dtau,g,h,ratio,u,xilim

REAL(KIND=rspec) :: &
  w(n,2)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
w(2,1)=tau(3)-tau(1)

IF(w(2,1) <= zero) THEN

  iflag=1
  message='FLUXAV_NUC(1)/ERROR:grid not monotonic'
  GOTO 9999

ENDIF

DO k=1,m

  vs(k,1,1)=gtau(1,k)

ENDDO

xilim=tau(1)

IF(n >= 5) THEN

  xilim=tau(n-2)

  DO i=2,n-3

    j=i+1
    w(j,1)=tau(i+2)-tau(j)

    IF(w(j,1) <= zero) THEN

      iflag=1
      message='FLUXAV_NUC(2)/ERROR:grid not monotonic'
      GOTO 9999

    ENDIF

    DO k=1,m

      vs(k,1,i)=gtau(j,k)

    ENDDO

  ENDDO

ENDIF

w(n-2,1)=tau(n)-xilim

IF(w(n-2,1) <= zero) THEN

  iflag=1
  message='FLUXAV_NUC(3)/ERROR:grid not monotonic'
  GOTO 9999

ENDIF

DO k=1,m

  vs(k,1,n-2)=gtau(n,k)

ENDDO

DO i=2,n-2

  DO k=1,m

    vs(k,2,i)=(vs(k,1,i)-vs(k,1,i-1))/w(i,1)

  ENDDO

ENDDO

dtau=tau(2)-tau(1)
ratio=dtau/w(2,1)
w(1,2)=(ratio-one)**2
w(1,1)=ratio*(ratio-one)
c1=ratio*(2*ratio-3.0_rspec)

DO k=1,m

  vs(k,2,1)=(gtau(2,k)-gtau(1,k))/dtau+vs(k,2,2)*c1

ENDDO

IF(n >= 5) THEN

  DO i=2,n-3

    j=i+1
    jj=i-1
    g=-w(j,1)/w(jj,2)
    c1=3*w(i,1)
    c2=3*w(j,1)

    DO k=1,m

      vs(k,2,i)=g*vs(k,2,jj)+c1*vs(k,2,j)+c2*vs(k,2,i)

    ENDDO

    w(i,2)=g*w(jj,1)+2*(w(i,1)+w(j,1))

  ENDDO

ENDIF

dtau=tau(n-1)-xilim
ratio=dtau/w(n-2,1)
g=-(ratio-one)**2/w(n-3,2)
w(n-2,2)=ratio*(ratio-one)
c1=ratio*(2*ratio-3.0_rspec)

DO k=1,m

  vs(k,2,n-2)=(gtau(n-1,k)-vs(k,1,n-3))/dtau+vs(k,2,n-2)*c1

ENDDO

w(n-2,2)=g*w(n-3,1)+w(n-2,2)

DO k=1,m

  vs(k,2,n-2)=(g*vs(k,2,n-3)+vs(k,2,n-2))/w(n-2,2)

ENDDO

DO j=n-3,1,-1

  DO k=1,m

    vs(k,2,j)=(vs(k,2,j)-w(j,1)*vs(k,2,j+1))/w(j,2)

  ENDDO

ENDDO

DO k=1,m

  DO jj=1,n

    j=n+1-jj

    IF(j == 1) THEN

      jm1=j

    ELSEIF(j == n) THEN

      jm1=j-2

    ELSE

      jm1=j-1

    ENDIF

    DO i=1,2

      vs(k,i,j)=vs(k,i,jm1)

    ENDDO

  ENDDO

  DO j=2,n-1,n-3

    jm1=j-1
    jp1=j+1
    IF(jm1 == 2) jm1=1
    IF(jp1 == n-1) jp1=n
    h=tau(jp1)-tau(jm1)
    u=tau(j)-tau(jm1)
    aa=vs(k,1,jm1)
    bb=vs(k,2,jm1)
    cc=(3*(vs(k,1,jp1)-vs(k,1,jm1))/h-(vs(k,2,jp1)+2*vs(k,2,jm1)))/h
    dd=(2*(vs(k,1,jm1)-vs(k,1,jp1))/h+(vs(k,2,jp1)+vs(k,2,jm1)))/h**2
    vs(k,1,j)=aa+u*(bb+u*(cc+dd*u))
    vs(k,2,j)=bb+u*(2*cc+3*dd*u)

  ENDDO

ENDDO

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FLUXAV_NUC

SUBROUTINE FLUXAV_EVAL(x,nx,y,ny,c,xl,yl,npds, &
                       pds,iflag,message)
!-------------------------------------------------------------------------------
!FLUXAV_EVAL evaluates bicubic splines
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  npds,                & !number of values and partial derivatives [-]
  nx,                  & !number of x elements [-]
  ny                     !number of y elements [-]

REAL(KIND=rspec), INTENT(IN) :: &
  xl,                  & !target x value [arb]
  yl                     !target y value [arb]

REAL(KIND=rspec), INTENT(IN) :: &
  c(:,:,:),            & !array of bicubic spline coefficients [arb]
  x(:),                & !array of x ordinates [arb]
  y(:)                   !array of y ordinates [arb]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag
                         !=-1 warning
                         !=0 no warnings or errors
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  pds(:)                 !array of values and partial derivatives [arb]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j,k,km1,kp1,kp2,lxl,lx,ly,l,lx1

REAL(KIND=rspec) :: &
  dy,hx,u,v

REAL(KIND=rspec) :: &
  su(2),sv(2),sux(2),suy(2),svx(2),sxy(2)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
lx=0
ly=0
iflag=0

!-------------------------------------------------------------------------------
!Find interval and evaluate
!-------------------------------------------------------------------------------
!Find lower bound
!Correlated table search for xl
CALL FLUXAV_SEARCH(x,nx,xl,lx)

IF(lx == 0) THEN

  iflag=1
  message='FLUXAV_EVAL(1)/ERROR:x<x(1)'
  GOTO 9999

ELSEIF(lx == nx) THEN

  iflag=1
  message='FLUXAV_EVAL(2)/ERROR:x>x(n)'
  GOTO 9999

ENDIF

!Correlated table search for yl
CALL FLUXAV_SEARCH(y,ny,yl,ly)

IF(ly == 0) THEN

  iflag=1
  message='FLUXAV_EVAL(3)/ERROR:y<y(1)'
  GOTO 9999

ELSEIF(ly == ny) THEN

  iflag=1
  message='FLUXAV_EVAL(4)/ERROR:y>y(n)'
  GOTO 9999

ENDIF

lx1=lx+1
hx=x(lx1)-x(lx)
dy=y(ly+1)-y(ly)
u=(xl-x(lx))/hx
v=(yl-y(ly))/dy
k=2*ly
kp1=k+1
kp2=k+2
km1=k-1

DO l=1,2

  lxl=lx-1+l
  i=2*(ly-1+l)
  j=i-1
  sv(l)=FLUXAV_SPLN(0,c(1,lxl,km1),c(1,lxl,kp1),c(1,lxl,k),c(1,lxl,kp2),dy,v)
  svx(l)=FLUXAV_SPLN(0,c(2,lxl,km1),c(2,lxl,kp1),c(2,lxl,k),c(2,lxl,kp2),dy,v)

  IF(npds >= 3) THEN

    su(l)=FLUXAV_SPLN(0,c(1,lx,j),c(1,lx1,j),c(2,lx,j),c(2,lx1,j),hx,u)
    suy(l)=FLUXAV_SPLN(0,c(1,lx,i),c(1,lx1,i),c(2,lx,i),c(2,lx1,i),hx,u)

    IF(npds >= 4) THEN

      sux(l)=FLUXAV_SPLN(1,c(1,lx,j),c(1,lx1,j),c(2,lx,j),c(2,lx1,j),hx,u)
      sxy(l)=FLUXAV_SPLN(1,c(1,lx,i),c(1,lx1,i),c(2,lx,i),c(2,lx1,i),hx,u)

    ENDIF

  ENDIF

ENDDO

pds(1)=FLUXAV_SPLN(0,sv(1),sv(2),svx(1),svx(2),hx,u)

IF(npds > 1) THEN

  pds(2)=FLUXAV_SPLN(1,sv(1),sv(2),svx(1),svx(2),hx,u)

  IF(npds > 2) THEN

    pds(3)=FLUXAV_SPLN(1,su(1),su(2),suy(1),suy(2),dy,v)

    IF(npds > 3) THEN

      pds(4)=FLUXAV_SPLN(1,sux(1),sux(2),sxy(1),sxy(2),dy,v)

      IF(npds > 4) THEN

        pds(5)=FLUXAV_SPLN(2,sv(1),sv(2),svx(1),svx(2),hx,u)

        IF(npds > 5) THEN

          pds(6)=FLUXAV_SPLN(2,su(1),su(2),suy(1),suy(2),dy,v)

        ENDIF

      ENDIF

    ENDIF

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FLUXAV_EVAL

SUBROUTINE FLUXAV_SEARCH(x,n,xl,jlo)
!-------------------------------------------------------------------------------
!FLUXAV_SEARCH is a correlated table search routine to find the indices of the
!  array x(j) that bound xl
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  This is a modified version of the Numerical Recipes routine HUNT
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n                      !number of elements in x

REAL(KIND=rspec), INTENT(IN) :: &
  xl,                  & !target value
  x(:)                   !monotonically increasing array

!Declaration of input/output variables
INTEGER, INTENT(INOUT) :: &
  jlo                    !input starting lower index
                         !<1     binary search
                         !=1,n-1 use value
                         !>n-1   binary search
                         !output starting lower index
                         !=0     xl     <  x(1) 
                         !=1     x(1)   <= xl   <= x(2)
                         !=2,n-1 x(jlo) <  xl   <= x(jlo+1)
                         !=n     x(jlo) >  x(n)

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  inc,jhi,jmid

!-------------------------------------------------------------------------------
!Search
!-------------------------------------------------------------------------------
!Check if xl below domain
IF(xl < x(1)) THEN

  !xl is out of range, below x(1)
  jlo=0                     

!Check if xl is above domain
ELSEIF(xl <= x(2)) THEN

  !x(1) <= xl <= x(2)
  jlo=1

ELSEIF(xl <= x(n)) THEN

  !Search for x(2) < xl <= x(n)
  !Check if jlo from previous call is usable
  IF((jlo < 1) .OR. &
     (jlo > n-1)) THEN

    jlo=2
    jhi=n

  ELSE

    !Bracket xl
    inc=1

    IF(xl > x(jlo)) THEN

      !Search up
      jhi=jlo+1

      DO WHILE(xl > x(jhi))

        inc=inc+inc
        jlo=jhi
        jhi=MIN0(jlo+inc,n)

      ENDDO

    ELSE

      !Search down
      jhi=jlo
      jlo=jlo-1

      DO WHILE(xl <= x(jlo))

        inc=inc+inc
        jhi=jlo
        jlo=MAX0(jlo-inc,1)

      ENDDO

    ENDIF

  ENDIF

  !Bisection
  DO WHILE(jhi-jlo > 1)

    jmid=(jhi+jlo)/2

    IF(xl > x(jmid)) THEN

      jlo=jmid

    ELSE

      jhi=jmid

    ENDIF

  ENDDO

ELSE

  !xl > x(n)
  jlo=n

ENDIF

END SUBROUTINE FLUXAV_SEARCH

FUNCTION FLUXAV_SPLN(k,s0,sh,sp0,sph,h,d)
!-------------------------------------------------------------------------------
!FLUXAV_SPLN calculates partial values of a bi-cubic spline
!
!References:
!  W.A.Houlberg, 8/2006
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k

REAL(KIND=rspec), INTENT(IN) :: &
  s0,sh,sp0,sph,h,d

!Declaration of output variables
REAL(KIND=rspec) :: &
  FLUXAV_SPLN

IF(k == 0) THEN

  FLUXAV_SPLN=s0+d*(h*sp0+d*(3*(sh-s0)-(sph+2*sp0)*h+d*(2*(s0-sh) &
              +(sph+sp0)*h)))

ELSEIF(k ==1) THEN

  FLUXAV_SPLN=sp0+d*(6*(sh-s0)/h-2*(sph+2*sp0)+3*d*(2*(s0-sh)/h+(sph+sp0)))

ELSEIF(k == 2) THEN

  FLUXAV_SPLN=6*(sh-s0)/h**2-2*(sph+2*sp0)/h+d*(2*(s0-sh)/h**2+(sph+sp0)/h)*6

ENDIF

END FUNCTION FLUXAV_SPLN

SUBROUTINE FLUXAV_SURF(iflag,message)
!-------------------------------------------------------------------------------
!FLUXAV_SURF calculates a set of flux surfaces from Psi(R,Z) data
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  Modified from routines in ONETWO created by H. StJohn
!-------------------------------------------------------------------------------

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 no warnings or errors
                         !=1 error

!-------------------------------------------------------------------------------
!Declaration of local variables
LOGICAL :: &
  l_auto                 !option for error recovery [-]
                         !=.TRUE.  automatic error recovery
                         !         may be set to 0 in if a failure occurs
                         !=.FALSE. no recovery

INTEGER :: &
  i,ii,iold,j,kflarc,k_conv,k_vopt(3),mp

REAL(KIND=rspec) :: &
  psifctr

!Temporary profiles
REAL(KIND=rspec) :: &
  psivlcpy(nr_f),value(3,np_f)

REAL(KIND=rspec) :: &
  pds(6),psil_xy(nx_f,ny_f)

!Number of points on a surface
INTEGER, PARAMETER :: &
  mxncon=1000

!Data on a contour
REAL(KIND=rspec), DIMENSION(mxncon) :: &
  arclen,bp,bp0,s0,xp,xp0,yp,yp0

REAL(KIND=rspec) :: &
  f(4,mxncon)

REAL(KIND=rspec) :: &
  xfun(np_f)

REAL(KIND=rspec) :: &
  a,arcl,bperr,capbr,capbz,dang,dbdr,dbdz,dbrdr,dbrdz,dbtdr,dbtdz,dbzdr,dbzdz, &
  dpsia,taxis,tlim,ydif,yref,psimtmp,psiltmp

REAL(KIND=rspec) :: &
  c1(2),c2(2)
   
!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Outer boundary at surface
psifctr=zero

!Extremes of limiter positions
xminlim_f=MINVAL(xlim_f(1:nlim_f))
xmaxlim_f=MAXVAL(xlim_f(1:nlim_f))
yminlim_f=MINVAL(ylim_f(1:nlim_f))
ymaxlim_f=MAXVAL(ylim_f(1:nlim_f))

!Check sign of poloidal flux
l_flip_f=.FALSE.
psil_xy(:,:)=psi_xy_f(:,:)

IF(psi_r_f(2) > psi_r_f(1)) THEN

  !Flip signs on 2-D poloidal flux array
  l_flip_f=.TRUE.
  psil_xy(:,:)=-psil_xy(:,:)

ENDIF

!Set up bicubic spline coefficient array
CALL FLUXAV_FIT(psil_xy,x_f,nx_f,y_f,ny_f, &
                cspln_xy_f,iflag,message)

IF(iflag /= 0) THEN

  message='FLUXAV_SURF(1)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

!Get a psimag value consistent with the spline fit
CALL FLUXAV_EVAL(x_f,nx_f,y_f,ny_f,cspln_xy_f,rmag_f,zmag_f,6, &
                 pds,iflag,message)

IF(iflag /= 0) THEN

  message='FLUXAV_SURF(2)/'//message
  IF(iflag == 1) GOTO 9999

ENDIF

!Make the sign consistent with the input value
SELECT CASE (l_flip_f)

CASE (.TRUE.)

  psimtmp=-pds(1)

CASE (.FALSE.)

  psimtmp=pds(1)

END SELECT  

!Get the original limiter value of psi
psiltmp = psi_r_f(nr_f)

!Set spline interpolation options
k_vopt(1)=1
k_vopt(2:3)=0

!Set normalized surface integral
xfun(1:np_f)=REAL((/ (i-1,i=1,np_f) /),rspec)/(np_f-1)

!Set initial values 
arcl=0.02_rspec
taxis=5.0_rspec
tlim=30.0_rspec

!Make copy of flux levels
!The target value of the flux will be overwritten with the actual value
10  SELECT CASE (l_flip_f)

    CASE (.TRUE.)

      psivlcpy(1:nr_f)=-psi_r_f(1:nr_f)

    CASE (.FALSE.)

      psivlcpy(1:nr_f)=psi_r_f(1:nr_f)

    END SELECT

    !Recast psivlcpy using the spline consistent value of psimag=psimtmp
    psivlcpy(1)=SIGN(psimtmp,psivlcpy(2))

    !From the plasma edge to the contour next to the magnetic axis
    !  - find the contour describing the surface (xp(i),yp(i)),i = 1...mp
    !  - form the required flux surface averages and get other surface info

    DO j=nr_f,2,-1 !Over flux surfaces

      kflarc=0

20    IF(j == nr_f) THEN

        l_auto=.TRUE.
        k_conv=1
        psivlcpy(nr_f)=psifctr*(psivlcpy(nr_f-1)-psivlcpy(nr_f))+psivlcpy(nr_f)
        a=(tlim-taxis)/(psivlcpy(nr_f)-psivlcpy(1))
        dpsia=psivlcpy(nr_f-1)-psivlcpy(nr_f)

      ELSE

        l_auto=.FALSE.
        k_conv=0

      ENDIF

      dang=(a*(psivlcpy(j)-psivlcpy(1))+taxis)*(kflarc+1)
      bperr=0.01_rspec

      CALL FLUXAV_EC(k_conv,mxncon,arcl,dpsia, &
                     l_auto,psivlcpy(j),bperr, &
                     xp,yp,mp,bp,iflag,message)

      IF(iflag == 1) THEN

        message='FLUXAV_SURF(3)/'//message
        GOTO 9999

      ENDIF

      IF(iflag == -1) THEN

        !Contour method failed
        iflag=1
        message='FLUXAV_SURF(4)/'//message
        GOTO 9999

      ELSEIF(iflag == 1) THEN

        !Nonrecoverable error in contour routines
        message='FLUXAV_SURF(5)/'//message
        GOTO 9999

      ELSEIF(iflag < -1) THEN

        !Too many points; reduce accuracy requirement and try again
        iflag=0
        kflarc=kflarc+1
        arcl=arcl*1.3
        IF(kflarc < 5) GOTO 20

      ENDIF

!      IF (cur0_f < zero) THEN
!
!        !Change sign of poloidal field
!        bp(1:mp)=-bp(1:mp)
!
!      ENDIF
      !Verify sign of poloidal field is same as current - LRB
      bp(1:mp)=SIGN(bp(1:mp), cur0_f)


      !Re-order contour arrays so that contour begins on the outside
      !Copy arrays
      xp0(1:mp)=xp(1:mp)
      yp0(1:mp)=yp(1:mp)
      bp0(1:mp)=bp(1:mp)

      !Find index of point nearest the plane of the axis on outside
      ii=0
      yref=1.0e6_rspec

      DO i=1,mp

        ydif=ABS(zmag_f-yp0(i))

        IF(ydif <= yref .AND. &
           xp0(i) >= rmag_f) THEN

          ii=i
          yref=ydif

        ENDIF

      ENDDO

      !Reorder arrays from new starting point
      DO i=1,mp-1

        iold=i+ii-1
        IF(iold > mp) iold=iold-mp+1
        xp(i)=xp0(iold)
        yp(i)=yp0(iold)
        bp(i)=bp0(iold)

      ENDDO

      xp(mp)=xp(1)
      yp(mp)=yp(1)
      bp(mp)=bp(1)

      !Compute arclength around the contour, and dl/Bp
      arclen(:)=zero
      s0(:)=zero

      DO i=1,mp-1

        arclen(i+1)=arclen(i)+SQRT((xp(i+1)-xp(i))**2 &
                        +(yp(i+1)-yp(i))**2)
        s0(i+1)=s0(i)+2*(arclen(i+1)-arclen(i))/(bp(i+1)+bp(i))

      ENDDO

      s0(:)=s0(:)/s0(mp)

      !(x,y) and Bp values evenly spaced in dl/Bp
      f(:,:)=zero
      f(1,:)=xp
      CALL SPLINE1_INTERP(k_vopt,mp,s0,f,np_f,xfun, &
                          value,iflag,message, &
                          -1)
      x_pr_f(:,j)=value(1,:)

      f(2:4,:)=zero
      f(1,:)=yp
      CALL SPLINE1_INTERP(k_vopt,mp,s0,f,np_f,xfun, &
                          value,iflag,message, &
                          -1)
      y_pr_f(:,j)=value(1,:)

      f(2:4,:)=zero
      f(1,:)=bp(:)
      CALL SPLINE1_INTERP(k_vopt,mp,s0,f,np_f,xfun, &
                          value,iflag,message, &
                          -1)
      bp_pr_f(:,j)=value(1,:)

      f(:,:)=zero
      f(1,:)=arclen(:)
      CALL SPLINE1_INTERP(k_vopt,mp,s0,f,np_f,xfun, &
                          value,iflag,message, &
                          -1)
      arclen_pr_f(:,j)=value(1,:)

      !(x,y) values evenly spaced in dl
      s0(1:mp)=arclen(1:mp)/arclen(mp)

      f(:,:)=zero
      f(1,:)=xp
      CALL SPLINE1_INTERP(k_vopt,mp,s0,f,np_f,xfun, &
                          value,iflag,message, &
                          -1)
      x_sr_f(:,j)=value(1,:)

      f(2:4,:)=zero
      f(1,:)=yp
      CALL SPLINE1_INTERP(k_vopt,mp,s0,f,np_f,xfun, &
                          value,iflag,message, &
                          -1)
      y_sr_f(:,j)=value(1,:)

      !Index of point nearest the plane of the axis on inside
      ii=0
      yref=1.0e6_rspec

      DO i=1,mp

        ydif=ABS(yp(i)-zmag_f)

        IF(xp(i) <= rmag_f .AND. &
           ydif <= yref) THEN

          ii=i
          yref=ydif

        ENDIF

        IF(yp(i) <= zmin_f(j)) THEN

          zmin_f(j)=yp(i)
          rzmin_f(j)=xp(i)

        ENDIF

        IF(yp(i) >= zmax_f(j)) THEN

          zmax_f(j)=yp(i)
          rzmax_f(j)=xp(i)

        ENDIF

      ENDDO

      rin_f(j)=xp(ii)

    ENDDO

    rin_f(1)=rmag_f

!------------------------------------------------------------------------------
!Volume enclosed
!------------------------------------------------------------------------------
vol_f(:)=zero

DO j=2,nr_f

  DO i=2,np_f

    vol_f(j)=vol_f(j)+(x_sr_f(i,j)-x_sr_f(i-1,j))*(x_sr_f(i,j)+x_sr_f(i-1,j)) &
             *(y_sr_f(i,j)+y_sr_f(i-1,j))/4

  ENDDO

  vol_f(j)=ABS(vol_f(j))*2*z_pi

ENDDO

!-----------------------------------------------------------------------------
!Integral dl/bpoloidal
!-----------------------------------------------------------------------------
DO j=2,nr_f

  xfun(:)=one/bp_pr_f(:,j)
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),xfun,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV_SURF(9)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  dlobp_f(j)=c2(2)

ENDDO

dlobp_f(1)=dlobp_f(2)

!-------------------------------------------------------------------------------
!Total magnetic field and drift parameter
!-------------------------------------------------------------------------------
DO j=nr_f,2,-1 !Over flux surfaces

  DO i=1,np_f

    !Get psi and derivatives
    CALL FLUXAV_EVAL(x_f,nx_f,y_f,ny_f,cspln_xy_f,x_pr_f(i,j),y_pr_f(i,j),6, &
                     pds,iflag,message)

    IF(iflag /= 0) THEN

      iflag=1
      message='FLUXAV_SURF(8)/'//message
      IF(iflag == 1) GOTO 9999

    ENDIF

    IF(l_flip_f) THEN

      !Real psi is opposite sign to what is returned
      pds(1)=-pds(1)
      pds(2)=-pds(2)
      pds(3)=-pds(3)
      pds(4)=-pds(4)
      pds(5)=-pds(5)
      pds(6)=-pds(6)

    ENDIF

    !B field components
    capbr=-pds(3)/x_pr_f(i,j)
    capbz=pds(2)/x_pr_f(i,j)
    btot_pr_f(i,j)=SQRT(capbr**2+capbz**2+(f_r_f(j)/x_pr_f(i,j))**2)

    !B field derivatives
    dbrdr=(pds(3)/x_pr_f(i,j)-pds(4))/x_pr_f(i,j)
    dbzdr=(-pds(2)/x_pr_f(i,j)+pds(5))/x_pr_f(i,j)
    dbtdr=(-f_r_f(j)/x_pr_f(i,j)+ffp_r_f(j)/f_r_f(j)*pds(2))/x_pr_f(i,j)
    dbdr=(capbr*dbrdr+capbz*dbzdr+f_r_f(j)/x_pr_f(i,j)*dbtdr)/btot_pr_f(i,j)
    dbrdz=-pds(6)/x_pr_f(i,j)
    dbzdz=pds(4)/x_pr_f(i,j)
    dbtdz=ffp_r_f(j)/f_r_f(j)*pds(3)/x_pr_f(i,j)
    dbdz=(capbr*dbrdz+capbz*dbzdz+f_r_f(j)/x_pr_f(i,j)*dbtdz)/btot_pr_f(i,j)

    !n dot grad(B)
    xndgrb_pr_f(i,j)=-(capbr*dbdr+capbz*dbdz)/btot_pr_f(i,j)

  ENDDO

ENDDO

!-----------------------------------------------------------------------------
!n.grad(theta)
!-----------------------------------------------------------------------------
!Use the identity Lp/2pi*B dot grad(theta) = bp
DO j=1,nr_f

  xnthi_pr_f(:,j)=ABS(btot_pr_f(:,j))/bp_pr_f(:,j)

ENDDO

!-----------------------------------------------------------------------------
!grth_r=n dot grad(CapTheta)
!-----------------------------------------------------------------------------
DO j=2,nr_f

  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),xnthi_pr_f(:,j),2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV_SURF(10)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  grth_f(j)=2*z_pi/c2(2)

ENDDO

grth_f(1)=grth_f(2)

!-------------------------------------------------------------------------------
!CapTheta
!-------------------------------------------------------------------------------
DO j=2,nr_f

  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),xnthi_pr_f(:,j),np_f, &
                     arclen_pr_f(:,j), &
                     ctheta_pr_f(:,j),iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV_SURF(11)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  ctheta_pr_f(:,j)=ctheta_pr_f(:,j)*grth_f(j)

ENDDO

!-------------------------------------------------------------------------------
!<1/R**2>
!-------------------------------------------------------------------------------
DO j=2,nr_f

  xfun(:)=one/bp_pr_f(:,j)/(x_pr_f(:,j)**2)
  c1(1)=arclen_pr_f(1,j)
  c1(2)=arclen_pr_f(np_f,j)
  CALL LINEAR1_INTEG(0,np_f,arclen_pr_f(:,j),xfun,2,c1, &
                     c2,iflag,message)

  IF(iflag /= 0) THEN

    iflag=1
    message='FLUXAV_SURF(12)/'//message
    IF(iflag == 1) GOTO 9999

  ENDIF

  rm2_f(j)=c2(2)/dlobp_f(j)

ENDDO

rm2_f(1)=one/rmag_f**2

!-------------------------------------------------------------------------------
!Toroidal flux at interior nodes, integrate dV/dPhi=2*pi/F/<R^-2>
!-------------------------------------------------------------------------------
phit_f(1)=zero

DO j=2,nr_f

  phit_f(j)=phit_f(j-1)+(f_r_f(j-1)*rm2_f(j-1)+f_r_f(j)*rm2_f(j)) &
                        /(4*z_pi)*(vol_f(j)-vol_f(j-1))

ENDDO

!-------------------------------------------------------------------------------
!Toroidal flux grid
!-------------------------------------------------------------------------------
rhot_f(1:nr_f)=SQRT(phit_f(1:nr_f)/phit_f(nr_f))

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FLUXAV_SURF

END MODULE FLUXAV_MOD