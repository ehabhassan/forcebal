MODULE NCLASS_MOD
!-------------------------------------------------------------------------------
!NCLASS-Calculates NeoCLASSical transport properties
!
!NCLASS_MOD is an F90 module to calculate neoclassical transport properties
!  for an axisymmetric toroidal plasma
!
!References:
!
!  S.P.Hirshman, D.J.Sigmar, Nucl Fusion 21 (1981) 1079
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!                analytic friction coefficients for all Ti/Tj 3/2006
!                recursion for friction and viscosity coefficients 8/2006
!
!Contains PUBLIC routines:
!
!    NCLASS            -neoclassical properties on a single surface
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  NCLASS_FLOW,         & !flow velocities within a surface
                         !  called from NCLASS
  NCLASS_FRICTION,     & !calculates friction coefficients to any order
                         !  called from NCLASS_MN
  NCLASS_INIT,         & !initialize arrays
                         !  called from NCLASS
  NCLASS_K,            & !velocity-dependent viscosities
                         !  called from NCLASS_MU
  NCLASS_LAGUERRE,     & !calculate Laguerre moments from velocity moments
                         !  called from NCLASS_FRICTION
  NCLASS_MN,           & !friction coefficients
                         !  called from NCLASS
  NCLASS_MU,           & !fluid viscosities
                         !  called from NCLASS
  NCLASS_NU,           & !collision frequencies
                         !  called from NCLASS_K
  NCLASS_TAU,          & !collision times
                         !  called from NCLASS
  NCLASS_VISCOSITY,    & !calculates viscosity coefficients to any order
                         !  called from ????
  NCLASS_BACKSUB,      & !LU matrix back substitution
                         !  called from NCLASS_FLOW
  NCLASS_DECOMP,       & !LU matrix decomposition
                         !  called from NCLASS_FLOW
  NCLASS_ERF             !error function
                         !  called from NCLASS_NU

!-------------------------------------------------------------------------------
! Private data
!-------------------------------------------------------------------------------
INTEGER, PRIVATE, PARAMETER :: &
  dpspec=SELECTED_REAL_KIND(15,300)

!Logical switches
LOGICAL, PRIVATE :: &
  l_banana_nc,         & !option to include banana viscosity [logical]
  l_classical_nc,      & !option to include classical transport [logical]
  l_fast_nc,           & !option to include fast ions [logical]
  l_pfirsch_nc,        & !option to include Pfirsch-Schluter transport [logical]
  l_potato_nc,         & !option to include potato orbits [logical]
  l_spitzer_nc           !option to solve Spitzer problem [logical]

!Options
INTEGER, PRIVATE :: &
  k_order_nc             !order of v moments to be solved [-]
                         !=2 u and q (default)
                         !=3 u, q, and u2
                         !...

!Constants
REAL(KIND=rspec), PRIVATE, SAVE :: &
  cden_nc,             & !density cutoff to ignore species (default 1.e10) [/m**3]
  cpotb_nc,            & !kappa(0)*Bt(0)/[2*q(0)**2] [T]
  cpotl_nc               !q(0)*R(0) [m]

!Terms summed over species
REAL(KIND=rspec), PRIVATE, SAVE :: &
  petap_nc,            & !parallel electrical resistivity [Ohm*m]
  pjbbs_nc,            & !<J_bs.B> [A*T/m**2]
  pjbf_nc,             & !<J_f.B> current response to fast ions [A*T/m**2]
  pjboh_nc,            & !<J_OH.B> Ohmic current [A*T/m**2]
  vc3_nc                 !total vcrit**3 for fast ions [m**3/s**3]

!Dimensions
INTEGER, PRIVATE, SAVE :: &
  mf_nc,               & !number of fast ion species [-]
  mfm_nc,              & !number of moments in Fm expansion [-]
  mi_nc,               & !number of isotopes [-]
  ms_nc,               & !number of thermal species (ms>1) [-]
  mz_nc                  !maximimum Z of any species [-]

!Electron isotope and species identification
INTEGER, PRIVATE, SAVE :: &
  imel_nc,             & !index of electron isotope [-]
  jsel_nc                !index of electron species [-]

!Moments arrays
REAL(KIND=dpspec), PRIVATE, SAVE, ALLOCATABLE :: &
  clag_nc(:)             !Laguerre coefficients [-]

!Isotope arrays
REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  vti_nc(:),           & !thermal velocities of isotopes [m/s]
  amntii_nc(:,:),      & !sum Mn/tau over charge states [kg/m**3/s]
  calmi_nc(:,:,:),     & !tp eff friction matrix [kg/m**3/s]
  calnii_nc(:,:,:,:),  & !fp eff friction matrix [kg/m**3/s]
  capmii_nc(:,:,:,:),  & !test particle (tp) friction matrix [-]
  capnii_nc(:,:,:,:)     !test particle (tp) friction matrix [-]

!Isotope and charge arrays
REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  grppiz_nc(:,:)         !p'+nePhi_E' for each species [keV/m**3/rho]

!Species arrays
INTEGER, PRIVATE, SAVE, ALLOCATABLE :: &
  jms_nc(:),           & !isotope number of s [-]
  jzs_nc(:)              !charge state of s [-]

REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  sqzs_nc(:),          & !orbit squeezing factor for species [-]
  vc3s_nc(:),          & !vcrit**3 for fast ions on species [m**3/s**3]
  xis_nc(:),           & !charge weighted density factor of species [-]
  ymus_nc(:,:,:),      & !normalized viscosity for species [kg/m**3/s]
  bsjbps_nc(:),        & !<J_bs.B> driven by unit p'/p of s [A*T*rho/m**2]
  bsjbts_nc(:),        & !<J_bs.B> driven by unit T'/T of s [A*T*rho/m**2]
  upars_nc(:,:,:),     & !(k,m,s) parallel flow of s from force m [T*m/s]
                         !m=1, p', T', Phi'
                         !m=2, <E.B>
                         !m=3, fast ion driven
  uthetas_nc(:,:,:),   & !poloidal flow of s from force m [m/s/T]
                         !m=1, p', T'
                         !m=2, <E.B>
                         !m=3, fast ion driven
  gfls_nc(:,:),        & !radial heat conduction flux comps of s [W*rho/m**3]
                         !m=1, banana-plateau, p' and T'
                         !m=2, Pfirsch-Schluter
                         !m=3, classical
                         !m=4, banana-plateau, <E.B>
                         !m=5, fast ion driven
  dpss_nc(:,:),        & !diffusion coefficient of s2 on p'/p of s1 [rho**2/s]
  dtss_nc(:,:),        & !diffusion coefficient of s2 on T'/T of s1 [rho**2/s]
  qfls_nc(:,:),        & !radial heat conduction flux comps of s [W*rho/m**3]
                         !m=1, banana-plateau, p' and T'
                         !m=2, Pfirsch-Schluter
                         !m=3, classical
                         !m=4, banana-plateau, <E.B>
                         !m=5, fast ion driven
  chipss_nc(:,:),      & !heat cond coefficient of s2 on p'/p of s1 [rho**2/s]
  chitss_nc(:,:),      & !heat cond coefficient of s2 on T'/T of s1 [rho**2/s]
  tauss_nc(:,:)          !90 degree scattering time [s]

!Fast ion arrays
REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &
  amuf_nc(:),          & !atomic mass of fast ion components [-]
  calmif_nc(:,:,:,:),  & !test part fast ion friction [kg/m**3/s
  denf_nc(:),          & !density of fast ion components [/m**3]
  ef_nc(:),            & !initial energy of fast ion components [keV]
  tausf_nc(:),         & !characteristic fast ion slowing down time on electrons [s]
  uf_nc(:),            & !<v.B> of fast ion components [T*m/s]
  zf_nc(:)               !charge number of fast ion components [-]

!Physical constants, mathematical constants, conversion factors
REAL(KIND=rspec), PRIVATE, PARAMETER :: &
  z_coulomb=1.6022e-19_rspec,    & !Coulomb charge [coul]
  z_epsilon0=8.8542e-12_rspec,   & !permittivity of free space [F/m]
  z_j7kv=1.6022e-16_rspec,       & !energy conversion factor [J/keV]
  z_pi=3.141592654_rspec,        & !pi [-]
  z_pmass=1.6726e-27_rspec         !proton mass [kg]

REAL(KIND=rspec), PRIVATE, PARAMETER :: &
  one=1.0_rspec,       & !REAL 1
  zero=0.0_rspec         !REAL 0

REAL(dpspec), PRIVATE, PARAMETER :: &
  onedp=1.0_dpspec,    & !REAL DP 1
  zerodp=0.0_dpspec      !REAL DP 0

!-------------------------------------------------------------------------------
! Procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE NCLASS(m_i,m_z,p_b2,p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi, &
                  p_gr2phi,p_ngrth,amu_i,grt_i,temp_i,den_iz,grp_iz, &
                  iflag,message, &
                  L_BANANA,L_PFIRSCH,L_CLASSICAL,L_POTATO,K_ORDER, &
                  C_DEN,C_POTB,C_POTL, &
                  AMU_F,Z_F,E_F,U_F,DEN_F, &
                  P_ETAP,P_JBBS,P_JBF,P_JBOH,P_FSHLD,&
                  M_S,JM_S,JZ_S, &
                  BSJBP_S,BSJBT_S, &
                  GFL_S,DN_S,VNNT_S,VNEB_S,VNF_S,DP_SS,DT_SS, &
                  UPAR_S,UTHETA_S, &
                  QFL_S,CHI_S,VQNT_S,VQEB_S,VQF_S, &
                  CHIP_SS,CHIT_SS, &
                  CALM_I,CALN_II,CAPM_II,CAPN_II,YMU_S, &
                  SQZ_S,XI_S,TAU_SS,TAUS_F,CALM_IF)
!-------------------------------------------------------------------------------
!NCLASS calculates the neoclassical transport properties of a multiple
!  species axisymmetric plasma using k_order parallel and radial force
!  balance equations for each species
!
!References:
!  S.P.Hirshman, D.J.Sigmar, Nucl Fusion 21 (1981) 1079
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!  W.A.Houlberg, named DO loop on EXIT 1/2005
!-------------------------------------------------------------------------------
!Declaration of input variables
INTEGER, INTENT(IN) :: &
  m_i,                 & !number of isotopes (> 1) [-]
  m_z                    !highest charge state [-]

REAL(KIND=rspec), INTENT(IN) :: &
  p_b2,                & !<B**2> [T**2]
  p_bm2,               & !<1/B**2> [/T**2]
  p_eb,                & !<E.B> [V*T/m]
  p_fhat,              & !mu_0*F/(dPsi/dr) [rho/m]
  p_fm(:),             & !poloidal moments of drift factor for PS [/m**2]
  p_ft,                & !trapped fraction [-]
  p_grbm2,             & !<grad(rho)**2/B**2> [rho**2/m**2/T**2]
  p_grphi,             & !potential gradient Phi' [V/rho]
  p_gr2phi,            & !second potential gradient Psi'(Phi'/Psi')' [V/rho**2]
  p_ngrth,             & !<n.grad(Theta)> [/m]
  amu_i(:),            & !atomic mass number [-]
  grt_i(:),            & !temperature gradient [keV/rho]
  temp_i(:),           & !temperature [keV]
  den_iz(:,:),         & !density [/m**3]
  grp_iz(:,:)            !pressure gradient [keV/m**3/rho]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional input variables
LOGICAL, INTENT(IN), OPTIONAL :: &
  L_BANANA,            & !option to include banana viscosity [logical]
  L_PFIRSCH,           & !option to include Pfirsch-Schluter transport [logical]
  L_CLASSICAL,         & !option to include classical transport [logical]
  L_POTATO               !option to include potato orbits [logical]

INTEGER, INTENT(IN), OPTIONAL :: &
  K_ORDER                !order of v moments to be solved [-]
                         !=2 u and q (default)
                         !=3 u, q, and u2
                         !... validated to 15

REAL(KIND=rspec), INTENT(IN), OPTIONAL :: &
  C_DEN,               & !density cutoff to ignore species (default 1.e10) [/m**3]
  C_POTB,              & !kappa(0)*Bt(0)/[2*q(0)**2] [T]
  C_POTL                 !q(0)*R(0) [m]

REAL(KIND=rspec), INTENT(IN), OPTIONAL :: &
  AMU_F(:),            & !atomic mass of fast ion components [-]
  Z_F(:),              & !charge number of fast ion components [-]
  E_F(:),              & !initial energy of fast ion components [keV]
  DEN_F(:),            & !density of fast ion components [/m**3]
  U_F(:)                 !<v.B> of fast ion components [T*m/s]

!Declaration of optional output variables
!Species mapping
INTEGER, INTENT(OUT), OPTIONAL ::&
  M_S,                 & !number of species (ms>1) [-]
  JM_S(:),             & !isotope number of s [-]
  JZ_S(:)                !charge state of s [-]

!Terms summed over species
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  P_ETAP,              & !parallel electrical resistivity [Ohm*m]
  P_FSHLD,             & !plasma shielding factor
  P_JBBS,              & !<J_bs.B> [A*T/m**2]
  P_JBF,               & !<J_f.B> current response to fast ions [A*T/m**2]
  P_JBOH                 !<J_OH.B> Ohmic current [A*T/m**2]

!Bootstrap current and electrical resistivity
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  BSJBP_S(:),          & !<J_bs.B> driven by unit p'/p of s [A*T*rho/m**2]
  BSJBT_S(:)             !<J_bs.B> driven by unit T'/T of s [A*T*rho/m**2]

!Continuity equation
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  DP_SS(:,:),          & !diffusion coefficient of s2 on p'/p of s1 [rho**2/s]
  DT_SS(:,:),          & !diffusion coefficient of s2 on T'/T of s1 [rho**2/s]
  GFL_S(:,:),          & !radial particle flux comps of s [rho/m**3/s]
                         !m=1, banana-plateau, p' and T'
                         !m=2, Pfirsch-Schluter
                         !m=3, classical
                         !m=4, banana-plateau, <E.B>
                         !m=5, fast ion driven
  DN_S(:),             & !diffusion coefficients (diag comp) [rho**2/s]
  VNEB_S(:),           & !<E.B> particle convection velocity [rho/s]
  VNF_S(:),            & !fast ion driven particle convection velocity [rho/s]
  VNNT_S(:)              !convection velocity (off diag p',T' comps) [rho/s]

!Momentum equation
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  UPAR_S(:,:,:),       & !parallel flow of s from force m [T*m/s]
                         !m=1, p', T', Phi'
                         !m=2, <E.B>
                         !m=3, fast ion driven
  UTHETA_S(:,:,:)        !poloidal flow of s from force m [m/s/T]
                         !m=1, p', T'
                         !m=2, <E.B>
                         !m=3, fast ion driven

!Energy equation
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  CHIP_SS(:,:),        & !heat cond coefficient of s2 on p'/p of s1 [rho**2/s]
  CHIT_SS(:,:),        & !heat cond coefficient of s2 on T'/T of s1 [rho**2/s]
  QFL_S(:,:),          & !radial heat conduction flux comps of s [W*rho/m**3]
                         !m=1, banana-plateau, p' and T'
                         !m=2, Pfirsch-Schluter
                         !m=3, classical
                         !m=4, banana-plateau, <E.B>
                         !m=5, fast ion driven
  CHI_S(:),            & !conduction coefficients (diag comp) [rho**2/s]
  VQEB_S(:),           & !<E.B> heat convection velocity [rho/s]
  VQF_S(:),            & !fast ion driven heat convection velocity [rho/s]
  VQNT_S(:)              !conduction velocity (off diag p',T' comps) [rho/s]

!Friction coefficients
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  CAPM_II(:,:,:,:),    & !test particle (tp) friction matrix [-]
  CAPN_II(:,:,:,:),    & !field particle (fp) friction matrix [-]
  CALM_I(:,:,:),       & !tp eff friction matrix [kg/m**3/s]
  CALN_II(:,:,:,:)       !fp eff friction matrix [kg/m**3/s]

!Viscosity coefficients
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  YMU_S(:,:,:)           !normalized viscosity for species [kg/m**3/s]

!Miscellaneous
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  SQZ_S(:),            & !orbit squeezing factor for species [-]
  XI_S(:),             & !charge weighted density factor of species [-]
  TAU_SS(:,:)            !90 degree scattering time [s]

!Fast ion terms
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  TAUS_F(:),           & !characteristic fast ion slowing down on electrons [s]
  CALM_IF(:,:,:,:)       !test part fast ion friction [kg/m**3/s

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,iza

REAL(KIND=rspec) :: &
  dent

REAL(KIND=rspec), ALLOCATABLE :: &
  denz2(:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Number of velocity moments
IF(PRESENT(K_ORDER)) THEN

  IF(K_ORDER < 2 .OR. &
     K_ORDER > 15) THEN

    iflag=1
    message='NCLASS(1)/ERROR: K_ORDER out of range'
    GOTO 9999

  ELSE

    k_order_nc=K_ORDER

  ENDIF

ELSE

  k_order_nc=2

ENDIF

!Banana contribution to viscosity
!On by default, unless user turns it off or specifies p_ft=0
l_banana_nc=.TRUE.
IF(ABS(p_ft) <= zero) l_banana_nc=.FALSE.
IF(PRESENT(L_BANANA) .AND. &
   (.NOT. L_BANANA)) l_banana_nc=.FALSE.

!Pfirsch-Schluter contribution to viscosity and transport
!On by default, unless user turns it off or specifies p_fm=0
IF(ABS(SUM(p_fm(:))) > zero) THEN

  l_pfirsch_nc=.TRUE.
  mfm_nc=SIZE(p_fm)

ELSE

  l_pfirsch_nc=.FALSE.
  mfm_nc=0

ENDIF

IF(PRESENT(L_PFIRSCH) .AND. &
   (.NOT. L_PFIRSCH)) l_pfirsch_nc=.FALSE.

!Classical contribution to transport
!On by default, unless user turns it off
l_classical_nc=.TRUE.
IF(PRESENT(L_CLASSICAL)) l_classical_nc=L_CLASSICAL

!Check whether this is a Spitzer problem
!If no neoclassical effects are included (ie the viscosity is zero because both
!  banana and Pfirsch-Schluter viscosities are off) this is a Spitzer problem
!In a Spitzer problem the bulk velocity is undetermined and the solution
!  requires setting a reference velocity
l_spitzer_nc=.TRUE.
IF(l_banana_nc .OR. &
   l_pfirsch_nc) l_spitzer_nc=.FALSE.

!Check potato orbit contribution to viscosity
!Off by default, unless user turns it on and provides constants
!No potato viscosity if either the banana or Pfirsch-Schluter are turned off
IF(.NOT. l_banana_nc .OR. &
   .NOT. l_pfirsch_nc) THEN

  l_potato_nc=.FALSE.
  cpotb_nc=zero
  cpotl_nc=zero

ELSEIF(PRESENT(L_POTATO) .AND. &
       PRESENT(C_POTB) .AND. &
       PRESENT(C_POTL)) THEN

  SELECT CASE (L_POTATO)

  CASE (.TRUE.)

    l_potato_nc=.TRUE.
    cpotb_nc=C_POTB
    cpotl_nc=C_POTL

  CASE (.FALSE.)

    l_potato_nc=.FALSE.
    cpotb_nc=zero
    cpotl_nc=zero

  END SELECT

ELSE

  l_potato_nc=.FALSE.
  cpotb_nc=zero
  cpotl_nc=zero

ENDIF

!Use input value of cutoff density, if present
cden_nc=1.0e10_rspec
IF(PRESENT(C_DEN)) cden_nc=C_DEN

!Trapped fraction between 0 and 1 inclusive
IF((p_ft < zero) .OR. &
   (p_ft > one)) THEN

  iflag=1
  message='NCLASS(2)/ERROR: must have 0 <= p_ft <= 1'
  GOTO 9999

ENDIF

!Fast ion contributions
l_fast_nc=.FALSE.

IF(PRESENT(AMU_F) .AND.&
   PRESENT(Z_F) .AND. &
   PRESENT(E_F) .AND. &
   PRESENT(U_F) .AND. &
   PRESENT(DEN_F)) THEN

  LOOP_I: DO i=1,SIZE(AMU_F) !Over fast components

    mf_nc=i
    IF(AMU_F(i) < 0.99_rspec .OR. &
       DEN_F(i) < cden_nc) THEN

      mf_nc=i-1
      EXIT LOOP_I

    ENDIF

  ENDDO LOOP_I !Over fast components

  IF(mf_nc > 0) l_fast_nc=.TRUE.

ENDIF

!Allocate private data
CALL NCLASS_INIT(m_i,m_z,amu_i,temp_i,den_iz, &
                 iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='NCLASS(3)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

!Set fast ion data
IF(l_fast_nc) THEN

  amuf_nc(1:mf_nc)=AMU_F(1:mf_nc)
  zf_nc(1:mf_nc)=Z_F(1:mf_nc)
  ef_nc(1:mf_nc)=E_F(1:mf_nc)
  uf_nc(1:mf_nc)=U_F(1:mf_nc)
  denf_nc(1:mf_nc)=DEN_F(1:mf_nc)

ENDIF

!-------------------------------------------------------------------------------
!Set various plasma properties
!-------------------------------------------------------------------------------
!Friction coefficients and collision times
CALL NCLASS_MN(amu_i,temp_i,den_iz)

!Species charge state density factor, total pressure, squeezing
dent=zero
ALLOCATE(denz2(mi_nc))
denz2(:)=zero

DO i=1,mi_nc !Over isotopes

  DO iza=1,mz_nc !Over charge states

    IF(den_iz(i,iza) > cden_nc) THEN

      denz2(i)=denz2(i)+den_iz(i,iza)*REAL(iza,rspec)**2
      dent=dent+den_iz(i,iza)*temp_i(i)

    ENDIF

  ENDDO !Over charge states

ENDDO !Over isotopes

DO i=1,ms_nc !Over species

  iza=ABS(jzs_nc(i))
  xis_nc(i)=den_iz(jms_nc(i),iza)*jzs_nc(i)**2/denz2(jms_nc(i))
  sqzs_nc(i)=one+p_fhat**2/p_b2*amu_i(jms_nc(i))*z_pmass*p_gr2phi &
                 /(z_coulomb*jzs_nc(i))
  IF(sqzs_nc(i) > 10.0_rspec) sqzs_nc(i)=10.0_rspec
  IF(sqzs_nc(i) < 0.5_rspec) sqzs_nc(i)=0.5_rspec

ENDDO !Over species

DEALLOCATE(denz2)

!Normalized viscosities
CALL NCLASS_MU(p_fm,p_ft,p_ngrth,amu_i,temp_i,den_iz)

!Add potential gradient to pressure gradient
DO i=1,ms_nc !Over species

  iza=ABS(jzs_nc(i))
  grppiz_nc(jms_nc(i),iza)=grp_iz(jms_nc(i),iza) &
                           +p_grphi*den_iz(jms_nc(i),iza)*jzs_nc(i) &
                           *z_coulomb/z_j7kv

ENDDO !Over species

!-------------------------------------------------------------------------------
!Solve for parallel and poloidal flows within a surface
!-------------------------------------------------------------------------------
iflag=0

CALL NCLASS_FLOW(p_b2,p_bm2,p_eb,p_fhat,p_grbm2,grt_i,temp_i,den_iz, &
                 iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='NCLASS(4)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Output
!-------------------------------------------------------------------------------
!Species mapping
IF(PRESENT(M_S)) M_S=ms_nc
IF(PRESENT(JM_S)) JM_S(1:ms_nc)=jms_nc(1:ms_nc)
IF(PRESENT(JZ_S)) JZ_S(1:ms_nc)=jzs_nc(1:ms_nc)

!Electrical resistivity and currents
IF(PRESENT(P_ETAP)) P_ETAP=petap_nc
IF(PRESENT(P_JBBS)) P_JBBS=pjbbs_nc
IF(PRESENT(P_JBOH)) P_JBOH=pjboh_nc
IF(PRESENT(BSJBP_S)) BSJBP_S(1:ms_nc)=bsjbps_nc(1:ms_nc)
IF(PRESENT(BSJBT_S)) BSJBT_S(1:ms_nc)=bsjbts_nc(1:ms_nc)

!Continuity equation
IF(PRESENT(DP_SS)) DP_SS(1:ms_nc,1:ms_nc)=dpss_nc(1:ms_nc,1:ms_nc)
IF(PRESENT(DT_SS)) DT_SS(1:ms_nc,1:ms_nc)=dtss_nc(1:ms_nc,1:ms_nc)
IF(PRESENT(GFL_S)) GFL_S(1:5,1:ms_nc)=gfls_nc(1:5,1:ms_nc)

IF(PRESENT(DN_S)) THEN

  DO i=1,ms_nc !Over species

    DN_S(i)=dpss_nc(i,i)

    IF(i /= jsel_nc .AND. &
       jsel_nc /= 0) THEN

      !Redistribute electron density gradient flux among ions
      iza=ABS(jzs_nc(i))
      DN_S(i)=DN_S(i)+dpss_nc(jsel_nc,i)*iza*den_iz(jms_nc(i),iza) &
                     /den_iz(jms_nc(jsel_nc),1)

    ENDIF

  ENDDO !Over species

ENDIF

IF(PRESENT(VNNT_S)) THEN

  DO i=1,ms_nc !Over species

    iza=ABS(jzs_nc(i))
    VNNT_S(i)=(SUM(gfls_nc(1:3,i))+dpss_nc(i,i)*(grp_iz(jms_nc(i),iza) &
              -den_iz(jms_nc(i),iza)*grt_i(jms_nc(i)))/temp_i(jms_nc(i))) &
              /den_iz(jms_nc(i),iza)

    IF(i /= jsel_nc .AND. &
       jsel_nc /= 0) THEN

      !Redistribute electron density gradient flux among ions
      VNNT_S(i)=VNNT_S(i)+dpss_nc(jsel_nc,i)*iza*(grp_iz(jms_nc(i),iza) &
                -den_iz(jms_nc(i),iza)*grt_i(jms_nc(i)))/temp_i(jms_nc(i)) &
                /den_iz(jms_nc(jsel_nc),1)

    ENDIF          

  ENDDO !Over species

ENDIF

IF(PRESENT(VNEB_S)) THEN

  DO i=1,ms_nc !Over species

    iza=ABS(jzs_nc(i))
    VNEB_S(i)=gfls_nc(4,i)/den_iz(jms_nc(i),iza)

  ENDDO !Over species

ENDIF

!Momentum equation

IF(PRESENT(UPAR_S)) UPAR_S(1:k_order_nc,1:3,1:ms_nc) &
                          =upars_nc(1:k_order_nc,1:3,1:ms_nc)
IF(PRESENT(UTHETA_S)) UTHETA_S(1:k_order_nc,1:3,1:ms_nc) &
                              =uthetas_nc(1:k_order_nc,1:3,1:ms_nc)

!Energy equation
IF(PRESENT(QFL_S)) QFL_S(1:5,1:ms_nc)=qfls_nc(1:5,1:ms_nc)

IF(PRESENT(CHI_S)) THEN

  DO i=1,ms_nc !Over species

    CHI_S(i)=chipss_nc(i,i)+chitss_nc(i,i)

  ENDDO !Over species

ENDIF

IF(PRESENT(VQNT_S)) THEN

  DO i=1,ms_nc !Over species

    iza=ABS(jzs_nc(i))
    VQNT_S(i)=(SUM(qfls_nc(1:3,i))/z_j7kv/den_iz(jms_nc(i),iza)&
              +chi_s(i)*grt_i(jms_nc(i)))/temp_i(jms_nc(i))

  ENDDO !Over species

ENDIF

IF(PRESENT(VQEB_S)) THEN

  DO i=1,ms_nc !Over species

    iza=ABS(jzs_nc(i))
    VQEB_S(i)=qfls_nc(4,i)/den_iz(jms_nc(i),iza)/temp_i(jms_nc(i))/z_j7kv

  ENDDO !Over species

ENDIF

IF(PRESENT(CHIP_SS)) CHIP_SS(1:ms_nc,1:ms_nc)=chipss_nc(1:ms_nc,1:ms_nc)
IF(PRESENT(CHIT_SS)) CHIT_SS(1:ms_nc,1:ms_nc)=chitss_nc(1:ms_nc,1:ms_nc)

!Friction coefficients
IF(PRESENT(CALM_I)) CALM_I(1:k_order_nc,1:k_order_nc,1:mi_nc) &
                      =calmi_nc(1:k_order_nc,1:k_order_nc,1:mi_nc)
IF(PRESENT(CALN_II)) CALN_II(1:k_order_nc,1:k_order_nc,1:mi_nc,1:mi_nc) &
                       =calnii_nc(1:k_order_nc,1:k_order_nc,1:mi_nc,1:mi_nc)
IF(PRESENT(CAPM_II)) CAPM_II(1:k_order_nc,1:k_order_nc,1:mi_nc,1:mi_nc) &
                       =capmii_nc(1:k_order_nc,1:k_order_nc,1:mi_nc,1:mi_nc)
IF(PRESENT(CAPN_II)) CAPN_II(1:k_order_nc,1:k_order_nc,1:mi_nc,1:mi_nc) &
                       =capnii_nc(1:k_order_nc,1:k_order_nc,1:mi_nc,1:mi_nc)

!Viscosity coefficients
IF(PRESENT(YMU_S)) YMU_S(1:k_order_nc,1:k_order_nc,1:ms_nc) &
                     =ymus_nc(1:k_order_nc,1:k_order_nc,1:ms_nc)

!Miscellaneous
IF(PRESENT(SQZ_S)) SQZ_S(1:ms_nc)=sqzs_nc(1:ms_nc)
IF(PRESENT(XI_S)) XI_S(1:ms_nc)=xis_nc(1:ms_nc)
IF(PRESENT(TAU_SS)) TAU_SS(1:ms_nc,1:ms_nc)=tauss_nc(1:ms_nc,1:ms_nc)

!Fast ions
IF(l_fast_nc) THEN

  IF(PRESENT(P_FSHLD)) P_FSHLD=pjbf_nc/z_coulomb &
                       /SUM(denf_nc(1:mf_nc)*zf_nc(1:mf_nc)*uf_nc(1:mf_nc))
  IF(PRESENT(P_JBF)) P_JBF=pjbf_nc
  IF(PRESENT(TAUS_F)) TAUS_F(1:mf_nc)=tausf_nc(1:mf_nc)

  IF(PRESENT(VNF_S)) THEN

    DO i=1,ms_nc !Over species

      iza=ABS(jzs_nc(i))
      VNF_S(i)=gfls_nc(5,i)/den_iz(jms_nc(i),iza)

    ENDDO !Over species

  ENDIF

  IF(PRESENT(VQF_S)) THEN

    DO i=1,ms_nc !Over species

      iza=ABS(jzs_nc(i))
      VQF_S(i)=qfls_nc(5,i)/den_iz(jms_nc(i),iza)/temp_i(jms_nc(i))/z_j7kv

    ENDDO !Over species

  ENDIF

  IF(PRESENT(CALM_IF)) CALM_IF(1:k_order_nc,1:k_order_nc,1:mi_nc,1:mf_nc)= &
                       calmif_nc(1:k_order_nc,1:k_order_nc,1:mi_nc,1:mf_nc)

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE NCLASS

SUBROUTINE NCLASS_FLOW(p_b2,p_bm2,p_eb,p_fhat,p_grbm2,grt_i,temp_i,den_iz, &
                       iflag,message)
!-------------------------------------------------------------------------------
!NCLASS_FLOW calculates the neoclassical parallel flows u0, u1=q/p (and
!  optionally u2), the poloidal flows, then other transport properties
!
!References:                                                     
!  S.P.Hirshman, D.J.Sigmar, Nucl Fusion 21 (1981) 1079
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  p_b2,                & !<B**2> [T**2]
  p_bm2,               & !<1/B**2> [/T**2]
  p_eb,                & !<E.B> [V*T/m]
  p_fhat,              & !mu_0*F/(dPsi/dr) [rho/m]
  p_grbm2,             & !<grad(rho)**2/B**2> [rho**2/m**2/T**2]
  grt_i(:),            & !temperature gradient [keV/rho]
  temp_i(:),           & !temperature [keV]
  den_iz(:,:)            !density [/m**3]

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
  i,im,iz,iza,j,jm,k,l,l1,m,m1

INTEGER, DIMENSION(k_order_nc) :: &
  indxk,indxkf

INTEGER, DIMENSION(k_order_nc*mi_nc) :: &
  indxki,indxkif

REAL(KIND=rspec) :: &
  c,cbp,cbpa,cbpaq,cc,ccl,ccla,cclaq,cclb,cclbq,cps,cpsa,cpsaq,cpsb,cpsbq,d, &
  denzc

REAL(KIND=rspec), DIMENSION(k_order_nc) :: &
  xk

REAL(KIND=rspec), DIMENSION(k_order_nc,k_order_nc) :: &
  akk,akkf

REAL(KIND=rspec), DIMENSION(k_order_nc*mi_nc) :: &
  xkif

REAL(KIND=rspec), DIMENSION(k_order_nc*mi_nc,k_order_nc*mi_nc) :: &
  akiki,akikif

REAL(KIND=rspec) :: &
  crkit(k_order_nc,k_order_nc+2,mi_nc), &
  crkif(k_order_nc,k_order_nc+1,mi_nc), &
  srcth(k_order_nc,ms_nc), &
  srcthp(ms_nc),srctht(ms_nc), &
  srcf(k_order_nc,mi_nc), &
  crhatp(k_order_nc,ms_nc,mi_nc),crhatt(k_order_nc,ms_nc,mi_nc), &
  uaip(k_order_nc,ms_nc,ms_nc),uait(k_order_nc,ms_nc,ms_nc)

 REAL(KIND=rspec), TARGET :: &
  rkst(k_order_nc,k_order_nc+2,ms_nc), &
  rksf(k_order_nc,k_order_nc+1,ms_nc), &
  rhatp(k_order_nc,ms_nc,ms_nc),rhatt(k_order_nc,ms_nc,ms_nc), &
  xkit(k_order_nc*mi_nc,2), &
  xabp(k_order_nc*mi_nc,ms_nc),xabt(k_order_nc*mi_nc,ms_nc)

 REAL(KIND=rspec), POINTER :: &
  pk(:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Local constants
cc=(p_fhat/z_coulomb)*z_j7kv           !Coulomb
cbp=p_fhat/p_b2/z_coulomb              !Banana-plateau
cps=(p_fhat/z_coulomb)*(1/p_b2-p_bm2)  !Pfirsch-Schluter
ccl=(p_grbm2/z_coulomb)/p_fhat         !Classical

!Local arrays
akk(:,:)=zero                          !Species coeff matrix - nc
rkst(:,:,:)=zero                       !Species response - nc
crkit(:,:,:)=zero                      !Isotopic response - nc
srcth(:,:)=zero                        !Source - nc
indxk(:)=0                             !Decomposition record - nc
srcthp(:)=zero                         !Source - unit p'
srctht(:)=zero                         !Source - unit T'
rhatp(:,:,:)=zero                      !Species response - unit p'
rhatt(:,:,:)=zero                      !Species response - unit T'
crhatp(:,:,:)=zero                     !Isotopic response - unit p'
crhatt(:,:,:)=zero                     !Isotopic response - unit T'

IF(l_fast_nc) THEN

  akkf(:,:)=zero                       !Coeff matrix - fast
  rksf(:,:,:)=zero                     !Species response - fast
  crkif(:,:,:)=zero                    !Isotopic response - fast
  srcf(:,:)=zero                       !Source - fast
  indxkf(:)=0                          !Decomposition record - fast

  DO k=1,k_order_nc !Over velocity moments

    DO i=1,mi_nc !Over isotopes

      srcf(k,i)=-SUM(calmif_nc(k,1,i,:)*uf_nc(:))

    ENDDO !Over isotopes

  ENDDO !Over velocity moments

ENDIF

!-------------------------------------------------------------------------------
!Calculate species and isotopic information for reduced charge state formalism
!Here we construct species and isotropic responses to various sources (forces)
!  so that we can construct a reduced matrix to solve for the isotopic flows
!-------------------------------------------------------------------------------
DO i=1,ms_nc !Over species

  !Isotopic and charge state indices
  im=jms_nc(i)
  iz=jzs_nc(i)
  iza=IABS(iz)

  !Response matrix for thermal forces
  akk(:,:)=xis_nc(i)*calmi_nc(:,:,im)-ymus_nc(:,:,i)

  !Response matrix for fast ion forces
  IF(l_fast_nc) akkf(:,:)=akk(:,:)+xis_nc(i)*SUM(calmif_nc(:,:,im,1:mf_nc))

  !LU decomposition of response matrix for thermal forces
  iflag=0
  message=''
  CALL NCLASS_DECOMP(akk,k_order_nc,indxk, &
                     d,c,iflag,message)

  !Check messages
  IF(iflag == 1) THEN

    message='NCLASS_FLOW(1)/'//message
    GOTO 9999

  ENDIF

  !LU decomposition of response matrix for fast ion forces
  IF(l_fast_nc) THEN

    iflag=0
    message=''
    CALL NCLASS_DECOMP(akkf,k_order_nc,indxkf, &
                       d,c,iflag,message)

          !Check messages
    IF(iflag == 1) THEN

      message='NCLASS_FLOW(2)/'//message
      GOTO 9999

    ENDIF

  ENDIF

  !Response to lambda terms involving isotopic flows
  DO k=1,k_order_nc !Over velocity moments

    !Thermal forces
    rkst(k,k,i)=xis_nc(i)
    pk => rkst(:,k,i)
    CALL NCLASS_BACKSUB(akk,k_order_nc,indxk,pk)
    crkit(:,k,im)=crkit(:,k,im)+xis_nc(i)*rkst(:,k,i)

    !Fast ion forces
    IF(l_fast_nc) THEN

      rksf(k,k,i)=xis_nc(i)
      pk => rksf(:,k,i)
      CALL NCLASS_BACKSUB(akkf,k_order_nc,indxkf,pk)
      crkif(:,k,im)=crkif(:,k,im)+xis_nc(i)*rksf(:,k,i)

    ENDIF

  ENDDO !Over velocity moments

  !Response to poloidal source (p' and T') terms 
  srcth(1,i)=(cc/iz)*grppiz_nc(im,iza)/den_iz(im,iza)
  srcth(2,i)=(cc/iz)*grt_i(im)

  DO k=1,k_order_nc !Over velocity moments

    rkst(k,k_order_nc+1,i)=SUM(srcth(:,i)*ymus_nc(k,:,i))

  ENDDO !Over velocity moments

  pk => rkst(:,k_order_nc+1,i)
  CALL NCLASS_BACKSUB(akk,k_order_nc,indxk,pk)
  crkit(:,k_order_nc+1,im)=crkit(:,k_order_nc+1,im) &
                           +xis_nc(i)*rkst(:,k_order_nc+1,i)

  !Response to unit p'/p and T'/T terms for transport coefficients
  srcthp(i)=-(cc/iz)*temp_i(im)
  rhatp(:,i,i)=srcthp(i)*ymus_nc(:,1,i)
  pk => rhatp(:,i,i)
  CALL NCLASS_BACKSUB(akk,k_order_nc,indxk,pk)
  crhatp(:,i,im)=crhatp(:,i,im)+xis_nc(i)*rhatp(:,i,i)
  srctht(i)=-(cc/iz)*temp_i(im)
  rhatt(:,i,i)=srctht(i)*ymus_nc(:,2,i)
  pk => rhatt(:,i,i)
  CALL NCLASS_BACKSUB(akk,k_order_nc,indxk,pk)
  crhatt(:,i,im)=crhatt(:,i,im)+xis_nc(i)*rhatt(:,i,i)

  !Response to parallel electric field term for resistivity 
  rkst(1,k_order_nc+2,i)=-REAL(iz,rspec)*z_coulomb*den_iz(im,iza)
  pk => rkst(:,k_order_nc+2,i)
  CALL NCLASS_BACKSUB(akk,k_order_nc,indxk,pk)
  crkit(:,k_order_nc+2,im)=crkit(:,k_order_nc+2,im) &
                           +xis_nc(i)*rkst(:,k_order_nc+2,i)

  !Response to external force terms
  IF(l_fast_nc) THEN

    rksf(:,k_order_nc+1,im)=-xis_nc(i)*srcf(:,im)
    pk => rksf(:,k_order_nc+1,i)
    CALL NCLASS_BACKSUB(akkf,k_order_nc,indxkf,pk)
    crkif(:,k_order_nc+1,im)=crkif(:,k_order_nc+1,im) &
                             +xis_nc(i)*rksf(:,k_order_nc+1,i)

  ENDIF

ENDDO !Over species

!-------------------------------------------------------------------------------
!Load coefficient matrix and source terms for isotopic flows
!-------------------------------------------------------------------------------
!Initialization
!  akiki(i,j) - i=(1,mi_nc) momentum equations
!                =mi_nc+(1,mi_nc) heat flux equations
!                =2*mi_nc+(1,mi_nc) u2 equations (k_order=3)
!               j=(1,mi_nc) flow velocities
!                =mi_nc+(1,mi_nc) heat flows
!                =2*mi_nc+(1,mi_nc) u2 flows (k_order=3)
akiki(:,:)=zero                        !Isotopic coeff matrix - nc
xkit(:,:)=zero                         !Isotopic response - nc
xabp(:,:)=zero                         !Isotopic response to p' - nc
xabt(:,:)=zero                         !Isotopic response to T' - nc

IF(l_fast_nc) THEN

  akikif(:,:)=zero                     !Isotopic coeff matrix - fast
  xkif(:)=zero                         !Isotopic response - fast

ENDIF

DO im=1,mi_nc !Over isotopes 1

  DO m=1,k_order_nc !Over velocity moments 1

    m1=im+(m-1)*mi_nc

    !Standard neoclassical
    !Diagonal coefficients
    akiki(m1,m1)=one

    !p' and T' force terms
    xkit(m1,1)=crkit(m,k_order_nc+1,im)

    !Unit p'/p and T'/T
    xabp(m1,:)=crhatp(m,:,im)
    xabt(m1,:)=crhatt(m,:,im)

    !<E.B> force terms
    xkit(m1,2)=crkit(m,k_order_nc+2,im)

    !Fast ion driven neoclassical
    IF(l_fast_nc) THEN

      !Diagonal coefficients
      akikif(m1,m1)=one

      !Fast ion force
      xkif(m1)=crkif(m,k_order_nc+1,im)

    ENDIF

    !Field particle friction       
    DO jm=1,mi_nc !Over isotopes 2

      DO l=1,k_order_nc !Over velocity moments 2

        !Standard neoclassical
        l1=jm+(l-1)*mi_nc
        ! EHAB MODIFIED LINES BELOW
       !WRITE(*,*) SHAPE(calnii_nc(:,l,im,jm)),SHAPE(crkit(m,:,im))
       !akiki(m1,l1)=akiki(m1,l1)+SUM(calnii_nc(:,l,im,jm)*crkit(m,:,im))
        akiki(m1,l1)=akiki(m1,l1)+SUM(calnii_nc(:,l,im,jm)*crkit(m,1:k_order_nc,im))
        ! EHAB MODIFIED LINES ABOVE

        !Fast ion driven neoclassical
        IF(l_fast_nc) THEN

          akikif(m1,l1)=akikif(m1,l1)+SUM(calnii_nc(:,l,im,jm)*crkif(m,:,im))

        ENDIF

      ENDDO !Over velocity moments 2

    ENDDO !Over isotopes 2

  ENDDO !Over velocity moments 1
     
ENDDO !Over isotopes 1

!Check if Spitzer problem, where bulk flow is indeterminate
IF(l_spitzer_nc) THEN

  !Force last species velocity
  akiki(mi_nc,:)=zero
  akiki(mi_nc,mi_nc)=one

  IF(l_fast_nc) THEN

    !Force last species velocity
    akikif(mi_nc,:)=zero
    akikif(mi_nc,mi_nc)=one

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Get lu decomposition of coefficient matrices for isotopic flows
!-------------------------------------------------------------------------------
!Initialization
indxki(:)=0                                !Decomposition record - nc
IF(l_fast_nc) indxkif(:)=0                 !Decomposition record - fast

!Standard neoclassical
iflag=0
message=''
CALL NCLASS_DECOMP(akiki,k_order_nc*mi_nc,indxki, &
                   d,c,iflag,message)

!Check messages
IF(iflag == 1) THEN

  message='NCLASS_FLOW(3)/'//message
  GOTO 9999

ENDIF

!Fast ion driven neoclassical
IF(l_fast_nc) THEN

  iflag=0
  message=''
  CALL NCLASS_DECOMP(akikif,k_order_nc*mi_nc,indxkif, &
                     d,c,iflag,message)

        !Check messages
  IF(iflag == 1) THEN

    message='NCLASS_FLOW(4)/'//message
    GOTO 9999

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Evaluate isotopic and species flows from back substitution for each source 
!xkit(1:mi_nc,m) are the isotopic velocities 
!xkit(mi_nc+1:2*mi_nc,m) are the isotopic heat flows 
!xkit(2*mi_nc+1:3*mi_nc,m) are the u2 flows
!m=1 for pressure and temperature gradient driven flows
!m=2 for <E.B> driven flows
!Similarly for xkif except force index not needed
!-------------------------------------------------------------------------------
xk(:)=zero
uaip(:,:,:)=zero
uait(:,:,:)=zero

!Evaluate species response for standard neoclassical
DO m=1,2 !Over forces

  pk => xkit(:,m)
  CALL NCLASS_BACKSUB(akiki,k_order_nc*mi_nc,indxki,pk)

ENDDO !Over forces

!Evaluate species response for fast ion driven neoclassical
IF(l_fast_nc) CALL NCLASS_BACKSUB(akikif,k_order_nc*mi_nc,indxkif,xkif)

!Evaluate unit response for standard neoclassical transport coefficients
DO i=1,ms_nc !Over species

  pk => xabp(:,i)
  CALL NCLASS_BACKSUB(akiki,k_order_nc*mi_nc,indxki,pk)
  pk => xabt(:,i)
  CALL NCLASS_BACKSUB(akiki,k_order_nc*mi_nc,indxki,pk)

ENDDO !Over species

!Evaluate parallel flows
DO i=1,ms_nc !Over species 1

  im=jms_nc(i)

  !Thermal force contributions 
  DO m=1,2 !Over forces

    upars_nc(:,m,i)=rkst(:,k_order_nc+m,i)

    !Thermal response contributions           
    DO jm=1,mi_nc !Over isotopes

      xk(:)=zero

      DO l=1,k_order_nc !Over velocity moments

        l1=jm+(l-1)*mi_nc
        xk(:)=xk(:)-calnii_nc(:,l,im,jm)*xkit(l1,m)

      ENDDO !Over velocity moments

      DO l=1,k_order_nc !Over velocity moments

        upars_nc(:,m,i)=upars_nc(:,m,i)+xk(l)*rkst(:,l,i)

      ENDDO !Over velocity moments

    ENDDO !Over isotopes

  ENDDO !Over forces

  !Fast ion flows
  IF(l_fast_nc) THEN

    !Thermal force contributions 
    upars_nc(:,3,i)=rksf(:,k_order_nc+1,i)

    !Thermal response contributions           
    DO jm=1,mi_nc !Over isotopes

      xk(:)=zero

      DO l=1,k_order_nc !Over velocity moments

        l1=jm+(l-1)*mi_nc
        xk(:)=xk(:)-calnii_nc(:,l,im,jm)*xkif(l1)

      ENDDO !Over velocity moments

      DO l=1,k_order_nc !Over velocity moments

        upars_nc(:,3,i)=upars_nc(:,3,i)+xk(l)*rksf(:,l,i)

      ENDDO !Over velocity moments

    ENDDO !Over isotopes

  ENDIF 

!Unit p'/p and T'/T
  DO j=1,ms_nc !Over species 2

    uaip(:,j,i)=rhatp(:,j,i)
    uait(:,j,i)=rhatt(:,j,i)

  !Response contributions
    DO jm=1,mi_nc !Over isotopes

      xk(:)=zero

      DO l=1,k_order_nc !Over velocity moments

        l1=jm+(l-1)*mi_nc
        xk(:)=xk(:)-calnii_nc(:,l,im,jm)*xabp(l1,j)

      ENDDO !Over velocity moments

      DO l=1,k_order_nc !Over velocity moments

        uaip(:,j,i)=uaip(:,j,i)+xk(l)*rkst(:,l,i)

      ENDDO !Over velocity moments

      xk(:)=zero

      DO l=1,k_order_nc !Over velocity moments

        l1=jm+(l-1)*mi_nc
        xk(:)=xk(:)-calnii_nc(:,l,im,jm)*xabt(l1,j)

      ENDDO !Over velocity moments

      DO l=1,k_order_nc !Over velocity moments

        uait(:,j,i)=uait(:,j,i)+xk(l)*rkst(:,l,i)

      ENDDO !Over velocity moments

    ENDDO !Over isotopes

  ENDDO !Over species 2

ENDDO !Over species 1

!-------------------------------------------------------------------------------
!Calculate poloidal flows from parallel flows
!-------------------------------------------------------------------------------
uthetas_nc(:,:,:)=upars_nc(:,:,:)/p_b2

DO i=1,ms_nc !Over species

  iza=IABS(jzs_nc(i))

  DO k=1,k_order_nc !Over moments

    IF(k == 1) THEN

      !Add p' and Phi'
      uthetas_nc(k,1,i)=uthetas_nc(k,1,i) &
                        +p_fhat*grppiz_nc(jms_nc(i),iza)*z_j7kv &
                        /(z_coulomb*REAL(jzs_nc(i),rspec) &
                        *den_iz(jms_nc(i),iza))/p_b2

    ELSEIF(k == 2) THEN

      !Add T'
      uthetas_nc(k,1,i)=uthetas_nc(k,1,i) &
                        +p_fhat*z_j7kv*grt_i(jms_nc(i)) &
                        /(REAL(jzs_nc(i),rspec)*z_coulomb*p_b2)

    ENDIF

  ENDDO !Over moments

ENDDO !Over species

!-------------------------------------------------------------------------------
!Evaluate species currents, transport coefficients and fluxes
!-------------------------------------------------------------------------------
!Initialization
petap_nc=zero
pjbbs_nc=zero
pjbf_nc=zero
pjboh_nc=zero
bsjbps_nc(:)=zero
bsjbts_nc(:)=zero
dpss_nc(:,:)=zero
dtss_nc(:,:)=zero
chipss_nc(:,:)=zero
chitss_nc(:,:)=zero
gfls_nc(:,:)=zero
qfls_nc(:,:)=zero

!Fast ion current from stacking
IF(l_fast_nc) pjbf_nc=z_coulomb*SUM(zf_nc(:)*denf_nc(:)*uf_nc(:))

!Loop over thermal species contributions
DO i=1,ms_nc !Over species 1

  im=jms_nc(i)
  iz=jzs_nc(i)
  iza=IABS(iz)

  !Currents
  denzc=den_iz(im,iza)*iz*z_coulomb

  !Bootstrap <J_bs.B> 
  pjbbs_nc=pjbbs_nc+denzc*upars_nc(1,1,i)

  !Ohmic <J_OH.B>
  pjboh_nc=pjboh_nc+denzc*upars_nc(1,2,i)

  !Response of thermal species to fast ions <J_f.B>
  pjbf_nc=pjbf_nc+denzc*upars_nc(1,3,i)

  !Unit p'/p and T'/T contributions to bootstrap
  bsjbps_nc(:)=bsjbps_nc(:)+denzc*uaip(1,:,i)
  bsjbts_nc(:)=bsjbts_nc(:)+denzc*uait(1,:,i)

  !Banana-plateau transport coefficients and fluxes
  cbpa=cbp/iz
  cbpaq=cbpa*(z_j7kv*temp_i(im))

  !Transport coefficients from unit p'/p and T'/T
  dpss_nc(i,i)=dpss_nc(i,i)-cbpa*ymus_nc(1,1,i)*srcthp(i)
  dtss_nc(i,i)=dtss_nc(i,i)-cbpa*ymus_nc(1,2,i)*srctht(i)
  chipss_nc(i,i)=chipss_nc(i,i)-cbpaq*ymus_nc(2,1,i)*srcthp(i)
  chitss_nc(i,i)=chitss_nc(i,i)-cbpaq*ymus_nc(2,2,i)*srctht(i)

  DO k=1,k_order_nc !Over velocity moments

    dpss_nc(:,i)=dpss_nc(:,i)-cbpa*ymus_nc(1,k,i)*uaip(k,:,i)
    dtss_nc(:,i)=dtss_nc(:,i)-cbpa*ymus_nc(1,k,i)*uait(k,:,i)
    chipss_nc(:,i)=chipss_nc(:,i)-cbpaq*ymus_nc(2,k,i)*uaip(k,:,i)
    chitss_nc(:,i)=chitss_nc(:,i)-cbpaq*ymus_nc(2,k,i)*uait(k,:,i)

  ENDDO !Over velocity moments

  !Transport fluxes from total p' and T'
  gfls_nc(1,i)=gfls_nc(1,i)-cbpa*p_b2*SUM(ymus_nc(1,:,i)*uthetas_nc(:,1,i))
  qfls_nc(1,i)=qfls_nc(1,i)-cbpaq*p_b2*SUM(ymus_nc(2,:,i)*uthetas_nc(:,1,i))

  !Transport fluxes from unit <E.B> 
  gfls_nc(4,i)=gfls_nc(4,i)-cbpa*p_b2*SUM(ymus_nc(1,:,i)*uthetas_nc(:,2,i))
  qfls_nc(4,i)=qfls_nc(4,i)-cbpaq*p_b2*SUM(ymus_nc(2,:,i)*uthetas_nc(:,2,i))

  !Transport fluxes from fast ions
  IF(l_fast_nc) THEN

    gfls_nc(5,i)=gfls_nc(5,i)-cbpa*p_b2*SUM(ymus_nc(1,:,i)*uthetas_nc(:,3,i))
    qfls_nc(5,i)=qfls_nc(5,i)-cbpaq*p_b2*SUM(ymus_nc(2,:,i)*uthetas_nc(:,3,i))

  ENDIF

  !Pfirsch-Schluter and classical fluxes
  !Test particle               
  cpsa=cps*(xis_nc(i)/iz)
  cpsaq=cpsa*(z_j7kv*temp_i(im))
  ccla=ccl*(xis_nc(i)/iz)
  cclaq=ccla*(z_j7kv*temp_i(im))

  IF(l_pfirsch_nc) THEN

    !Pfirsch-Schluter, thermal forces
    gfls_nc(2,i)=gfls_nc(2,i)-cpsa*SUM(calmi_nc(1,:,im)*srcth(:,i))
    qfls_nc(2,i)=qfls_nc(2,i)-cpsaq*SUM(calmi_nc(2,:,im)*srcth(:,i))

    !Pfirsch-Schluter, fast ion forces
    IF(l_fast_nc) THEN

      gfls_nc(5,i)=gfls_nc(5,i)+cpsa*srcf(1,im)
      qfls_nc(5,i)=qfls_nc(5,i)+cpsaq*srcf(2,im)

      DO k=1,k_order_nc

        gfls_nc(5,i)=gfls_nc(5,i)-cpsa*srcth(k,i)*SUM(calmif_nc(1,k,im,:))
        qfls_nc(5,i)=qfls_nc(5,i)-cpsaq*srcth(k,i)*SUM(calmif_nc(2,k,im,:))

      ENDDO

    ENDIF

    !Unit p'/p and T'/T 
    dpss_nc(i,i)=dpss_nc(i,i)-cpsa*calmi_nc(1,1,im)*srcthp(i)
    dtss_nc(i,i)=dtss_nc(i,i)-cpsa*calmi_nc(1,2,im)*srctht(i)
    chipss_nc(i,i)=chipss_nc(i,i)-cpsaq*calmi_nc(2,1,im)*srcthp(i)
    chitss_nc(i,i)=chitss_nc(i,i)-cpsaq*calmi_nc(2,2,im)*srctht(i)

    !Field particle 
    DO j=1,ms_nc !Over species 2

      jm=jms_nc(j)
      cpsb=cpsa*xis_nc(j)
      cpsbq=cpsb*(z_j7kv*temp_i(im))
      gfls_nc(2,i)=gfls_nc(2,i)-cpsb*SUM(calnii_nc(1,:,im,jm)*srcth(:,j))
      qfls_nc(2,i)=qfls_nc(2,i)-cpsbq*SUM(calnii_nc(2,:,im,jm)*srcth(:,j))

      !Unit p'/p and T'/T 
      dpss_nc(j,i)=dpss_nc(j,i)-cpsb*calnii_nc(1,1,im,jm)*srcthp(j)
      dtss_nc(j,i)=dtss_nc(j,i)-cpsb*calnii_nc(1,2,im,jm)*srctht(j)
      chipss_nc(j,i)=chipss_nc(j,i)-cpsbq*calnii_nc(2,1,im,jm)*srcthp(j)
      chitss_nc(j,i)=chitss_nc(j,i)-cpsbq*calnii_nc(2,2,im,jm)*srctht(j)

    ENDDO !Over species 2

  ENDIF

  IF(l_classical_nc) THEN

    !Classical 
    gfls_nc(3,i)=gfls_nc(3,i)+ccla*SUM(calmi_nc(1,:,im)*srcth(:,i))
    qfls_nc(3,i)=qfls_nc(3,i)+cclaq*SUM(calmi_nc(2,:,im)*srcth(:,i))

    !Unit p'/p and T'/T 
    dpss_nc(i,i)=dpss_nc(i,i)+ccla*calmi_nc(1,1,im)*srcthp(i)
    dtss_nc(i,i)=dtss_nc(i,i)+ccla*calmi_nc(1,2,im)*srctht(i)
    chipss_nc(i,i)=chipss_nc(i,i)+cclaq*calmi_nc(2,1,im)*srcthp(i)
    chitss_nc(i,i)=chitss_nc(i,i)+cclaq*calmi_nc(2,2,im)*srctht(i)

   !Field particle 
    DO j=1,ms_nc !Over species 2

      jm=jms_nc(j)
      cclb=ccla*xis_nc(j)
      cclbq=cclb*(z_j7kv*temp_i(im))
      gfls_nc(3,i)=gfls_nc(3,i)+cclb*SUM(calnii_nc(1,:,im,jm)*srcth(:,j))
      qfls_nc(3,i)=qfls_nc(3,i)+cclbq*SUM(calnii_nc(2,:,im,jm)*srcth(:,j))

      !Unit p'/p and T'/T 
      dpss_nc(j,i)=dpss_nc(j,i)+cclb*calnii_nc(1,1,im,jm)*srcthp(j)
      dtss_nc(j,i)=dtss_nc(j,i)+cclb*calnii_nc(1,2,im,jm)*srctht(j)
      chipss_nc(j,i)=chipss_nc(j,i)+cclbq*calnii_nc(2,1,im,jm)*srcthp(j)
      chitss_nc(j,i)=chitss_nc(j,i)+cclbq*calnii_nc(2,2,im,jm)*srctht(j)

    ENDDO !Over species 2

  ENDIF

ENDDO !Over species 1

!Electrical resistivity
petap_nc=one/pjboh_nc

!Add <E.B> normalization
pjboh_nc=p_eb*pjboh_nc
upars_nc(:,2,:)=p_eb*upars_nc(:,2,:)
uthetas_nc(:,2,:)=p_eb*uthetas_nc(:,2,:)
gfls_nc(4,:)=p_eb*gfls_nc(4,:)
qfls_nc(4,:)=p_eb*qfls_nc(4,:)

!-------------------------------------------------------------------------------
!Convert full coefficient matrices to diffusivities and conductivities
!-------------------------------------------------------------------------------
DO i=1,ms_nc !Over species

  im=jms_nc(i)
  iza=IABS(jzs_nc(i))
  dpss_nc(:,i)=dpss_nc(:,i)/den_iz(im,iza)
  dtss_nc(:,i)=dtss_nc(:,i)/den_iz(im,iza)
  chipss_nc(:,i)=chipss_nc(:,i)/den_iz(im,iza)/temp_i(im)/z_j7kv
  chitss_nc(:,i)=chitss_nc(:,i)/den_iz(im,iza)/temp_i(im)/z_j7kv

ENDDO !Over species

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE NCLASS_FLOW

SUBROUTINE NCLASS_FRICTION(testm,fieldn,xab,tab_in)
!--------------------------------------------------------------------------------
!NCLASS_FRICTION computes the Laguerre order 3/2 matrix elements of the
!  Coulomb collision operator for the l=1 spherical harmonic (friction) to all
!  orders in the Laguerre polynomial energy expansion
!
!References:
!  S.P.Hirshman, original code 9/2004
!  W.A.Houlberg, S.P.Hirshman, adapted for NCLASS 8/2006
!
!Comments:
!  Hij(L)=4*Gamma[1/2(L+4)]/sqrt(pi)    [Ref: Eq.(18)]
!         We never need to compute this explicitly, computed recursively from:
!         Hij(L+2)=(1/2)(L+4)*Hij(L), Hij(L=1)=3
!
!  ilag32  Laguerre "i" index (=2*i1, i1=raw index in memo)
!  jlag32  Laguerre "j" index (=2*j1, j1=raw index in memo)
!  yab2    1/xab**2
!  n_tilda memo Eq. 22d
!  delta_n memo Eq. 22e
!  source  x**(j+3)/(1+ x**2)**(-(L+2)/2)
!
!  gmat0(L,xab)=Gij(L,0,xab), for L odd, memo Eq. 21
!  gmat(i,j,xab) for i= *i1+1, =Gij(i,  j,  xab)
!                for i=2*i1  , =Gij(i-1,j,1/xab)
!  onepx(L)=(1+xab*xab)**(-L/2), L=1,2,....
!  onepxi(L)=[1+(xab*xab)-1]**(-L/2)
!  Store gmat(i,j,xab) at odd i locations      (i=i1+1)
!  Store gmat(i,j,xab**-1) at even i locations (i=i1)
!----------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=dpspec), INTENT(IN) :: &
  xab,                 & !va/vb, ratio of thermal velocities of particles a,b
  tab_in                 !Ta/Tb, ratio of temperatures of particles a,b

!Declaration of output variables
REAL(KIND=dpspec), DIMENSION(:,:), INTENT(OUT) :: &
  testm,               & !2D array of test-particle matrix elements to be computed
  fieldn                 !2D array of field-particle matrix elements to be computed

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i1,j1,ilag32,jlag32,llag32

REAL(KIND=dpspec) :: &
  onepx(1:4*k_order_nc+1), &
  onepxi(1:4*k_order_nc+1), &
  gmat0(1:4*k_order_nc-1), &
  gmat(0:2*k_order_nc+1,0:2*k_order_nc+2)

REAL(KIND=dpspec) :: &
  rawm(0:k_order_nc-1,0:k_order_nc), &
  rawn(0:k_order_nc-1,0:k_order_nc)

REAL(KIND=dpspec) :: &
  delta_n,gmat2,L00,n_tilda,source,t1,tab,yab2
         
!--------------------------------------------------------------------------------
!Local variables
!--------------------------------------------------------------------------------
onepx(1)=onedp/SQRT(onedp+xab*xab)
onepxi(1)=xab*onepx(1)

tab=tab_in
IF(tab == zerodp) tab=EPSILON(tab)

gmat(:,:)=zerodp
gmat(1,0)=xab*onepx(1)      !g(1,0,x)
gmat(0,0)=onepx(1)          !g(1,0,x**-1)
gmat0(1)=gmat(1,0)

t1=TINY(t1)/onepx(1)

DO llag32=2,4*k_order_nc+1

  IF(onepx(llag32-1) > t1) onepx(llag32)=onepx(llag32-1)*onepx(1)
  onepxi(llag32)=onepxi(llag32-1)*onepxi(1)

END DO

!--------------------------------------------------------------------------------
!Initialization loop
!--------------------------------------------------------------------------------
!Eq (21a) of memo, for j=0, L-odd only
DO llag32=1,4*k_order_nc-3,2

  t1=xab*onepx(llag32+2)
  gmat0(llag32+2)=((llag32+1)*gmat0(llag32)+t1)/(llag32+2)

ENDDO

!Eq (21a,b) of memo
!To prevent overflow, write x^(j+1)(1+x**2)^[(i+j+2)/2] in terms of onepx, onepxi
DO i1=0,2*(k_order_nc-1),2

  ilag32=i1+1

  DO j1=0,2*(k_order_nc-1),2

    llag32=ilag32+j1+2
    t1=onepx(ilag32+1)*onepxi(j1+1)

    !gmat(i,j,x) stored at ilag32
    gmat(ilag32+2,j1)=((ilag32+1)*gmat(ilag32,j1)+t1)/llag32
    gmat(ilag32,j1+2) = ((j1+1)*gmat(ilag32,j1)-t1)/llag32

    !gmat(i,j,x**-1) stored at i1
    t1=onepxi(ilag32+1)*onepx(j1+1)
    gmat(i1+2,j1)=((ilag32+1)*gmat(i1,j1)+t1)/llag32
    gmat(i1,j1+2)=((j1+1)*gmat(i1,j1)-t1)/llag32

  ENDDO

ENDDO

!Now, llag32=ilag32+jlag32+1
IF(xab == zerodp) THEN

  yab2=onedp/EPSILON(xab)**2

ELSE

  yab2=onedp/xab**2

ENDIF

!--------------------------------------------------------------------------------
!Raw test and field particle coefficients
!--------------------------------------------------------------------------------
DO i1=0,k_order_nc-1

  ilag32=2*i1

  DO j1=0,k_order_nc-1

    jlag32=2*j1
    llag32=ilag32+jlag32+1

    !Test particles, Eq. (22a) of memo
    source=xab*onepx(llag32+2)
    rawm(i1,j1)=gmat0(llag32)
    IF(llag32 /= 1) rawm(i1,j1)=rawm(i1,j1)+(ilag32+jlag32+ilag32*jlag32) &
                                            *(gmat0(llag32)-xab*onepx(llag32)) &
                                            *yab2/(llag32-1)
    gmat2=gmat0(llag32)-source
    rawm(i1,j1)=rawm(i1,j1)+(tab-onedp)*(ilag32+1)*yab2*gmat2
    rawm(i1,j1)=rawm(i1,j1)/(llag32+2) !*Hij(llag32)
        
    !Field particles, Eq. (22b) of memo
    source=xab*onepx(ilag32+1)*onepxi(jlag32+2)            
    n_tilda=(jlag32+3)*gmat(ilag32+1,jlag32+2)-source
    rawn(i1,j1)=(onedp+jlag32*onedp/3)*gmat(ilag32+1,jlag32+2) &
                +yab2*(onedp/3+jlag32*onedp/5/tab)*n_tilda
    delta_n=n_tilda*(onedp+3*yab2)
    !i <-> j, x <-> xinv part, gmat(xinv) was stored at EVEN first index values
    IF(jlag32 == 0) THEN

      t1=onepx(ilag32+3)

    ELSE

      t1=onepx(ilag32+3)*onepxi(jlag32)

    ENDIF

    n_tilda=(ilag32+3)*gmat(jlag32,ilag32+2)-t1
    rawn(i1,j1)=rawn(i1,j1)+xab*((onedp+ilag32*onedp/3)*gmat(jlag32,ilag32+2) &
                           +xab**2*(onedp/3+ilag32*onedp/5/tab)*n_tilda)
    delta_n=delta_n+xab**3*n_tilda+source
    rawn(i1,j1)=(tab*rawn(i1,j1)+delta_n*(onedp-tab)/3) !*Hij(llag32)

      WRITE(36,200) i1,j1,rawm(i1,j1),i1,j1,rawn(i1,j1)
 200  FORMAT(' MRAW(',i3,',',i3,') = ', 1pe12.3, &
             ' NRAW(',i3,',',i3,') = ',1e12.3)

  ENDDO

ENDDO

!--------------------------------------------------------------------------------
!Laguerre test and field particle coefficients
!--------------------------------------------------------------------------------
!Note that indexing on testm and fieldn starts at 1, not 0!
testm(:,:)=zerodp
fieldn(:,:)=zerodp
rawm(:,k_order_nc)=zerodp
rawn(:,k_order_nc)=zerodp

DO i1=0,k_order_nc-1

  DO j1=0,k_order_nc-1

    L00=3*clag_nc(i1+1)*clag_nc(j1+1)
    testm(i1+1,j1+1)=L00*NCLASS_LAGUERRE(1,rawm,i1,j1)
    fieldn(i1+1,j1+1)=L00*NCLASS_LAGUERRE(1,rawn,i1,j1)

      WRITE(37,210) i1,j1,testm(i1+1,j1+1),i1,j1,fieldn(i1+1,j1+1),L00
 210  FORMAT(' TESTM(',i3,',',i3,') = ', 1pe12.3, &
             ' FIELDN(',i3,',',i3,') = ',1e12.3, &
             ' L00 = ',1e12.3)

  ENDDO

ENDDO

END SUBROUTINE NCLASS_FRICTION

SUBROUTINE NCLASS_INIT(m_i,m_z,amu_i,temp_i,den_iz, &
                       iflag,message)
!-------------------------------------------------------------------------------
!NCLASS_INIT initializes the species information and allocates arrays
!
!References:                                                     
!  W.A.Houlberg, F90 free format 8/2004
!                arbitrary order of velocity moments 7/2006
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  m_i,                 & !number of isotopes (> 1) [-]
  m_z                    !highest charge state [-]

REAL(KIND=rspec), INTENT(IN) :: &
  amu_i(:),            & !atomic mass number [-]
  temp_i(:),           & !temperature [keV]
  den_iz(:,:)            !density [/m**3]

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
  i,im,iza,j,k,nset1,nset2,nset3

!-------------------------------------------------------------------------------
!Find number of significant isotopes and charge states, and max charge
!-------------------------------------------------------------------------------
mi_nc=0
mz_nc=0
ms_nc=0

DO i=1,m_i !Over isotopes

  DO j=1,m_z !Over charge states

    IF(den_iz(i,j) >= cden_nc) THEN

      ms_nc=ms_nc+1
      IF(i > mi_nc) mi_nc=i
      IF(j > mz_nc) mz_nc=j

    ENDIF

  ENDDO !Over charge states

ENDDO !Over isotopes

!At least two species
IF(mi_nc < 2) THEN

  iflag=1
  message='NCLASS_INIT(1)/ERROR:m_i must be >= 2'
  GOTO 9999

ENDIF

!Highest charge state at least 1
IF(mz_nc < 1) THEN

  iflag=1
  message='NCLASS_INIT(2)/ERROR:m_z must be >= 1'
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Allocate or reallocate moments arrays
!-------------------------------------------------------------------------------
!Note we start the Laguerre coefficients at index 1 instead of 0
IF(ALLOCATED(clag_nc)) THEN

  nset1=SIZE(clag_nc)

  !If storage requirements have increased, reallocate
  IF(nset1 /= k_order_nc) THEN

    !Reallocate
    DEALLOCATE(clag_nc)
    ALLOCATE(clag_nc(k_order_nc))

    clag_nc(1)=onedp
    clag_nc(2:k_order_nc)=zerodp

    DO i=2,k_order_nc

      clag_nc(i)=REAL(2*i+1,dpspec)*clag_nc(i-1)/REAL(2*(i-1),dpspec)

    ENDDO

  ENDIF

ELSE

  !Allocate
  ALLOCATE(clag_nc(k_order_nc))

  clag_nc(1)=onedp
  clag_nc(2:k_order_nc)=zerodp

  DO i=2,k_order_nc

    clag_nc(i)=REAL(2*i+1,dpspec)*clag_nc(i-1)/REAL(2*(i-1),dpspec)

  ENDDO

ENDIF

!-------------------------------------------------------------------------------
!Allocate or reallocate isotope arrays
!-------------------------------------------------------------------------------
IF(ALLOCATED(calmi_nc)) THEN

  nset1=SIZE(calmi_nc,1)
  nset2=SIZE(calmi_nc,3)

  !If storage requirements have increased, reallocate
  IF(nset1 /= k_order_nc .OR. &
     nset2 /= mi_nc) THEN

    !Reallocate
    DEALLOCATE(vti_nc, &
               amntii_nc, &
               calmi_nc, &
               calnii_nc, &
               capmii_nc, &
               capnii_nc)
    ALLOCATE(vti_nc(mi_nc), &
             amntii_nc(mi_nc,mi_nc), &
             calmi_nc(k_order_nc,k_order_nc,mi_nc), &
             calnii_nc(k_order_nc,k_order_nc,mi_nc,mi_nc), &
             capmii_nc(k_order_nc,k_order_nc,mi_nc,mi_nc), &
             capnii_nc(k_order_nc,k_order_nc,mi_nc,mi_nc))

  ENDIF

ELSE

  !Allocate
  ALLOCATE(vti_nc(mi_nc), &
           amntii_nc(mi_nc,mi_nc), &
           calmi_nc(k_order_nc,k_order_nc,mi_nc), &
           calnii_nc(k_order_nc,k_order_nc,mi_nc,mi_nc), &
           capmii_nc(k_order_nc,k_order_nc,mi_nc,mi_nc), &
           capnii_nc(k_order_nc,k_order_nc,mi_nc,mi_nc))

ENDIF

vti_nc(:)=zero
amntii_nc(:,:)=zero
calmi_nc(:,:,:)=zero
calnii_nc(:,:,:,:)=zero
capmii_nc(:,:,:,:)=zero
capnii_nc(:,:,:,:)=zero

!-------------------------------------------------------------------------------
!Allocate or reallocate isotope/charge state arrays
!-------------------------------------------------------------------------------
IF(ALLOCATED(grppiz_nc)) THEN

  nset1=SIZE(grppiz_nc,1)
  nset2=SIZE(grppiz_nc,2)

  !If storage requirements have increased, reallocate
  IF(nset1 /= mi_nc .OR. &
     nset2 /= mz_nc) THEN

    !Reallocate
    DEALLOCATE(grppiz_nc)
    ALLOCATE(grppiz_nc(mi_nc,mz_nc))

  ENDIF

ELSE

  !Allocate
  ALLOCATE(grppiz_nc(mi_nc,mz_nc))

ENDIF

grppiz_nc(:,:)=zero

!-------------------------------------------------------------------------------
!Allocate or reallocate species arrays
!-------------------------------------------------------------------------------
IF(ALLOCATED(jms_nc)) THEN

  nset1=SIZE(jms_nc)

  !If storage requirements have increased, reallocate
  IF(nset1 /= ms_nc) THEN

    !Reallocate
    DEALLOCATE(jms_nc, &
               jzs_nc, &
               sqzs_nc,&
               xis_nc, &
               bsjbps_nc, &
               bsjbts_nc, &
               gfls_nc,&
               dpss_nc,&
               dtss_nc,&
               qfls_nc,&
               chipss_nc, &
               chitss_nc, &
               tauss_nc)

    ALLOCATE(jms_nc(ms_nc), &
             jzs_nc(ms_nc), &
             sqzs_nc(ms_nc), &
             xis_nc(ms_nc), &
             bsjbps_nc(ms_nc), &
             bsjbts_nc(ms_nc), &
             gfls_nc(5,ms_nc), &
             dpss_nc(ms_nc,ms_nc), &
             dtss_nc(ms_nc,ms_nc), &
             qfls_nc(5,ms_nc), &
             chipss_nc(ms_nc,ms_nc), &
             chitss_nc(ms_nc,ms_nc), &
             tauss_nc(ms_nc,ms_nc))

    IF(l_fast_nc) THEN

      DEALLOCATE(vc3s_nc)
      ALLOCATE(vc3s_nc(ms_nc))

        vc3s_nc(:)=zero

    ENDIF

  ENDIF

ELSE

  !Allocate
  ALLOCATE(jms_nc(ms_nc), &
           jzs_nc(ms_nc), &
           sqzs_nc(ms_nc), &
           xis_nc(ms_nc), &
           bsjbps_nc(ms_nc), &
           bsjbts_nc(ms_nc), &
           gfls_nc(5,ms_nc), &
           dpss_nc(ms_nc,ms_nc), &
           dtss_nc(ms_nc,ms_nc), &
           qfls_nc(5,ms_nc), &
           chipss_nc(ms_nc,ms_nc), &
           chitss_nc(ms_nc,ms_nc), &
           tauss_nc(ms_nc,ms_nc))

  IF(l_fast_nc) THEN

    ALLOCATE(vc3s_nc(ms_nc))

      vc3s_nc(:)=zero

  ENDIF

ENDIF

jms_nc(:)=0
jzs_nc(:)=0
sqzs_nc(:)=zero
xis_nc(:)=zero
bsjbps_nc(:)=zero
bsjbts_nc(:)=zero
gfls_nc(:,:)=zero
dpss_nc(:,:)=zero
dtss_nc(:,:)=zero
qfls_nc(:,:)=zero
chipss_nc(:,:)=zero
chitss_nc(:,:)=zero
tauss_nc(:,:)=zero

!-------------------------------------------------------------------------------
!Allocate or reallocate species/moments arrays
!-------------------------------------------------------------------------------
IF(ALLOCATED(ymus_nc)) THEN

  nset1=SIZE(ymus_nc,1)
  nset2=SIZE(ymus_nc,3)

  !If storage requirements have increased, reallocate
  IF(nset1 /= ms_nc .OR. &
     nset2 /= k_order_nc) THEN

    !Reallocate
    DEALLOCATE(ymus_nc,&
               upars_nc, &
               uthetas_nc)

    ALLOCATE(ymus_nc(k_order_nc,k_order_nc,ms_nc), &
             upars_nc(k_order_nc,3,ms_nc), &
             uthetas_nc(k_order_nc,3,ms_nc))

  ENDIF

ELSE

  !Allocate
  ALLOCATE(ymus_nc(k_order_nc,k_order_nc,ms_nc), &
           upars_nc(k_order_nc,3,ms_nc), &
           uthetas_nc(k_order_nc,3,ms_nc))

ENDIF

ymus_nc(:,:,:)=zero
upars_nc(:,:,:)=zero
uthetas_nc(:,:,:)=zero

!-------------------------------------------------------------------------------
!Allocate or reallocate fast ion arrays
!-------------------------------------------------------------------------------
IF(l_fast_nc) THEN

  IF(ALLOCATED(ef_nc)) THEN

    nset1=SIZE(ef_nc,1)

    !If storage requirements have increased, reallocate
    IF(nset1 /= mf_nc) THEN

      !Reallocate
      DEALLOCATE(amuf_nc, &
                 zf_nc,&
                 ef_nc,&
                 uf_nc, &
                 denf_nc, &
                 tausf_nc)

      ALLOCATE(amuf_nc(mf_nc), &
               zf_nc(mf_nc), &
               ef_nc(mf_nc), &
               uf_nc(mf_nc), &
               denf_nc(mf_nc), &
               tausf_nc(mf_nc))

    ENDIF

  ELSE

    !Allocate
    ALLOCATE(amuf_nc(mf_nc), &
             zf_nc(mf_nc), &
             ef_nc(mf_nc), &
             uf_nc(mf_nc), &
             denf_nc(mf_nc), &
             tausf_nc(mf_nc))

  ENDIF

  amuf_nc(:)=zero
  zf_nc(:)=zero
  ef_nc(:)=zero
  uf_nc(:)=zero
  denf_nc(:)=zero
  tausf_nc(:)=zero

  IF(ALLOCATED(calmif_nc)) THEN

    nset1=SIZE(calmif_nc,1)
    nset2=SIZE(calmif_nc,3)
    nset3=SIZE(calmif_nc,4)

    !If storage requirements have increased, reallocate
    IF(nset1 /= k_order_nc .OR. &
       nset2 /= mi_nc .OR. &
       nset3 /= mf_nc) THEN

      !Reallocate
      DEALLOCATE(calmif_nc)

      ALLOCATE(calmif_nc(k_order_nc,k_order_nc,mi_nc,mf_nc))

    ENDIF

  ELSE

    !Allocate
    ALLOCATE(calmif_nc(k_order_nc,k_order_nc,mi_nc,mf_nc))

  ENDIF

  calmif_nc(:,:,:,:)=zero

ENDIF

!-------------------------------------------------------------------------------
!Map internal species to external species and identify electrons
!-------------------------------------------------------------------------------
k=0
imel_nc=0
jsel_nc=0

DO i=1,mi_nc !Over isotopes

  DO j=1,mz_nc !Over charge states

    IF(den_iz(i,j) >= cden_nc) THEN

      !Set isotope number and charge state for this species
      k=k+1
      jms_nc(k)=i

      IF(amu_i(i) < 0.5_rspec) THEN

        !Electrons
        jzs_nc(k)=-j
        imel_nc=i
        jsel_nc=k

      ELSE

        !Ions
        jzs_nc(k)=j

      ENDIF

    ENDIF

  ENDDO !Over charge states

ENDDO !Over isotopes

!In case there are no electrons, use the first ion species
IF(imel_nc == 0) imel_nc=1

!-------------------------------------------------------------------------------
!Thermal velocities
!-------------------------------------------------------------------------------
vti_nc(:)=SQRT(2*z_j7kv*temp_i(:)/amu_i(:)/z_pmass)

!-------------------------------------------------------------------------------
!Fast ion effects
!-------------------------------------------------------------------------------
IF(l_fast_nc) THEN

!Critical velocities
  vc3s_nc(:)=zero

  DO j=1,ms_nc
 
    im=jms_nc(j)
    iza=jzs_nc(j)

    IF(j /= jsel_nc) THEN

      vc3s_nc(j)=3*SQRT(z_pi)/4*(amu_i(imel_nc)/amu_i(im)) &
                 *(den_iz(im,iza)*iza**2/den_iz(imel_nc,1))*vti_nc(imel_nc)**3

    ENDIF

  ENDDO

  vc3_nc=SUM(vc3s_nc(:))

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE NCLASS_INIT

SUBROUTINE NCLASS_K(p_fm,p_ft,p_ngrth,x,amu_i,temp_i, &
                    ykb_s,ykp_s,ykpo_s,ykpop_s)
!-------------------------------------------------------------------------------
!NCLASS_K calculates the velocity-dependent neoclassical viscosity coefficients
!
!References:                                                      
!  S.P.Hirshman, D.J.Sigmar, Nucl Fusion 21 (1981) 1079
!  C.E.Kessel, Nucl Fusion 34 (1994) 1221                              
!  K.C.Shaing, M.Yokoyama, M.Wakatani, C.T.Hsu, Phys Plasmas 3 (1996) 965
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  p_fm(:),             & !poloidal moments of drift factor for PS [/m**2]
  p_ft,                & !trapped fraction [-]
  p_ngrth,             & !<n.grad(Theta)> [/m]
  x,                   & !normalized velocity v/(2kT/m)**0.5 [-]
  amu_i(:),            & !atomic mass number [-]
  temp_i(:)              !temperature [keV]

!Declaration of output variables
REAL(KIND=rspec) :: &
  ykb_s(:),            & !banana viscosity for s [/s]
  ykp_s(:),            & !Pfirsch-Schluter viscosity for s [/s]
  ykpo_s(:),           & !potato viscosity for s [/s]
  ykpop_s(:)             !potato plateau viscosity for s [/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,im,iz

REAL(KIND=rspec) :: &
  c1,c2,c3,c4

REAL(KIND=rspec) :: &
  ynud_s(ms_nc),ynut_s(ms_nc),ynuti_s(mfm_nc,ms_nc)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
ykb_s(:)=zero
ykp_s(:)=zero
ykpo_s(:)=zero
ykpop_s(:)=zero

!Get collisional frequencies
CALL NCLASS_NU(p_ngrth,x,temp_i,ynud_s,ynut_s,ynuti_s)

!-------------------------------------------------------------------------------
!Set velocity dependent viscosities (K's)
!-------------------------------------------------------------------------------
c1=3*x**2/2
c2=3*2.19_rspec/(2.0**1.5)*x**0.33333_rspec

IF(l_potato_nc) c3=3*z_pi/(64*2.0**0.33333_rspec)/ABS(cpotl_nc)

DO i=1,ms_nc !Over species

  im=jms_nc(i)
  iz=jzs_nc(i)

  !Banana
  IF(l_banana_nc) THEN

    !Provide cutoff to eliminate failure at unity trapped fraction
    !At A=>1 viscosity will go over to Pfirsch-Schluter value
    c4=one-p_ft
    IF(c4 < 1.0e-3_rspec) c4=1.0e-3_rspec
    ykb_s(i)=p_ft/c4/sqzs_nc(i)**1.5*ynud_s(i)

  ENDIF

  !Pfirsch-Schuter
  IF(l_pfirsch_nc) THEN

    ykp_s(i)=ykp_s(i)+c1*vti_nc(im)**2*SUM(p_fm(:)*ynuti_s(:,i))/ynut_s(i)

  ENDIF

  !Potato
  IF(l_potato_nc) THEN

    c4=ABS(amu_i(im)*z_pmass*vti_nc(im)/(iz*z_coulomb*cpotb_nc*cpotl_nc))
    ykpo_s(i)=c2*c4**(0.333333_rspec)*ynud_s(i)/sqzs_nc(i)**(1.666667_rspec)
    ykpop_s(i)=c3*vti_nc(im)*c4**(1.333333_rspec)

  ENDIF

ENDDO !Over species

END SUBROUTINE NCLASS_K

FUNCTION NCLASS_LAGUERRE(k_sum,x,i1,j1)
!-------------------------------------------------------------------------------
!NCLASS_LAGUERRE calculates the Laguerre moments from the velocity moments
!
!References:
!  S.P.Hirshman 9/2004
!  W.A.Houlberg, S.P.Hirshman adapted for NCLASS module 7/2006
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_sum,               & !option for reduced coefficients
                         !=0 off, for viscosity use Lagurre sum
                         !=1 on, for friction use weight0d Laguerre sum
  i1,                  & !first index of Laguerre matrix element
  j1                     !second index of Laguerre matrix element

REAL(KIND=dpspec), INTENT(IN) :: &
  x(:,:)                 !raw velocity moments

!Declaration of output variables
REAL(KIND=dpspec) :: &
  NCLASS_LAGUERRE

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j

REAL(KIND=dpspec) :: &
  partsumi,ratio

REAL(KIND=dpspec) :: &
  partsum(0:i1)

!-------------------------------------------------------------------------------
!Compute i partial sums
!-------------------------------------------------------------------------------
DO i=0,i1

  partsumi=x(i+1,j1+1)

  DO j=j1,1,-1

    ratio=2*REAL(j-j1-1,dpspec)/REAL(j*(2*j+3),dpspec)
    IF(k_sum == 1) ratio=ratio*(REAL(i+j,dpspec)+1.5_dpspec)
    partsumi=ratio*partsumi+x(i+1,j)

  ENDDO

  partsum(i)=partsumi

ENDDO

!Sum over partial sums, using a partial sum method

partsumi=partsum(i1)

DO i=i1,1,-1

  ratio=REAL(i-i1-1,dpspec)/REAL(i,dpspec)
  IF(k_sum == 0) ratio=ratio/(REAL(i,dpspec)+1.5_dpspec)
  partsumi=ratio*partsumi+partsum(i-1)

ENDDO

NCLASS_LAGUERRE=partsumi

END FUNCTION NCLASS_LAGUERRE

SUBROUTINE NCLASS_MN(amu_i,temp_i,den_iz)
!-------------------------------------------------------------------------------
!NCLASS_MN calculates neoclassical friction coefficients
!
!References:                                         
!  S.P.Hirshman, D.J.Sigmar, Nucl Fusion 21 (1981) 1079
!  S.P.Hirshman, Phys Fluids 20 (1977) 589
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!                analytic friction coefficients for all Ta/Tb 3/2006
!                recursive generation of thermal friction matrix 9/2006
!
!Comments:                                                     
!  The k_order*korder matrix of test particle (M) and field particle (N)
!    coefficients of the collision operator use the Laguerre polynomials of
!    order 3/2 as basis functions for each isotopic combination
!  The indices on the M and N matrices are one greater than the notation in the
!    H&S review article so as to avoid 0 as an index
!New expressions implemented for inclusion of all values of Ta/Tb.
!  For comparison see:
!  Eqn 4.11 (HS81) for M00
!  Eqn 4.12 (HS81) for M01
!  Eqn 4.8 (HS81) for M10
!  Eqn 4.13 (HS81) for M11
!  Eqn 4.15 (HS81) for M02
!  Eqn 4.16 (HS81) for M12
!  Eqn 4.8 (HS81) for M20
!  Eqn 4.8 (HS81) for M21
!  Eqn 5.21 (HS81) for M22
!  Eqn 4.11 (HS81) for N00
!  Eqn 4.9 and 4.12 (HS81) for N01
!  Eqn 4.12 (HS81) for N10
!  Eqn 4.14 (HS81) for N11
!  Eqn 4.15 for N02 (HS81)
!  Eqn 4.17 for N12 (HS81)
!  Momentum conservation for N20 (HS81)
!  Eqn 4.9 and 4.17 for N21 (HS81)
!  Eqn 5.22 for N22 (HS81) 
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  amu_i(:),            & !atomic mass number [-]
  temp_i(:),           & !temperature [keV]
  den_iz(:,:)            !density [/m**3]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j,k,im,jm

REAL(KIND=rspec) :: &
  denzi2m,tauef,t,tpv2,v,v2,v3,v5,y32,y52,y72,y92,y112

REAL(KIND=dpspec) :: &
  tdp,vdp,mdp(1:k_order_nc,1:k_order_nc),ndp(1:k_order_nc,1:k_order_nc)

!-------------------------------------------------------------------------------
!Use recursion for k_order > 3, and analytic for lower values
!-------------------------------------------------------------------------------
!Recursive thermal friction matrix - good up to ~15 moments in double precision
IF(k_order_nc > 3) THEN

  DO im=1,mi_nc !Over isotopes a

    DO jm=1,mi_nc !Over isotopes b

      !initialize matrix elements and set parameters for this interaction pair
      mdp(:,:)=zerodp                  !Test particle matrix
      ndp(:,:)=zerodp                  !Field particle matrix
      tdp=temp_i(im)/temp_i(jm)        !Temperature ratio, Ta/Tb
      vdp=vti_nc(im)/vti_nc(jm)        !Thermal velocity ratio, vta/vtb

      !Get matrix elements for this interaction pair
      CALL NCLASS_FRICTION(mdp,ndp,vdp,tdp)

      !Change signs to adhere to convention in NCLASS
      DO i=1,k_order_nc !Over rows

        DO j=1,k_order_nc !Over columns

          IF(MOD(i+j,2) > 0) THEN
            
            ndp(i,j)=-ndp(i,j) !Even

          ELSE

            mdp(i,j)=-mdp(i,j) !Odd

          ENDIF

        ENDDO !Over columns

      ENDDO !Over rows

      !Load this interaction pair into full thermal friction matrices
      capmii_nc(:,:,im,jm)=mdp(:,:)
      capnii_nc(:,:,im,jm)=ndp(:,:)

    ENDDO !Over isotopes b

  ENDDO !Over isotopes a

ELSE

  !Analytic thermal friction matrix - up to 3 moments
  DO im=1,mi_nc !Over isotopes a

    DO jm=1,mi_nc !Over isotopes b

      !Set parameters
      t=temp_i(im)/temp_i(jm)          !Temperature ratio, Ta/Tb
      v=vti_nc(im)/vti_nc(jm)          !Thermal velocity ratio, vta/vtb
      v2=v**2
      v3=v2*v
      v5=v2*v3
      tpv2=t+v2
      y32=(one+v2)*SQRT(one+v2)
      y52=(one+v2)*y32
      y72=(one+v2)*y52
      y92=(one+v2)*y72
      y112=(one+v2)*y92

      !Test particle coefficients, M
      capmii_nc(1,1,im,jm)=-v/y32*tpv2
      capmii_nc(1,2,im,jm)=3*v3/2/y52*tpv2
      capmii_nc(2,1,im,jm)=capmii_nc(1,2,im,jm)-2*(t-one)*v/y52*(5*one/2+v2)
      capmii_nc(2,2,im,jm)=-v/y72*(tpv2*(15*one/2+4*v2+13*v2**2/4) &
                           -(t-one)*v2*(11*one/2+v2))

      IF(k_order_nc == 3) THEN

        capmii_nc(1,3,im,jm)=-15*v5/8/y72*tpv2
        capmii_nc(2,3,im,jm)=3*v3/y92*(tpv2*(21*one/4+2*v2+23*v2**2/16) &
                             -v2/4*(t-one)*(17*one/2+v2))
        capmii_nc(3,1,im,jm)=capmii_nc(1,3,im,jm)+3*v3/y72*(t-one)*(7*one/2+v2)
        capmii_nc(3,2,im,jm)=capmii_nc(2,3,im,jm)-v/y92*(t-one) &
                             *(35*one/2+21*v2+171*v2**2/8+19*v2**3/4)
        capmii_nc(3,3,im,jm)=-v/y112*(tpv2*(175*one/8+28*v2+459*v2**2/8 &
                             +17*v2**3+433*v2**4/64)-(t-one)*v2 &
                             *(119*one/4+39*v2/2+411*v2**2/16+25*v2**3/8))

      ENDIF 
    
      !Field particle coefficients, N
       capnii_nc(1,1,im,jm)=-capmii_nc(1,1,im,jm)
       capnii_nc(1,2,im,jm)=-one/v2*capmii_nc(1,2,im,jm)
       capnii_nc(2,1,im,jm)=-v3/y52*(1.5*tpv2+3*(t-one))
       capnii_nc(2,2,im,jm)=v3/y72*(27*tpv2/4+9*(t-one)/2)

      IF(k_order_nc == 3) THEN

        capnii_nc(1,3,im,jm)=-capmii_nc(1,3,im,jm)/v2**2
        capnii_nc(2,3,im,jm)=-v3/y92*(225*tpv2/16+45*(t-one)/8)
        capnii_nc(3,1,im,jm)=v2**2*capnii_nc(1,3,im,jm)+15*v5/2/y72*(t-one)
        capnii_nc(3,2,im,jm)=v2*capnii_nc(2,3,im,jm)-105*v5/8/y92*(t-one)
        capnii_nc(3,3,im,jm)=525*v5/16/y112*(5*tpv2/4+(t-one))

      ENDIF

    ENDDO !Over isotopes b
   
  ENDDO !Over isotopes a

ENDIF

!-------------------------------------------------------------------------------
!Collision times
!-------------------------------------------------------------------------------
CALL NCLASS_TAU(amu_i,temp_i,den_iz)

!-------------------------------------------------------------------------------
!Reduced friction coefficients for thermal species
!-------------------------------------------------------------------------------
calmi_nc(:,:,:)=zero
calnii_nc(:,:,:,:)=zero

DO im=1,mi_nc !Over isotopes

  DO j=1,k_order_nc !Over rows

    DO k=1,k_order_nc !Over columns

      !Test particle component
      calmi_nc(j,k,im)=SUM(amntii_nc(im,1:mi_nc)*capmii_nc(j,k,im,1:mi_nc))

      !Field particle component
      calnii_nc(j,k,im,1:mi_nc)=amntii_nc(im,1:mi_nc)*capnii_nc(j,k,im,1:mi_nc)

    ENDDO !Over columns

  ENDDO !Over rows

ENDDO !Over isotopes

!-------------------------------------------------------------------------------
!Reduced friction coefficients for fast ions
!-------------------------------------------------------------------------------
IF(l_fast_nc) THEN

  !Sum n*Z^2/m over all charge states of all ions
  denzi2m=zero

  DO im=1,mi_nc !Over isotopes

    IF(im /= imel_nc) THEN

      DO j=1,mz_nc !Over charge states

        denzi2m=denzi2m+j**2*den_iz(im,j)/amu_i(im)

      ENDDO !Over charge states

    ENDIF

  ENDDO !Over isotopes

  !Loop over fast ion energy components
  DO k=1,mf_nc !Over energy components

    tauef=tauss_nc(jsel_nc,jsel_nc)*den_iz(imel_nc,1)/(denf_nc(k)*zf_nc(k)**2)

    DO im=1,mi_nc !Over isotopes

      IF(im == imel_nc) THEN

        !Electrons
        calmif_nc(1:2,1:2,im,k)=z_pmass*amu_i(imel_nc)*den_iz(imel_nc,1)/tauef
        calmif_nc(1,1,im,k)=-calmif_nc(1,1,im,k)
        calmif_nc(1,2,im,k)=-1.5*calmif_nc(1,2,im,k)
        calmif_nc(2,1,im,k)=1.5*calmif_nc(2,1,im,k)
        calmif_nc(2,2,im,k)=3.25*calmif_nc(2,2,im,k)

      ELSE

        !Ions
!        gca=zero
!
!        DO j=1,mz_nc !Over charge states
!
!          gca=gca+REAL(j,rspec)**2*den_iz(im,j)
!
!        ENDDO !Over charge states
!
!        gca=gca/denzi2m/amu_i(im)
!        giv2=G_i*v_0^2/v_th^2
!        giv2=(2*ef_nc(k)*z_j7kv/z_pmass/amuf_nc(k))*gif_nc(k)/vti_nc(im)**2
!        calmif_nc(1:2,1,im,k)=gca*z_pmass*amu_i(imel_nc)/denf_nc(k)/3 &
!                              *den_iz(imel_nc,1)/tauef*dendotf_nc(k)*tausf_nc(k)
!        calmif_nc(1,1,im,k)=-(one+amu_i(im)/amuf_nc(k))*calmif_nc(1,1,im,k)
!        calmif_nc(2,1,im,k)=10*giv2/3*calmif_nc(2,1,im,k)
!        calnif_nc(1:2,im,k)=gca*z_pmass*amu_i(imel_nc)*den_iz(imel_nc,1) &
!                            /tauef
!        calnif_nc(1,im,k)=(one+amu_i(im)/amuf_nc(k))*calnif_nc(1,im,k)
!        calnif_nc(2,im,k)=-calnif_nc(2,im,k)

      ENDIF

    ENDDO !Over isotopes

  ENDDO !Over energy components

ENDIF

END SUBROUTINE NCLASS_MN

SUBROUTINE NCLASS_MU(p_fm,p_ft,p_ngrth,amu_i,temp_i,den_iz)
!-------------------------------------------------------------------------------
!NCLASS_MU calculates the matrix of neoclassical fluid viscosities
!
!References:                                             
!  K.C.Shaing, M.Yokoyama, M.Wakatani, C.T.Hsu, Phys Plasmas 3 (1996) 965
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!                recursive generation of viscosity matrix 9/2006
!
!Comments:
!  Integrates the velocity-dependent banana and Pfirsch-Schluter contributions
!    over velocity space
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  p_fm(:),             & !poloidal moments of drift factor for PS [/m**2]
  p_ft,                & !trapped fraction [-]
  p_ngrth,             & !<n.grad(Theta)> [/m]
  amu_i(:),            & !atomic mass number [-]
  temp_i(:),           & !temperature [keV]
  den_iz(:,:)            !density [/m**3]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  mpnts=13

REAL(KIND=rspec), PARAMETER :: &
  bmax=3.2_rspec

!Saved
INTEGER, SAVE :: &
  init=0

REAL(KIND=rspec), SAVE :: &
  c1,h

REAL(KIND=rspec), ALLOCATABLE, SAVE :: &
  x(:),w(:,:)

!Other
INTEGER :: &
  i,im,iza,k,l,m

REAL(KIND=rspec) :: &
  c2,x2,xk,xx

REAL(KIND=dpspec) :: &
  L00

REAL(KIND=rspec), DIMENSION(ms_nc) :: &
  ykb_s,ykp_s,ykpo_s,ykpop_s

REAL(KIND=dpspec), DIMENSION(k_order_nc,k_order_nc,ms_nc) :: &
  rymubp_s

REAL(KIND=rspec), DIMENSION(k_order_nc,k_order_nc,ms_nc) :: &
  ymubp_s

REAL(KIND=dpspec), ALLOCATABLE :: &
  rymub_s(:,:,:),rymupo_s(:,:,:),rymupopp_s(:,:,:)

REAL(KIND=rspec), ALLOCATABLE :: &
  ymub_s(:,:,:),ymupo_s(:,:,:),ymupopp_s(:,:,:)

LOGICAL :: &
  l_recursive=.FALSE.

!-------------------------------------------------------------------------------
!Recursive viscosity matrix
!-------------------------------------------------------------------------------
IF(l_recursive) THEN

  CALL NCLASS_VISCOSITY(p_ft)
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Initialization for numerical integration
!-------------------------------------------------------------------------------
IF(init == 0) THEN

  !Set integration points and weights
  ALLOCATE(x(mpnts), &
           w(mpnts,2*k_order_nc-1))
  x(:)=zero
  w(:,:)=zero
  h=bmax/(mpnts-1)

  DO m=2,mpnts !Over velocity nodes for integration

    x(m)=h*(m-1)
    x2=x(m)*x(m)
    w(m,1)=x2**2*EXP(-x2)

    DO k=2,2*k_order_nc-1 !Over polynomial moments

      w(m,k)=x2*w(m,k-1)

    ENDDO !Over polynomial moments

  ENDDO !Over velocity nodes for integration

  !Half weight on end point for trapezoidal integration
  w(mpnts,:)=w(mpnts,:)/2
  c1=8*h/3/SQRT(z_pi)

  !Allocate arrays for potato viscosities
  IF(l_potato_nc) THEN

    ALLOCATE(rymub_s(k_order_nc,k_order_nc,ms_nc), &
             rymupo_s(k_order_nc,k_order_nc,ms_nc), &
             rymupopp_s(k_order_nc,k_order_nc,ms_nc), &
             ymub_s(k_order_nc,k_order_nc,ms_nc), &
             ymupo_s(k_order_nc,k_order_nc,ms_nc), &
             ymupopp_s(k_order_nc,k_order_nc,ms_nc))

  ENDIF

  init=1

ENDIF

!Null out local viscosity arrays
rymubp_s(:,:,:)=zero
ymubp_s(:,:,:)=zero

IF(l_potato_nc) THEN

  ymub_s(:,:,:)=zero
  ymupo_s(:,:,:)=zero
  ymupopp_s(:,:,:)=zero
  rymub_s(:,:,:)=zero
  rymupo_s(:,:,:)=zero
  rymupopp_s(:,:,:)=zero

ENDIF

!-------------------------------------------------------------------------------
!Integrate over velocity space for fluid viscosities
!-------------------------------------------------------------------------------
IF(l_banana_nc .OR. &
   l_pfirsch_nc) THEN

  DO m=2,mpnts !Over velocity grid

    xx=x(m)
    !Get velocity-dependent k values        
    CALL NCLASS_K(p_fm,p_ft,p_ngrth,xx,amu_i,temp_i, &
                  ykb_s,ykp_s,ykpo_s,ykpop_s)

    DO i=1,ms_nc !Over species

      IF(.NOT. l_banana_nc) THEN

        !Only Pfirsch-Schluter
        xk=ykp_s(i)

      ELSEIF(.NOT. l_pfirsch_nc) THEN

        !Only banana
        xk=ykb_s(i)

      ELSE

        !Both banana and Pfirsch-Schluter
        xk=ykb_s(i)*ykp_s(i)/(ykb_s(i)+ykp_s(i))

      ENDIF

      DO k=1,k_order_nc

        DO l=1,k_order_nc

          rymubp_s(k,l,i)=rymubp_s(k,l,i)+xk*w(m,k+l-1)

        ENDDO

      ENDDO

      IF(l_potato_nc) THEN

        DO k=1,k_order_nc

          DO l=1,k_order_nc

            rymub_s(k,l,i)=rymub_s(k,l,i)+ykb_s(i)*w(m,k+l-1)
            rymupo_s(k,l,i)=ymupo_s(k,l,i)+ykpo_s(i)*w(m,k+l-1)
            xk=ykpo_s(i)*ykpop_s(i)/(ykpo_s(i)+ykpop_s(i))
            rymupopp_s(k,l,i)=ymupopp_s(k,l,i)+xk*w(m,k+l-1)

          ENDDO

        ENDDO

      ENDIF

    ENDDO  !Over species

  ENDDO !Over velocity grid

  !Convert from polynomial moments to Laguerre moments and scale
  DO i=1,ms_nc !Over species

    DO k=1,k_order_nc !Over rows

      DO l=1,k_order_nc !Over columns

        L00=(-1)**(k+l)*clag_nc(k)*clag_nc(l)
        ymubp_s(k,l,i)=L00*NCLASS_LAGUERRE(0,rymubp_s(:,:,i),k-1,l-1)

        IF(l_potato_nc) THEN

          ymub_s(k,l,i)=L00*NCLASS_LAGUERRE(0,rymub_s(:,:,i),k-1,l-1)
          ymupo_s(k,l,i)=L00*NCLASS_LAGUERRE(0,rymupo_s(:,:,i),k-1,l-1)
          ymupopp_s(k,l,i)=L00*NCLASS_LAGUERRE(0,rymupopp_s(:,:,i),k-1,l-1)

        ENDIF

      ENDDO

    ENDDO

  ENDDO !Over species

  SELECT CASE (l_potato_nc)

  CASE (.TRUE.)

    !Banana Pfirsch-Schluter plus potato potato-plateau
    ymus_nc(:,:,:)=(ymupo_s(:,:,:)**3*ymupopp_s(:,:,:) &
                   +ymub_s(:,:,:)**3*ymubp_s(:,:,:)) &
                   /(ymupo_s(:,:,:)**3+ymub_s(:,:,:)**3)

  CASE (.FALSE.)

    !Banana Pfirsch-Schluter
    ymus_nc(:,:,:)=ymubp_s(:,:,:)

  END SELECT

  !Add constants and mass density scale factor
  DO i=1,ms_nc

    im=jms_nc(i)
    iza=IABS(jzs_nc(i))
    c2=c1*den_iz(im,iza)*amu_i(im)*z_pmass
    ymus_nc(:,:,i)=c2*ymus_nc(:,:,i)

  ENDDO

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE NCLASS_MU

SUBROUTINE NCLASS_NU(p_ngrth,x,temp_i, &
                     ynud_s,ynut_s,ynuti_s)
!-------------------------------------------------------------------------------
!NCLASS_NU calculates the velocity dependent pitch angle diffusion and
!  anisotropy relaxation rates, nu_D, nu_T, and nu_T*I_Rm
!
!References:                                                     
!  S.P.Hirshman, D.J.Sigmar, Phys Fluids 19 (1976) 1532
!  S.P.Hirshman, D.J.Sigmar, Nucl Fusion 21 (1981) 1079
!  K.C.Shaing, M.Yokoyama, M.Wakatani, C.T.Hsu, Phys Plasmas 3 (1996) 965
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  p_ngrth,             & !<n.grad(Theta)> [/m]
  x,                   & !normalized velocity v/(2kT/m)**0.5 [-]
  temp_i(:)              !temperature [keV]

!Declaration of output variables
REAL(KIND=rspec) :: &
  ynud_s(:),           & !pitch angle diffusion rate for s [/s]
  ynut_s(:),           & !anisotropy relaxation rate for s [/s]
  ynuti_s(:,:)           !PS anisotropy relaxation rates for s [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,im,j,jm,m

REAL(KIND=rspec) :: &
  c,c1,c2,g,phi,wma

!-------------------------------------------------------------------------------
!Initializaton
!-------------------------------------------------------------------------------
ynud_s(:)=zero
ynut_s(:)=zero
ynuti_s(:,:)=zero

!-------------------------------------------------------------------------------
!Collision rates
!-------------------------------------------------------------------------------
DO i=1,ms_nc !Over species a

  im=jms_nc(i)

  !nu_D and nu_T
  DO j=1,ms_nc !Over species b

    jm=jms_nc(j)
    c1=vti_nc(jm)/vti_nc(im)
    c2=x/c1
    phi=NCLASS_ERF(c2)
    g=(phi-c2*(2/SQRT(z_pi))*EXP(-c2**2))/(2*c2**2)
    ynud_s(i)=ynud_s(i)+(3*SQRT(z_pi)/4)*(phi-g)/x**3/tauss_nc(i,j)
    ynut_s(i)=ynut_s(i)+((3*SQRT(z_pi)/4)*((phi-3*g)/x**3 &
                       +4*(temp_i(im)/temp_i(jm)+1/c1**2)*g/x))/tauss_nc(i,j)

  ENDDO !Over species b
   
  !nu_T*I_m
  DO m=1,mfm_nc !Over poloidal moments

    IF(ABS(p_ngrth) > zero) THEN

      wma=x*vti_nc(im)*m*p_ngrth
      c1=ynut_s(i)/wma
      c2=c1**2

      IF(c2 > 9.0_rspec) THEN

        !Use asymptotic limit for efficiency
        ynuti_s(m,i)=0.4_rspec

      ELSE

        !Use full calculation
        c=c1*ATAN(one/c1)
        ynuti_s(m,i)=c/2+c2*(3*(c-0.5_rspec)+c2*4.5*(c-one))

      ENDIF

    ELSE

      ynuti_s(m,i)=0.4_rspec

    ENDIF

  ENDDO !Over poloidal moments

ENDDO !Over species a

END SUBROUTINE NCLASS_NU

SUBROUTINE NCLASS_TAU(amu_i,temp_i,den_iz)
!-------------------------------------------------------------------------------
!NCLASS_TAU calculates various characteristic collision times
!
!References:                                      
!  S.P.Hirshman, D.J.Sigmar, Nucl Fusion 21 (1981) 1079
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  amu_i(:),            & !atomic mass number of i [-]
  temp_i(:),           & !temperature of i [keV]
  den_iz(:,:)            !density of i,z [/m**3]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,im,iz,iza,j,jm,jz,jza

REAL(KIND=rspec) :: &
  c1,c2,cln

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
amntii_nc(:,:)=zero
tauss_nc(:,:)=zero

!Use same Coulomb logarithm for all species
!Electrons
cln=37.8_rspec-LOG(SQRT(den_iz(imel_nc,1))/temp_i(imel_nc))
!Use constant for all species
cln=17.0_rspec

!-------------------------------------------------------------------------------
!Collision times
!-------------------------------------------------------------------------------
c1=4/SQRT(z_pi)/3*(4*z_pi)*cln*(z_coulomb/(4*z_pi*z_epsilon0))**2 &
   *(z_coulomb/z_pmass)**2

DO i=1,ms_nc

  im=jms_nc(i)
  iz=jzs_nc(i)
  iza=IABS(iz)
  c2=(vti_nc(im)**3)*amu_i(im)**2/c1

  DO j=1,ms_nc

    jm=jms_nc(j)
    jz=jzs_nc(j)
    jza=ABS(jz)
    tauss_nc(i,j)=c2/iz**2/(den_iz(jm,jza)*jz**2)
    amntii_nc(im,jm)=amntii_nc(im,jm)+amu_i(im)*z_pmass*den_iz(im,iza) &
                     /tauss_nc(i,j)

  ENDDO

ENDDO

!-------------------------------------------------------------------------------
!Characteristic fast ion slowing down time on electrons
!-------------------------------------------------------------------------------
IF(l_fast_nc) tausf_nc(1:mf_nc)=amuf_nc(1:mf_nc)/amu_i(imel_nc) &
                                /zf_nc(1:mf_nc)**2*tauss_nc(jsel_nc,jsel_nc)

END SUBROUTINE NCLASS_TAU

SUBROUTINE NCLASS_VISCOSITY(p_ft)
!--------------------------------------------------------------------------------
!NCLASS_VISCOSITY computes the Laguerre order 3/2 matrix elements of the
!  Coulomb collision operator for the viscosity to all orders in the Laguerre
!  polynomial energy expansion
!
!References:
!  W.A.Houlberg, S.P.Hirshman, 9/2006
!--------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  p_ft                   !trapped fraction [-]

!Declaration of local variables
INTEGER :: &
  i,im,j,jm

REAL(KIND=dpspec), DIMENSION(0:2*k_order_nc-1) :: &
  p,q

REAL(KIND=dpspec), DIMENSION(k_order_nc,k_order_nc) :: &
  testm

REAL(KIND=dpspec), DIMENSION(k_order_nc,k_order_nc,mi_nc) :: &
  ymu

REAL(KIND=dpspec) :: &
  beta,L00,onepv2,v,v2,v3

!-------------------------------------------------------------------------------
!Isotopic viscosities in banana regime using recursion
!-------------------------------------------------------------------------------
ymu(:,:,:)=zero

DO im=1,mi_nc !Over im isotopes

  DO jm=1,mi_nc !Over jm isotopes

    !Set lowest order terms for recursion relationships
    v=vti_nc(im)/vti_nc(jm)
    v2=v**2
    v3=v*v2
    onepv2=onedp+v2
    beta=1/SQRT(onepv2)
    p(0)=v/SQRT(onepv2)
    q(0)=-v/SQRT(onepv2)+LOG(v+SQRT(onepv2))

    !Use recursion for higher order recursive elements
    DO i=1,2*k_order_nc-1

      beta=(2*i-1)*beta/2/onepv2
      p(i)=i*p(i-1)+v*beta
      q(i)=(i-1)*q(i-1)+v3*beta

    ENDDO

    !Construct test matrix from recursive elements
    DO i=1,k_order_nc !Over rows

      DO j=1,k_order_nc !Over columns

        testm(i,j)=p(i+j-2)-q(i+j-2)/v2

      ENDDO !Over columns

    ENDDO !Over rows

    !Convert from polynomial moments to Laguerre moments
    DO i=0,k_order_nc-1

      DO j=0,k_order_nc-1

        L00=(-1)**(i+j)*clag_nc(i+1)*clag_nc(j+1)*amntii_nc(im,jm)
        ymu(i+1,j+1,im)=ymu(i+1,j+1,im)+L00*NCLASS_LAGUERRE(0,testm,i,j)

      ENDDO

    ENDDO

  ENDDO !Over jm isotopes

ENDDO !Over im isotopes

!-------------------------------------------------------------------------------
!Scaling in banana regime
!-------------------------------------------------------------------------------
!Trapped fraction
ymu(:,:,:)=ymu(:,:,:)*p_ft/(one-p_ft)

!Charge state wieghting and orbit squeezing
DO i=1,ms_nc
  
  ymus_nc(:,:,i)=ymu(:,:,jms_nc(i))*xis_nc(i)/sqzs_nc(i)

ENDDO

END SUBROUTINE NCLASS_VISCOSITY

SUBROUTINE NCLASS_BACKSUB(a,n,indx,b)
!-------------------------------------------------------------------------------
!NCLASS_BACKSUB solves the matrix equation a x = b, where a is in LU form as
!  generated by a prior call to NCLASS_DECOMP
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments: 
!  The complete procedure is, given the equation a*x = b,
!    CALL NCLASS_DECOMP( )
!    CALL NCLASS_BACKSUB( ) 
!    and b is now the solution vector
!  To solve with a different right hand side, just reload b as  desired
!    and use the same LU decomposition
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n,                   & !number of equations to be solved [-]
  indx(:)                !row permutations due to partial pivoting [-]

REAL(KIND=rspec), INTENT(IN) :: &
  a(:,:)                 !coefficient matrix in lu decomposed form [-]

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &
  b(:)                   !(input) right hand side of equation [-]
                         !(output) solution vector [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,ii,j,k

REAL(KIND=rspec) :: &
  sum

!-------------------------------------------------------------------------------
!Find the index of the first nonzero element of b
!-------------------------------------------------------------------------------
ii=0

DO i=1,n

  k=indx(i)
  sum=b(k)
  b(k)=b(i)

  IF(ii /= 0) THEN

    DO j=ii,i-1

      sum=sum-a(i,j)*b(j)

    ENDDO
 
  ELSEIF (sum /= zero) THEN

    ii=i

  ENDIF

  b(i)=sum

ENDDO
    
!-------------------------------------------------------------------------------
!Back substitution
!-------------------------------------------------------------------------------
DO i=n,1,-1

  sum=b(i)

  IF(i < n) THEN

    DO j=i+1,n

      sum=sum-a(i,j)*b(j)

    ENDDO

  ENDIF

  b(i)=sum/a(i,i)

ENDDO   

END SUBROUTINE NCLASS_BACKSUB

SUBROUTINE NCLASS_DECOMP(a,n,indx, &
                         d,c,iflag,message)
!-------------------------------------------------------------------------------
!NCLASS_DECOMP performs an LU decomposition of the matrix a and is called
!  prior to NCLASS_BACKSUB to solve linear equations or to invert a matrix
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n                      !number of equations to be solved [-]

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &
  a(:,:)                 !coefficient matrix, overwritten on return [-]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag,               & !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error
  indx(:)                !row permutations due to partial pivoting [-]

REAL(KIND=rspec), INTENT(OUT) :: &
  c,                   & !condition of matrix = Min|a(i,i)|/Max|a(i,i)|
  d                      !flag for number of row exchanges [-]
                         !=1.0 even number
                         !=-1.0 odd number

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,imax,j,k

REAL(KIND=rspec) :: &
  aamax,al,as,dum,sum, &
  vv(n)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
d=one
c=one

!-------------------------------------------------------------------------------
!Loop over rows to get the implicit scaling information
!-------------------------------------------------------------------------------
DO i=1,n

  aamax=zero

  DO j=1,n

    IF(ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))

  ENDDO

  IF(aamax == zero) THEN

    iflag=1
    message='NCLASS_DECOMP/ERROR:singular matrix(1)'
    GOTO 9999

  ENDIF

  vv(i)=one/aamax

ENDDO
   
!-------------------------------------------------------------------------------
!Use Crout's method for decomposition
!-------------------------------------------------------------------------------
!Loop over columns
DO j=1,n

  DO i=1,j-1

    sum=a(i,j)

    DO k=1,i-1

      sum=sum-a(i,k)*a(k,j)

    ENDDO

    a(i,j)=sum

  ENDDO

  !Search for largest pivot element using dum as a figure of merit
  aamax=zero

  DO i=j,n

    sum=a(i,j)

    DO k=1,j-1

      sum=sum-a(i,k)*a(k,j)

    ENDDO

    a(i,j)=sum
    dum=vv(i)*ABS(sum)

    IF(dum >= aamax) THEN

      imax=i
      aamax=dum

    ENDIF

  ENDDO

  IF(j /= imax) THEN

    !Interchange rows
    DO k=1,n

      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum

    ENDDO

    d=-d
    vv(imax)=vv(j)

  ENDIF

  indx(j)=imax

  IF(a(j,j) == zero) THEN

    iflag=1
    message='NCLASS_DECOMP/ERROR:singular matrix(2)'
    GOTO 9999

  ENDIF

  IF(j /= n) THEN

    !Divide by pivot element
    dum=one/a(j,j)

    DO i=j+1,n

      a(i,j)=a(i,j)*dum

    ENDDO

  ENDIF

ENDDO

!Check condition of decomposed matrix
al=zero
as=1.0e20_rspec

DO i=1,n

  IF(ABS(a(i,i)) > al) al=a(i,i)
  IF(ABS(a(i,i)) < as) as=a(i,i)

ENDDO

c=as/al

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE NCLASS_DECOMP

FUNCTION NCLASS_ERF(x)
!-------------------------------------------------------------------------------
!NCLASS_ERF evaluates the error function erf(x)
!
!References:
!  M.Abramowitz, L.A.Stegun, Handbook of Math. Functions, p. 299
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The error function can consume a lot of time deep in the multiple species
!    loop so a very efficient calculation is called for
!  Time consumption is much more critical than accuracy as suggested by T.Amano
!  A three term expansion from Abramowitz is not sufficiently accurate because
!    it generates viscosities with singularities in the vicinity of the BP-PS
!    transition for ions
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  x                      !argument of error function [-]

!Declaration of output variables
REAL(KIND=rspec) :: &
  NCLASS_ERF             !value of error function [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec) :: &
  t

REAL(KIND=rspec), PARAMETER :: &
  c = 0.3275911_rspec,   &
  a1= 0.254829592_rspec, &
  a2=-0.284496736_rspec, &
  a3= 1.421413741_rspec, &
  a4=-1.453152027_rspec, &
  a5= 1.061405429_rspec

!-------------------------------------------------------------------------------
!Apply fit
!-------------------------------------------------------------------------------
t=one/(one+c*x)
NCLASS_ERF=one-(a1+(a2+(a3+(a4+a5*t)*t)*t)*t)*t*EXP(-x**2)

END FUNCTION NCLASS_ERF

END MODULE NCLASS_MOD
