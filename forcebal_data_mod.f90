MODULE FORCEBAL_DATA_MOD
!-------------------------------------------------------------------------------
!FORCEBAL_DATA_MOD is an F90 module that contains FORCEBAL profile data
!
!References:
!  W.A.Houlberg, F90 free format 8/2006
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD

!Options
LOGICAL :: &
  l_banana,            & !option to include banana viscosity [-]
  l_pfirsch,           & !option to include Pfirsch-Schluter transport [-]
  l_classical,         & !option to include classical transport [-]
  l_potato,            & !option to include potato orbit effects [-]
  l_squeeze,           & !option to include orbit squeezing [-]
  l_reduce_out           !reduce output for mult charge impurities [-]
                         !=.TRUE. only highest charge state
                         !=.FALSE. all charge states

!Array dimensions
!An 'isotope' is any species with the same nuclear mass and charge
!  -examples are e, D, H, T, He, ...
!  _i in array names
!A 'species' is an isotope with a given atomic charge state
!  -examples are e-, D+, H+, T+, He+, He++, ...
!  _s in array names
!The charge state can be identified by the position arrays
!  _z in array names

INTEGER :: &
  k_order,             & !order of v moments to be solved [-]
                         !=2 u and q (default)
                         !=3 u, q, and u2
                         !=else error
  m_i,                 & !number of isotopes (> 1) [-]
  m_s,                 & !number of species [ms>1]
  m_z,                 & !highest charge state [-]
  nr_r                   !number of radial nodes analysis [-]

!Following are JET-specific for rerieving data from PPFs
CHARACTER(len=4) :: &
  cid_dda(25),         & !DDA name for profiles [-]
                         !1 - T_e profile
                         !2 - T_i profile
                         !3 - n_e profile
                         !4-13 - single charge state ion densities
                         !14-23 - multiple charge state ion densities
                         !24 - diagnostic impurity for toroidal rotation
                         !25 - toroial rotation velocity of idagnostic impurity
  cid_dda_efit,        & !DDA name for EFIT solution [-]
  cid_dtype(25)          !dtype name for profiles [-]
                         !see above

CHARACTER(len=7) :: &
  cid_uid(25)            !uid name for profiles [-]
                         !see above
CHARACTER(len=8) :: &
  cid_uid_efit           !uid name for EFIT solution [-]

INTEGER :: &
  id_seqp(25),         & !sequence number for profiles [-]
                         !see above
  id_seqp_efit           !sequence number for EFIT solution [-]

!Species identification
LOGICAL :: &
  l_diag                 !option for including diagnostic impurity [-]

CHARACTER(len=3), ALLOCATABLE :: &
  cid_i(:)               !isotope identification [-]

INTEGER :: &
  js_diag,             & !species index of diagnostic ion [-]
  js_el,               & !species index of electrons [-]
  n_scsion,            & !number of single charge state ions [-]
  n_mcsion               !number of multiple charge state ions [-]

!Arrays for mapping between species and isotopic variables
INTEGER, ALLOCATABLE :: &
  izmax_i(:),          & !nuclear charge of each isotope [-]
  jif_s(:),            & !isotope number for each species [-]
  jzf_s(:)               !charge number for each species [-]

REAL(KIND=rspec), ALLOCATABLE :: &
  amu_i(:)               !atomic mass number of isotopes [-]

!Radial grid and geometry
REAL(KIND=rspec), ALLOCATABLE :: &
  b2_r(:),             & !!<B**2> [T**2]
  bm2_r(:),            & !<1/B**2> [/T**2]
  btout_r(:),          & !toroidal field at rout_r(i) [T]
  bpout_r(:),          & !poloidal field at rout_r(i) [T]
  cur_r(:),            & !total curent [A]
  cur_bs_r(:),         & !calculated bootstrap curent [A]
  cur_ex_r(:),         & !experimental bootstrap current [A]
  cur_nb_r(:),         & !neutral beam current [A]
  elong_r(:),          & !plasma elongation [-]
  f_r(:),              & !poloidal current=2*pi*R*B_t/mu0 [A]
  fhat_r(:),           & !RB_t/(dpsi_r/drho) [-]
  fm_r(:,:),           & !geometric factor [-]
  ftrap_r(:),          & !trapped particle fraction [-]
  gph_r(:),            & !<1/R**2>V'/(2*pi)**2 [-]
  gr2bm2_r(:),         & !a0**2*<|grad(rhot_r)|**2/B**2> [/T**2]
  grho1_r(:),          & !a0*<|grad(rhot_r)|> [-]
  grho2_r(:),          & !a0**2*<|grad(rhot_r)|**2> [-]
  grth_r(:),           & !<n.grad(theta)> [/m]
  gth_r(:),            & !<gtt/sqrt(g)>-theta average
  phit_r(:),           & !toroidal magnetic flux [Wb]
  psi_r(:),            & !poloidal magnetic flux [Wb/rad]
  psip_r(:),           & !Psi' [Wb/m]
  q_r(:),              & !safety factor [-]
  rhop_r(:),           & !normalized poloidal flux grid proportional to psi [-]
  rhot_r(:),           & !normalized toroidal flux grid proportional to sqrt(tor flux) [-]
  rin_r(:),            & !major radius grid on inside of torus in axis plane [m]
  rm2_r(:),            & !<1/R**2> [/m**2]
  rout_r(:),           & !major radius grid on outside of torus in axis plane [m]
  vol_r(:),            & !volume enclosed [m**3]
  vp_r(:)                !dV/d rhot_r/a0 [m**2]

!Plasma profiles
REAL(KIND=rspec), ALLOCATABLE :: &
  den_riz(:,:,:),      & !density profiles [/m**3]
  grp_riz(:,:,:),      & !pressure gradient profiles [keV/m**4]
  grt_ri(:,:),         & !temperature gradient profiles [keV/m]
  temp_ri(:,:),        & !temperature profiles [keV]
  te_r(:),             & !e temperature profile [keV]
  ti_r(:),             & !i temperature profile [keV]
  zeff_r(:),           & !calculated Zeff profile [-]
  zeff_ex_r(:)           !expt Zeff profile [-]

!Particle transport
REAL(KIND=rspec), ALLOCATABLE :: &
  d_eff_rs(:,:),       & !effective particle diffusivity [m**2/s]
  d_eff_g_rs(:,:),     & !effective particle diffusivity/grho2 [m**2/s]
  d_n_rs(:,:),         & !(diagonal) particle diffusivity [m**2/s]
  d_n_g_rs(:,:),       & !particle diffusivity/grho2 [m**2/s]
  g_n_rss(:,:,:),      & !ns2(r)~ns1(r)**g [-]
  g_te_rs(:,:),        & !ns(r)~te(r)**g [-]
  g_ti_rs(:,:),        & !ns(r)~ti(r)**g [-]
  gam_rs(:,:),         & !radial particle flux [/m**2/s]
  v_eb_rs(:,:),        & !Ware conv velocity [m/s]
  v_eb_g_rs(:,:),      & !Ware conv velocity/grho1 [m/s]
  v_nt_rs(:,:),        & !off-diagonal pinch velocity [m/s]
  v_nt_g_rs(:,:)         !off-diagonal pinch velocity/grho1 [m/s]

!Heat transport
REAL(KIND=rspec), ALLOCATABLE :: &
  chi_eff_r(:,:),      & !effective (e,i) thermal diffusivity [m**/2]
  chi_eff_g_r(:,:),    & !effective (e,i) thermal diffusivity/grho2 [m**/2]
  q_con_r(:,:)           !

!Current density, electrical resistivity, parallel electric field
REAL(KIND=rspec), ALLOCATABLE :: &
  eta_par_r(:),        & !parallel electrical resistivity [Ohm*m]
  xj_r(:),             & !total parallel current density, <J.B>/Bt0 [A/m**2]
  xj_bs_r(:),          & !total bootstrap current density, <J.B>/Bt0 [A/m**2]
  xj_ex_r(:),          & !expt total parallel current density, <J.B>/Bt0 [A/m**2]
  xj_nb_r(:),          & !total NBI current density, <J.B>/Bt0 [A/m**2]
  edotb_r(:),          & !<E.B> [V*T/m]
  e_par_r(:),          & !parallel electric field, <E.B>/Bt0 [V/m]
  e_par_ex_r(:)          !expt parallel electric field, <E.B>/Bt0 [V/m]

!Rotation velocities
REAL(KIND=rspec), ALLOCATABLE :: &
  u_par_rs(:,:),       & !parallel flow function <u.B> [m*T/s]
  u_p_rs(:,:),         & !pol flow function <u.gr th>/<B.gr th> [m/T/s]
  v_p_o_rs(:,:),       & !pol flow velocity on outside [m/s]
  v_t_o_rs(:,:),       & !tor flow velocity on outside [m/s]
  vt_im_ex_r(:),       & !expt toroidal flow velocity on outside [m/s]
  vp_im_ex_r(:),       & !expt poloidal flow velocity on outside [m/s]
  kk_im_ex_r(:),       & !Flux surf ave poloidal flow velocity on outside [m/s/T]
  v_par_o_rs(:,:),     & !parallel flow velocity on outside [m/s]
  v_per_o_rs(:,:),     & !perpendicular flow velocity on outside [m/s]
  xm_p_o_rs(:,:),      & !pol Mach no on outside [-]
  xm_t_o_rs(:,:),      & !tor Mach no on outside [-]
  xm_t_im_ex_o_r(:),   & !expt tor Mach no on outside [-]
  omega_rs(:,:),       & !tor rotation frequency [rad/s]
  omega_im_ex_r(:)       !expt tor rotation frequency [rad/s]

!Radial electric field
REAL(KIND=rspec), ALLOCATABLE :: &
  e_rad_r(:,:),        & !radial electric field from diag ion=Phi' [V/m]
                         !1-total
                         !2-tor rotation contribution
                         !3-pol rotation contribution
                         !4-pressure gradient contribution
  e_rad_o_r(:,:),      & !Er on outside midplane from diag ion [V/m]
                         !same contributions as above
  e_rad_rs(:,:,:)        !radial electric field=Phi' [V/m]
                         !same contributions as above

!Rotational shear damping
REAL(KIND=rspec), ALLOCATABLE :: &
  omexb_o_r(:),        & !Hahm-Burrell ExB shear damping on outside [/s]
  sqz_rs(:,:)            !orbit squeezing factor [-]

!Fast ions from NBI
LOGICAL :: &
  l_beam                 !option for including beam ions [logical]

REAL(KIND=rspec) :: &
  amu_b(3),            & !atomic mass number of beam ions [-]
  pmw_b(3),            & !beam injection power [MW]
  e_b(3),              & !beam energy [keV]
  v0_b(3),             & !beam ion birth energy
  z_b(3)                 !beam ion charge [-]

REAL(KIND=rspec), ALLOCATABLE :: &
  fshld_r(:),          & !beam ion shielding factor [-]
  den_rb(:,:),         & !beam ion density [/m**3]
  dendot_rb(:,:),      & !beam ion source rate [/m**3/s]
  taus_rb(:,:),        & !beam slowing downing down rate on electrons [s]
  u_rb(:,:)              !beam ion parallel velocity [m/s]

END MODULE FORCEBAL_DATA_MOD