!COMFBL-radial variables for FORCEBAL
!Electrons, main ions and single charge state impurity
!  d_eff_r-eff part diffusivity=-gam_r/(dn/(dr_t)) (m**2/s)
!  d_eff_g_r-eff part diffusivity=d_eff_r/<grad(r_t)^2> (m**2/s)
!  d_n_r-diagonal diffusivity (m**2/s)
!  d_n_g_r-d_n_r/<grad(r_t)^2> (m**2/s)
!  gam_r-net nc radial flux=<Gam.grad(r_t)> (/m**2/s)
!  sqz_r-orbit squeezing factor
!  omega_r-toroidal rotation frequency (rad/s)
!  u_par_r-<u.B> (T*m/s)
!  u_p_r-<u.grad(theta)>/<B.grad(theta)> (m/s/T)
!  v_n_r-nondiffusive (non-WP) velocity (m/s)
!  v_n_g_r-v_n_r/(<grad(r_t)^2>)**0.5 (m/s)
!  v_p_o_r-poloidal velocity on outside of torus (m/s)
!  v_t_o_r-toroidal velocity on outside of torus (m/s)
!  xm_p_o_r-poloidal Mach number on outside of torus
!  xm_p_ex_im_o_r-expt poloidal Mach number of im on outside of torus
!  xm_t_o_r-toroidal Mach number on outside of torus
!  xm_t_ex_im_o_r-expt toroidal Mach number of im on outside of torus
!  v_wp_r-Ware pinch velocity (m/s)
!  v_wp_g_r-v_wp_r/(<grad(r_t)^2>)**0.5 (m/s)
!  chi_eff_r-eff cond=-q_con_r/(ne*dkTe/(drt)) (m**2/s)
!  chi_eff_g_r-chi_eff_r/<grad(r_t)^2> (m**2/s)
!  eta_par_r-parallel electrical resistivity (Ohm*m)
!  q_con_r-net nc cond flux=<q_con.grad(r_t)> (w/m**2)
!  e_rad_r-E_r and components=-dPhi/(dr_t) (V/m)
!  e_rad_ex_r-E_r using expt p', vp and vt (V/m)
!  e_rad_o_r-E_r and components on outside of torus (V/m)
!  e_rad_ex_o_r-E_r on outside using expt p', vp and vt (V/m)
!  edotb_r-<E.B> (V*T/m)
!  omexb_o_r-ExB shear damping rate on outside of torus (rad/s)
!  xj_bs_r-bootstrap current density=<Jbc.B>/Bt0 (A/m**2)
!Multiple charge state impurity
!  gam_im2_r-net nc radial flux=<Gam.grad(r_t)> (/m**2/s)
!  d_n_g_im2_r-d_n_r/<grad(r_t)^2> (m**2/s)
!  v_n_g_im2_r-v_n_r/(<grad(r_t)^2>)**0.5 (m/s)
!  v_wp_im2_r-Ware pinch velocity (m/s)
      REAL           d_eff_r,                 d_eff_g_r,
     #               d_n_r,                   d_n_g_r,
     #               gam_r,
     #               sqz_r,                   omega_r,
     #               u_par_r,                 u_p_r,
     #               v_n_r,                   v_n_g_r,
     #               v_p_o_r,                 v_t_o_r,
     #               xm_p_o_r,                xm_p_im_ex_o_r,
     #               xm_t_o_r,                xm_t_im_ex_o_r,
     #               v_wp_r,                  v_wp_g_r,
     #               chi_eff_r,               chi_eff_g_r,
     #               eta_par_r,               q_con_r,
     #               e_rad_r,                 e_rad_ex_r,
     #               e_rad_o_r,               e_rad_ex_o_r,
     #               e_par_r,                 omexb_o_r,
     #               xj_bs_r,                 gam_im2_r,
     #               d_n_g_im2_r,             v_n_g_im2_r,
     #               v_wp_g_im2_r
      COMMON/COMFBL/ d_eff_r(mxnr_r,mxns+2),  d_eff_g_r(mxnr_r,mxns+2),
     #               d_n_r(mxnr_r,mxns+2),    d_n_g_r(mxnr_r,mxns+2),
     #               gam_r(mxnr_r,mxns+2),    sqz_r(mxnr_r,mxns+2),
     #               omega_r(mxnr_r,mxns+2),
     #               u_par_r(mxnr_r,mxns+2),  u_p_r(mxnr_r,mxns+2),
     #               v_n_r(mxnr_r,mxns+2),    v_n_g_r(mxnr_r,mxns+2),
     #               v_p_o_r(mxnr_r,mxns+2),  v_t_o_r(mxnr_r,mxns+2),
     #               xm_p_o_r(mxnr_r,mxns+2), xm_p_im_ex_o_r(mxnr_r),
     #               xm_t_o_r(mxnr_r,mxns+2), xm_t_im_ex_o_r(mxnr_r),
     #               v_wp_r(mxnr_r,mxns+2),   v_wp_g_r(mxnr_r,mxns+2),
     #               chi_eff_r(mxnr_r,2),     chi_eff_g_r(mxnr_r,2),
     #               eta_par_r(mxnr_r),       q_con_r(mxnr_r,2),
     #               e_rad_r(mxnr_r,4),       e_rad_ex_r(mxnr_r),
     #               e_rad_o_r(mxnr_r,4),     e_rad_ex_o_r(mxnr_r),
     #               e_par_r(mxnr_r),         omexb_o_r(mxnr_r),
     #               xj_bs_r(mxnr_r),         gam_im2_r(mxnr_r,mxmz),
     #               d_n_g_im2_r(mxnr_r,mxmz),v_n_g_im2_r(mxnr_r,mxmz),
     #               v_wp_g_im2_r(mxnr_r,mxmz)
