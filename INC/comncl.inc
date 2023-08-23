!COMNCL-internal variables for NCLASS
!  dencut-density cutoff below which species is ignored
!  etap-parallel electrical resistivity
!  xb2-<B^2>
!  xbm2-<1/B^2>
!  xedotb-<E.B>
!  xfhat-mu_0*F/(dPsi/dr)
!  xft-trapped fraction
!  xgrbm2-<grad(rho)^2/B^2>
!  xjdotb-<J_bs.B>
!  xngrth-<n.grad(Theta)>
!  xbpot-kappa(0)*bt0/(2*q(0)**2)
!  fm-poloidal moments of geometric factor for PS viscosity
!  amuai-isotopic atomic masses
!  tai-isotopic temperatures
!  vtai-isotopic thermal velocities
!  xgrt-isotopic temperature gradients
!  xompot-vt(0)/(R(0)q(0))
!  amntau-effective isotopic momentum relaxation rate
!  tauab-90 degree velocity-dependent collisional scattering times
!  denai-species densities
!  xi-species charge weighted density factors
!  sqz-species squeezing factors
!  xgrp-species pressure gradients
!  xkban-species velocity-dependent banana viscosity
!  xkps-species velocity-dependent Pfirsch-Schluter viscosity
!  xkpot-species velocity-dependent potato viscosity
!  xkpp-species velocity-dependent potato-plateau viscosity
!  xnud-species velocity-dependent pitch angle diffusion rate
!  xnut-species velocity-dependent rate of anisotropy relaxation
!  capm-test particle friction matrix
!  capn-field particle friction matrix
!  calm-test particle charge-weighted friction matrix
!  caln-field particle charge-weighted friction matrix
!  uai-normalized parallel species flows
!  xmu-species viscosity matrix
!  xmuban-species banana viscosity
!  xmubps-species banana Pfirsch-Schluter viscosity
!  xmupot-species potato viscosity
!  xmuppp=species potato potato-plateau viscosity
!  xfxt-external source terms
!  xnuti-species velocity-dependent rate of anisotropy relaxation for p-PS
!  gfl-species radial particle flux by regime
!  qfl-species radial heat conduction flux by regime
!  dnai-species diagonal particle diffusion coefficient
!  vnai-species particle convection velocity from off-diagonal p' and T'
!  vwpai-species <E.B> driven Ware pinch velocity
!  xjdotbp-<J_bs.B> driven by unit species p'/p
!  xjdotbt-<J_bs.B> driven by unit species T'/T
!  dpai-species diffusion coefficient array for p'/p terms
!  dtai-species diffusion coefficient array for T'/T terms
!  xpai-species heat conduction coefficient array for p'/p terms
!  xtai-species heat conduction coefficient array for T'/T terms
!  kboot-bootstrap current model option switch
!  kncord-order of equations solved (number of v moments)
!  kpotato-option to include potato orbits (=0 is off)
!  mi-number of isotopes
!  mz-higher charge state of all species
      REAL           dencut,                  etap,
     #               xb2,                     xbm2,
     #               xedotb,                  xfhat,
     #               xft,                     xgrbm2,
     #               xjdotb,                  xngrth,
     #               xbpot,                    
     #               fm,
     #               amuai,                   tai,
     #               vtai,                    xgrt,
     #               xompot,
     #               amntau,
     #               tauab,
     #               denai,                   xi,
     #               sqz,                     xgrp,
     #               xkban,                   xkps,
     #               xkpot,                   xkpp,
     #               xnud,                    xnut,
     #               capm,                    capn,
     #               calm,                    caln,
     #               uai,                     xmu,
     #               xmuban,                  xmubps,
     #               xmupot,                  xmuppp,
     #               xfxt,                    xnuti,
     #               gfl,                     qfl,
     #               dnai,                    vnai,
     #               vwpai,                   xjdotbp,
     #               xjdotbt,
     #               dpai,                    dtai,
     #               xpai,                    xtai
      INTEGER        kboot,                   kncord,
     #               kpotato,
     #               mi,                      mz
      COMMON/COMNCL/ dencut,                  etap,
     #               xb2,                     xbm2,
     #               xedotb,                  xfhat,
     #               xft,                     xgrbm2,
     #               xjdotb,                  xngrth,
     #               xbpot,
     #               fm(3),
     #               amuai(mxmi),             tai(mxmi),
     #               vtai(mxmi),              xgrt(mxmi),
     #               xompot(mxmi),
     #               amntau(mxmi,mxmi),
     #               tauab(mxmi,mxmi,mxmz,mxmz),
     #               denai(mxmi,mxmz),        xi(mxmi,mxmz),
     #               sqz(mxmi,mxmz),          xgrp(mxmi,mxmz),
     #               xkban(mxmi,mxmz),        xkps(mxmi,mxmz),
     #               xkpot(mxmi,mxmz),        xkpp(mxmi,mxmz),
     #               xnud(mxmi,mxmz),         xnut(mxmi,mxmz),
     #               capm(3,3,mxmi,mxmi),     capn(3,3,mxmi,mxmi),
     #               calm(3,3,mxmi),          caln(3,3,mxmi,mxmi),
     #               uai(3,3,mxmi,mxmz),      xmu(3,3,mxmi,mxmz),
     #               xmuban(3,3,mxmi,mxmz),   xmubps(3,3,mxmi,mxmz),
     #               xmupot(3,3,mxmi,mxmz),   xmuppp(3,3,mxmi,mxmz),
     #               xfxt(3,mxmi,mxmz),       xnuti(3,mxmi,mxmz),
     #               gfl(5,mxmi,mxmz),        qfl(5,mxmi,mxmz),
     #               dnai(mxsp),              vnai(mxsp),
     #               vwpai(mxsp),             xjdotbp(mxsp),
     #               xjdotbt(mxsp),
     #               dpai(mxsp,mxsp),         dtai(mxsp,mxsp),
     #               xpai(mxsp,mxsp),         xtai(mxsp,mxsp),
     #               kboot,                   kncord,
     #               kpotato,
     #               mi,                      mz
