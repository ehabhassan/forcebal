PROGRAM FORCEBAL
!-------------------------------------------------------------------------------
!FORCEBAL is a program that calculates neoclassical transport properties,
!  poloidal rotation velocities and the radial electric field using the parallel
!  and radial force balance equations.
!General sequence:
!  1) read MHD equilibrium data
!  2) read profile data
!  3) call NCLASS
!  4) calculate additional dependent variables
!  5) output to various files
!
!References:
!  W.A.Houlberg, K.C.Shaing, S.P.Hirshman, M.C.Zarnstorff, Phys Plasmas 4 (1997)
!    3230
!  W.A.Houlberg, F90 free format 8/2004
!  W.A.Houlberg, repaired toroidal rotation constraint 8/2006
!
!Comments:
!  Units are contained in [] in the description of variables
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE FORCEBAL_DATA_MOD
USE NETCDF
USE READ_PRO_MOD       !USE READ_PRO_JET_MOD !for use at JET
USE WRITE_MOD
USE WRITE_NETCDF_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
!Declaration of namelist input variables
!See FORCEBAL_DATA_MOD for declarations of:
! l_banana,            & !option to include banana viscosity [logical]
! l_pfirsch,           & !option to include Pfirsch-Schluter transport [logical]
! l_classical,         & !option to include classical transport [logical]
! l_potato,            & !option to include potato orbit effects [logical]
! l_squeeze,           & !option to include orbit squeezing [logical]
! l_reduce_out           !reduce output for mult charge impurities [logical]
!                        !=.TRUE. only highest charge state
!                        !=.FALSE. all charge states
! k_order                !order of v moments to be solved [integer]
!                        !=2 u and q (default)
!                        !=3 u, q, and u2
!                        !=else error
! nr_r                   !number of radial nodes analysis [integer]
!Following are JET-specific for rerieving data from PPFs:
! cid_uid_efit           !UID name for EFIT solution [character]
! cid_dda_efit           !DDA name for EFIT solution [character]
! id_seqp_efit           !sequence number for EFIT solution [integer]
! cid_uid(25)            !UID names for profiles [character]
! cid_dda(25)            !DDA names for profiles [character]
! cid_dtype(25)          !dtype names for profiles [character]
! id_seqp(25)            !sequence numbers for profiles [integer]
!Following are for fast ions from NBI:
! amu_b(3),            & !atomic mass number of beam ions [-]
! pmw_b(3)               !beam injection power [MW] 
! z_b(3)                 !beam ion charge [-]

INTEGER, PARAMETER :: &
  mx_ni=10               !Maximum number of scs plus mcs ions in NAMELIST input [-]
                         !scs denotes single charge state ions
                         !mcs denotes multiple charge state ions

CHARACTER(len=15) :: &
  cid_device,          & !name of device being analyzed [character]
  cid_mhdeq              !name of equilibrium code data used [character]
                         !='efit' for EFIT equilibria
                         !=else failure

CHARACTER(len=3) :: &
  cid_run,             & !run identification [character]
  cid_diag,            & !name of diagnositc ion [character]
  cid_mcsion(mx_ni),   & !names of multiple charge state ions [character]
  cid_scsion(mx_ni)      !names of single charge state ions [character]

INTEGER :: &
  k_edotb,             & !option to use experimental or calculated <E.B> [-]
                         !=1    use experimental
                         !=else use calculated
  id_shot                !shot number [-]

REAL(KIND=rspec) :: &
  time                   !analysis time [s]

!-------------------------------------------------------------------------------
!Declaration of other variables
CHARACTER(len=25) :: &
  cn_cdf,cn_msg,cn_sum,cn_tmp

CHARACTER(len=5) :: &
  ctimems

CHARACTER(len=9) :: &
  cidshot

CHARACTER(len=120) :: &
  label,message

INTEGER :: &
  n_cdf,n_msg,n_sum,n_tmp

INTEGER :: &
  i,ierr,iflag,itimems,j,k

REAL(KIND=rspec) :: &
  a0,a1,bt0,c_den,ecrit,e0_b,r0,tau,vcrit,vol,denz2

!Physical constants, mathematical constants, conversion factors
REAL(KIND=rspec) :: &
  z_coulomb=1.6022e-19_rspec,    & !Coulomb charge [C]
  z_emass=9.1095e-31_rspec,      & !electron mass [kg]
  z_epsilon0=8.8542e-12_rspec,   & !permittivity of free space [F/m]
  z_j7kv=1.6022e-16_rspec,       & !energy conversion factor [J/keV]
  z_mu0=1.2566e-06_rspec,        & !permeability of free space [H/m]
  z_pi=3.141592654_rspec,        & !pi [-]
  z_pmass=1.6726e-27_rspec         !proton mass [kg]

REAL(KIND=rspec), PARAMETER :: &
  one=1.0_rspec,       & !REAL 1
  zero=0.0_rspec         !REAL 0

!Output arrays
INTEGER, PARAMETER :: &
  mx_npro=400, &
  mx_nsca=100

CHARACTER(len=15) :: &
  namepro(mx_npro),unitpro(mx_npro),namesca(mx_nsca),unitsca(mx_nsca)

CHARACTER(len=95) :: &
  descpro(mx_npro),descsca(mx_nsca)

INTEGER :: &
  npro,nsca

REAL(KIND=rspec) :: &
  valsca(mx_nsca)

REAL(KIND=rspec), ALLOCATABLE :: &
  valpro(:,:)

!-------------------------------------------------------------------------------
! EHAB MODIFIED LINES BELOW
!NAMELIST/inforce/l_banana,l_pfirsch,l_classical,l_potato,l_squeeze, &
NAMELIST/inforcebal/l_banana,l_pfirsch,l_classical,l_potato,l_squeeze, &
                 l_reduce_out, &
                 cid_device,cid_mhdeq,cid_run,cid_diag,cid_scsion,cid_mcsion, &
                 k_edotb,k_order, &
                 id_shot,nr_r, &
                 time, &
                 cid_uid_efit,cid_dda_efit,id_seqp_efit, &
                 cid_uid,cid_dda,cid_dtype,id_seqp, &
                 amu_b,z_b,e0_b,pmw_b
! EHAB MODIFIED LINES ABOVE

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Warning/error flag and message
iflag=0
message=' '

!Namelist input - logical
l_banana=.TRUE.
l_pfirsch=.TRUE.
l_classical=.TRUE.
l_potato=.FALSE.
l_squeeze=.FALSE.
l_reduce_out=.TRUE.

!Namelist input - character
cid_device='diiid'
cid_mhdeq='efit'
cid_run='A01'
cid_diag=''
cid_scsion(:)=''
cid_mcsion(:)=''
cid_uid_efit=''
cid_dda_efit=''
cid_uid(:)=''
cid_dda(:)=''
cid_dtype(:)=''

!Namelist input - integer
k_edotb=0
k_order=2
id_shot=10000
nr_r=41
id_seqp_efit=0
id_seqp(:)=0

!Namelist input - real
time=zero
amu_b(:)=zero
e0_b=zero
pmw_b(:)=zero

!Other
c_den=1.0e10_rspec       !density cutoff below which species is ignored [/m**3]

!-------------------------------------------------------------------------------
!Set the input namelist unit, open, read and close file
!-------------------------------------------------------------------------------
n_tmp=20
! EHAB MODIFIED LINES BELOW
!cn_tmp='nml_forcebal.dat'
cn_tmp='inforcebal'
! EHAB MODIFIED LINES ABOVE
OPEN(UNIT=n_tmp, &                                                    
     FILE=cn_tmp, &
     STATUS='old', &
     ACCESS='sequential')

! EHAB MODIFIED LINES BELOW
!READ(n_tmp,inforce)
READ(n_tmp,inforcebal)
! EHAB MODIFIED LINES ABOVE

CLOSE(UNIT=n_tmp)

!-------------------------------------------------------------------------------
!Consistency checks
!-------------------------------------------------------------------------------
!NBI fast ions
IF(amu_b(1) < 0.99_rspec .OR. &
   e0_b < 0.99_rspec .OR. &
   z_b(1) < 0.99_rspec .OR. &
   SUM(pmw_b(:)) < 0.1_rspec) THEN

  l_beam=.FALSE.
  e0_b=zero
  amu_b(:)=zero
  pmw_b(:)=zero
  e_b(:)=zero
  z_b(:)=zero

ELSE

  l_beam=.TRUE.

  !Full, half, third energies
  amu_b(2:3)=amu_b(1)
  e_b(1:3)=e0_b/(/ (i,i=1,3) /)
  z_b(2:3)=z_b(1)

ENDIF

!-------------------------------------------------------------------------------
!Open output files and write namelist to summary file
!-------------------------------------------------------------------------------
!Set file names and device dependent options
itimems=INT(1.0e3*(time+1.0e-5_rspec))

WRITE(ctimems,'(i5)') itimems

DO i=1,5

  IF(ctimems(i:i) == ' ') ctimems(i:i)='0'

ENDDO

WRITE(cidshot,'(i9)') id_shot

!Summary
n_sum=10
cn_sum='sum_'//TRIM(ADJUSTL(cidshot))//'_'//TRIM(ADJUSTL(ctimems))//'_' &
       //TRIM(ADJUSTL(cid_run))//'.dat'
OPEN(UNIT=n_sum, &
     FILE=cn_sum, &  
     STATUS='unknown', &    
     FORM='formatted')

!netCDF
cn_cdf='cdf_'//TRIM(ADJUSTL(cidshot))//'_'//TRIM(ADJUSTL(ctimems))//'_' &
       //TRIM(ADJUSTL(cid_run))//'.nc'
ierr=NF90_CREATE(PATH=TRIM(cn_cdf), &
                 CMODE=NF90_CLOBBER, &
                 NCID=n_cdf)

IF(ierr /= nf90_noerr) THEN

  iflag=1
  message='NF90_CREATE/ERROR: '//TRIM(nf90_strerror(ierr))
  CALL WRITE_LINE(n_msg,message,1,1)
  GOTO 9999

ELSE

  ierr=NF90_CLOSE(n_cdf)

  IF(ierr /= nf90_noerr) THEN

    iflag=1
    message='NF90_CLOSE/ERROR: '//TRIM(nf90_strerror(ierr))
    CALL WRITE_LINE(n_msg,message,1,1)
    GOTO 9999

  ENDIF

ENDIF

!Messages
n_msg=11
cn_msg='msg_'//TRIM(ADJUSTL(cidshot))//'_'//TRIM(ADJUSTL(ctimems))//'_' &
        //TRIM(ADJUSTL(cid_run))//'.dat'
OPEN(UNIT=n_msg, &
     FILE=cn_msg, &
     STATUS='unknown', &
     FORM='formatted')

!Listing of namelist variables to summary file
label='*** Namelist Data ***'
CALL WRITE_LINE(n_sum,label,1,1)
!EHAB MODIFIED LINES BELOW
!WRITE(n_sum,inforce)
WRITE(n_sum,inforcebal)
!EHAB MODIFIED LINES ABOVE

!-------------------------------------------------------------------------------
!Allocate and initialize arrays
!-------------------------------------------------------------------------------
!FORCEBAL data module
CALL FORCEBAL_INIT(mx_ni,cid_scsion,cid_diag,cid_mcsion, &
                   iflag,message)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

!-------------------------------------------------------------------------------
!Get geometry information
!-------------------------------------------------------------------------------
rhot_r(1:nr_r)=REAL((/ (i-1,i=1,nr_r) /),rspec)/(nr_r-1)

IF(cid_mhdeq == 'efit') THEN

  !EFIT type, Psi(R,Z)
  CALL FORCEBAL_EFIT(cid_device,id_shot,time,cid_uid_efit,cid_dda_efit, &
                id_seqp_efit,nr_r,rhot_r, &
                a0,bt0,r0, &
                b2_r,bm2_r,bpout_r,btout_r,elong_r,f_r,fhat_r,fm_r, &
                ftrap_r,gph_r,grho1_r,grho2_r,gr2bm2_r,grth_r,gth_r, &
                phit_r,psi_r,q_r,rhop_r,rin_r,rout_r,vol_r,vp_r, &
                iflag,message)

ELSE

  !Not an implemented option
  iflag=1
  message='FORCEBAL/ERROR:cid_mhdeq option not implemented'

ENDIF

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

IF(cid_device == 'jet') THEN

  !Renormalize metrics to JETTO minor radius
  a1=SQRT(ABS(phit_r(nr_r)/z_pi/bt0))
  grho1_r(1:nr_r)=grho1_r(1:nr_r)*(a1/a0)
  grho2_r(1:nr_r)=grho2_r(1:nr_r)*(a1/a0)**2
  gr2bm2_r(1:nr_r)=gr2bm2_r(1:nr_r)*(a1/a0)**2
  fhat_r(1:nr_r)=fhat_r(1:nr_r)*(a1/a0)
  vp_r(1:nr_r)=vp_r(1:nr_r)*(a0/a1)

ELSE

  a1=a0

ENDIF

rm2_r(2:nr_r)=(2*z_pi)**2*gph_r(2:nr_r)/vp_r(2:nr_r)
rm2_r(1)=rm2_r(2)

!-------------------------------------------------------------------------------
!Get plasma profile information
!-------------------------------------------------------------------------------
CALL READ_PRO(id_shot,time,c_den,iflag,message)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

!Rescale <J.B>--><J.B>/bt0, <E.B>--><E.B>/bt0
e_par_ex_r(:)=e_par_ex_r/bt0
xj_ex_r(:)=xj_ex_r(:)/bt0

!Set Psi'
psip_r(1)=zero
psip_r(2:nr_r)=z_mu0*f_r(2:nr_r)/fhat_r(2:nr_r)

!Set integrated currents
cur_r(1)=zero
cur_ex_r(1)=zero

DO i=2,nr_r

  cur_r(i)=gth_r(i)*gph_r(i)/q_r(i)*f_r(i)
  cur_ex_r(i)=f_r(i)*(cur_ex_r(i-1)/f_r(i-1) &
              +a1*bt0*(rhot_r(i)-rhot_r(i-1))/z_mu0 &
              *(vp_r(i)*xj_ex_r(i)/f_r(i)**2 &
               +vp_r(i-1)*xj_ex_r(i-1)/f_r(i-1)**2)/2)

ENDDO

!Reset total current using <J.B> from second to last node
cur_r(nr_r)=f_r(nr_r)*(cur_r(nr_r-1)/f_r(nr_r-1) &
            +(f_r(nr_r-1)/f_r(nr_r))**2*vp_r(nr_r)/vp_r(nr_r-1) &
            *(rhot_r(nr_r)-rhot_r(nr_r-1))/(rhot_r(nr_r-1)-rhot_r(nr_r-2)) &
            *(cur_r(nr_r-1)/f_r(nr_r-1)-cur_r(nr_r-2)/f_r(nr_r-2)))

!Total current density from equilibrium
DO i=2,nr_r-1

  xj_r(i)=f_r(i)**2*z_mu0/vp_r(i)/bt0/a1/(rhot_r(i+1)-rhot_r(i-1)) &
          *(cur_r(i+1)/f_r(i+1)-cur_r(i-1)/f_r(i-1))

ENDDO

xj_r(1)=xj_r(2)
xj_r(nr_r)=xj_r(nr_r-1)

!Set beam ion arrays
IF(l_beam) THEN

  vol=SUM(vol_r(1:nr_r-1))

  DO k=1,nr_r !Over radial grid

    !Particle source density assumes H(r)=1
    dendot_rb(k,1:3)=1.0e6*pmw_b(1:3)/e_b(1:3)/z_j7kv/vol

    !Characteristic slowing down on electrons
    tau=(3*SQRT(z_pi)/4)/(4*z_pi*17.0)/(z_coulomb/(4*z_pi*z_epsilon0))**2 &
        /(z_coulomb/z_pmass)**2*amu_i(1)*amu_b(1)/den_riz(k,1,1) &
        *(2*te_r(k)*z_j7kv/z_emass)**1.5

    !v_crit
    vcrit=zero

    DO i=2,m_i !Over ion isotopes

      DO j=1,m_z !Over charge states

        IF(den_riz(k,i,j) > c_den) THEN

          vcrit=vcrit+REAL(j,rspec)**2*den_riz(k,i,j)/amu_i(i)

        ENDIF

      ENDDO !Over charge states

    ENDDO !Over ion isotopes

    vcrit=(vcrit*3*SQRT(z_pi)/4*amu_i(1)/den_riz(k,1,1))**(1.0/3.0) &
          *(2*te_r(k)*z_j7kv/z_emass)**0.5
    ecrit=amu_b(1)*z_pmass*vcrit**2/2/z_j7kv
    den_rb(k,1:3)=dendot_rb(k,1:3)*tau/3 &
                  *LOG((e_b(1:3)**1.5+ecrit**1.5) &
                      /(temp_ri(k,2)**1.5+ecrit**1.5))
    v0_b(1:3)=SQRT(2*e_b(1:3)*z_j7kv/amu_b(1)/z_pmass)

    !Need better estimates of this
    u_rb(k,1:3)=0.2*bt0*v0_b(1:3)

    !Readjust main ion density for charge neutrality
    den_riz(k,2,1)=den_riz(k,2,1)-SUM(den_rb(k,1:3))

    !Recalculate Zeff
    denz2=SUM(den_rb(k,1:3))

    DO i=2,m_i !Over ion isotopes

      DO j=1,m_z !Over charge states

        denz2=denz2+j**2*den_riz(k,i,j)

      ENDDO !Over charge states

    ENDDO !Over ion isotopes

    zeff_r(k)=denz2/den_riz(k,1,1)

  ENDDO  !Over radial grid

ENDIF

!-------------------------------------------------------------------------------
!Call NCLASS interface to get rotation and transport properties from NCLASS
!-------------------------------------------------------------------------------
CALL FORCEBAL_NCLASS(k_edotb,c_den,a1,bt0,q_r(1), &
                     iflag,message)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

!Integrated bootstrap current
cur_bs_r(1)=zero

DO i=2,nr_r

  cur_bs_r(i)=f_r(i)*(cur_bs_r(i-1)/f_r(i-1) &
              +bt0/z_mu0*a1*(rhot_r(i)-rhot_r(i-1)) &
              *(vp_r(i)*xj_bs_r(i)/f_r(i)**2 &
               +vp_r(i-1)*xj_bs_r(i-1)/f_r(i-1)**2)/2)
ENDDO

!Integrated neutral beam current
cur_nb_r(1)=zero

DO i=2,nr_r

  cur_nb_r(i)=f_r(i)*(cur_nb_r(i-1)/f_r(i-1) &
              +bt0/z_mu0*a1*(rhot_r(i)-rhot_r(i-1)) &
              *(vp_r(i)*xj_nb_r(i)/f_r(i)**2 &
               +vp_r(i-1)*xj_nb_r(i-1)/f_r(i-1)**2)/2)
ENDDO

!-------------------------------------------------------------------------------
!Scalar data
!-------------------------------------------------------------------------------
namesca(:)=''
unitsca(:)='-'
descsca(:)=''
valsca(:)=zero
nsca=0

nsca=nsca+1
namesca(nsca)='l_banana'
descsca(nsca)='Include banana viscosity (1=yes)'
IF(l_banana) valsca(nsca)=one

nsca=nsca+1
namesca(nsca)='l_pfirsch'
descsca(nsca)='Include Pfirsch-Schlueter viscosity (1=yes)'
IF(l_pfirsch) valsca(nsca)=one

nsca=nsca+1
namesca(nsca)='l_classical'
descsca(nsca)='Include classical transport (1=yes)'
IF(l_classical) valsca(nsca)=one

nsca=nsca+1
namesca(nsca)='l_potato'
descsca(nsca)='Include potato viscosity (1=yes)'
IF(l_potato) valsca(nsca)=one

nsca=nsca+1
namesca(nsca)='l_squeeze'
descsca(nsca)='Include orbit squeezing (1=yes)'
IF(l_squeeze) valsca(nsca)=one

nsca=nsca+1
namesca(nsca)=cid_device
descsca(nsca)='Device'
valsca(nsca)=one

nsca=nsca+1
namesca(nsca)='id_shot'
descsca(nsca)='Shot number'
valsca(nsca)=id_shot

nsca=nsca+1
namesca(nsca)='time'
unitsca(nsca)='s'
descsca(nsca)='Time'
valsca(nsca)=time

IF(cid_device == 'jet') THEN

  nsca=nsca+1
  namesca(nsca)=cid_uid_efit
  descsca(nsca)='UID for EFIT solution'

  nsca=nsca+1
  namesca(nsca)=cid_dda_efit
  descsca(nsca)='DDA for EFIT solution'

  nsca=nsca+1
  namesca(nsca)='SEQP_EFIT'
  descsca(nsca)='SEQP for EFIT solution'
  valsca(nsca)=id_seqp_efit

  nsca=nsca+1
  namesca(nsca)=cid_uid(1)
  descsca(nsca)='UID for Te profile'

  nsca=nsca+1
  namesca(nsca)=cid_dda(1)
  descsca(nsca)='DDA for Te profile'

  nsca=nsca+1
  namesca(nsca)=cid_dtype(1)
  descsca(nsca)='dtype for Te profile'

  nsca=nsca+1
  namesca(nsca)='SEQP_Te'
  descsca(nsca)='SEQP for Te profile'
  valsca(nsca)=id_seqp(1)

  nsca=nsca+1
  namesca(nsca)=cid_uid(2)
  descsca(nsca)='UID for Ti profile'

  nsca=nsca+1
  namesca(nsca)=cid_dda(2)
  descsca(nsca)='DDA for Ti profile'

  nsca=nsca+1
  namesca(nsca)=cid_dtype(2)
  descsca(nsca)='dtype for Ti profile'

  nsca=nsca+1
  namesca(nsca)='SEQP_Ti'
  descsca(nsca)='SEQP for Ti profile'
  valsca(nsca)=id_seqp(2)

  nsca=nsca+1
  namesca(nsca)=cid_uid(3)
  descsca(nsca)='UID for ne profile'

  nsca=nsca+1
  namesca(nsca)=cid_dda(3)
  descsca(nsca)='DDA for ne profile'

  nsca=nsca+1
  namesca(nsca)=cid_dtype(3)
  descsca(nsca)='dtype for ne profile'

  nsca=nsca+1
  namesca(nsca)='SEQP_ne'
  descsca(nsca)='SEQP for ne profile'
  valsca(nsca)=id_seqp(3)

  DO i=4,3+n_scsion

    nsca=nsca+1
    namesca(nsca)=cid_uid(i)
    descsca(nsca)='UID for single charge state of '//cid_scsion(i-3)

    nsca=nsca+1
    namesca(nsca)=cid_dda(i)
    descsca(nsca)='DDA for single charge state of '//cid_scsion(i-3)

    nsca=nsca+1
    namesca(nsca)=cid_dtype(i)
    descsca(nsca)='dtype for single charge state of '//cid_scsion(i-3)

    nsca=nsca+1
    namesca(nsca)='SEQP_'//cid_scsion(i-3)
    descsca(nsca)='SEQP for single charge state of '//cid_scsion(i-3)
    valsca(nsca)=id_seqp(i)
    
  ENDDO

  DO i=14,13+n_mcsion

    nsca=nsca+1
    namesca(nsca)=cid_uid(i)
    descsca(nsca)='UID for multiple charge states of '//cid_mcsion(i-13)

    nsca=nsca+1
    namesca(nsca)=cid_dda(i)
    descsca(nsca)='DDA for multiple charge states of '//cid_mcsion(i-13)

    nsca=nsca+1
    namesca(nsca)=cid_dtype(i)
    descsca(nsca)='dtype for multiple charge states of '//cid_mcsion(i-13)

    nsca=nsca+1
    namesca(nsca)='SEQP_'//cid_mcsion(i-13)
    descsca(nsca)='SEQP for multiple charge states of '//cid_mcsion(i-13)
    valsca(nsca)=id_seqp(i)
    
  ENDDO

  IF(TRIM(ADJUSTL(cid_uid(24))) /= '') THEN

    nsca=nsca+1
    namesca(nsca)=cid_uid(24)
    descsca(nsca)='UID for diagnsotic impurity density '//cid_diag

    nsca=nsca+1
    namesca(nsca)=cid_dda(24)
    descsca(nsca)='DDA for diagnsotic impurity density '//cid_diag

    nsca=nsca+1
    namesca(nsca)=cid_dtype(24)
    descsca(nsca)='dtype for diagnsotic impurity density '//cid_diag

    nsca=nsca+1
    namesca(nsca)='SEQP_'//cid_diag
    descsca(nsca)='SEQP for diagnsotic impurity density '//cid_diag
    valsca(nsca)=id_seqp(24)

  ENDIF

  IF(TRIM(ADJUSTL(cid_uid(25))) /= '') THEN

    nsca=nsca+1
    namesca(nsca)=cid_uid(25)
    descsca(nsca)='UID for diagnsotic impurity tor rotation '//cid_diag

    nsca=nsca+1
    namesca(nsca)=cid_dda(25)
    descsca(nsca)='DDA for diagnsotic impurity tor rotation '//cid_diag

    nsca=nsca+1
    namesca(nsca)=cid_dtype(25)
    descsca(nsca)='dtype for diagnsotic impurity tor rotation '//cid_diag

    nsca=nsca+1
    namesca(nsca)='SEQP_'//cid_diag//'_vel'
    descsca(nsca)='SEQP for diagnsotic impurity tor rotation '//cid_diag
    valsca(nsca)=id_seqp(25)

  ENDIF

ENDIF

nsca=nsca+1
namesca(nsca)='nr_r'
descsca(nsca)='Number of radial points'
valsca(nsca)=nr_r

nsca=nsca+1
namesca(nsca)='a0'
unitsca(nsca)='m'
descsca(nsca)='Plasma minor radius in midplane'
valsca(nsca)=a0

nsca=nsca+1
namesca(nsca)='a1'
unitsca(nsca)='m'
descsca(nsca)='Reference plasma minor radius'
valsca(nsca)=a1

nsca=nsca+1
namesca(nsca)='R0'
unitsca(nsca)='m'
descsca(nsca)='Reference plasma major radius'
valsca(nsca)=r0

nsca=nsca+1
namesca(nsca)='Bt0'
unitsca(nsca)='T'
descsca(nsca)='Toroidal field at R0'
valsca(nsca)=bt0

nsca=nsca+1
namesca(nsca)='elongation'
descsca(nsca)='Plasma elongation at edge'
valsca(nsca)=elong_r(nr_r)

!Print scalar data to summary file
label='*** Scalar Data ***'
! EHAB ADDED LINES BELOW
!nsca = SIZE(namesca)
! EHAB ADDED LINES ABOVE
CALL WRITE_LINE(n_sum,label,1,1)
CALL WRITE_OUT0(n_sum,nsca,namesca,unitsca,descsca,valsca)

!Print scalar data to netCDF file
CALL WRITE_NETCDF_OUT0(cn_cdf,nsca,namesca,unitsca,descsca,valsca, &
                       iflag,message)
!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

!-------------------------------------------------------------------------------
!Profile data
!-------------------------------------------------------------------------------
!Allocate and initialize output arrays
ALLOCATE(valpro(nr_r,mx_npro))

  namepro(:)=''
  unitpro(:)='-'
  descpro(:)=''
  valpro(:,:)=zero
  npro=0

CALL FORCEBAL_OUT(mx_npro,a1,r0,npro,valpro,namepro,unitpro,descpro)

!Print profile data to summary data file
label='*** Profiles ***'
CALL WRITE_LINE(n_sum,label,1,1)
CALL WRITE_OUT1(n_sum,'sum',nr_r,npro,namepro,unitpro,descpro,valpro,1)

!Print profile data to 1D data file
n_tmp=20
cn_tmp='1d_'//TRIM(ADJUSTL(cidshot))//'_'//TRIM(ADJUSTL(ctimems))//'_' &
       //TRIM(ADJUSTL(cid_run))//'.dat'
OPEN(UNIT=n_tmp, &
     FILE=cn_tmp, &
     STATUS='unknown', &
     FORM='formatted')

CALL WRITE_OUT1(n_tmp,'1d',nr_r,npro,namepro,unitpro,descpro,valpro,-1)

CLOSE(UNIT=n_tmp)

!Print profile data to netCDF data file
CALL WRITE_NETCDF_OUT1(cn_cdf,1,nr_r,npro,namepro,unitpro,descpro,valpro, &
                  iflag,message)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END PROGRAM FORCEBAL

SUBROUTINE FORCEBAL_EFIT(device,idshot,time,cid_uid_efit,cid_dda_efit, &
              id_seqp_efit,nr_r,rhot_r, &
              a0,bt0,r0, &
              b2_r,bm2_r,bpout_r,btout_r,elong_r,f_r,fhat_r,fm_r, &
              ftrap_r,gph_r,grho1_r,grho2_r,gr2bm2_r,grth_r,gth_r, &
              phit_r,psi_r,q_r,rhop_r,rin_r,rout_r,vol_r,vp_r, &
              iflag,message)
!-------------------------------------------------------------------------------
!FORCEBAL_EFIT generates plasma geometry information for FORCEBAL by calling
!  getting EQDSK information about the MHD equilibrium and then calling
!  FLUXAV_LOAD to load the data in the FLUXAV module, and then calling FLUXAV to 
!  generate needed flux surface quantities
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!                Updated call to READ_EFIT_JET 8/2006
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE FLUXAV_MOD
IMPLICIT NONE

!Declaration of input variables
CHARACTER(len=*), INTENT(IN) :: &
  device,              & !experimental device
  cid_uid_efit,        & !UID for EFIT solution in PPFs at JET
  cid_dda_efit           !DDA FOR EFIT solution in PPFs at JET

INTEGER, INTENT(IN) :: &
  idshot,              & !shot number/id [-]
  id_seqp_efit,        & !sequence number for EFIT solution in PPFs at JET [-]
  nr_r                   !number of radial points [-]

REAL(KIND=rspec), INTENT(IN) :: &
  time,                & !analysis time [s]
  rhot_r(nr_r)           !normalized tor flux grid proportional to (Phi)**0.5 [-]

!Declaration of output variables
INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

REAL(KIND=rspec), INTENT(OUT) :: &
  a0,                  & !minor radius, half diameter of boundary flux surface [m]
  bt0,                 & !toroidal field at r0 [T]
  r0                     !major radius, center of boundary flux suface [m]

REAL(KIND=rspec), INTENT(OUT) :: &
  b2_r(nr_r),          & !<B**2> [T**2]
  bm2_r(nr_r),         & !<1/B**2> [/T**2]
  bpout_r(nr_r),       & !poloidal field at rout_r(i) [T]
  btout_r(nr_r),       & !toroidal field at rout_r(i) [T]
  elong_r(nr_r),       & !plasma elongation [-]
  f_r(nr_r),           & !poloidal current=2*pi*R*B_t/mu0 [A]
  fhat_r(nr_r),        & !RB_t/(dpsi_r/drho) [-]
  fm_r(3,nr_r),        & !geometric factor [-]
  ftrap_r(nr_r),       & !trapped particle fraction [-]
  gph_r(nr_r),         & !<1/R**2>V'/(2*pi)**2 [-]
  grho1_r(nr_r),       & !a0*<|grad(rhot_r)|> [-]
  grho2_r(nr_r),       & !a0**2*<|grad(rhot_r)|**2> [-]
  gr2bm2_r(nr_r),      & !a0**2*<|grad(rhot_r)|**2/B**2> [/T**2]
  grth_r(nr_r),        & !<n.grad(theta)> [/m]
  gth_r(nr_r),         & !<gtt/sqrt(g)>-theta average [-]
  phit_r(nr_r),        & !toroidal magnetic flux [Wb]
  psi_r(nr_r),         & !poloidal magnetic flux [Wb/rad]
  q_r(nr_r),           & !safety factor [-]
  rhop_r(nr_r),        & !normalized poloidal flux grid proportional to psi [-]
  rin_r(nr_r),         & !major radius grid on inside of torus in axis plane [m]
  rout_r(nr_r),        & !major radius grid on outside of torus in axis plane [m]
  vol_r(nr_r),         & !volume enclosed [m**3]
  vp_r(nr_r)             !dV/d rhot_r/a0 [m**2]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER, PARAMETER :: &
  mxnx_xy=130, &
  mxny_xy=130,  &
  mxn_lim=200, &
  mxn_bdry=1500

REAL(KIND=rspec), PARAMETER :: &
  z_mu0=1.2566e-06_rspec,        & !permeability of free space [H/m]
  z_pi=3.141592654_rspec           !pi [-]

REAL(KIND=rspec), PARAMETER :: &
  one=1.0_rspec,       & !REAL 1
  zero=0.0_rspec         !REAL 0

CHARACTER :: &
  cidshot*9,ctimems*5,cnin*25

INTEGER :: &
  i,itimems,k_grid,nin

INTEGER :: &
  n_lim,nx_xy,ny_xy

REAL(KIND=rspec) :: &
  cur,psimag,psilim,rmag,zmag

REAL(KIND=rspec) :: &
  f_x(mxnx_xy),ffp_x(mxnx_xy),psi_x(mxnx_xy),q_x(mxnx_xy),rhop_x(mxnx_xy), &
  r2_r(nr_r),rm2_r(nr_r),triang_r(nr_r), &
  x_xy(mxnx_xy),y_xy(mxny_xy),psi_xy(mxnx_xy,mxny_xy), &
  x_lim(mxn_lim),y_lim(mxn_lim)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Local arrays
x_xy(:)=zero
y_xy(:)=zero
psi_xy(:,:)=zero
f_x(:)=zero
ffp_x(:)=zero
rhop_x(:)=zero
psi_x(:)=zero
q_x(:)=zero
rhop_x(:)=zero
x_lim(:)=zero
y_lim(:)=zero

r2_r(:)=zero
rm2_r(:)=zero
triang_r(:)=zero

!-------------------------------------------------------------------------------
!Get EFIT data
!-------------------------------------------------------------------------------
IF(device == 'jet') THEN

  !Read PPF files for JET EFIT MHD equilibrium information
!  CALL READ_EFIT_JET(idshot,time,cid_uid_efit,cid_dda_efit,id_seqp_efit, &
!                     nx_xy,ny_xy,r0,bt0,cur,psimag,psilim,rmag,zmag, &
!		                   x_xy,y_xy,psi_xy,psi_x,rhop_x,f_x,ffp_x,q_x, &
!                     n_lim,x_lim,y_lim,iflag,message)

ELSE

  !Read EQDSK file to get EFIT MHD equilibrium information
  itimems=1.0e3*(time+.00001_rspec)
  WRITE(ctimems,'(i5)') itimems

  DO i=1,5

    IF(ctimems(i:i) == ' ') ctimems(i:i)='0'

  ENDDO

  WRITE(cidshot,'(i9)') idshot

  DO i=1,9

    IF(cidshot(i:i) == ' ') cidshot(i:i)='0'

  ENDDO

  nin=20
  cnin='g'//cidshot(4:9)//'.'//ctimems
  CALL READ_EFIT_EQDSK(nin,cnin,mxnx_xy,mxny_xy,mxn_lim, &
                       bt0,cur,psimag,psilim,r0,rmag,zmag,nx_xy,ny_xy,x_xy, &
                       y_xy,psi_xy,f_x,ffp_x,psi_x,q_x,rhop_x,n_lim,x_lim, &
                       y_lim,iflag,message)
ENDIF

!Check messages
IF(iflag /= 0) THEN

  message='FORCEBAL_EFIT(1)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Call FLUXAV to generate metrics from EFIT MHD equilibrium
!-------------------------------------------------------------------------------
CALL FLUXAV_LOAD(cur,r0,rmag,zmag,psimag,psilim, &
                 nx_xy,ny_xy,x_xy,y_xy,psi_xy,f_x,ffp_x,rhop_x,q_x, &
                 n_lim,x_lim,y_lim, &
                 iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='FORCEBAL_EFIT(2)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

k_grid=0

CALL FLUXAV(k_grid,nr_r,rhot_r, &
            a0,b2_r,bm2_r,bpout_r,btout_r,elong_r,triang_r,f_r,fhat_r,fm_r, &
            ftrap_r,gph_r,gr2bm2_r,grho1_r,grho2_r,grth_r,gth_r,phit_r, &
            psi_r,q_r,r2_r,rin_r,rm2_r,rout_r,vol_r,vp_r, &
            iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='FORCEBAL_EFIT(3)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

!Change grid normalization to a0
fhat_r(:)=fhat_r(:)*a0
gph_r(:)=gph_r(:)/a0
gr2bm2_r(:)=gr2bm2_r(:)*a0**2
grho1_r(:)=grho1_r(:)*a0
grho2_r(:)=grho2_r(:)*a0**2
vp_r(:)=vp_r(:)/a0

!Define rhop_r
rhop_r(1:nr_r)=(psi_r(1:nr_r)-psi_r(1))/(psi_r(nr_r)-psi_r(1))

!At this point q_r is always positive, correct sign for coordinate consistency
q_r(1:nr_r)=q_r(1:nr_r)*SIGN(one,bpout_r(nr_r)*btout_r(nr_r))

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FORCEBAL_EFIT

SUBROUTINE FORCEBAL_INIT(mx_ni,cid_scsion,cid_diag,cid_mcsion, &
                         iflag,message)
!-------------------------------------------------------------------------------
!FORCEBAL_INIT initializes the arrays in the FORCEBAL data module
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE FORCEBAL_DATA_MOD
IMPLICIT NONE

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  mx_ni                 !dimension of ion species arrays in input [-]

CHARACTER(len=*), INTENT(IN) :: &
  cid_diag,            & !name of diagnositc ion [character]
  cid_scsion(mx_ni),   & !names of single charge state ions [character]
  cid_mcsion(mx_ni)      !names of multiple charge state ions [character]

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
   i

REAL(KIND=rspec), PARAMETER :: &
  zero=0.0_rspec         !REAL 0

!-------------------------------------------------------------------------------
!Determine number of isotopes
!-------------------------------------------------------------------------------
!Electrons
m_i=1

!Check single charge state (scs) ions until first null entry is encountered
i=1
n_scsion=0

DO WHILE(cid_scsion(i) /= '' .AND. &
         i <= mx_ni)  !Over scs ions

  n_scsion=n_scsion+1
  m_i=m_i+1
  i=i+1

ENDDO !Over single charge state ions

!Check for diagnostic impurity
l_diag=.FALSE.

IF(cid_diag /= '') THEN

  l_diag=.TRUE.
  m_i=m_i+1

ENDIF

!Check multiple charge state (mcs) ions until first null entry is encountered
i=1
n_mcsion=0

DO WHILE(cid_mcsion(i) /= '' .AND. &
         i <= mx_ni)  !Over scs ions

  n_mcsion=n_mcsion+1
  m_i=m_i+1
  i=i+1

ENDDO !Over single charge state ions

!-------------------------------------------------------------------------------
!Determine charge states and masses of each isotope, then total species
!-------------------------------------------------------------------------------
!Allocate arrays for isotope identification, maximum charge and atomic mass
ALLOCATE(cid_i(m_i), &
         izmax_i(m_i), &
         amu_i(m_i))

  cid_i(:)=''
  izmax_i(:)=0
  amu_i(:)=zero

!Store isotope names in single stacked array: e, scs, diag, mcs
cid_i(1)='e'
js_el=1
IF(n_scsion > 0) cid_i(2:1+n_scsion)=cid_scsion(1:n_scsion)

IF(l_diag) THEN

  cid_i(2+n_scsion)=cid_diag
  js_diag=2+n_scsion
  IF(n_mcsion > 0) cid_i(3+n_scsion:2+n_scsion+n_mcsion)=cid_mcsion(1:n_mcsion)

ELSE

  js_diag=0
  IF(n_mcsion > 0) cid_i(2+n_scsion:1+n_scsion+n_mcsion)=cid_mcsion(1:n_mcsion)

ENDIF

!Get maximum charge and mass for each isotope, and number of species
m_s=0

DO i=1,m_i

  CALL FORCEBAL_SPECIES(cid_i(i), &
                        izmax_i(i),amu_i(i),iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='FORCEBAL_INIT/'//message
    IF(iflag > 0) GOTO 9999

  ELSE

    IF(l_diag) THEN

      IF(i <= js_diag) THEN

        m_s=m_s+1

      ELSE

        m_s=m_s+izmax_i(i)

      ENDIF

    ELSE

      IF(i <= 1+n_scsion) THEN

        m_s=m_s+1

      ELSE

        m_s=m_s+izmax_i(i)

      ENDIF

    ENDIF

  ENDIF

ENDDO

m_z=MAXVAL(izmax_i)

!-------------------------------------------------------------------------------
!Allocate and initialize arrays
!-------------------------------------------------------------------------------
ALLOCATE(rhop_r(nr_r), &
         rhot_r(nr_r), &
         rin_r(nr_r), &
         rout_r(nr_r), &
         elong_r(nr_r), &
         vol_r(nr_r), &
         vp_r(nr_r), &
         cur_r(nr_r), &
         cur_bs_r(nr_r), &
         cur_ex_r(nr_r), &
         cur_nb_r(nr_r), &
         f_r(nr_r), &
         gph_r(nr_r), &
         gth_r(nr_r), &
         btout_r(nr_r), &
         bpout_r(nr_r), &
         phit_r(nr_r), &
         psi_r(nr_r), &
         psip_r(nr_r), &
         q_r(nr_r), &
         b2_r(nr_r), &
         bm2_r(nr_r), &
         fhat_r(nr_r), &
         ftrap_r(nr_r), &
         fm_r(3,nr_r), &
         grth_r(nr_r), &
         grho1_r(nr_r), &
         grho2_r(nr_r), &
         gr2bm2_r(nr_r), &
         rm2_r(nr_r))

  rhop_r(:)=zero
  rhot_r(:)=zero
  rin_r(:)=zero
  rout_r(:)=zero
  elong_r(:)=zero
  vol_r(:)=zero
  vp_r(:)=zero
  cur_r(:)=zero
  cur_bs_r(:)=zero
  cur_ex_r(:)=zero
  cur_nb_r(:)=zero
  f_r(:)=zero
  gph_r(:)=zero
  gth_r(:)=zero
  btout_r(:)=zero
  bpout_r(:)=zero
  phit_r(:)=zero
  psi_r(:)=zero
  psip_r(:)=zero
  q_r(:)=zero
  b2_r(:)=zero
  bm2_r(:)=zero
  fhat_r(:)=zero
  ftrap_r(:)=zero
  fm_r(:,:)=zero
  grth_r(:)=zero
  grho1_r(:)=zero
  grho2_r(:)=zero
  gr2bm2_r(:)=zero
  rm2_r(:)=zero

ALLOCATE(te_r(nr_r), &
         ti_r(nr_r), &
         den_riz(nr_r,m_i,m_z), &
         grp_riz(nr_r,m_i,m_z), &
         temp_ri(nr_r,m_i), &
         grt_ri(nr_r,m_i), &
         zeff_r(nr_r), &
         zeff_ex_r(nr_r))

  te_r(:)=zero
  ti_r(:)=zero
  den_riz(:,:,:)=zero
  grp_riz(:,:,:)=zero
  temp_ri(:,:)=zero
  grt_ri(:,:)=zero
  zeff_r(:)=zero
  zeff_ex_r(:)=zero

ALLOCATE(jif_s(m_s), &
         jzf_s(m_s), &
         d_eff_rs(nr_r,m_s), &
         d_eff_g_rs(nr_r,m_s), &
         d_n_rs(nr_r,m_s), &
         d_n_g_rs(nr_r,m_s), &
         g_n_rss(nr_r,m_s,m_s), &
         g_te_rs(nr_r,m_s), &
         g_ti_rs(nr_r,m_s), &
         gam_rs(nr_r,m_s), &
         v_eb_rs(nr_r,m_s), &
         v_eb_g_rs(nr_r,m_s), &
         v_nt_rs(nr_r,m_s), &
         v_nt_g_rs(nr_r,m_s))

  jif_s(:)=0
  jzf_s(:)=0
  d_eff_rs(:,:)=zero
  d_eff_g_rs(:,:)=zero
  d_n_rs(:,:)=zero
  d_n_g_rs(:,:)=zero
  g_n_rss(:,:,:)=zero
  g_te_rs(:,:)=zero
  g_ti_rs(:,:)=zero
  gam_rs(:,:)=zero
  v_eb_rs(:,:)=zero
  v_eb_g_rs(:,:)=zero
  v_nt_rs(:,:)=zero
  v_nt_g_rs(:,:)=zero

ALLOCATE(chi_eff_r(nr_r,2), &
         chi_eff_g_r(nr_r,2), &
         q_con_r(nr_r,2))

  chi_eff_r(:,:)=zero
  chi_eff_g_r(:,:)=zero
  q_con_r(:,:)=zero

ALLOCATE(eta_par_r(nr_r), &
         xj_r(nr_r), &
         xj_bs_r(nr_r), &
         xj_ex_r(nr_r), &
         xj_nb_r(nr_r), &
         edotb_r(nr_r), &
         e_par_r(nr_r), &
         e_par_ex_r(nr_r))

  eta_par_r(:)=zero
  xj_r(:)=zero
  xj_bs_r(:)=zero
  xj_ex_r(:)=zero
  xj_nb_r(:)=zero
  edotb_r(:)=zero
  e_par_r(:)=zero
  e_par_ex_r(:)=zero

ALLOCATE(u_par_rs(nr_r,m_s), &
         u_p_rs(nr_r,m_s), &
         v_p_o_rs(nr_r,m_s), &
         v_t_o_rs(nr_r,m_s), &
         vt_im_ex_r(nr_r), &
         ! EHAB ADDED TWO LINES BELOW
         vp_im_ex_r(nr_r), &
         kk_im_ex_r(nr_r), &
         ! EHAB ADDED TWO LINES ABOVE
         v_par_o_rs(nr_r,m_s), &
         v_per_o_rs(nr_r,m_s), &
         xm_p_o_rs(nr_r,m_s), &
         xm_t_o_rs(nr_r,m_s), &
         xm_t_im_ex_o_r(nr_r), &
         omega_rs(nr_r,m_s), &
         omega_im_ex_r(nr_r))

  u_par_rs(:,:)=zero
  u_p_rs(:,:)=zero
  v_p_o_rs(:,:)=zero
  v_t_o_rs(:,:)=zero
  vt_im_ex_r(:)=zero
  ! EHAB ADDED TWO LINES BELOW
  vp_im_ex_r(:)=zero
  kk_im_ex_r(:)=zero
  ! EHAB ADDED TWO LINES ABOVE
  v_par_o_rs(:,:)=zero
  v_per_o_rs(:,:)=zero
  xm_p_o_rs(:,:)=zero
  xm_t_o_rs(1:,:)=zero
  xm_t_im_ex_o_r(:)=zero
  omega_rs(:,:)=zero
  omega_im_ex_r(:)=zero

ALLOCATE(e_rad_r(nr_r,4), &
         e_rad_o_r(nr_r,4), &
         e_rad_rs(nr_r,m_s,4))

  e_rad_r(:,:)=zero
  e_rad_o_r(:,:)=zero
  e_rad_rs(:,:,:)=zero

ALLOCATE(omexb_o_r(nr_r), &
         sqz_rs(nr_r,m_s))

  omexb_o_r(:)=zero
  sqz_rs(:,:)=zero

IF(l_beam) THEN

  ALLOCATE(den_rb(nr_r,3), &
           dendot_rb(nr_r,3), &
           taus_rb(nr_r,3), &
           u_rb(nr_r,3), &
           fshld_r(nr_r))

    den_rb(:,:)=zero
    dendot_rb(:,:)=zero
    taus_rb(:,:)=zero
    u_rb(:,:)=zero
!EHAB ADDED LINES BELOW
ELSE

  ALLOCATE(den_rb(nr_r,3), &
           dendot_rb(nr_r,3), &
           taus_rb(nr_r,3), &
           u_rb(nr_r,3), &
           fshld_r(nr_r))

    den_rb(:,:)=zero
    dendot_rb(:,:)=zero
    taus_rb(:,:)=zero
    u_rb(:,:)=zero
!EHAB ADDED LINES ABOVE

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FORCEBAL_INIT

SUBROUTINE FORCEBAL_NCLASS(k_edotb,c_den,a1,bt0,q0, &
                           iflag,message)
!-------------------------------------------------------------------------------
!FORCEBAL_NCLASS generates profiles of NCLASS quantities using the MHD
!  equilibrium information along with plasma density, temperature and toroidal
!  rotation profiles
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE FORCEBAL_DATA_MOD
USE NCLASS_MOD
IMPLICIT NONE

!Declaration of input variables
INTEGER :: &
  k_edotb                !option to use experimental or calculated <E.B> [-]
                         !=1 use experimental
                         !=else use calculated
                    
REAL(KIND=rspec), INTENT(IN) :: &
  a1,                  & !reference minor radius [m]
  bt0,                 & !toroidal field at r0 [T]
  c_den,               & !density cutoff below which species is ignored [/m**3]
  q0                     !safety factor at magnetic axis [-]

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
!Array size declarations      
INTEGER :: &
  i,im,ima,iz,iza,j,jiter,jitermx,k,kk

INTEGER, ALLOCATABLE :: &
  jm_s(:),jz_s(:)

REAL(KIND=rspec)  :: &
  b,c_potb,c_potl,p_b2,p_bm2,p_bsjb,p_eb,p_etap,p_fhat,p_fm(1:3),p_fshld,p_ft, &
  p_grbm2,p_grphi,p_gr2phi,p_nbjb,p_ngrth

REAL(KIND=rspec), ALLOCATABLE :: &
  den_iz(:,:),dp_ss(:,:),dt_ss(:,:),gfl_s(:,:),grp_iz(:,:),grt_i(:),temp_i(:), &
  qfl_s(:,:),upar_s(:,:,:),utheta_s(:,:,:)

REAL(KIND=rspec) :: &
  dr,drout,dent,erref,grad,ppr,toler,vtherm,rdum(5)

REAL(KIND=rspec), ALLOCATABLE :: &
  erold1_r(:),erold2_r(:),ertemp_r(:),omegas_rs(:,:)

!Physical constants, mathematical constants, conversion factors
REAL(KIND=rspec) :: &
  z_coulomb=1.6022e-19_rspec,    & !Coulomb charge [C]
  z_j7kv=1.6022e-16_rspec,       & !energy conversion factor [J/keV]
  z_pi=3.141592654_rspec,        & !pi [-]
  z_pmass=1.6726e-27_rspec         !proton mass [kg]

REAL(KIND=rspec), PARAMETER :: &
  one=1.0_rspec,       & !REAL 1
  zero=0.0_rspec         !REAL 0

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Options and constants
jitermx=50
toler=0.01_rspec
c_potb=elong_r(1)*bt0/2/q0**2
c_potl=rout_r(1)*q0
p_gr2phi=zero

!Allocate NCLASS arrays
ALLOCATE(jm_s(m_s), &
         jz_s(m_s), &
         grt_i(m_i), &
         temp_i(m_i), &
         den_iz(m_i,m_z), &
         grp_iz(m_i,m_z), &
         gfl_s(5,m_s), &
         qfl_s(5,m_s), &
         upar_s(k_order,3,m_s), &
         utheta_s(k_order,3,m_s), &
         dp_ss(m_s,m_s), &
         dt_ss(m_s,m_s))

  jm_s(:)=0
  jz_s(:)=0
  grt_i(:)=zero
  temp_i(:)=zero
  den_iz(:,:)=zero
  grp_iz(:,:)=zero
  gfl_s(:,:)=zero
  qfl_s(:,:)=zero
  upar_s(:,:,:)=zero
  utheta_s(:,:,:)=zero
  dp_ss(:,:)=zero
  dt_ss(:,:)=zero

!Allocate local arrays
ALLOCATE(erold1_r(nr_r), &
         erold2_r(nr_r), &
         ertemp_r(nr_r), &
         omegas_rs(nr_r,m_s))

  erold1_r(:)=zero
  erold2_r(:)=zero
  ertemp_r(:)=zero
  omegas_rs(:,:)=zero

!-------------------------------------------------------------------------------
!Plasma profiles and gradients
!-------------------------------------------------------------------------------
!Electron temperatures
temp_ri(1:nr_r,1)=te_r(1:nr_r)

!Ion temperatures
DO k=2,m_i !Over ion species

  temp_ri(1:nr_r,k)=ti_r(1:nr_r)

ENDDO !Over species

!Gradients
DO k=1,m_s !Over species

  im=jif_s(k)
  iza=jzf_s(k)

  DO i=1,nr_r !Over radial grid

    IF(i == 1) THEN

      grt_ri(i,im)=zero
      grp_riz(i,im,iza)=zero

    ELSEIF(i == nr_r) THEN

      grt_ri(i,im)=grt_ri(nr_r-1,im)
      grp_riz(i,im,iza)=grp_riz(nr_r-1,im,iza)

    ELSE

      dr=(rhot_r(i+1)-rhot_r(i-1))*a1
      grt_ri(i,im)=(temp_ri(i+1,im)-temp_ri(i-1,im))/dr
      grp_riz(i,im,iza)=(den_riz(i+1,im,iza)*temp_ri(i+1,im) &
                        -den_riz(i-1,im,iza)*temp_ri(i-1,im))/dr

    ENDIF

  ENDDO !Over radial grid

ENDDO !Over species

!-------------------------------------------------------------------------------
!Initialize Er
!-------------------------------------------------------------------------------
DO i=2,nr_r-1

  !Use dr for flux surface averaged balances
  dr=(rhot_r(i+1)-rhot_r(i-1))*a1

  !Use drout for outside balances
  drout=rout_r(i+1)-rout_r(i-1)

  !Toroidal rotation contribution
  e_rad_r(i,2)=btout_r(i)/fhat_r(i)*vt_im_ex_r(i)
!  e_rad_r(i,2)=bpout_r(i)*vt_im_ex_r(i)*(drout/dr)
  e_rad_o_r(i,2)=e_rad_r(i,2)*(dr/drout)

  !Poloidal rotation contribution
  e_rad_r(i,3)=zero
  e_rad_o_r(i,3)=zero

  !Pressure gradient contribution
  IF(l_diag) THEN

    im=jif_s(js_diag)
    iza=jzf_s(js_diag)
    e_rad_r(i,4)=grp_riz(i,im,iza)*z_j7kv/(iza*z_coulomb*den_riz(i,im,iza))
    e_rad_o_r(i,4)=e_rad_r(i,4)*(dr/drout)

  ENDIF

  !Totals
  e_rad_r(i,1)=SUM(e_rad_r(i,2:4))
  e_rad_o_r(i,1)=SUM(e_rad_o_r(i,2:4))

ENDDO

!Edge values
e_rad_r(nr_r,1:4)=e_rad_r(nr_r-1,1:4)
e_rad_o_r(nr_r,1:4)=e_rad_o_r(nr_r-1,1:4)
erold1_r(:)=e_rad_r(:,1)
erold2_r(:)=e_rad_r(:,1)
ertemp_r(:)=e_rad_r(:,1)

!-------------------------------------------------------------------------------
!Iteration for Er
!-------------------------------------------------------------------------------
DO j=1,jitermx !Over iteration

  !Initialize arrays
  e_par_r(:)=zero
  eta_par_r(:)=zero
  omexb_o_r(:)=zero
  xm_t_im_ex_o_r(:)=zero
  xj_bs_r(:)=zero
  chi_eff_r(:,:)=zero
  chi_eff_g_r(:,:)=zero
  q_con_r(:,:)=zero
  omegas_rs(:,:)=zero
  gam_rs(:,:)=zero
  d_eff_rs(:,:)=zero
  d_n_rs(:,:)=zero
  v_nt_rs(:,:)=zero
  v_eb_rs(:,:)=zero
  d_eff_g_rs(:,:)=zero
  d_n_g_rs(:,:)=zero
  g_n_rss(:,:,:)=zero
  g_te_rs(:,:)=zero
  g_ti_rs(:,:)=zero
  v_nt_g_rs(:,:)=zero
  v_eb_g_rs(:,:)=zero
  u_par_rs(:,:)=zero
  u_p_rs(:,:)=zero
  v_t_o_rs(:,:)=zero
  xm_t_o_rs(:,:)=zero
  xm_p_o_rs(:,:)=zero
  v_p_o_rs(:,:)=zero
  v_par_o_rs(:,:)=zero
  v_per_o_rs(:,:)=zero
  omega_rs(:,:)=zero
  e_rad_rs(:,:,:)=zero

  !-----------------------------------------------------------------------------
  !Loop over radial grid, except for end points where BCs are applied
  !-----------------------------------------------------------------------------
  DO i=2,nr_r-1 !Over radial grid

    !Radial step for gradients
    dr=(rhot_r(i+1)-rhot_r(i-1))*a1
    drout=rout_r(i+1)-rout_r(i-1)

    !Radial potential gradient
    p_grphi=-e_rad_r(i,1)

    !Squeezing factor
    IF(l_squeeze) THEN

      p_gr2phi=(ertemp_r(i)*(psip_r(i+1)-psip_r(i-1))/psip_r(i) &
               -(ertemp_r(i+1)-ertemp_r(i-1)))/dr

    ENDIF

    !Parallel electric field
    p_eb=1

    !Geometric quantities
    p_b2=b2_r(i)
    p_bm2=bm2_r(i)
    p_fhat=fhat_r(i)
    p_ft=ftrap_r(i)
    p_grbm2=gr2bm2_r(i)
    p_ngrth=grth_r(i)
    p_fm(:)=fm_r(:,i)

    !Temperatures, pressures and gradients
    temp_i(:)=zero
    grt_i(:)=zero
    den_iz(:,:)=zero
    grp_iz(:,:)=zero

    temp_i(1:m_i)=temp_ri(i,1:m_i)
    grt_i(1:m_i)=grt_ri(i,1:m_i)

    DO k=1,m_s !Over species

      im=jif_s(k)
      iza=jzf_s(k)
      grp_iz(im,iza)=grp_riz(i,im,iza)
      den_iz(im,iza)=den_riz(i,im,iza)
      IF(den_iz(im,iza) < c_den) den_iz(im,iza)=1.01*c_den
      omegas_rs(i,k)=-grp_iz(im,iza)*z_j7kv*fhat_r(i)/z_coulomb &
                     /iza/den_riz(i,im,iza)/rout_r(i)/btout_r(i)
      IF(k == js_el) omegas_rs(i,k)=-omegas_rs(i,k)

    ENDDO !Over species

    !---------------------------------------------------------------------------
    !Call NCLASS
    !---------------------------------------------------------------------------
    CALL NCLASS(m_i,m_z, &
                p_b2,p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2, &
                p_grphi,p_gr2phi,p_ngrth,amu_i,grt_i,temp_i, &
                den_iz,grp_iz, &
                iflag,message, &
                L_BANANA=l_banana, &
                L_PFIRSCH=l_pfirsch, &
                L_CLASSICAL=l_classical, &
                L_POTATO=l_potato, &
                K_ORDER=k_order, &
                C_DEN=c_den, &
                C_POTB=c_potb, &
                C_POTL=c_potl, &
                AMU_F=amu_b, &
                Z_F=z_b, &
                E_F=e_b, &
                U_F=u_rb(i,:), &
                DEN_F=den_rb(i,:), &
                P_ETAP=p_etap, &
                P_JBBS=p_bsjb, &
                P_JBF=p_nbjb, &
                P_FSHLD=p_fshld, &
                M_S=m_s, &
                JM_S=jm_s, &
                JZ_S=jz_s, &
                GFL_S=gfl_s, &
                DN_S=d_n_rs(i,:), &
                VNNT_S=v_nt_rs(i,:), &
                VNEB_S=v_eb_rs(i,:), &
                DP_SS=dp_ss, &
                DT_SS=dt_ss, &
                UPAR_S=upar_s, &
                UTHETA_S=utheta_s, &
                QFL_S=qfl_s, &
                SQZ_S=sqz_rs(i,:))

    IF(iflag == 1) GOTO 9999

    !---------------------------------------------------------------------------
    !Store values
    !---------------------------------------------------------------------------
    !Parallel electrical resistivity
    eta_par_r(i)=p_etap

    !Bootstrap current density
    xj_bs_r(i)=p_bsjb/bt0

    !Beam ion current density and shielding factor
    IF(l_beam) THEN

      xj_nb_r(i)=p_nbjb/bt0
      fshld_r(i)=p_fshld

    ENDIF      

    !Parallel electric field from calculated resistivity and bootstrap
    e_par_r(i)=p_etap*(xj_r(i)-xj_nb_r(i)-p_bsjb/bt0)

    IF(k_edotb == 1) THEN

      edotb_r(i)=e_par_ex_r(i)*bt0

    ELSE

      edotb_r(i)=e_par_r(i)*bt0

    ENDIF

    !Reset <E.B> dependent results
    !Neoclassical pinch
    v_eb_rs(i,:)=edotb_r(i)*v_eb_rs(i,:)

    !Radial fluxes
    gfl_s(4,:)=edotb_r(i)*gfl_s(4,:)
    qfl_s(4,:)=edotb_r(i)*qfl_s(4,:)

    !Flow velocities
    upar_s(:,2,:)=edotb_r(i)*upar_s(:,2,:)
    utheta_s(:,2,:)=edotb_r(i)*utheta_s(:,2,:)

    !Poloidal flow contribution to radial electric field
    IF(l_diag) THEN

      !Increment parallel flow by the Er from the poloidal flow on first loop
      IF(j == 1) upar_s(1,1,:)=upar_s(1,1,:)-btout_r(i)**2*SUM(utheta_s(1,:,js_diag))

      im=jm_s(js_diag)
      iza=jz_s(js_diag)
      e_rad_r(i,3)=-btout_r(i)**2*SUM(utheta_s(1,1:3,js_diag))/fhat_r(i)
      e_rad_o_r(i,3)=e_rad_r(i,3)*(dr/drout)

    ENDIF

    !Total radial electric field
    e_rad_r(i,1)=SUM(e_rad_r(i,2:4))
    e_rad_o_r(i,1)=SUM(e_rad_o_r(i,2:4))

    !Mass-density-weighted radial electric field
    rdum(1:5)=zero

    DO k=1,m_s !Over ion species

      IF(k /= js_el) THEN

        im=jm_s(k)
        iza=jz_s(k)
        e_rad_rs(i,k,1)=e_rad_r(i,1)
        e_rad_rs(i,k,3)=-fhat_r(i)*rm2_r(i)*(psip_r(i)/2/z_pi)**2 &
                        *SUM(utheta_s(1,1:3,k))
        e_rad_rs(i,k,4)=grp_iz(im,iza)*z_j7kv/iza/z_coulomb/den_iz(im,iza)
        e_rad_rs(i,k,2)=e_rad_rs(i,k,1)-e_rad_rs(i,k,3)-e_rad_rs(i,k,4)
        rdum(1:4)=rdum(1:4)+amu_i(im)*den_iz(im,iza)*e_rad_rs(i,k,1:4)
        rdum(5)=rdum(5)+amu_i(im)*den_iz(im,iza)

      ENDIF

    ENDDO !Over ion species

    e_rad_rs(i,1,1:4)=rdum(1:4)/rdum(5)

    !Particle fluxes
    d_n_g_rs(i,:)=d_n_rs(i,:)/grho2_r(i)
    v_nt_g_rs(i,:)=v_nt_rs(i,:)/grho1_r(i)
    v_eb_g_rs(i,:)=v_eb_rs(i,:)/grho1_r(i)

    DO k=1,m_s !Over species

      im=jm_s(k)
      iz=jz_s(k)
      iza=IABS(iz)
      gam_rs(i,k)=SUM(gfl_s(:,k))
      grad=(den_riz(i+1,im,iza)-den_riz(i-1,im,iza))/dr

      !Protect against small density gradient
      IF(ABS(grad*a1/den_riz(i,im,iza)) > 1.0e-5) THEN

        d_eff_rs(i,k)=-gam_rs(i,k)/grad

      ELSE

        d_eff_rs(i,k)=zero

      ENDIF

      d_eff_g_rs(i,k)=d_eff_rs(i,k)/grho2_r(i)

      !Average parallel flow
      u_par_rs(i,k)=SUM(upar_s(1,1:3,k))
      ppr=grp_iz(im,iza)*z_j7kv/(z_coulomb*iz*den_iz(im,iza))

      !Average poloidal flow
      u_p_rs(i,k)=SUM(utheta_s(1,1:3,k))

      !Toroidal flow velocity on outside
      v_t_o_rs(i,k)=u_p_rs(i,k)*btout_r(i) &
                    -fhat_r(i)*(ppr-e_rad_r(i,1))/btout_r(i)

      !Toroidal Mach number on outside
      vtherm=SQRT(2*temp_i(im)*z_j7kv/amu_i(im)/z_pmass)
      xm_t_o_rs(i,k)=v_t_o_rs(i,k)/vtherm

      !Poloidal flow velocity on outside
      v_p_o_rs(i,k)=u_p_rs(i,k)*bpout_r(i)

      !Poloidal Mach number on outside
      xm_p_o_rs(i,k)=xm_t_o_rs(i,k)-e_rad_o_r(i,1)/bpout_r(i)/vtherm

      !Parallel flow velocity on outside
      b=btout_r(i)*SQRT(one+(bpout_r(i)/btout_r(i))**2)
      v_par_o_rs(i,k)=u_p_rs(i,k)*b-(ppr-e_rad_r(i,1))*fhat_r(i)/b

      !Perpendicular flow velocity on outside
      v_per_o_rs(i,k)=(ppr-e_rad_r(i,1))*fhat_r(i)*bpout_r(i)/btout_r(i)/b

      !Exponent factors for particle flux
      g_te_rs(i,k)=-(dp_ss(1,k)+dt_ss(1,k))/dp_ss(k,k)

      DO kk=1,m_s !Over driving species

        g_ti_rs(i,k)=g_ti_rs(i,k)-(dp_ss(kk,k)+dt_ss(kk,k))/dp_ss(k,k)
        g_n_rss(i,kk,k)=-dp_ss(kk,k)/dp_ss(k,k)

      ENDDO !Over driving species

    ENDDO !Over isotopes

    IF(l_diag) THEN

      vtherm=SQRT(2*temp_i(m_i)*z_j7kv/amu_i(m_i)/z_pmass)
      xm_t_im_ex_o_r(i)=vt_im_ex_r(i)/vtherm 

    ENDIF

    !---------------------------------------------------------------------------
    !Conduction fluxes and conductivities
    !---------------------------------------------------------------------------
    !Electrons
    q_con_r(i,1)=SUM(qfl_s(:,1))

    !Protect against small ion temperature gradient
    IF(ABS(grt_i(1)*a1/temp_i(1)) > 1.0e-5_rspec) THEN

      chi_eff_r(i,1)=-q_con_r(i,1)/den_riz(i,1,1)/z_j7kv/grt_i(1)

    ELSE

      chi_eff_r(i,1)=zero

    ENDIF

    chi_eff_g_r(i,1)=chi_eff_r(i,1)/grho2_r(i)

    !Ions - sum over all species
    dent=zero

    DO k=2,m_s !Over ion species

      ima=jm_s(k)
      iza=jz_s(k)
      q_con_r(i,2)=q_con_r(i,2)+SUM(qfl_s(:,k))
      dent=dent+den_iz(ima,iza)

    ENDDO !Over ion species

    !Protect against small ion temperature gradient
    IF(ABS(grt_i(2)*a1/temp_i(2)) > 1.0e-5_rspec) THEN

      chi_eff_r(i,2)=-q_con_r(i,2)/dent/z_j7kv/grt_i(2)

    ELSE

      chi_eff_r(i,2)=zero

    ENDIF

    chi_eff_g_r(i,2)=chi_eff_r(i,2)/grho2_r(i)

  ENDDO !Over radial grid

  !Without orbit squeezing, the linear correction to Er above is OK
  jiter=0

  IF(.NOT. l_squeeze) THEN

    erref=1.0e3*te_r(1)/a1

    !Check convergence
    DO i=2,nr_r-1 !Over radial grid

      IF(ABS((e_rad_r(i,1)-erold1_r(i))/erref) > toler) jiter=1
      ertemp_r(i)=(e_rad_r(i,1)+2*erold1_r(i)+erold2_r(i))/4
      erold2_r(i)=erold1_r(i)
      erold1_r(i)=e_rad_r(i,1)

    ENDDO !Over radial grid

  ENDIF

  !Toroidal rotation constant
  DO i=2,nr_r-1 !Over radial grid

    omega_rs(i,1:m_s)=omegas_rs(i,1:m_s)+e_rad_r(i,1)*fhat_r(i)/rout_r(i) &
                     /btout_r(i)

  ENDDO !Over radial grid

  !-------------------------------------------------------------------------------
  !Axial and edge values
  !-------------------------------------------------------------------------------
  !Fluxes, transport coefficients and flow velocities
  !Axis
  d_eff_rs(1,1:m_s)=d_eff_rs(2,1:m_s)
  d_n_rs(1,1:m_s)=d_n_rs(2,1:m_s)
  v_nt_rs(1,1:m_s)=zero
  v_eb_rs(1,1:m_s)=zero
  d_eff_g_rs(1,1:m_s)=d_eff_rs(1,1:m_s)/grho2_r(1)
  d_n_g_rs(1,1:m_s)=d_n_rs(1,1:m_s)/grho2_r(1)
  v_nt_g_rs(1,1:m_s)=zero
  v_eb_g_rs(1,1:m_s)=zero
  sqz_rs(1,:)=sqz_rs(2,:)
  u_par_rs(1,1:m_s)=u_par_rs(2,1:m_s)
  u_p_rs(1,1:m_s)=u_p_rs(2,1:m_s)
  v_t_o_rs(1,1:m_s)=v_t_o_rs(2,1:m_s)
  xm_t_o_rs(1,1:m_s)=xm_t_o_rs(2,1:m_s)
  v_p_o_rs(1,1:m_s)=zero
  xm_p_o_rs(1,1:m_s)=xm_p_o_rs(2,1:m_s)
  v_par_o_rs(1,1:m_s)=v_par_o_rs(2,1:m_s)
  v_per_o_rs(1,1:m_s)=zero
  g_te_rs(1,1:m_s)=g_te_rs(2,1:m_s)
  g_ti_rs(1,1:m_s)=g_ti_rs(2,1:m_s)
  g_n_rss(1,1:m_s,1:m_s)=g_n_rss(2,1:m_s,1:m_s)
  omega_rs(1,1:m_s)=omega_rs(2,1:m_s)

  IF(l_beam) THEN

    fshld_r(1)=one-z_b(1)/zeff_r(1)
    xj_nb_r(1)=xj_nb_r(2)*fshld_r(1)/fshld_r(2)

  ENDIF

  !Edge
  gam_rs(nr_r,1:m_s)=gam_rs(nr_r-1,1:m_s)
  d_eff_rs(nr_r,1:m_s)=d_eff_rs(nr_r-1,1:m_s)
  d_n_rs(nr_r,1:m_s)=d_n_rs(nr_r-1,1:m_s)
  v_nt_rs(nr_r,1:m_s)=v_nt_rs(nr_r-1,1:m_s)
  v_eb_rs(nr_r,1:m_s)=v_eb_rs(nr_r-1,1:m_s)
  d_eff_g_rs(nr_r,1:m_s)=d_eff_rs(nr_r,1:m_s)/grho2_r(nr_r)
  d_n_g_rs(nr_r,1:m_s)=d_n_rs(nr_r,1:m_s)/grho2_r(nr_r)
  v_nt_g_rs(nr_r,1:m_s)=v_nt_rs(nr_r,1:m_s)/grho1_r(nr_r)
  v_eb_g_rs(nr_r,1:m_s)=v_eb_rs(nr_r,1:m_s)/grho1_r(nr_r)
  sqz_rs(nr_r,:)=sqz_rs(nr_r-1,:)
  u_par_rs(nr_r,1:m_s)=u_par_rs(nr_r-1,1:m_s)
  u_p_rs(nr_r,1:m_s)=u_p_rs(nr_r-1,1:m_s)
  v_t_o_rs(nr_r,1:m_s)=v_t_o_rs(nr_r-1,1:m_s)
  xm_t_o_rs(nr_r,1:m_s)=xm_t_o_rs(nr_r-1,1:m_s)
  v_p_o_rs(nr_r,1:m_s)=v_p_o_rs(nr_r-1,1:m_s)
  v_par_o_rs(nr_r,1:m_s)=v_par_o_rs(nr_r-1,1:m_s)
  v_per_o_rs(nr_r,1:m_s)=v_per_o_rs(nr_r-1,1:m_s)
  xm_p_o_rs(nr_r,1:m_s)=xm_p_o_rs(nr_r-1,1:m_s)
  g_te_rs(nr_r,1:m_s)=g_te_rs(nr_r-1,1:m_s)
  g_ti_rs(nr_r,1:m_s)=g_ti_rs(nr_r-1,1:m_s)
  g_n_rss(nr_r,1:m_s,1:m_s)=g_n_rss(nr_r-1,1:m_s,1:m_s)
  omega_rs(nr_r,1:m_s)=omega_rs(nr_r-1,1:m_s)

  IF(l_beam) THEN

    fshld_r(nr_r)=fshld_r(nr_r-1)
    xj_nb_r(nr_r)=xj_nb_r(nr_r-1)

  ENDIF

  !Experimental flow velocities for diagnostic impurity
  !Axis
  xm_t_im_ex_o_r(1)=xm_t_im_ex_o_r(2)

  !Edge
  xm_t_im_ex_o_r(nr_r)=xm_t_im_ex_o_r(nr_r-1)

  !Electron conduction
  !Axis
  q_con_r(1,1)=zero
  chi_eff_r(1,1)=chi_eff_r(2,1)
  chi_eff_g_r(1,1)=chi_eff_r(1,1)/grho2_r(1)

  !Edge
  q_con_r(nr_r,1)=q_con_r(nr_r-1,1)
  chi_eff_r(nr_r,1)=chi_eff_r(nr_r-1,1)
  chi_eff_g_r(nr_r,1)=chi_eff_r(nr_r,1)/grho2_r(nr_r)

  !Ion conduction
  !Axis
  q_con_r(1,1)=zero
  chi_eff_r(1,2)=chi_eff_r(2,2)
  chi_eff_g_r(1,2)=chi_eff_r(1,2)/grho2_r(1)

  !Edge
  q_con_r(nr_r,2)=q_con_r(nr_r-1,2)
  chi_eff_r(nr_r,2)=chi_eff_r(nr_r-1,2)
  chi_eff_g_r(nr_r,2)=chi_eff_r(nr_r,2)/grho2_r(nr_r)

  !Parallel electrical resistivity
  !Axis
  eta_par_r(1)=eta_par_r(2)

  !Edge
  eta_par_r(nr_r)=eta_par_r(nr_r-1)+(eta_par_r(nr_r-1)-eta_par_r(nr_r-2)) &
                  *(rhot_r(nr_r)-rhot_r(nr_r-1))/(rhot_r(nr_r-1)-rhot_r(nr_r-2))

  !Bootstrap current density
  !Axis
  xj_bs_r(1)=zero

  !Edge
  xj_bs_r(nr_r)=xj_bs_r(nr_r-1)+(xj_bs_r(nr_r-1)-xj_bs_r(nr_r-2)) &
                *(rhot_r(nr_r)-rhot_r(nr_r-1))/(rhot_r(nr_r-1)-rhot_r(nr_r-2))

  !Parallel electric field
  !Axis
  e_par_r(1)=e_par_r(2)

  !Edge
  e_par_r(nr_r)=e_par_r(nr_r-1)+(e_par_r(nr_r-1)-e_par_r(nr_r-2)) &
                *(rhot_r(nr_r)-rhot_r(nr_r-1))/(rhot_r(nr_r-1)-rhot_r(nr_r-2))

  !Radial electric field
  !Axis
  e_rad_r(1,1:4)=zero
  e_rad_rs(1,1:m_s,1:4)=zero
  e_rad_o_r(1,1:4)=zero

  !Edge
  e_rad_r(nr_r,1:4)=e_rad_r(nr_r-1,1:4)
  e_rad_rs(nr_r,1:m_s,1:4)=e_rad_rs(nr_r-1,1:m_s,1:4)
  e_rad_o_r(nr_r,1:4)=e_rad_o_r(nr_r-1,1:4)

  !ExB shear damping on outside
  DO i=3,nr_r-1 !Over radial grid

    omexb_o_r(i)=(rout_r(i)*bpout_r(i))**2 &
                 /SQRT(btout_r(i)**2+bpout_r(i)**2) &
                 *(e_rad_o_r(i+1,1)/rout_r(i+1)/bpout_r(i+1) &
                 -e_rad_o_r(i-1,1)/rout_r(i-1)/bpout_r(i-1)) &
                 /(psi_r(i+1)-psi_r(i-1))

  ENDDO !Over radial grid

  omexb_o_r(2)=omexb_o_r(3)
  omexb_o_r(1)=zero
  omexb_o_r(nr_r)=omexb_o_r(nr_r-1)
  IF(jiter == 0) GOTO 9999

ENDDO !Over iteration

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FORCEBAL_NCLASS

SUBROUTINE FORCEBAL_OUT(mx_npro,a1,r0,npro,valpro,namepro,unitpro,descpro)
!-------------------------------------------------------------------------------
!FORCEBAL_OUT loads the various radial output arrays into a two-dimensional
!  array and also loads label, description, and unit arrays
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE FORCEBAL_DATA_MOD
IMPLICIT NONE

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  mx_npro                !maximum number of radial arrays [-]

REAL(KIND=rspec), INTENT(IN) :: &
  a1,                   & !reference minor radius [m]
  r0                      !geometric center [m]

!Declaration of output variables
INTEGER, INTENT(OUT) :: &
  npro                   !number of radial arrays [-]

CHARACTER(len=*), INTENT(OUT) :: &
  namepro(mx_npro),    & !variable names [character]
  unitpro(mx_npro),    & !variable units [character]
  descpro(mx_npro)       !variable descriptions [character]

REAL(KIND=rspec), INTENT(OUT) :: &
  valpro(nr_r,mx_npro)   !radial values of variables [see unitpro]

!-------------------------------------------------------------------------------
!Declaration of local variables
CHARACTER(len=5) :: &
  cidz_s(m_s),ciz

INTEGER :: &
   j,j1,k,k1,kscs

REAL(KIND=rspec), PARAMETER :: &
  one=1.0_rspec,       & !REAL 1
  zero=0.0_rspec         !REAL 0

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output arrays
namepro(1:mx_npro)=''
unitpro(1:mx_npro)=''
descpro(1:mx_npro)=''
valpro(1:nr_r,1:mx_npro)=zero
npro=0

!Set species (isotope + charge state) identifier
cidz_s(:)=''

DO j=1,m_s !Over species

  WRITE(ciz,'(i5)') jzf_s(j)
  cidz_s(j)=TRIM(ADJUSTL(cid_i(jif_s(j))))//TRIM(ADJUSTL(ciz))

ENDDO !Over species

!Species identifiers - set max index of scs species
IF(l_diag) THEN

  kscs=2+n_scsion

ELSE

  kscs=1+n_scsion

ENDIF

!-------------------------------------------------------------------------------
!Radial grids
!-------------------------------------------------------------------------------
npro=npro+1
namepro(npro)='rho_p'
unitpro(npro)='-'
descpro(npro)='Normalized poloidal flux grid - proportional to poloidal flux'
valpro(1:nr_r,npro)=rhop_r(1:nr_r)

npro=npro+1
namepro(npro)='rho_t'
unitpro(npro)='-'
descpro(npro)='Normalized toroidal flux grid - ' &
              //'proportional to square root toroidal flux'
valpro(1:nr_r,npro)=rhot_r(1:nr_r)

npro=npro+1
namepro(npro)='r_p'
unitpro(npro)='m'
descpro(npro)='a1 * normalized poloidal flux grid'
valpro(1:nr_r,npro)=a1*rhop_r(1:nr_r)

npro=npro+1
namepro(npro)='r_t'
unitpro(npro)='m'
descpro(npro)='a1 * normalized toroidal flux grid'
valpro(1:nr_r,npro)=a1*rhot_r(1:nr_r)

npro=npro+1
namepro(npro)='R_in'
unitpro(npro)='m'
descpro(npro)='Major radius inboard midplane'
valpro(1:nr_r,npro)=rin_r(1:nr_r)

npro=npro+1
namepro(npro)='R_o'
unitpro(npro)='m'
descpro(npro)='Major radius outboard midplane'
valpro(1:nr_r,npro)=rout_r(1:nr_r)

!-------------------------------------------------------------------------------
!Magnetic fields and fluxes
!-------------------------------------------------------------------------------
npro=npro+1
namepro(npro)='B_p_o'
unitpro(npro)='T'
descpro(npro)='Poloidal magnetic field outboard midplane'
valpro(1:nr_r,npro)=bpout_r(1:nr_r)

npro=npro+1
namepro(npro)='B_t_o'
unitpro(npro)='T'
descpro(npro)='Toroidal magnetic field outboard midplane'
valpro(1:nr_r,npro)=btout_r(1:nr_r)

npro=npro+1
namepro(npro)='Phi_t'
unitpro(npro)='Wb'
descpro(npro)='Toroidal magnetic flux'
valpro(1:nr_r,npro)=phit_r(1:nr_r)

npro=npro+1
namepro(npro)='Psi'
unitpro(npro)='Wb/rad'
descpro(npro)='Poloidal magnetic flux divided by 2 pi'
valpro(1:nr_r,npro)=psi_r(1:nr_r)

npro=npro+1
namepro(npro)='q'
unitpro(npro)='-'
descpro(npro)='Safety factor'
valpro(1:nr_r,npro)=ABS(q_r(1:nr_r))

!-------------------------------------------------------------------------------
!Flux functions and metrics
!-------------------------------------------------------------------------------
npro=npro+1
namepro(npro)='F'
unitpro(npro)='A'
descpro(npro)='Poloidal current external to surface'
valpro(1:nr_r,npro)=f_r(1:nr_r)

npro=npro+1
namepro(npro)='f_trap'
unitpro(npro)='-'
descpro(npro)='Trapped particle fraction'
valpro(1:nr_r,npro)=ftrap_r(1:nr_r)

npro=npro+1
namepro(npro)='grad_rho_sq'
unitpro(npro)='-'
descpro(npro)='Metric, <grad(rho)**2>'
valpro(1:nr_r,npro)=grho2_r(1:nr_r)

npro=npro+1
namepro(npro)='grad_rho'
unitpro(npro)='-'
descpro(npro)='Metric, <|grad(rho)|>'
valpro(1:nr_r,npro)=grho1_r(1:nr_r)

npro=npro+1
namepro(npro)='dV7dr_t'
unitpro(npro)='m**2'
descpro(npro)='Radial derivative of volume, dV/drho/a1'
valpro(1:nr_r,npro)=vp_r(1:nr_r)

npro=npro+1
namepro(npro)='Vol'
unitpro(npro)='m**3'
descpro(npro)='Enclosed plasma volume'
valpro(1:nr_r,npro)=vol_r(1:nr_r)

!-------------------------------------------------------------------------------
!Plasma profiles
!-------------------------------------------------------------------------------
DO j1=1,m_s !Over species

  IF(l_reduce_out) THEN

    IF(j1 > m_i) CYCLE

    IF(j1 <= kscs) THEN

      j=j1

    ELSE

      j=j+izmax_i(j1)

    ENDIF

  ELSE

    j=j1

  ENDIF    


  npro=npro+1
  namepro(npro)='den_'//TRIM(cidz_s(j))
  unitpro(npro)='/m**3'
  descpro(npro)='Density of '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=den_riz(1:nr_r,jif_s(j),jzf_s(j))

  npro=npro+1
  namepro(npro)='T_'//TRIM(cidz_s(j))
  unitpro(npro)='keV'
  descpro(npro)='Temperature of '//TRIM(cidz_s(j))
  !EHAB MODIFIED LINES BELOW
 !valpro(1:nr_r,npro)=temp_ri(1:npro,jif_s(j))
  valpro(1:nr_r,npro)=temp_ri(1:nr_r,jif_s(j))
  !EHAB MODIFIED LINES ABOVE

  npro=npro+1
  namepro(npro)='prs_'//TRIM(cidz_s(j))
  unitpro(npro)='keV/m**3'
  descpro(npro)='Pressure of '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=den_riz(1:nr_r,jif_s(j),jzf_s(j))*temp_ri(1:nr_r,jif_s(j))

ENDDO !Over species

IF(MAXVAL(zeff_ex_r) >= one) THEN

  npro=npro+1
  namepro(npro)='Z_eff_ex'
  unitpro(npro)='-'
  descpro(npro)='Effective charge, Zeff - experimental data'
  valpro(1:nr_r,npro)=zeff_ex_r(1:nr_r)

ENDIF

npro=npro+1
namepro(npro)='Z_eff'
unitpro(npro)='-'
descpro(npro)='Effective charge, Zeff'
valpro(1:nr_r,npro)=zeff_r(1:nr_r)

!-------------------------------------------------------------------------------
!Particle fluxes, diffusivities and pinches
!-------------------------------------------------------------------------------
DO j1=1,m_s !Over species

  IF(l_reduce_out) THEN

    IF(j1 > m_i) CYCLE

    IF(j1 <= kscs) THEN

      j=j1

    ELSE

      j=j+izmax_i(j1)

    ENDIF

  ELSE

    j=j1

  ENDIF    

  npro=npro+1
  namepro(npro)='Sqz_'//TRIM(cidz_s(j))
  unitpro(npro)='-'
  descpro(npro)='Orbit squeezing factor for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=sqz_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='Gam_'//TRIM(cidz_s(j))
  unitpro(npro)='#/m**2/s'
  descpro(npro)='Particle flux of '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=gam_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='D_eff_'//TRIM(cidz_s(j))
  unitpro(npro)='m**2/s'
  descpro(npro)='Effective particle diffusivity for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=d_eff_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='D_eff_g_'//TRIM(cidz_s(j))
  unitpro(npro)='m**2/s'
  descpro(npro)='Effective particle diffusivity/g for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=d_eff_g_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='D_n_'//TRIM(cidz_s(j))
  unitpro(npro)='m**2/s'
  descpro(npro)='Diagonal particle diffusivity for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=d_n_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='v_n_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Particle pinch from grad(T,n) for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_nt_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='v_wp_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Ware pinch for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_eb_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='v_nwp_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Total pinch for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_eb_rs(1:nr_r,j)+v_nt_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='D_n_g_'//TRIM(cidz_s(j))
  unitpro(npro)='m**2/s'
  descpro(npro)='Diagonal particle diffusivity/g for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=d_n_g_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='v_n_g_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Particle pinch from grad(T,n)/g for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_nt_g_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='v_wp_g_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Ware pinch/g for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_eb_g_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='v_nwp_g_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Total pinch/g for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_eb_g_rs(1:nr_r,j)+v_nt_g_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='g_Te_'//TRIM(cidz_s(j))
  unitpro(npro)='-'
  descpro(npro)='Te exponent in the density profile of '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=g_te_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='g_Ti_'//TRIM(cidz_s(j))
  unitpro(npro)='-'
  descpro(npro)='Ti exponent in the density profile of '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=g_ti_rs(1:nr_r,j)

  DO k1=1,m_s !Over target species

    IF(l_reduce_out) THEN

      IF(k1 > m_i) CYCLE

      IF(k1 <= kscs) THEN

        k=k1

      ELSE

        k=k+izmax_i(k1)

      ENDIF

    ELSE

      k=k1

    ENDIF    

    !g factor - first index is gradient and second is flux
    npro=npro+1
    namepro(npro)='g_n'//TRIM(cidz_s(k))//'_'//TRIM(cidz_s(j))
    unitpro(npro)='-'
    descpro(npro)=TRIM(cidz_s(k))//' exponent in the density profile of ' &
                  //TRIM(cidz_s(j))
    valpro(1:nr_r,npro)=g_n_rss(1:nr_r,k,j)

  ENDDO !Over target species

ENDDO !Over species

!-------------------------------------------------------------------------------
!Conduction heat fluxes and conductivities
!-------------------------------------------------------------------------------
npro=npro+1
namepro(npro)='q_con_e'
unitpro(npro)='w/m**2'
descpro(npro)='Electron thermal conduction flux'
valpro(1:nr_r,npro)=q_con_r(1:nr_r,1)

npro=npro+1
namepro(npro)='q_con_i'
unitpro(npro)='w/m**2'
descpro(npro)='Ion thermal conduction flux'
valpro(1:nr_r,npro)=q_con_r(1:nr_r,2)

npro=npro+1
namepro(npro)='chi_eff_e'
unitpro(npro)='m**2/s'
descpro(npro)='Effective electron thermal conductivity'
valpro(1:nr_r,npro)=chi_eff_r(1:nr_r,1)

npro=npro+1
namepro(npro)='chi_eff_i'
unitpro(npro)='m**2/s'
descpro(npro)='Effective ion thermal conductivity'
valpro(1:nr_r,npro)=chi_eff_r(1:nr_r,2)

npro=npro+1
namepro(npro)='chi_eff_g_e'
unitpro(npro)='m**2/s'
descpro(npro)='Effective electron thermal conductivity/g'
valpro(1:nr_r,npro)=chi_eff_g_r(1:nr_r,1)

npro=npro+1
namepro(npro)='chi_eff_g_i'
unitpro(npro)='m**2/s'
descpro(npro)='Effective ion thermal conductivity/g'
valpro(1:nr_r,npro)=chi_eff_g_r(1:nr_r,2)

!-------------------------------------------------------------------------------
!Electrical resistivity
!-------------------------------------------------------------------------------
npro=npro+1
namepro(npro)='eta_par_r'
unitpro(npro)='Ohm*m'
descpro(npro)='Parallel electrical resistivity'
valpro(1:nr_r,npro)=eta_par_r(1:nr_r)

!-------------------------------------------------------------------------------
!Parallel currents
!-------------------------------------------------------------------------------
npro=npro+1
namepro(npro)='J'
unitpro(npro)='A/m**2'
descpro(npro)='Total parallel current density, <J.B>/Bt0'
valpro(1:nr_r,npro)=xj_r(1:nr_r)

npro=npro+1
namepro(npro)='I'
unitpro(npro)='A'
descpro(npro)='Enclosed total toroidal current'
valpro(1:nr_r,npro)=cur_r(1:nr_r)

npro=npro+1
namepro(npro)='J_bs'
unitpro(npro)='A/m**2'
descpro(npro)='Bootstrap parallel current density, <Jbs.B>/Bt0'
valpro(1:nr_r,npro)=xj_bs_r(1:nr_r)

npro=npro+1
namepro(npro)='I_bs'
unitpro(npro)='A'
descpro(npro)='Enclosed bootstrap toroidal current'
valpro(1:nr_r,npro)=cur_bs_r(1:nr_r)

IF(ABS(xj_ex_r(1)) > zero) THEN

  npro=npro+1
  namepro(npro)='J_ex'
  unitpro(npro)='A/m**2'
  descpro(npro)='Expt total parallel current density, <Jex.B>/Bt0'
  valpro(1:nr_r,npro)=xj_ex_r(1:nr_r)

  npro=npro+1
  namepro(npro)='I_ex'
  unitpro(npro)='A'
  descpro(npro)='Expt enclosed total toroidal current'
  valpro(1:nr_r,npro)=cur_ex_r(1:nr_r)

ENDIF

IF(l_beam) THEN

  npro=npro+1
  namepro(npro)='fspitz'
  unitpro(npro)='-'
  descpro(npro)='Spitzer beam current shielding factor'
  valpro(1:nr_r,npro)=one-z_b(1)/zeff_r(1:nr_r)

  npro=npro+1
  namepro(npro)='fms'
  unitpro(npro)='-'
  descpro(npro)='Mikkelsen-Singer beam current shielding factor'
  valpro(1:nr_r,npro)=one-z_b(1)/zeff_r(1:nr_r)*(one &
                      -(1.55_rspec+0.85/zeff_r(1:nr_r)) &
                       *SQRT(rhot_r(1:nr_r)*a1/r0) &
                      +(0.20_rspec+1.55/zeff_r(1:nr_r)) &
                       *rhot_r(1:nr_r)*a1/r0)

  npro=npro+1
  namepro(npro)='ftad'
  unitpro(npro)='-'
  descpro(npro)='Tani-Azumi-Devoto beam current shielding factor'
  valpro(1:nr_r,npro)=one-z_b(1)/zeff_r(1:nr_r)*(one-ftrap_r(1:nr_r))

  npro=npro+1
  namepro(npro)='fshld'
  unitpro(npro)='-'
  descpro(npro)='Beam current shielding factor'
  valpro(1:nr_r,npro)=fshld_r(1:nr_r)

  npro=npro+1
  namepro(npro)='J_nb'
  unitpro(npro)='A/m**2'
  descpro(npro)='Beam parallel current density, <Jnb.B>/Bt0'
  valpro(1:nr_r,npro)=xj_nb_r(1:nr_r)

  npro=npro+1
  namepro(npro)='I_nb'
  unitpro(npro)='A'
  descpro(npro)='Beam enclosed total toroidal current'
  valpro(1:nr_r,npro)=cur_nb_r(1:nr_r)

ENDIF

!-------------------------------------------------------------------------------
!Parallel electric field
!-------------------------------------------------------------------------------
npro=npro+1
namepro(npro)='E_par'
unitpro(npro)='V/m'
descpro(npro)='Parallel electric field, <E.B/Bt0>'
valpro(1:nr_r,npro)=e_par_r(1:nr_r)

IF(ABS(e_par_ex_r(1)) > zero) THEN

  npro=npro+1
  namepro(npro)='E_par_ex'
  unitpro(npro)='V/m'
  descpro(npro)='Expt parallel electric field, <E.B/Bt0>'
  valpro(1:nr_r,npro)=e_par_ex_r(1:nr_r)

ENDIF

!-------------------------------------------------------------------------------
!Radial electric field
!-------------------------------------------------------------------------------
!Diagnostic impurity
IF(l_diag) THEN

  npro=npro+1
  namepro(npro)='E_rad_tot'
  unitpro(npro)='V/m'
  descpro(npro)='Radial electric field, Erho=-d(Phi)/drho/a1, from ' &
                //TRIM(cidz_s(js_diag))//' force balance'
  valpro(1:nr_r,npro)=e_rad_r(1:nr_r,1)

  npro=npro+1
  namepro(npro)='E_rad_t'
  unitpro(npro)='V/m'
  descpro(npro)='Vtor contribution to Erho from ' &
                //TRIM(cidz_s(js_diag))//' force balance'
  valpro(1:nr_r,npro)=e_rad_r(1:nr_r,2)

  npro=npro+1
  namepro(npro)='E_rad_p'
  unitpro(npro)='V/m'
  descpro(npro)='Vpol contribution to Erho from ' &
                //TRIM(cidz_s(js_diag))//' force balance'
  valpro(1:nr_r,npro)=e_rad_r(1:nr_r,3)

  npro=npro+1
  namepro(npro)='E_rad_prs'
  unitpro(npro)='V/m'
  descpro(npro)='Diamagnetic contribution to Erho from ' &
                //TRIM(cidz_s(js_diag))//' force balance'
  valpro(1:nr_r,npro)=e_rad_r(1:nr_r,4)

  npro=npro+1
  namepro(npro)='E_rad_tot_o'
  unitpro(npro)='V/m'
  descpro(npro)='Er on the outside midplane for '//TRIM(cidz_s(js_diag))
  valpro(1:nr_r,npro)=e_rad_o_r(1:nr_r,1)

  npro=npro+1
  namepro(npro)='E_rad_t_o'
  unitpro(npro)='V/m'
  descpro(npro)='Vtor contribution to Er on the outside midplane for ' &
                //TRIM(cidz_s(js_diag))
  valpro(1:nr_r,npro)=e_rad_o_r(1:nr_r,2)

  npro=npro+1
  namepro(npro)='E_rad_p_o'
  unitpro(npro)='V/m'
  descpro(npro)='Vpol contribution to Er on the outside midplane for ' &
                //TRIM(cidz_s(js_diag))
  valpro(1:nr_r,npro)=e_rad_o_r(1:nr_r,3)

  npro=npro+1
  namepro(npro)='E_rad_prs_o'
  unitpro(npro)='V/m'
  descpro(npro)='Diamagnetic contribution to Er on the outside midplane for ' &
                //TRIM(cidz_s(js_diag))
  valpro(1:nr_r,npro)=e_rad_o_r(1:nr_r,4)

  npro=npro+1
  namepro(npro)='omega_exb_o'
  unitpro(npro)='/s'
  descpro(npro)='Hahm-Burell ExB shear damping for '//TRIM(cidz_s(js_diag))
  valpro(1:nr_r,npro)=omexb_o_r(1:nr_r)

ENDIF

!Mass density weighted
npro=npro+1
namepro(npro)='E_rad_m_tot'
unitpro(npro)='V/m'
descpro(npro)='Mass-density weighted E_rho'
valpro(1:nr_r,npro)=e_rad_rs(1:nr_r,1,1)

npro=npro+1
namepro(npro)='E_rad_m_t'
unitpro(npro)='V/m'
descpro(npro)='Vtor contribution to mass-density weighted E_rho'
valpro(1:nr_r,npro)=e_rad_rs(1:nr_r,1,2)

npro=npro+1
namepro(npro)='E_rad_m_p'
unitpro(npro)='V/m'
descpro(npro)='Vpol contribution to mass-density weighted E_rho'
valpro(1:nr_r,npro)=e_rad_rs(1:nr_r,1,3)

npro=npro+1
namepro(npro)='E_rad_m_prs'
unitpro(npro)='V/m'
descpro(npro)='Diamagnetic contribution to mass-density weighted '//' E_rho'
valpro(1:nr_r,npro)=e_rad_rs(1:nr_r,1,4)

!Each ion species
DO j1=2,m_s !Over ion species

  IF(l_reduce_out) THEN

    IF(j1 > m_i) CYCLE

    IF(j1 <= kscs) THEN

      j=j1

    ELSE

      j=j+izmax_i(j1)

    ENDIF

  ELSE

    j=j1

  ENDIF    

  npro=npro+1
  namepro(npro)='E_rad_'//TRIM(cidz_s(j))//'_t'
  unitpro(npro)='V/m'
  descpro(npro)='Vtor contribution to Erho in the '//TRIM(cidz_s(j)) &
                //' force balance'
  valpro(1:nr_r,npro)=e_rad_rs(1:nr_r,j,2)

  npro=npro+1
  namepro(npro)='E_rad_'//TRIM(cidz_s(j))//'_p'
  unitpro(npro)='V/m'
  descpro(npro)='Vpol contribution to Erho in the '//TRIM(cidz_s(j)) &
                //' force balance'
  valpro(1:nr_r,npro)=e_rad_rs(1:nr_r,j,3)

  npro=npro+1
  namepro(npro)='E_rad_'//TRIM(cidz_s(j))//'_prs'
  unitpro(npro)='V/m'
  descpro(npro)='Diamagnetic contribution to Erho in the '//TRIM(cidz_s(j)) &
                //' force balance'
  valpro(1:nr_r,npro)=e_rad_rs(1:nr_r,j,4)

ENDDO !Over ion species

!-------------------------------------------------------------------------------
!Particle toroidal and poloidal flow velocities on the outside
!-------------------------------------------------------------------------------
DO j1=1,m_s !Over species

  IF(l_reduce_out) THEN

    IF(j1 > m_i) CYCLE

    IF(j1 <= kscs) THEN

      j=j1

    ELSE

      j=j+izmax_i(j1)

    ENDIF

  ELSE

    j=j1

  ENDIF
  
  npro=npro+1
  namepro(npro)='v_t_o_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Vtor on the outside midplane for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_t_o_rs(1:nr_r,j)

  IF(j == js_diag .AND. &
     l_diag) THEN

    npro=npro+1
    namepro(npro)='v_t_o_'//TRIM(cidz_s(j))//'_ex'
    unitpro(npro)='m/s'
    descpro(npro)='Experimental Vtor on the outside midplane for ' &
                  //TRIM(cidz_s(j))
    valpro(1:nr_r,npro)=vt_im_ex_r(1:nr_r)

   npro=npro+1
    namepro(npro)='v_p_o_'//TRIM(cidz_s(j))//'_ex'
    unitpro(npro)='m/s'
    descpro(npro)='Experimental Vtor on the outside midplane for ' &
                  //TRIM(cidz_s(j))
    valpro(1:nr_r,npro)=vp_im_ex_r(1:nr_r)

   npro=npro+1
    namepro(npro)='k_p_o_'//TRIM(cidz_s(j))//'_ex'
    unitpro(npro)='m/s'
    descpro(npro)='k flux surf ave for ' &
                  //TRIM(cidz_s(j))
    valpro(1:nr_r,npro)=kk_im_ex_r(1:nr_r)

  ENDIF

  npro=npro+1
  namepro(npro)='M_t_o_'//TRIM(cidz_s(j))
  unitpro(npro)='-'
  descpro(npro)='Tor Mach no on the outside midplane for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=xm_t_o_rs(1:nr_r,j)

  IF(j == js_diag .AND. &
     l_diag) THEN

    npro=npro+1
    namepro(npro)='M_t_o_'//TRIM(cidz_s(j))//'_ex'
    unitpro(npro)='-'
    descpro(npro)='Experimental tor Mach no on the outside midplane for ' &
                  //TRIM(cidz_s(j))
    valpro(1:nr_r,npro)=xm_t_im_ex_o_r(1:nr_r)

  ENDIF

  npro=npro+1
  namepro(npro)='v_p_o_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Vpol on the outside midplane for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_p_o_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='M_p_o_'//TRIM(cidz_s(j))
  unitpro(npro)='-'
  descpro(npro)='Pol Mach no on the outside midplane for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=xm_p_o_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='v_par_o_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Vpar on the outside midplane for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_par_o_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='v_per_o_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s'
  descpro(npro)='Vperp on the outside midplane for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=v_per_o_rs(1:nr_r,j)

ENDDO !Over species

!-------------------------------------------------------------------------------
!Flux surface particle flow velocities
!-------------------------------------------------------------------------------
DO j1=1,m_s !Over species

  IF(l_reduce_out) THEN

    IF(j1 > m_i) CYCLE

    IF(j1 <= kscs) THEN

      j=j1

    ELSE

      j=j+izmax_i(j1)

    ENDIF

  ELSE

    j=j1

  ENDIF    

  npro=npro+1
  namepro(npro)='Omega_'//TRIM(cidz_s(j))
  unitpro(npro)='rad/s'
  descpro(npro)='Toroidal rotation frequency of '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=omega_rs(1:nr_r,j)

  IF(j == js_diag .AND. &
     l_diag) THEN

    npro=npro+1
    namepro(npro)='Omega_'//TRIM(cidz_s(j))//'_ex'
    unitpro(npro)='rad/s'
    descpro(npro)='Experimental toroidal rotation frequency of ' &
                  //TRIM(cidz_s(j))
    valpro(1:nr_r,npro)=omega_im_ex_r(1:nr_r)

  ENDIF

  npro=npro+1
  namepro(npro)='u_par_'//TRIM(cidz_s(j))
  unitpro(npro)='T*m/s'
  descpro(npro)='Parallel velocity parameter <v.B> for '//TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=u_par_rs(1:nr_r,j)

  npro=npro+1
  namepro(npro)='u_p_'//TRIM(cidz_s(j))
  unitpro(npro)='m/s/T'
  descpro(npro)='Poloidal velocity parameter <v.theta>/<B.theta> for ' &
                //TRIM(cidz_s(j))
  valpro(1:nr_r,npro)=u_p_rs(1:nr_r,j)

ENDDO !Over species

END SUBROUTINE FORCEBAL_OUT

SUBROUTINE FORCEBAL_SPECIES(cs, &
                            izs,amus,iflag,message)
!------------------------------------------------------------------------------
!FORCEBAL_SPECIES sets the mass number and atomic charge of a species
!  identified by its name
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!------------------------------------------------------------------------------
USE SPEC_KIND_MOD
IMPLICIT NONE

!Declaration of input variables
CHARACTER(len=*), INTENT(IN) :: &
  cs                     !species name [character]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag,               & !error and warning flag [-]
                         !=-1 warning
                         !=0 no warnings or errors
                         !=1 error
  izs                    !nuclear charge [-]

REAL(KIND=rspec), INTENT(OUT) :: &
  amus                   !atomic mass number [-]

!------------------------------------------------------------------------------
!Declaration of local variables
INTEGER, PARAMETER :: &
  mx_d=99

CHARACTER(len=3) :: &
  cd(mx_d)=(/'H  ','D  ','T  ','He3','He4', &
             'Li ','Be ','B  ','C  ','N  ', &
             'O  ','F  ','Ne ','Na ','Mg ', &
             'Al ','Si ','P  ','S  ','Cl ', &
             'A  ','K  ','Ca ','Sc ','Ti ', &
             'V  ','Cr ','Mn ','Fe ','Co ', &
             'Ni ','Cu ','Zn ','Ga ','Ge ', &
             'As ','Se ','Br ','Kr ','Rb ', &
             'Sr ','Y  ','Zr ','Nb ','Mo ', &
             'Tc ','Ru ','Rh ','Pd ','Ag ', &
             'Cd ','In ','Sn ','Sb ','Te ', &
             'I  ','Xe ','Cs ','Ba ','La ', &
             'Ce ','Pr ','Nd ','Pm ','Sm ', &
             'Eu ','Gd ','Tb ','Dy ','Ho ', &
             'Er ','Tm ','Yb ','Lu ','Hf ', &
             'Ta ','W  ','Re ','Os ','Ir ', &
             'Pt ','Au ','Hg ','Tl ','Pb ', &
             'Bi ','Po ','At ','Rn ','Fr ', &
             'Ra ','Ac ','Th ','Pa ','U  ', &
             'Np ','Pu ','Ar ','e  '/)

INTEGER :: &
  i,idone

INTEGER :: &
  izd(mx_d)=(/  1,  1,  1,  2,  2,  3,  4,  5,  6,  7, &
                8,  9, 10, 11, 12, 13, 14, 15, 16, 17, &
               18, 19, 20, 21, 22, 23, 24, 25, 26, 27, &
               28, 29, 30, 31, 32, 33, 34, 35, 36, 37, &
               38, 39, 40, 41, 42, 43, 44, 45, 46, 47, &
               48, 49, 50, 51, 52, 53, 54, 55, 56, 57, &
               58, 59, 60, 61, 62, 63, 64, 65, 66, 67, &
               68, 69, 70, 71, 72, 73, 74, 75, 76, 77, &
               78, 79, 80, 81, 82, 83, 84, 85, 86, 87, &
               88, 89, 90, 91, 92, 93, 94, 18,  1/)

REAL(KIND=rspec) :: &
  amud(mx_d)=(/ 1.000,    2.000,    3.000,    3.000,    4.000, &
                6.939,    9.012,   10.811,   12.011,   14.007, &
               15.999,   18.998,   20.183,   22.990,   24.312, &
               26.982,   28.086,   30.974,   32.064,   35.453, &
               39.948,   39.102,   40.080,   44.956,   47.900, &
               50.942,   51.996,   54.938,   55.847,   58.933, &
               58.710,   63.540,   65.370,   69.720,   72.590, &
               74.922,   78.960,   79.909,   83.800,   85.470, &
               87.620,   88.905,   91.220,   92.906,   95.940, &
               99.000,  101.070,  102.905,  106.400,  107.870, &
              112.400,  114.820,  118.690,  121.750,  127.600, &
              126.904,  131.300,  132.905,  137.340,  138.910, &
              140.120,  140.907,  144.240,  145.000,  150.350, &
              151.960,  157.250,  158.924,  162.500,  164.930, &
              167.260,  168.934,  173.040,  174.970,  178.490, &
              180.948,  183.850,  186.200,  190.200,  192.200, &
              195.090,  196.967,  200.590,  204.370,  207.190, &
              208.980,  210.000,  210.000,  222.000,  223.000, &
              226.000,  227.000,  232.038,  231.000,  238.030, &
              237.000,  242.000,   39.948,  5.44631e-4/)

REAL(KIND=rspec), PARAMETER :: &
  zero=0.0_rspec

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''
izs=0
amus=zero

!-------------------------------------------------------------------------------
!Identify species from name and get characteristics
!-------------------------------------------------------------------------------
!Check for null input
IF(LEN_TRIM(cs) < 1) THEN

  iflag=-1
  message='FORCEBAL_SPECIES/WARNING:null species name'
  GOTO 9999

ENDIF

i=1
idone=0

DO WHILE(idone == 0 .AND. &
         i <= mx_d)

  IF(TRIM(ADJUSTL(cs)) == TRIM(ADJUSTL(cd(i)))) THEN

    izs=izd(i)
    amus=amud(i)
    idone=1

  ENDIF

  i=i+1

ENDDO

!Check for invalid name
IF(idone == 0) THEN

  iflag=1
  message='FORCEBAL_SPECIES/ERROR:invalid species name'
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE FORCEBAL_SPECIES

SUBROUTINE READ_EFIT_EQDSK(nin,cnin,mxnx_xy,mxny_xy,mxn_lim, &
                           bt0,cur,psimag,psilim,r0,rmag,zmag, &
                           nx_xy,ny_xy,x_xy,y_xy,psi_xy, &
                           f_x,ffp_x,psi_x,q_x,rhop_x, &
                           n_lim,x_lim,y_lim, &
                           iflag,message)
!-------------------------------------------------------------------------------
!READ_EFIT_EQDSK reads an EQDSK file from EFIT
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  This routine reads and discards some of the stored EIFT data (e.g., the
!    boundary points that are not presently used with this application).
!  It also constructs the relevant 1-D and 2-D grids that are implicit in the
!    stored data.
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
IMPLICIT NONE

!Declaration of input variables
CHARACTER(len=*), INTENT(IN) :: &
  cnin                   !input file name [character]

INTEGER,INTENT(IN) :: &
  mxnx_xy,             & !maximum number of x points on psi(x,y) grid [-]
  mxny_xy,             & !maximum number of y points on psi(x,y) grid [-]
  mxn_lim,             & !maximum number of points on limiter surface [-]
  nin                    !input unit number [-]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag,               & !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error
  nx_xy,               & !number of x points on psi(x,y) grid [-]
  ny_xy,               & !number of y points on psi(x,y) grid [-]
  n_lim                  !number of points on limiter [-]

REAL(KIND=rspec), INTENT(OUT) :: &
  bt0,                 & !toroidal field at r0 [T]
  cur,                 & !toroidal plasma current [A]
  psimag,              & !poloidal flux/(2*pi) at axis [Wb/rad]
  psilim,              & !poloidal flux/(2*pi) at limiter/separatrix [Wb/rad]
  r0,                  & !reference major radius, center of limiter surface [m]
  rmag,                & !horizontal position of magnetic axis [m]
  zmag                   !vertical position of magnetic axis [m]

REAL(KIND=rspec), INTENT(OUT) :: &
  x_xy(mxnx_xy),             & !vertical grid for 2-D poloidal flux [m]
  y_xy(mxny_xy),             & !horizontal grid for 2-D poloidal flux [m]
  psi_xy(mxnx_xy,mxny_xy),   & !poloidal flux/(2*pi) on 2-D grid [Wb/rad]
  f_x(mxnx_xy),              & !F=R*B_t on equilibrium psi grid [m*T]
  ffp_x(mxnx_xy),            & !F*dF/dpsi on equilibrium psi grid [rad*T]
  psi_x(mxnx_xy),            & !poloidal flux/(2*pi) = equilibrium psi grid [Wb/rad]
  q_x(mxnx_xy),              & !safety factor on equilibrium psi grid [-]
  rhop_x(mxnx_xy),           & !normalized poloidal flux grid proportional to psi [-]
  x_lim(mxn_lim),            & !horizontal positions of limiter points [m]
  y_lim(mxn_lim)               !vertical positions of limiter points [m]

!-------------------------------------------------------------------------------
!Declaration of local variables
!Input from EQDSK file that is not retained
INTEGER :: &
  n_bdry                 !number of points on plasma boundary [-]


REAL(KIND=rspec) :: &
  rmin,                & !horizontal inside of computational domain [m]
  zmid,                & !vertical center of comoputational domain [m]
  rdim,                & !width of computational domain [m]
  zdim,                & !height of computational domain [m]
  x_bdry,              & !horizontal positions of boundary points [m]
  y_bdry,              & !vertical positions of boundary points [m]
  p_x(mxnx_xy),        & !plasma kinetic pressure [N/m**2]
  pp_x(mxnx_xy)          !dp/dpsi on equilibrium psi grid [rad*N/m**2/Wb]

!Other
INTEGER :: &
  i,j

REAL(KIND=rspec) :: &
  dum

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Open the EQDSK file
OPEN(UNIT=nin, &
     STATUS='old', & 
     FILE=cnin, &
     FORM='formatted')

!-------------------------------------------------------------------------------
!Read the EQDSK file
!-------------------------------------------------------------------------------
!Point data - dum values are duplicate information or not used
READ(nin,'(52x,2i4)') nx_xy,ny_xy

!Check if x dimension is exceeded
IF(nx_xy > mxnx_xy) THEN

  !Horizontal grid points exceed dimensions set by parameters
  iflag=1
  message='READ_EFIT_EQDSK(1)/ERROR:x grid dimension exceeded'
  GOTO 9999

ENDIF

!Check if y dimension is exceeded
IF(ny_xy > mxny_xy) THEN

  !Vertical grid points exceed dimensions set by parameters
  iflag=1
  message='READ_EFIT_EQDSK(2)/ERROR:y grid dimension exceeded'
  GOTO 9999

ENDIF

READ(nin,'(5e16.9)') rdim,zdim,r0,rmin,zmid
READ(nin,'(5e16.9)') rmag,zmag,psimag,psilim,bt0
READ(nin,'(5e16.9)') cur
READ(nin,'(5e16.9)') dum

!Read 1-D and 2-D data, radial grid is equally spaced in poloidal flux (1:nx_xy)
READ(nin,'(5e16.9)') (f_x(i),i=1,nx_xy)
READ(nin,'(5e16.9)') (p_x(i),i=1,nx_xy)
READ(nin,'(5e16.9)') (ffp_x(i),i=1,nx_xy)
READ(nin,'(5e16.9)') (pp_x(i),i=1,nx_xy)
READ(nin,'(5e16.9)') ((psi_xy(i,j),i=1,nx_xy),j=1,ny_xy)
READ(nin,'(5e16.9)') (q_x(i),i=1,nx_xy)

!Boundary and limiter data
READ(nin,'(2i5)') n_bdry,n_lim

!Check if boundary dimension is exceeded
!IF(n_bdry > mxn_bdry) THEN
!
!  !Boundary points exceed dimensions set by parameters
!  iflag=1
!  message='READ_EFIT_EQDSK(3)/ERROR:bdry grid dim exceeded'
!  GOTO 9999
!
!ENDIF

!Check if limiter dimension is exceeded
IF(n_lim > mxn_lim) THEN

  !Limiter points exceed dimensions set by parameters
  iflag=1
  message='READ_EFIT_EQDSK(4)/ERROR:lim grid dim exceeded'
  GOTO 9999

ENDIF

READ(nin,'(5e16.9)') (x_bdry,y_bdry,i=1,n_bdry)
READ(nin,'(5e16.9)') (x_lim(i),y_lim(i),i=1,n_lim)

!Construct implied grids
!2D grid
x_xy(1:nx_xy)=rmin+rdim*(/ (i-1,i=1,nx_xy) /)/(nx_xy-1)
y_xy(1:ny_xy)=zmid-zdim/2+zdim*(/ (i-1,i=1,ny_xy) /)/(ny_xy-1)

!1D radial grid and and poloidal flux
psi_x(1:nx_xy)=psimag+(psilim-psimag)*(/ (i-1,i=1,nx_xy) /)/(nx_xy-1)
rhop_x(1:nx_xy)=(psi_x(1:nx_xy)-psi_x(1))/(psi_x(nx_xy)-psi_x(1))

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

!Close the EQDSK file
CLOSE(unit=nin)

END SUBROUTINE READ_EFIT_EQDSK
