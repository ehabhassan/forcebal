MODULE READ_PRO_MOD
!-------------------------------------------------------------------------------
!READ_PRO_MOD is an F90 module to read profile data from files
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  Other modules exist to obtain profile data from storage systems for given
!    experiments at JET, PPPL, GA, ...
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE FORCEBAL_DATA_MOD
USE LINEAR1_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private data
!-------------------------------------------------------------------------------
!Constants
REAL(KIND=rspec), PRIVATE, PARAMETER :: &
  one=1.0_rspec,       & !REAL 1
  zero=0.0_rspec         !REAL 0

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  READ_REC               !reads a data file to determine the number of records
                         !  called from READ_PRO

!-------------------------------------------------------------------------------
! Procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE READ_PRO(idshot,time,c_den, &
                    iflag,message)
!-------------------------------------------------------------------------------
!READ_PRO reads ASCII data files then interpolates to the working grid
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------
! te: electron temperature
! ti: ion temperature
! ne: electron density
! nc: impurity density
! zf: Zeff
! om: impurity flux surface average toroidal rotation frequency
! vt: impurity toroidal rotation velocity on outside
! vp: impurity poloidal rotation velocity on outside
! kk: impurity flux surface average poloidal rotation  [m/s/T]
! eb: <E.B>
! jb: <J.B>
! nb: n_NBI
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  idshot                 !shot number/id [-]

REAL(KIND=rspec), INTENT(IN) :: &
  time,                & !analysis time [s]
  c_den                  !density cutoff below which species is ignored [/m**3]

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
CHARACTER :: &
  cidshot*9,cnin*25,ctimems*5

INTEGER :: &
  i,ierr,itimems,j,jpro,k,kscs,nin
     
INTEGER :: &
  nr_4d,nset1

REAL(KIND=rspec), ALLOCATABLE :: &
  r_4d(:),y_4d(:),denb_ex_r(:) 

REAL(KIND=rspec) :: &
  fscale,yim,yim2

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
nin=20
te_r(:)=zero
ti_r(:)=zero
zeff_ex_r(:)=zero
vt_im_ex_r(:)=zero
vp_im_ex_r(:)=zero
kk_im_ex_r(:)=zero
e_par_ex_r(:)=zero
xj_ex_r(:)=zero
omega_im_ex_r(:)=zero

ALLOCATE(denb_ex_r(1:nr_r))
  denb_ex_r(:)=zero

ALLOCATE(r_4d(1:100))   !LRB
  r_4d(:)=zero
ALLOCATE(y_4d(1:100))
  y_4d(:)=zero

!Species identifiers - set max index of scs species
kscs=1+n_scsion
IF(l_diag) kscs=kscs+1

!Species info for scs species - e + scs ions + diag ion
jif_s(1:kscs)=(/ (i,i=1,kscs) /)
jzf_s(1:kscs)=izmax_i(1:kscs)

!Species info for mcs species
IF(n_mcsion > 0) THEN

  k=kscs

  DO i=1,n_mcsion

    DO j=1,izmax_i(kscs+i)

      k=k+1
      jif_s(k)=kscs+i
      jzf_s(k)=j

    ENDDO

  ENDDO

ENDIF

!-------------------------------------------------------------------------------
!Set up file naming convention
!-------------------------------------------------------------------------------
itimems=1.0e3*(time+1.0e-5)
WRITE(ctimems,'(i5)') itimems

DO i=1,5

  IF(ctimems(i:i) == ' ') ctimems(i:i)='0'

ENDDO

WRITE(cidshot,'(i9)') idshot

DO i=1,9

  IF(cidshot(i:i) == ' ') cidshot(i:i)='0'

ENDDO

!-------------------------------------------------------------------------------
!Get profile data
!-------------------------------------------------------------------------------
!Set name and scale
DO jpro=1,12

  IF(jpro == 1) THEN

    !te: electron temperature
    cnin='te'//cidshot(4:9)//'.'//ctimems
    ierr=1
    fscale=one

  ELSEIF(jpro == 2) THEN

    !Ion temperature
    cnin='ti'//cidshot(4:9)//'.'//ctimems
    ierr=1
    fscale=one

  ELSEIF(jpro == 3) THEN

    !Electron density
    cnin='ne'//cidshot(4:9)//'.'//ctimems
    ierr=1
    fscale=1.0e19_rspec

  ELSEIF(jpro == 4) THEN

    !Zeff
    cnin='zf'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=one

  ELSEIF(jpro == 5) THEN

    !Impurity flux surface average toroidal rotation frequency
    cnin='om'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=one

  ELSEIF(jpro == 6) THEN

    !Impurity density
    cnin='nc'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=1.0e19_rspec

  ELSEIF(jpro == 7) THEN

    !Impurity toroidal rotation velocity on outside
    cnin='vt'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=one

  ELSEIF(jpro == 8) THEN

    !Impurity poloidal rotation velocity on outside
    cnin='vp'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=one

  ELSEIF(jpro == 9) THEN

    !Impurity flux surface average poloidal rotation  [m/s/T]
    cnin='kk'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=one

  ELSEIF(jpro == 10) THEN

    !<E.B>
    cnin='eb'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=one

  ELSEIF(jpro == 11) THEN

    !<J.B>
    cnin='jb'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=one

  ELSEIF(jpro == 12) THEN

    !n_NBI
    cnin='nb'//cidshot(4:9)//'.'//ctimems
    ierr=0
    fscale=1.0e19_rspec

  ENDIF

  !Read number of records in data
  CALL READ_REC(nin,cnin,nr_4d)

  !Check for errors
  IF(nr_4d == 0 .AND. &
     ierr == 1) THEN

    !Critical data missing
    iflag=1
    message='READ_PRO_DIIID(1)/ERROR:missing file '//cnin
    GOTO 9999

  ELSEIF(nr_4d == 0 .AND. &
         ierr == 0) THEN

    !Non-critical data missing
    iflag=0

  ELSE

    !Open data file
    OPEN(UNIT=nin, &
    STATUS='old', &
    FILE=cnin, &
    FORM='formatted')

    !Allocate arrays for reading data
    nset1=SIZE(r_4d)

    IF(ALLOCATED(r_4d)) THEN

      IF(nr_4d > nset1) THEN

        DEALLOCATE(r_4d,y_4d)
        ALLOCATE(r_4d(nr_4d), &
                 y_4d(nr_4d))

      ENDIF

    ELSE

      ALLOCATE(r_4d(1:nr_4d), &
               y_4d(1:nr_4d))

    ENDIF

    !Read data, scale it, and close file
    DO i=1,nr_4d

      READ(nin,*) r_4d(i),y_4d(i)

    ENDDO

    y_4d(:)=fscale*y_4d(:)

    CLOSE(unit=nin)

    !Interpolate from data grid to computational grid
    IF(jpro == 1) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          te_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(te_r)/'//message
        GOTO 9999

      ENDIF

    ELSEIF(jpro == 2) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          ti_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(ti_r)/'//message
        GOTO 9999

      ENDIF

    ELSEIF(jpro == 3) THEN

      !Electron density
      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          den_riz(:,1,1),iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(electron density)/'//message
        GOTO 9999

      ENDIF

      !Insure that density is above cutoff
      DO i=1,nr_r !Over radial nodes

        IF(den_riz(i,1,1) <= c_den) den_riz(i,1,1)=1.01_rspec*c_den

      ENDDO !Over radial nodes


    ELSEIF(jpro == 4) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,  &
                          zeff_ex_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(zeff_ex_r)/'//message
        GOTO 9999

      ENDIF

    ELSEIF(jpro == 5) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &   
                          omega_im_ex_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(omega_ex_r)/'//message
        GOTO 9999

      ENDIF

    ELSEIF(jpro == 6 .AND. &
           l_diag) THEN

      !Diagnostic impurity density
      j=jif_s(js_diag)
      k=jzf_s(js_diag)
      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          den_riz(:,j,k),iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(diag impurity density)/'//message
        GOTO 9999

      ENDIF

      !Insure that density is above cutoff
      DO i=1,nr_r !Over radial nodes

        IF(den_riz(i,j,k) <= c_den) den_riz(i,j,k)=1.01*c_den

      ENDDO !Over radial nodes

    ELSEIF(jpro == 7) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          vt_im_ex_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(vt_im_ex_r)/'//message
        GOTO 9999

      ENDIF
    ELSEIF(jpro == 8) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          vp_im_ex_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(vp_im_ex_r)/'//message
        GOTO 9999

      ENDIF
    ELSEIF(jpro == 9) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          kk_im_ex_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(kk_im_ex_r)/'//message
        GOTO 9999

      ENDIF

    ELSEIF(jpro == 10) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          e_par_ex_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(e_par_ex_r)/'//message
        GOTO 9999

      ENDIF

    ELSEIF(jpro == 11) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          xj_ex_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(xj_ex_r)/'//message
        GOTO 9999

      ENDIF

    ELSEIF(jpro == 12) THEN

      CALL LINEAR1_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r, &
                          denb_ex_r,iflag,message)

      IF(iflag == 1) THEN

        message='READ_PRO_DIIID(denb_ex_r)/'//message
        GOTO 9999

      ENDIF

    ENDIF

  ENDIF

ENDDO

!-------------------------------------------------------------------------------
!For other 'test' species, use a small percentage of electron density
!-------------------------------------------------------------------------------
!Single charge state ions
IF(n_scsion > 1) THEN

  DO k=3,1+n_scsion

    den_riz(1:nr_r,k,izmax_i(k))=1.0e-5*den_riz(1:nr_r,1,1)

  ENDDO

ENDIF

!Multiple charge state ions
IF(n_mcsion > 0) THEN

  DO k=kscs+1,m_s

    den_riz(1:nr_r,jif_s(k),jzf_s(k))=1.0e-5*den_riz(1:nr_r,1,1)

  ENDDO

ENDIF
    
!-------------------------------------------------------------------------------
!Set main ion density using charge neutrality
!-------------------------------------------------------------------------------
DO i=1,nr_r !Over radial nodes

  !Start with beam ion contributions
  yim=denb_ex_r(i)
  yim2=denb_ex_r(i)

  !Contributions to charge and Zeff from other than e and main ions
  IF(m_i > 2) THEN

    DO k=3,m_s

      yim=yim+jzf_s(k)*den_riz(i,jif_s(k),jzf_s(k))
      yim2=yim2+jzf_s(k)**2*den_riz(i,jif_s(k),jzf_s(k))

    ENDDO !Over remaining scs ions

  ENDIF

  den_riz(i,2,izmax_i(2))=(den_riz(i,1,1)-yim)/izmax_i(2)
  zeff_r(i)=(izmax_i(2)**2*den_riz(i,2,izmax_i(2))+yim2)/den_riz(i,1,1)

ENDDO !Over radial nodes

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

IF(ALLOCATED(r_4d)) DEALLOCATE(r_4d,y_4d)
DEALLOCATE(denb_ex_r)

END SUBROUTINE READ_PRO

SUBROUTINE READ_REC(nin,cnin,nrec)
!-------------------------------------------------------------------------------
!READ_REC reads a data file to determine the number of records
!
!References: 
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
CHARACTER(len=*), INTENT(IN) :: &
  cnin                   !input file name [-]

INTEGER, INTENT(IN) :: &
  nin                    !input file unit number [-]

!Declaration of output variables
INTEGER, INTENT(OUT) :: &
  nrec                   !number of records [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  iflag

REAL(KIND=rspec) :: &
  x

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
nrec=0
iflag=0

!-------------------------------------------------------------------------------
!Open and read number of records in file
!-------------------------------------------------------------------------------
OPEN(UNIT=nin, &
     STATUS='old', &
     FILE=cnin, &
     FORM='formatted', &
     IOSTAT=iflag)

IF(iflag /= 0) GOTO 9999

DO WHILE(iflag == 0)

  READ(nin,*,IOSTAT=iflag) x
  nrec=nrec+1

ENDDO

nrec=nrec-1

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

CLOSE(unit=nin)

END SUBROUTINE READ_REC

END MODULE READ_PRO_MOD

