MODULE WRITE_NETCDF_MOD
!-------------------------------------------------------------------------------
!WRITE_NETCDF-WRITEs output to a netCDF file
!
!WRITE_NETCDF_MOD is an F90 module of standarized output routines for netCDF
!
!References:
!
!  P.I.Strand, 5/2003
!  W.A.Houlberg, 7/2006
!
!Contains PUBLIC routines:
!
!  WRITE_NETCDF_OUT0 -writes 0D data
!  WRITE_NETCDF_OUT1 -writes 1D data
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE NETCDF
IMPLICIT NONE

!-------------------------------------------------------------------------------
! Procedures
!-------------------------------------------------------------------------------
CONTAINS 

SUBROUTINE WRITE_NETCDF_OUT0(cnetcdf,nv,pname,units,description,value, &
                             iflag,message)
!-------------------------------------------------------------------------------      
!WRITE_NETCDF_OUT0 writes scalar data to a netcdf file 'cnetcdf'
!
!References:
!  P.I.Strand, 5/2003
!  W.A.Houlberg, 7/2006
!
!Comments:
!  This subroutine requires netCDF version 3.5 including the f90 interface
!-------------------------------------------------------------------------------

!Declaration of input variables
CHARACTER(len=*), INTENT(IN) :: &
  cnetcdf,             & !existing output file name [character]
  pname(:),            & !plot names of variables [character]
  units(:),            & !units of variables [character]
  description(:)         !descriptions of variables [character]

INTEGER, INTENT(IN) :: &
  nv                     !no. of variables [-]

REAL(KIND=rspec), INTENT(IN) :: &
  value(:)               !data values (nx,nv) [see unitpro]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

CHARACTER(len=LEN(pname)) :: &
  name(nv)

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,id_var(nv),ierr,j,ncid

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
ierr=0
id_var(:)=0
ncid=0

!Replace invalid characters in plot name to get valid variable name
! EHAB MODIFIED LINES BELOW
!name(:)=pname(:)
name(:)=pname(1:nv)
! EHAB MODIFIED LINES ABOVE

DO i=1,nv

  DO j=1,LEN(name(i))

    IF(name(i)(j:j) == '(') name(i)(j:j)='_'
    IF(name(i)(j:j) == ')') name(i)(j:j)='_'
    IF(name(i)(j:j) == '+') name(i)(j:j)='p'
    IF(name(i)(j:j) == '-') name(i)(j:j)='m'
    IF(name(i)(j:j) == '/') name(i)(j:j)='O'
    IF(name(i)(j:j) == '*') name(i)(j:j)='T'

  ENDDO

ENDDO

!-------------------------------------------------------------------------------
!Open the existing NETCDF file
!-------------------------------------------------------------------------------
ierr=NF90_OPEN(TRIM(cnetcdf),NF90_WRITE,ncid)
IF(ierr /= nf90_noerr) GOTO 9999

ierr=NF90_REDEF(ncid)
IF(ierr /= nf90_noerr) GOTO 9999

!-------------------------------------------------------------------------------
!Define the variables
!-------------------------------------------------------------------------------

DO i=1,nv

    !Variable id 
    ierr=NF90_DEF_VAR(ncid,TRIM(name(i)),NF90_FLOAT,id_var(i))
    IF(ierr /= nf90_noerr) GOTO 9999

    !Plot name
    ierr=NF90_PUT_ATT(ncid,id_var(i),'plotname',TRIM(pname(i)))
    IF(ierr /= nf90_noerr) GOTO 9999

    !Description
    ierr=NF90_PUT_ATT(ncid,id_var(i),'description',TRIM(description(i)))
    IF(ierr /= nf90_noerr) GOTO 9999

    !Units
    ierr=NF90_PUT_ATT(ncid,id_var(i),'units',TRIM(units(i)))
    IF(ierr /= nf90_noerr) GOTO 9999
    
ENDDO

!-------------------------------------------------------------------------------
!Close the definition section of the file 
!-------------------------------------------------------------------------------
ierr=NF90_ENDDEF(ncid)
IF(ierr /= nf90_noerr) GOTO 9999

!-------------------------------------------------------------------------------
!Variable values
!-------------------------------------------------------------------------------
DO i=1,nv

  ierr=NF90_PUT_VAR(ncid,id_var(i),value(i))
  IF(ierr /= nf90_noerr) GOTO 9999

ENDDO  

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

IF(ierr/=nf90_noerr) THEN

  iflag=1
  message='WRITE_NETCDF_OUT0:ERROR:'//TRIM(nf90_strerror(ierr))

ENDIF

!Close netCDF file
ierr=NF90_CLOSE(ncid)

IF(ierr /= nf90_noerr) THEN

  iflag=1
  message='WRITE_NETCDF_OUT0:ERROR:'//TRIM(nf90_strerror(ierr))

ENDIF

END SUBROUTINE WRITE_NETCDF_OUT0

SUBROUTINE WRITE_NETCDF_OUT1(cnetcdf,ix,nx,nv,pname,units,description,value, &
                             iflag,message)
!-------------------------------------------------------------------------------      
!WRITE_NETCDF_OUT1 writes 1D data to a netcdf file 'cnetcdf'
!
!References:
!  P.I.Strand, 5/2003
!  W.A.Houlberg, 7/2006
!
!Comments:
!  This subroutine requires netCDF version 3.5 including the f90 interface
!-------------------------------------------------------------------------------

!Declaration of input variables
CHARACTER(len=*), INTENT(IN) :: &
  cnetcdf,             & !existing output file name [character]
  pname(:),            & !names of variables [character]
  units(:),            & !units of variables [character]
  description(:)         !descriptions of variables [character]

INTEGER, INTENT(IN) :: &
  ix,                  & !index of variable for coordinate [-]
  nx,                  & !no. of points for coordinate [-]
  nv                     !no. of variables [-]

REAL(KIND=rspec), INTENT(IN) :: &
  value(:,:)             !data values (nx,nv) [see unitpro]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

CHARACTER(len=LEN(pname)) :: &
  name(nv)

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,id_dim_x,id_var_x,id_var(nv),ierr,j,ncid

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
ierr=0
id_dim_x=0
id_var_x=0
id_var(:)=0
ncid=0

!Replace invalid characters in plot name to get valid variable name
! EHAB MODIFIED LINES BELOW
!name(:)=pname(:)
name(:)=pname(1:nv)
! EHAB MODIFIED LINES ABOVE

DO i=1,nv

  DO j=1,LEN(name(i))

    IF(name(i)(j:j) == '(') name(i)(j:j)='_'
    IF(name(i)(j:j) == ')') name(i)(j:j)='_'
    IF(name(i)(j:j) == '+') name(i)(j:j)='p'
    IF(name(i)(j:j) == '-') name(i)(j:j)='m'
    IF(name(i)(j:j) == '/') name(i)(j:j)='O'
    IF(name(i)(j:j) == '*') name(i)(j:j)='T'

  ENDDO

ENDDO

!-------------------------------------------------------------------------------
!Open the existing NETCDF file
!-------------------------------------------------------------------------------
ierr=NF90_OPEN(TRIM(cnetcdf),NF90_WRITE,ncid)
IF(ierr /= nf90_noerr) GOTO 9999

ierr=NF90_REDEF(ncid)
IF(ierr /= nf90_noerr) GOTO 9999

!-------------------------------------------------------------------------------      
!Define the grid
!-------------------------------------------------------------------------------
!Dimension id
ierr=NF90_DEF_DIM(ncid,TRIM(name(ix)),nx,id_dim_x)
IF(ierr /= nf90_noerr) GOTO 9999

!Variable id
ierr=NF90_DEF_VAR(ncid,TRIM(name(ix)),NF90_FLOAT,(/id_dim_x/),id_var_x)
IF(ierr /= nf90_noerr) GOTO 9999

!Description (long name)
ierr=NF90_PUT_ATT(ncid,id_var_x,'description',TRIM(description(ix)))
IF(ierr /= nf90_noerr) GOTO 9999
      
!Units
ierr=NF90_PUT_ATT(ncid,id_var_x,'units',TRIM(units(ix)))
IF(ierr /= nf90_noerr) GOTO 9999

!-------------------------------------------------------------------------------
!Define the variables
!-------------------------------------------------------------------------------
DO i=1,nv

  IF(i /= ix) THEN

    !Variable id 
    ierr=NF90_DEF_VAR(ncid,TRIM(name(i)),NF90_FLOAT,(/id_dim_x/),id_var(i))
    IF(ierr /= nf90_noerr) GOTO 9999

    !Plot name
    ierr=NF90_PUT_ATT(ncid,id_var(i),'plotname',TRIM(pname(i)))
    IF(ierr /= nf90_noerr) GOTO 9999

    !Description
    ierr=NF90_PUT_ATT(ncid,id_var(i),'description',TRIM(description(i)))
    IF(ierr /= nf90_noerr) GOTO 9999

    !Units
    ierr=NF90_PUT_ATT(ncid,id_var(i),'units',TRIM(units(i)))
    IF(ierr /= nf90_noerr) GOTO 9999
    
  END IF

ENDDO

!-------------------------------------------------------------------------------
!Close the definition section of the file 
!-------------------------------------------------------------------------------
ierr=NF90_ENDDEF(ncid)
IF(ierr /= nf90_noerr) GOTO 9999

!-------------------------------------------------------------------------------
!Grid values 
!-------------------------------------------------------------------------------
ierr=NF90_PUT_VAR(ncid,id_var_x,value(1:nx,ix))
IF(ierr /= nf90_noerr) GOTO 9999

!-------------------------------------------------------------------------------
!Variable values
!-------------------------------------------------------------------------------
DO i=1,nv

  IF(i /= ix) THEN

    ierr=NF90_PUT_VAR(ncid,id_var(i),value(1:nx,i))
    IF(ierr /= nf90_noerr) GOTO 9999

  ENDIF

ENDDO  

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

IF(ierr/=nf90_noerr) THEN

  iflag=1
  message='WRITE_NETCDF_OUT1:ERROR:'//TRIM(nf90_strerror(ierr))

ENDIF

!Close netCDF file
ierr=NF90_CLOSE(ncid)

IF(ierr /= nf90_noerr) THEN

  iflag=1
  message='WRITE_NETCDF_OUT1:ERROR:'//TRIM(nf90_strerror(ierr))

ENDIF

END SUBROUTINE WRITE_NETCDF_OUT1

END MODULE WRITE_NETCDF_MOD

      
