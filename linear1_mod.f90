MODULE LINEAR1_MOD
!-------------------------------------------------------------------------------
!LINEAR1-LINEAR interpolation in 1d
!
!LINEAR1_MOD is an F90 module of linear interpolating routines in 1d
!
!References:
!
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routines:
!
!  LINEAR1_INTERP      -interpolate from one grid to another
!  LINEAR1_INTEG       -integrate the linear fit
!
!Comments:
!
!  Linear interpolation routines are C0 (only f is continuous)
!
!  The modernization of the code structure into an F90 module takes advantage of
!    some of the more attractive features of F90:
!    -use of KIND for precision declarations
!    -optional arguments for I/O
!    -generic names for all intrinsic functions
!    -compilation using either free or fixed form
!    -no common blocks or other deprecated Fortran features
!    -dynamic and automatic alocation of variables
!    -array syntax for vector operations
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  LINEAR1_SEARCH         !find the indices that bracket an abscissa value
                         !  called from SPLINE1_EVAL

!-------------------------------------------------------------------------------
! Private data
!-------------------------------------------------------------------------------
!Constants
REAL(KIND=rspec), PRIVATE, PARAMETER :: &
  one=1.0_rspec,       & !REAL 1
  zero=0.0_rspec         !REAL 0

!-------------------------------------------------------------------------------
! Procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE LINEAR1_INTERP(n0,x0,y0,n1,x1, &
                          y1,iflag,message)
!-------------------------------------------------------------------------------
!W_LIN_INTERP performs a linear interpolation from the x0 mesh to the x1 mesh
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  If the target mesh is outside the domain of the input mesh, the end data
!    value is returned
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n0,                  & !number of source abscissas [-]
  n1                     !number of target abscissas [-]

REAL(KIND=rspec), INTENT(IN) :: &
  x0(:),               & !source abscissas (in increasing order) [arb]
  x1(:),               & !target abscissas [arb]
  y0(:)                  !source values [arb]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  y1(:)                  !target values [arb]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,il

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
IF(n1 <= 0) THEN

  iflag=-1
  message='LINEAR1_INTERP/WARNING(1):no points in output array'
  GOTO 9999

ENDIF

!Make sure there are at least two nodes in the input grid
IF(n0 < 2) THEN

  iflag=1
  message='LINEAR1_INTERP/ERROR(2):<2 points in source array'
  GOTO 9999

ENDIF

!Set starting index of source grid
il=1

!-------------------------------------------------------------------------------
!Interpolate from x0 to x1
!-------------------------------------------------------------------------------
DO i=1,n1 !Over index of target grid

  10 IF(x1(i) < x0(1)) THEN

    !Target is below data range, use innermost data value
    y1(i)=y0(1)
    iflag=-1
    message='LINEAR1_INTERP(3)/WARNING:x<x(1), use end point'

  ELSEIF(x1(i) == x0(1)) THEN

    !Target and source nodes coincide
    y1(i)=y0(1)

  ELSEIF(x1(i) > x0(il+1)) THEN

    !Beyond next source node
    !Step x0 grid forward and loop
    IF(il < n0-1) THEN

      il=il+1
      GOTO 10

    ELSE

      !Target is above data range, set to last value
      y1(i)=y0(n0)
      iflag=-1
      message='LINEAR1_INTERP(4)/WARNING:x>x(n0), use end point'

    ENDIF

  ELSE

    !Between the proper set of source nodes, interpolate
    y1(i)=y0(il)+(y0(il+1)-y0(il))*(x1(i)-x0(il))/(x0(il+1)-x0(il))

  ENDIF

ENDDO !Over index of target grid

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE
      
END SUBROUTINE LINEAR1_INTERP

SUBROUTINE LINEAR1_INTEG(k,n0,x0,f,n1,x1, &
                         value,iflag,message)
!-------------------------------------------------------------------------------
!LINEAR1_INTEG evaluates the integral f(x)*x**k, where f(x) is a linear function
!  and k >= 0, over the domain of x1
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!  W.A.Houlberg, generalized to arbitrary k, arbitrary domains 8/2006
!
!Comments:
!  The domain of x1 may be larger or smaller than the domain of x0
!  For values of x1 outside the domain of x0 the function value is assumed 0
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k,                   & !exponent of weighting factor (s(x)*x**k) [-]
                         !>=0
  n0,                  & !number of source abcissas [-]
  n1                     !number of output abcissas (may be 1) [-]

REAL(KIND=rspec), INTENT(IN) :: &
  f(:),                & !source ordinates [arb]
  x0(:),               & !source abcissas [arb]
  x1(:)                  !output abcissas [arb]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  value(:)               !integral of f(x)*x**k_order from x0(1) to x1(i) [arb]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,ilo,jlo,ido,jdo

REAL(KIND=rspec) :: &
  add,f2,sum,xnew,xold

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
value(1:n1)=zero

!If no target values are requested, return
IF(n1 <= 0) THEN

  iflag=-1
  message='LINEAR1_INTEG(1)/WARNING:no target points'
  GOTO 9999

ENDIF

!Find lowest index of source abscissas
CALL LINEAR1_SEARCH(x0,n0,x1(1),jlo)

IF(jlo == n0) THEN

  !Starting abscissa of target is above domain of source
  iflag=-1
  message='LINEAR1_INTEG(2)/WARNING:target above range'
  GO TO 9999

ELSEIF (jlo== 0) THEN

  !Find lowest valid index of target abscissas
  CALL LINEAR1_SEARCH(x1,n1,x0(1),ilo)

  IF(ilo == n1) THEN

    !End abscissa of target is below domain of source
    iflag=-1
    message='LINEAR1_INTEG(3)/WARNING:target below range'
    GO TO 9999

  ELSE

    !Start calculation at first node of source
    jlo=1
    xold=x0(1)
    !Integrate toward next target node
    ilo=ilo+1

  ENDIF

ELSE

  !Start calculation at first node of target and integrate toward next target
  xold=x1(1)
  ilo=2

ENDIF

f2=(f(jlo+1)-f(jlo))/(x0(jlo+1)-x0(jlo))
sum=zero

!-------------------------------------------------------------------------------
!Integrate over target nodes
!-------------------------------------------------------------------------------
LOOP_I: DO i=ilo,n1 !Over target nodes

  !Find dx to next source or target nodes
  ido=0

  DO WHILE(ido == 0) !Over source nodes

    IF(x1(i) < x0(jlo+1)) THEN

      !Hit target node
      xnew=x1(i)
      ido=1
      jdo=0

    ELSEIF(x1(i) == x0(jlo+1)) THEN

      !Hit both source and target nodes
      xnew=x1(i)
      ido=1
      jdo=1

    ELSEIF(x1(i) > x0(jlo+1)) THEN

      !Hit source node
      xnew=x0(jlo+1)
      ido=0
      jdo=1

    ENDIF

!Integrate over dx
    add= (xnew**(k+1)-xold**(k+1))/(k+1)*(f(jlo)-xold*f2) &
        +(xnew**(k+2)-xold**(k+2))/(k+2)*f2

!Add increment and update endpoint
    sum=sum+add
    xold=xnew

    IF(jdo == 1) THEN

      !Increment source node and derivative
      jlo=jlo+1

      IF(jlo == n0) THEN

        !Completed calculation
        value(i:n1)=sum
        EXIT LOOP_I

      ELSE

        f2=(f(jlo+1)-f(jlo))/(x0(jlo+1)-x0(jlo))

      ENDIF

    ENDIF

!Set integral value
    value(i)=sum

  ENDDO !Over source nodes

ENDDO LOOP_I !Over target nodes
   
!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE LINEAR1_INTEG

SUBROUTINE LINEAR1_SEARCH(x,n,xl, &
                          jlo)
!-------------------------------------------------------------------------------
!LINEAR1_SEARCH is a correlated table search routine to find the indices of the
!  array x that bound xl
!
!References:
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!
!Comments:
!  This is similar to the Numerical Recipes routine HUNT
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n                      !number of abscissas [-]

REAL(KIND=rspec), INTENT(IN) :: &
  xl,                  & !target value [arb]
  x(:)                   !monotonically increasing array of abscissas [arb]

!Declaration of input/output variables
INTEGER, INTENT(INOUT) :: &
  jlo                    !input starting lower index [-]
                         !<1     binary search
                         !=1,n-1 use value
                         !>n-1   binary search
                         !output starting lower index [-]
                         !=0     xl < x(1) 
                         !=1     x(1) <= xl <= x(2)
                         !=2,n-1 x(jlo) < xl <= x(jlo+1)
                         !=n     x(jlo) > x(n)

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  inc,jhi,jmid

!-------------------------------------------------------------------------------
!Check lower end of array, first two points
!-------------------------------------------------------------------------------
IF(xl < x(1)) THEN

  !Target is below node 1
  jlo=0

ELSEIF(xl <= x(2)) THEN

  !Target is between nodes 1 and 2 inclusive
  jlo=1
   
!-------------------------------------------------------------------------------
!Check middle range
!-------------------------------------------------------------------------------
ELSEIF(xl <= x(n)) THEN

  !Target is between nodes 2 and n
  IF(jlo < 1 .OR. &
     jlo > (n-1)) THEN

    !jlo from previous call is unusable
    jlo=2
    jhi=n

  !Bracket target value
  ELSE

    !Start with jlo from previous call
    inc=1

    IF(xl > x(jlo)) THEN

      !Search up
      jhi=jlo+1

      DO WHILE(xl > x(jhi))

        inc=2*inc
        jlo=jhi
        jhi=MIN(jlo+inc,n)

      ENDDO

    ELSE

      !Search down
      jhi=jlo
      jlo=jlo-1

      DO WHILE(xl <= x(jlo))

        inc=inc+inc
        jhi=jlo
        jlo=MAX(jlo-inc,1)

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

!-------------------------------------------------------------------------------
!Target is above node n
!-------------------------------------------------------------------------------
ELSE

  jlo=n

ENDIF

END SUBROUTINE LINEAR1_SEARCH

END MODULE LINEAR1_MOD
