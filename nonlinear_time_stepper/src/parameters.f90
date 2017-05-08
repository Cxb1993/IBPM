MODULE parameters

  ! c Kunihiko (Sam) Taira and Tim Colonius
  ! December 2005
  ! calling program for ibfs v1.0

  ! Updated: March 2006 
  !   ibfs v1.2 - now includes actuation
  !
  ! ibfs v2.1 modified by Hsieh-Chen Tsai
  ! January 2013 
  ! - now solves flows with arbitrary body motion in body-fixed frame

  IMPLICIT NONE
  
  ! parameters
  INTEGER :: istart                ! initial time index
  INTEGER :: istop                ! last time index
  INTEGER :: isave                 ! save a restart every isave steps
  INTEGER :: m                    ! cells in x
  INTEGER :: n                   ! cells in y
  REAL(KIND(0.0D0)) :: dt     ! time step
  REAL(KIND(0.0D0)) :: Re      ! Reynolds number
  REAL(KIND(0.0D0)) :: cgtol   ! tol. for cg convergence (poission eq)
  INTEGER :: cg_max_iter           ! max. iterations for any cg iteration
  INTEGER :: n_actuator                 ! number of actuators 
  INTEGER :: n_body          ! number of moving bodies
                                        ! (will be counted in grid.f90) 
  REAL(KIND(0.D0)) :: len   ! length scale for grid
  REAL(KIND(0.D0)) :: offsetx   ! offset for grid in x
  REAL(KIND(0.D0)) :: offsety    ! offset for grid in y
  INTEGER :: mgridlev
  REAL(KIND(0.D0)) :: pi

  REAL(KIND(0.D0)) :: rot_angle  ! rotating angle of the grid
  REAL(KIND(0.D0)) :: rox  ! x-coord. of center of rotation rotation
  REAL(KIND(0.D0)) :: roy  ! y-coord. of center of rotation rotation

  LOGICAL :: stationary       ! all stationary bodies w.r.t the grid
  LOGICAL :: compute_pressure  ! whether to output pressure

  real(kind(0.d0)) :: R_rho
  real(kind(0.d0)) :: R_E
  real(kind(0.d0)) :: R_sh
  real(kind(0.d0)) :: R_th

  LOGICAL :: PINNED
  LOGICAL :: STANDARD

CONTAINS
  
  SUBROUTINE input

    LOGICAL :: readinput

    NAMELIST /read_parameters/ istart,istop,isave,m,n,dt,Re,cgtol, &
                                cg_max_iter,len,offsetx,offsety, &
                                mgridlev, compute_pressure, STANDARD, &
                                    PINNED, R_rho, &
                                R_E, R_sh, R_th

    pi  = 4.0d0*atan(1.0d0)
                                 
    ! read input
    INQUIRE(file='input/ib.inp',exist=readinput)

    IF (readinput) THEN
       OPEN(unit=3,file='input/ib.inp',form='formatted',status='old')
       READ(unit=3,nml=read_parameters)
       CLOSE(3)
      ELSE
       STOP 'cannot find input file'
    END IF

    write(*,*) 'read input file'

  END SUBROUTINE input

END MODULE parameters


