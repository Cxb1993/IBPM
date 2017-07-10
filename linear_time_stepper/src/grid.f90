MODULE grid

  ! c Kunihiko (Sam) Taira and Tim Colonius
  ! December 2005
  ! calling program for ibfs v1.0
  ! ibfs v2.1 modified by Hsieh-Chen Tsai
  ! January 2013
  ! - now solves flows with arbitrary body motion in body-fixed frame
  !
  !
  ! ibfs v.2.2 modified by Thibault Flinois
  ! February 2013
  ! Now includes forward linearized simulation, based on Won Tae Joe's implementation
  ! Now includes adjoint equations optimization, also based on Won Tae Joe's implementation
  !


  USE parameters
  IMPLICIT NONE

  ! Variables for immersed boundary geometry and its motion (if any)
  REAL(KIND(0.0D0)) :: support = 6.0D0 ! support for smearing delta functions
  REAL(KIND(0.D0)) :: delta ! near field grid spacing

  ! coordinates and velocity on body
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: xb,vb, fb , ab, f_rdst
  INTEGER, DIMENSION(:), ALLOCATABLE :: codeb
  
  ! arrays for smearing coefficients
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: smear
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ismear
  INTEGER :: nsmear
  
  ! Numbers of cells, edges, etc.
  INTEGER :: nq, nb, nf, ns

  ! Integer pointer for field variables
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: f ! f(i,k) gives kth-comp. of immersed body force (or position or vel) at ith point on body
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: u ! u(i,j) gives point in the q array where u at face of cell i,j lives
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: v ! v(i,j) gives point in the q array where v at face of cell i,j lives

  ! For adjoint optimization
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: fb_base_old, fb_base, fb_base_new, rhs_adj, f_rdst_base

real(kind(0.d0)), dimension(:), allocatable :: vb0, ab0, xbp
integer :: nel
real(kind(0.d0)), dimension(:), allocatable :: u_ib, ud_ib, udd_ib
real(kind(0.d0)), dimension(:), allocatable :: u_ib_base, ud_ib_base, udd_ib_base
real(kind(0.d0)), dimension(:), allocatable :: xb0_base, xb_base
real(kind(0.d0)), dimension(:,:), allocatable :: sol_mat, K_sol, M_sol
real(kind(0.d0)), dimension(:), allocatable :: h_el, h0_el

real(kind(0.d0)), dimension(:,:), allocatable :: JHmat, JWmat, JEmat

  ! for bcs
  INTEGER :: top,bottom,left,right
  INTEGER :: top_phi,bottom_phi,left_phi,right_phi

  ! a special type for bodies
  TYPE body
     LOGICAL :: moving
     INTEGER :: npts
     INTEGER :: ipos  ! position in overall array
     REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: x,y
  END TYPE body
  TYPE actuator
     LOGICAL :: slaved
     INTEGER :: slavebody, slavept
     INTEGER :: npts
     INTEGER :: ipos ! position in overall array
     REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: x,y,s ! s is the strength of actuator for multipoint actuators
  END TYPE actuator

  INTEGER, PARAMETER :: maxbodies = 999 ! a large number
  TYPE(body), DIMENSION(maxbodies) :: bdy
  TYPE(actuator), DIMENSION(maxbodies) :: act

CONTAINS

  !======================================================================================!

  SUBROUTINE setup

    INTEGER :: i,j,next
 
    ! firt order of business is to setup the grid
    delta = len/REAL(m) 
 
    ! a few things to setup for multigrid laplace solution
    ! check that grid is divisible by 4
    IF ((MOD(m,4).ne.0).or.(MOD(n,4).ne.0)) THEN
       STOP 'grid must be divisible by 4'
    END IF
    ! for multigrid boundary conditions
    left = 0
    right = n+1
    bottom = 2*(n+1)
    top = 2*(n+1) + m+1
    
    left_phi = 0
    right_phi = n
    bottom_phi = 2*n
    top_phi = 2*n+m
    
    ! indexing for streamfunction
    ns = (m-1)*(n-1)

 
    ! indexing for velocities (flux components)
    nq = (m+1)*n + (n+1)*m
    ALLOCATE( u(1:m+1,1:n),  v(1:m,1:n+1) )
    next = 0
    DO j=1,n    
       DO i=1,m+1
          next = next+1
          u(i,j) = next
       END DO
    END DO
    DO j=1,n+1
       DO i=1,m
          next = next+1
          v(i,j) = next
       END DO
    END DO
    IF (next.ne.nq) STOP "error in setup_parms - u"

    ! now set up a body
    CALL setup_geometry
   


  END SUBROUTINE setup

!======================================================================================!

SUBROUTINE setup_reg(xx)

INTEGER :: ijsmear
REAL(KIND(0.D0)), DIMENSION(:) :: xx
REAL(KIND(0.D0)) :: x,y,d2, sup
real(kind(0.d0)) :: dist1, dist2
INTEGER :: i,j,k

! check to see if ismear and smear are already allocated
IF (ALLOCATED(ismear)) THEN
DEALLOCATE( ismear, smear )
END IF

d2 = 0.5d0*delta
sup = support*d2

! first count up how many boundary points smear over a given interior point
Nsmear = 0
DO i=1,m+1
    DO j=1,n
        x = delta*REAL(i-1)-offsetx
        y = delta*REAL(j-1)-offsety+d2
        ijsmear = 0
            DO k=1,Nb
                ! x-directional vector
                dist1 = ABS(x-xx(f(k,1)))
                dist2 = ABS(y-xx(f(k,2)))
                IF ( (dist1 < sup).AND.( (dist2) < sup)) THEN
                    ijsmear = ijsmear + 1
                END IF
            END DO
            Nsmear = MAX(Nsmear,ijsmear)
    END DO
END DO

DO i=1,m
    DO j=1,n+1
        x = delta*REAL(i-1)-offsetx+d2
        y = delta*REAL(j-1)-offsety
        ijsmear = 0
        DO k=1,Nb
            ! y-directional vector
            dist1 = ABS(x-xx(f(k,1)))
            dist2 = ABS(y-xx(f(k,2)))
            IF ((dist1 < sup).AND.( dist2 < sup) ) THEN
                ijsmear = ijsmear + 1
            END IF
        END DO
        Nsmear = MAX(Nsmear,ijsmear)
    END DO
END DO

ALLOCATE( ismear(Nq,Nsmear), smear(Nq,Nsmear) )

! now get the points over which the force is smeared
! default is to set the coefficient to 0 and the influence point to 1.  This avoids having to have
! any logic in the routines reg and regT
ismear = 1
smear = 0.d0

DO i=1,m+1
    DO j=1,n
        x = delta*REAL(i-1)-offsetx
        y = delta*REAL(j-1)-offsety+d2
        ijsmear = 0
        DO k=1,Nb
            ! x-directional vector
            dist1 = ABS(x-xx(f(k,1)))
            dist2 = ABS(y-xx(f(k,2)))
            IF ((dist1 < sup).AND.(dist2 < sup)) THEN
                ijsmear = ijsmear + 1
                ismear(u(i,j),ijsmear) = f(k,1)
                smear(u(i,j),ijsmear) =  (delta*delta)* deltafnc(dist1, delta) &
                * deltafnc(dist2, delta)
            END IF
        END DO
    END DO
END DO
DO i=1,m
    DO j=1,n+1
        x = delta*REAL(i-1)-offsetx + d2
        y = delta*REAL(j-1)-offsety
        ijsmear = 0
        DO k=1,Nb
            ! y-directional vector
            dist1 = ABS(x-xx(f(k,1)))
            dist2 = ABS(y-xx(f(k,2)))
            IF ((dist1< sup).AND.(dist2< sup)) THEN
                ijsmear = ijsmear+1
                ismear(v(i,j),ijsmear) = f(k,2)
                smear(v(i,j),ijsmear) = (delta*delta)* deltafnc(dist1, delta) &
                * deltafnc(dist2, delta)
            END IF
        END DO
    END DO
END DO


END SUBROUTINE setup_reg

!-----------------------------------------------------------------------------

FUNCTION deltafnc( r, dr )

REAL(KIND(0.D0)) :: r,dr,deltafnc
real(kind(0.d0)) :: r1, r2, r3, r4, a1, a5, a6, a7, a8, a9

r1 = r/dr
r2 = r1*r1
r3 = r2*r1
r4 = r3*r1

if (r1 .le. 1.d0) then
a5 = asin((1.d0/2.d0)*sqrt(3.d0)*(2.d0*r1-1.d0))
a8 = sqrt(1.d0-12.d0*r2+12.d0*r1)

deltafnc = 0.4166666667d-1*r4+(-.1388888889+0.3472222222d-1*a8)*r3+ &
(-0.7121664902d-1-0.5208333333d-1*a8+0.2405626122*a5)*r2+&
(-.2405626122*a5-.3792313933+.1012731481*a8)*r1+0.8018753741d-1*a5 &
-0.4195601852d-1*a8+.6485698427

elseif (r1 .le. 2.d0) then
a6 = asin((1.d0/2.d0)*sqrt(3.d0)*(-3.d0+2.d0*r1))
a9 = sqrt(-23.d0+36.d0*r1-12.d0*r2)

deltafnc = -0.6250000000d-1*r4+(.4861111111-0.1736111111d-1*a9)*r3 + &
(-1.143175026+0.7812500000d-1*a9-.1202813061*a6)*r2 + &
(.8751991178+.3608439183*a6-.1548032407*a9)*r1-.2806563809*a6 + &
0.822848104d-2+.1150173611*a9
elseif (r1 .le. 3.d0 ) then
a1 = asin((1.d0/2.d0*(2.d0*r1-5.d0))*sqrt(3.d0))
a7 = sqrt(-71.d0-12.d0*r2+60.d0*r1)

deltafnc = 0.2083333333d-1*r4+(0.3472222222d-2*a7-.2638888889)*r3+ &
(1.214391675-0.2604166667d-1*a7+0.2405626122d-1*a1)*r2+ &
(-.1202813061*a1-2.449273192+0.7262731481d-1*a7)*r1 +.1523563211*a1 &
+1.843201677-0.7306134259d-1*a7
!print *, deltafnc

else
deltafnc = 0.d0
end if

deltafnc = deltafnc / dr

END FUNCTION deltafnc

!----------------------------------------------------------------------------

  SUBROUTINE setup_geometry

    LOGICAL :: readinput
    INTEGER :: i,i_bdy,i_act, next
    CHARACTER(3) :: file_num

    ! look for bodies in input directory
    readinput = .TRUE.
    n_body = 0
    DO WHILE (readinput)
       WRITE(file_num,"(I3.3)") n_body+1
       INQUIRE(file="input/body."//file_num//".inp",exist=readinput)
       IF (readinput) THEN
          n_body=n_body+1
          OPEN(unit=8,file="input/body."//file_num//".inp",form='formatted',status='old')
          READ(8,*) bdy(n_body)%npts
          READ(8,*) bdy(n_body)%moving
          ALLOCATE( bdy(n_body)%x(bdy(n_body)%npts), bdy(n_body)%y(bdy(n_body)%npts) )
          DO i=1,bdy(n_body)%npts
             READ(8,*) bdy(n_body)%x(i), bdy(n_body)%y(i)
          END DO
          CLOSE(8)
       END IF
    END DO

    ! look for actuators in input directory
    readinput = .TRUE.
    n_actuator = 0
    DO WHILE (readinput)
       WRITE(file_num,"(I3.3)") n_actuator+1
       INQUIRE(file="input/actuator."//file_num//".inp",exist=readinput)
       IF (readinput) THEN
          n_actuator=n_actuator+1
          OPEN(unit=8,file="input/actuator."//file_num//".inp",form='formatted',status='old')
          READ(8,*) act(n_actuator)%npts
          READ(8,*) act(n_actuator)%slaved
          IF (act(n_actuator)%slaved) THEN
             READ(8,*) act(n_actuator)%slavebody, act(n_actuator)%slavept
          END IF

          ALLOCATE( act(n_actuator)%x(act(n_actuator)%npts), act(n_actuator)%y(act(n_actuator)%npts) )
          ALLOCATE( act(n_actuator)%s(act(n_actuator)%npts) )
          DO i=1,act(n_actuator)%npts
             READ(8,*) act(n_actuator)%x(i), act(n_actuator)%y(i)
          END DO
          CLOSE(8)
       END IF 
    END DO
    IF(ADJOINT)THEN !check nact=n_actuator
      IF(n_actuator.ne.nact) STOP 'nact from ib.inp is not equal to number of actuator files'
    ENDIF
    
    write(*,*) 'read all bodies and actuators'
    ! accumulate all bodies and actuators into global vector xb
    nb = 0
    DO i=1,n_body
       write(*,*) 'body no.',i,'has',bdy(i)%npts,'points.  Moving?',bdy(i)%moving
       nb = nb + bdy(i)%npts
    END DO
    DO i=1,n_actuator
       write(*,*) 'act. no.',i,'has',act(i)%npts,'points.  Slaved?',act(i)%slaved
       nb = nb + act(i)%npts
    END DO
    write(*,*) 'there are',nb,'lagrangian points'

    IF (nb==0) THEN
       STOP 'you must supply at least one body or actuator'
    END IF
   ! indexing for immersed body forces, positions, and velocities
    nf  = 2*nb                     ! number of immersed body forces
    ALLOCATE( f(1:nb,1:2)  )      
    next = 0
    DO i=1,nb
       next = next + 1
       f(i,1) = next
    END DO
    DO i = 1,nb
       next = next + 1
       f(i,2) = next
    END DO
    IF (next.ne.nf) STOP "error in setup_parms - f"

    ALLOCATE( xb(nf), vb(nf), ab(nf), fb(nf), codeb(nb), f_rdst(nf) ) ! for all simulations
    allocate( ab0(nf), vb0(nf), xbp(nf) )
    allocate( u_ib(3*nb), ud_ib(3*nb), udd_ib(3*nb) )
    allocate( Sol_mat(3*nb, 3*nb), K_sol(3*nb, 3*nb), M_sol(3*nb, 3*nb) )
    allocate( f_rdst_base(nf ), u_ib_base(3*nb), ud_ib_base(3*nb), udd_ib_base(3*nb) )
    allocate( xb0_base(nf), xb_base(nf) )
    allocate( h_el(nb-1), h0_el(nb-1))
    allocate( JHmat(nq, nf), JWmat(3*nb,nf), JEmat(nf, nf))

    Sol_mat = 0.d0
    K_sol = 0.d0
    M_sol = 0.d0

    ab0 = 0.d0
    vb0 = 0.d0
    xbp = 0.d0


    nel = nb - 1

    u_ib = 0.d0
    ud_ib = 0.d0
    udd_ib = 0.d0


    h_el = 0.d0
    h0_el = 0.d0

    JHmat = 0.d0
    JWmat = 0.d0
    JEmat = 0.d0


    ab = 0.d0
    vb = 0.d0
    fb = 0.d0

    IF(ADJOINT.OR.LINEAR) ALLOCATE(fb_base_old(nf),fb_base(nf),fb_base_new(nf)) !for adjoint and linear simulations only
    
    IF(ADJOINT) THEN !for adjoint simulations only
      ALLOCATE(rhs_adj(nf))
      fb_base_old = 0.d0
      fb_base = 0.d0
      fb_base_new = 0.d0
      rhs_adj = 0.d0
    ENDIF
    


    ! we initialize xb to the positions read from the files.  If this is a restart, this will be overwritten later
    
    CALL collect_bodies( xb,codeb)
    
    stationary = .TRUE.
    DO i=1,n_body
       IF ( bdy(i)%moving ) THEN
          stationary = .FALSE.
       END IF
    END DO

    write(*,*) 'setup global positions, velocities, and immersed body forces'
    CALL setup_reg(xb)
    write(*,*) 'setup regularization of initial geometry'
    
  END SUBROUTINE setup_geometry
  !======================================================================================!
  
    SUBROUTINE check_lin_files
    

    CHARACTER(7) :: charit
    LOGICAL anyactuators, foundallfiles
    INTEGER icheck

    ! Check for actuators as they may not be appropriate for linearized simulation
    anyactuators = .TRUE.
    INQUIRE(file="input/actuator.001.inp",exist=anyactuators)
    IF (anyactuators) THEN
      WRITE(*,*) '******* WARNING: Found actuator(s): Not recommended for linear simulation *******'
    END IF 
    
    
    foundallfiles=.TRUE.
    ! Check all required files exist in the output directory
    IF(STEADY)THEN
      WRITE(charit,"(I7.7)") istart
      INQUIRE(file="output/ib_lin"//charit//".var",exist=foundallfiles)
    ELSE
      icheck=istart
      DO WHILE (foundallfiles.AND.icheck.LE.istop+ibase)
        WRITE(charit,"(I7.7)") icheck
        INQUIRE(file="output/ib"//charit//".var",exist=foundallfiles)
        icheck=icheck+ibase
      END DO
    ENDIF

!    IF (foundallfiles) THEN
!      WRITE(*,*) 'Found all files required for linear simulation'
!    ELSE
!      STOP 'Could not find all files required for linear simulation'
!    ENDIF
    END SUBROUTINE check_lin_files
  
  !======================================================================================!
    SUBROUTINE check_adj_files
    

    CHARACTER(7) :: charit
    LOGICAL foundallfiles
    INTEGER icheck
    
    
    foundallfiles=.TRUE.
    ! Check all required restart files exist in the output directory
    
    icheck=istop
    DO WHILE (foundallfiles.AND.icheck.GE.istart+ibase)
      WRITE(charit,"(I7.7)") icheck
      INQUIRE(file="output/ib"//charit//".var",exist=foundallfiles)
      icheck=icheck-ibase
    END DO
    
    IF (foundallfiles) THEN
      WRITE(*,*) 'Found all files required for adjoint simulation'
    ELSE
      STOP 'Could not find all files required for adjoint simulation'
    ENDIF
    END SUBROUTINE check_adj_files
  
  !======================================================================================!

  SUBROUTINE collect_bodies( xx, code )

    REAL(KIND(0.D0)), DIMENSION(:) :: xx
    INTEGER, DIMENSION(:) :: code
    INTEGER :: i,i_bdy,i_act,next
    

    next = 1
    DO i_bdy=1,n_body
       bdy(i_bdy)%ipos = next
       DO i=1,bdy(i_bdy)%npts
          xx(f(next,1)) = bdy(i_bdy)%x(i)
          xx(f(next,2)) = bdy(i_bdy)%y(i)
          code(next) = i_bdy
          next = next+1
       END DO
    END DO
    DO i_act=1,n_actuator
       act(i_act)%ipos = next
       DO i=1,act(i_act)%npts
          xx(f(next,1)) = act(i_act)%x(i)
          xx(f(next,2)) = act(i_act)%y(i)
          code(next) = -i_act
          next = next+1
       END DO
    END DO


  END SUBROUTINE collect_bodies
  !======================================================================================!
  FUNCTION getbody( b, dir, bdyno, npts )
    
    INTEGER :: bdyno, dir, npts
    REAL(KIND(0.D0)), DIMENSION(:) :: b
    REAL(KIND(0.D0)), DIMENSION(npts) :: getbody

    getbody = b(f ( bdy(bdyno)%ipos:bdy(bdyno)%ipos+npts-1, dir ))
 
  END FUNCTION getbody
  !======================================================================================!
  FUNCTION getact( b, dir, actno, npts )
    
    INTEGER :: actno, dir, npts
    REAL(KIND(0.D0)), DIMENSION(:) :: b
    REAL(KIND(0.D0)), DIMENSION(npts) :: getact

    getact = b(f ( act(actno)%ipos:act(actno)%ipos+npts-1, dir ))
 
  END FUNCTION getact

END MODULE grid

