MODULE variables

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
  
  USE grid
  USE user
  IMPLICIT NONE

  ! in what follows, the last index refers to the grid level, first 1 (or 2) indices to point in space
  
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: omega, s, sbase     ! vorticity
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: rhs_old   ! nonlinear terms at previous time-step (for AB2 integration)
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: q, q0, q0p, q0r  
  
  ! necessary for the linearized simulations
REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: q_base_old,q_base,q_base_new

  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: q0_base_old,q0_base,q0_base_new, q0r_base_new, q0p_base_new
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: omega_base_old,omega_base,omega_base_new
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: rhs_old_base
  INTEGER :: counter_interp
  
  ! necessary for adjoint optimization
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: control,control_old,gradient
  REAL(KIND(0.D0)) :: tempcost
  LOGICAL :: FORWARD_SIM
  
  ! fluxes and additional potential flow and solid body rotation be added

  ! variables for Cholesky (stationary body)
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: cholmat
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: cholvec
real(kind(0.d0)), dimension(:), allocatable :: vb_base
CONTAINS
  
  !======================================================================================!

  SUBROUTINE setup_variables

    USE parameters
    USE myfft
    INTEGER :: k

! Allocate variables and set constants
    
    CALL setup_fft

    ALLOCATE( omega(2:m,2:n,mgridlev), rhs_old(2:m,2:n,mgridlev), s(2:m,2:n,mgridlev) , &
        sbase(2:m,2:n,mgridlev))
    ALLOCATE( q(Nq,mgridlev), q0(Nq,mgridlev), q0p(Nq,mgridlev), q0r(Nq,mgridlev) )
    
    IF(LINEAR.OR.ADJOINT) THEN

       ALLOCATE( omega_base_old(2:m,2:n,mgridlev), omega_base(2:m,2:n,mgridlev), omega_base_new(2:m,2:n,mgridlev) )
        ALLOCATE( rhs_old_base(2:m,2:n,mgridlev) )
       ALLOCATE( q_base_old(Nq,mgridlev),  q_base(Nq,mgridlev),  q_base_new(Nq,mgridlev) )
       ALLOCATE( q0_base_old(Nq,mgridlev),  q0_base(Nq,mgridlev),  q0_base_new(Nq,mgridlev) )
       ALLOCATE (q0r_base_new(Nq,mgridlev),q0p_base_new(Nq,mgridlev))
        allocate( vb_base( nf ) )
    ENDIF
    
    IF (stationary) THEN
       ALLOCATE( cholmat(Nf,Nf), cholvec(Nf) )
    END IF

! Initialize variables for the simulation  
 
    IF(FORWARD_SIM)THEN !forward simulations

      IF(istart.eq.0)THEN !forward simulations with istart=0
        CALL initial_condition
      ELSE ! forward simulations restarting from istart>0
        CALL read_variables(istart)
      ENDIF
      
      
    ELSE !backward adjoint simulations
      CALL initial_condition
    ENDIF
    CALL setup_reg(xb)
  END SUBROUTINE setup_variables
  
  !======================================================================================!
 
  SUBROUTINE setup_adj_variables

    USE parameters
    
    ALLOCATE(control(istart:istop,2*nbodyf+2*nact+3) )
    ALLOCATE(control_old(istart:istop,2*nbodyf+2*nact+3) )
    ALLOCATE(gradient(istart:istop,2*nbodyf+2*nact+3) )

    control = 0.d0
    control_old = 0.d0
    gradient = 0.d0  
    tempcost = 0.d0

  END SUBROUTINE setup_adj_variables
  
  !======================================================================================!
 
   SUBROUTINE destroy_variables

    USE parameters
    USE myfft

    CALL destroy_fft
    DEALLOCATE( omega, rhs_old, q, q0, s, q0p, q0r)
    IF (stationary) THEN
       DEALLOCATE( cholmat, cholvec )
    END IF
    IF(LINEAR.OR.ADJOINT) THEN

       DEALLOCATE( omega_base_old, omega_base, omega_base_new)
       DEALLOCATE( rhs_old_base)
       DEALLOCATE( q_base_old,q_base,q_base_new)
       DEALLOCATE( q0_base_old,q0_base,q0_base_new,q0r_base_new,q0p_base_new)
       
  
    ENDIF

  END SUBROUTINE destroy_variables
  
  !======================================================================================!

  SUBROUTINE destroy_adj_variables

    USE parameters
    USE myfft
       DEALLOCATE(control,control_old,gradient)
       
  END SUBROUTINE destroy_adj_variables

  !======================================================================================!
 
  SUBROUTINE destroy_grid

    USE parameters

    DEALLOCATE( f, u, v )
    DEALLOCATE( xb,vb,fb, codeb, smear, ismear )
    

  END SUBROUTINE destroy_grid

  !================================================================
 
  SUBROUTINE write_cholesky

    USE parameters

    OPEN(unit=100,file="output/ib.chd",form="unformatted",status="unknown")
    WRITE(100) cholmat, cholvec
    CLOSE(100)

  END SUBROUTINE write_cholesky
  
  !================================================================


!================================================================

SUBROUTINE write_bmat

USE parameters
use grid

OPEN(unit=100,file="output/bmat.chd",form="unformatted",status="unknown")
WRITE(100) M_sol
WRITE(100) K_sol
WRITE(100) sol_mat
WRITE(100) cholmat
write(100) h_el, h0_el
CLOSE(100)

END SUBROUTINE write_bmat

!================================================================

!================================================================

SUBROUTINE write_Jmat

USE parameters
use grid

OPEN(unit=100,file="output/Jmat.chd",form="unformatted",status="unknown")
WRITE(100) JHmat
WRITE(100) JWmat
WRITE(100) JEmat
CLOSE(100)

END SUBROUTINE write_Jmat

!================================================================


!================================================================

SUBROUTINE read_bmat

USE parameters
use grid

OPEN(unit=100,file="output/bmat.chd",form="unformatted",status="unknown")
read(100) M_sol
read(100) K_sol
read(100) sol_mat
read(100) cholmat
read(100) h_el, h0_el
CLOSE(100)

END SUBROUTINE read_bmat

!================================================================
!================================================================

SUBROUTINE read_Jmat

USE parameters
use grid

OPEN(unit=100,file="output/Jmat.chd",form="unformatted",status="unknown")
read(100) JHmat
read(100) JWmat
read(100) JEmat
CLOSE(100)

END SUBROUTINE read_Jmat

!================================================================




  SUBROUTINE read_cholesky

    USE parameters
    USE grid

    OPEN(unit=100,file="output/ib.chd",form="unformatted",status="unknown")
    
    READ(100) cholmat, cholvec
    CLOSE(100)
 
  END SUBROUTINE read_cholesky
  
 !======================================================================================!
  
  SUBROUTINE write_variables(it)

    USE parameters
    CHARACTER(7) :: charit
    INTEGER :: it

    write(*,*) 'writing variables at it=',it
    WRITE(charit,"(I7.7)") it
    IF(LINEAR)THEN
       OPEN(unit=100,file="output/ib_lin"//charit//".var",form="unformatted",status="unknown")
    ELSEIF(ADJOINT.and..not.FORWARD_SIM)THEN
       OPEN(unit=100,file="output/ib_adj"//charit//".var",form="unformatted",status="unknown")
    ELSE
      OPEN(unit=100,file="output/ib"//charit//".var",form="unformatted",status="unknown")
    ENDIF
    
    WRITE(100) m,n,mgridlev,nb
    WRITE(100) re,dt,len,offsetx,offsety
    WRITE(100) omega,xb,vb,fb,codeb,rhs_old, q, q0p, q0r
    WRITE(100) s
    WRITE(100) rot_angle, rox, roy
    write(100) u_ib, ud_ib, udd_ib
    write(100) f_rdst, xbp
    CLOSE(100)


  END SUBROUTINE write_variables
 !======================================================================================!

  SUBROUTINE write_control(ioptimize)

    INTEGER :: ioptimize, i,j
    CHARACTER(3) :: char_iopt

   WRITE(char_iopt,"(I3.3)") ioptimize
   OPEN(unit=203, file="output/control"//char_iopt//".dat", form='formatted',status='replace')
   
   DO i=istart+1,istop
     WRITE(203,*) i, (control(i,j),j=1,nact*2+nbodyf*2+3)
   END DO
   
   CLOSE(203)
   WRITE(*,*) 'Wrote Control to file at iopt = ',ioptimize
  END SUBROUTINE write_control
  
   !======================================================================================!

  SUBROUTINE write_gradient(ioptimize)

    INTEGER :: ioptimize, i,j
    CHARACTER(3) :: char_iopt

   WRITE(char_iopt,"(I3.3)") ioptimize
   OPEN(unit=203, file="output/gradient"//char_iopt//".dat", form='formatted',status='replace')
   
   DO i=istart+1,istop
     WRITE(203,*) i, (gradient(i,j),j=1,nact*2+nbodyf*2+3)
   END DO
   
   CLOSE(203)
   WRITE(*,*) 'Wrote Gradient to file at iopt=',ioptimize
  END SUBROUTINE write_gradient
 
   !======================================================================================!  
  
  SUBROUTINE read_control(ioptimize)
        
    INTEGER :: ioptimize, i,j
    REAL(KIND(0.d0)) :: temp
    CHARACTER(3) :: char_iopt
    LOGICAL :: readcontrol
    IF (ioptimize .eq. 1) THEN ! 1st optimization horizon: need to have an initial control guess
      INQUIRE(file="input/control001.dat",exist=readcontrol)
      
      IF(readcontrol)THEN ! Initial guess either from precomputed file...
        WRITE(*,*) 'Found pre-existing control file in input directory'
        OPEN(unit=206, file='input/control001.dat', form='formatted',status='old')
        DO i=(istart+1),(istop)
          READ(206,*) temp, (control(i,j),j=1,2*nbodyf+2*nact+3)
        END DO
        WRITE(*,*) 'Read Control file from input directory at iopt=',ioptimize
      ELSE ! ... OR default initial guess = constant value set in opt.inp
        DO j=1,2*nbodyf+2*nact+3
            control(:,j) = CINIT(j)
        ENDDO
           WRITE(*,*) 'Used default controls at iopt=',ioptimize
      END IF
    ELSE ! Not the first optimization iteration: MUST have an old file with controls

      WRITE(char_iopt,"(I3.3)") ioptimize
      OPEN(unit=206, file="output/control"//char_iopt//".dat", form='formatted',status='old')
        DO i=istart+1,istop
          READ(206,*) temp, (control(i,j),j=1,2*nbodyf+2*nact+3)
        END DO
        WRITE(*,*) 'Read old Control file at iopt=',ioptimize
    END IF
    CLOSE(206)
   END SUBROUTINE read_control
 
 !======================================================================================!
 
  SUBROUTINE read_variables(it)

    USE parameters
    CHARACTER(7) :: charit
    INTEGER :: it

    write(*,*) 'reading variables at it=',it
    WRITE(charit,"(I7.7)") it
    OPEN(unit=100,file="output/ib_lin"//charit//".var",form="unformatted",status="unknown")

    READ(100) m,n,mgridlev,nb
    READ(100) re,dt,len,offsetx,offsety
    READ(100) omega,xb,vb,fb,codeb,rhs_old, q, q0p, q0r
    READ(100) s
    READ(100) rot_angle, rox, roy
    read(100) u_ib, ud_ib, udd_ib
    read(100) f_rdst, xbp
    CLOSE(100)


  END SUBROUTINE read_variables
  
!======================================================================================!
     SUBROUTINE update_base(itime)

    !************************************************************************!
    !*    Update current base flow field state by reading from a            *!
    !*    new file or by interpolating between two sets of values           *!
    !*    base = current linearised value                                   *!
    !*    base_old = previous base flow field state relative to itime       *!
    !*    base_new = next base flow field state relative to itime           *!   
    !************************************************************************!

    
    INTEGER :: i,j,k, itime
    
    ! Check if nonlinear simulation base flow field state needs to be updated 
    
    IF(STEADY) THEN
      CALL read_base(0)
      omega_base = omega_base_new
      q_base  = q_base_new
      q0_base = q0_base_new
    ELSE
    
       IF (MOD(itime-istart, ibase).eq.0) THEN
         
         !For the first time step, need to load initial values
         IF(LINEAR.and.itime==istart) CALL read_base(itime)
         IF(ADJOINT.and..not.FORWARD_SIM.and.itime==istop) CALL read_base(itime)
         counter_interp = 0 ! reset interpolation counter to 0
       
         ! Update 'base' and 'base_old' variables
         omega_base = omega_base_new
         q_base  = q_base_new
         q0_base = q0_base_new
         fb_base = fb_base_new

         omega_base_old = omega_base_new
         q_base_old  = q_base_new
         q0_base_old = q0_base_new
         fb_base_old = fb_base_new
         
         ! Update 'base_new' variables
         IF(LINEAR) CALL read_base(itime+ibase)
         IF(ADJOINT.AND..NOT.FORWARD_SIM) THEN
           IF(itime.eq.istart+ibase.and.ibase.gt.1)THEN
             CALL read_base(istart+1)
             IBASE=IBASE-1
           ELSE
             CALL read_base(itime-ibase)
           ENDIF
         ENDIF
       ELSE    ! If base state does not need to be updated, interpolate variables
        
         counter_interp = counter_interp + 1   ! update interpolation counter

         ! Update interpolation of state and immersed body forces
         DO k=1,mgridlev
           DO j=2,n
             DO i=2,m
                omega_base(i,j,k) = omega_base_old(i,j,k) + &
                         REAL(counter_interp)/REAL(ibase)*(omega_base_new(i,j,k)-omega_base_old(i,j,k))
             END DO
           END DO

           DO i=1,Nq
             q_base(i,k) = q_base_old(i,k) + REAL(counter_interp)/REAL(ibase)*(q_base_new(i,k)-q_base_old(i,k))
             q0_base(i,k) = q0_base_old(i,k) + REAL(counter_interp)/REAL(ibase)*(q0_base_new(i,k)-q0_base_old(i,k))
           END DO
         END DO
       DO i=1,Nf
          fb_base(i) = fb_base_old(i) + REAL(counter_interp)/REAL(ibase)*(fb_base_new(i)-fb_base_old(i))
       END DO 
      END IF
    
    END IF
    
    END SUBROUTINE update_base
  
  
  !======================================================================================!
  
  SUBROUTINE read_base(it)
  
    USE parameters
    USE grid
    CHARACTER(7) :: charit
    INTEGER :: it, kk

    
    write(*,*) 'Reading state of base flow field at itime=',it
    WRITE(charit,"(I7.7)") it
    OPEN(unit=100,file="output/ib"//charit//".var",form="unformatted",status="unknown")
    READ(100) m,n,mgridlev,nb
    READ(100) re,dt,len,offsetx,offsety
    READ(100) omega_base_new, xb, vb_base, fb_base, codeb, rhs_old_base, q_base_new, q0p_base_new, q0r_base_new
    read(100) sbase
    read(100) rot_angle, rox, roy
    read(100) xb0_base, u_ib_base, ud_ib_base, udd_ib_base
    read(100) f_rdst_base
    CLOSE(100)
 
    q0_base_new=q0p_base_new+q0r_base_new


!position of beam in base configuration...
do kk = 1, nb

xb_base(kk) = xb0_base(kk) + u_ib_base(3*(kk-1) + 1)
xb_base(nb + kk) =  xb0_base(nb + kk) + u_ib_base(3*(kk-1) + 2)

end do

q_base = q_base_new


!In linear approximation, take beam position to be that of xb_base...
xb = xb_base

END SUBROUTINE read_base

  !======================================================================================!

  SUBROUTINE write_cost(ioptimize,totalcost)
  INTEGER :: ioptimize
  REAL(KIND(0.D0)) :: totalcost

   IF(ioptimize.EQ.1)THEN
     OPEN(unit=204, file="output/cost.dat", form='formatted',status='replace')
   ELSE
     OPEN(unit=204, file="output/cost.dat", form='formatted',position='append')
   ENDIF

   WRITE(204,*) ioptimize, totalcost
   CLOSE(204)
   WRITE(*,*) 'Wrote cost file'
  END SUBROUTINE write_cost
    
  !======================================================================================!
  
  SUBROUTINE initial_condition

    USE parameters
    USE user
    CHARACTER(7) :: charit
    REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: xold,vold,fold
    INTEGER, DIMENSION(:), ALLOCATABLE :: codeold
    INTEGER :: m_in,n_in,mgridlev_in,nb_in, i, j
    REAL(KIND(0.D0)) :: re_in,dt_in,len_in,offsetx_in,offsety_in ,xx,yy,gamma0,rc
    REAL(KIND(0.D0)), DIMENSION(5) :: uv
    LOGICAL :: readic
    
    INQUIRE(file="input/initial.var",exist=readic)
    
    IF (readic.and.FORWARD_SIM) THEN !found initial.var and forward simulation
       WRITE(*,*) 'Inputing initial condition files'
       OPEN(unit=100,file="input/initial.var",form="unformatted",status="unknown")
       READ(100) m_in,n_in,mgridlev_in,nb_in

       IF ((m_in.ne.m).or.(n_in.ne.n).or.(mgridlev.ne.mgridlev_in)) THEN
          WRITE(*,*) 'Initial condition file used different numbers of grid points'
          STOP 'Unable to continue'
       END IF
       READ(100) re_in,dt_in,len_in,offsetx_in,offsety_in
    
       IF (re_in.ne.re) THEN
          WRITE(*,*) 'Reynolds number has changed from initial condition'
          WRITE(*,*) 'Disregarding old value'
       END IF
       
       IF ((len_in.ne.len).or.(offsetx.ne.offsetx_in).or.(offsety.ne.offsety_in)) THEN
          WRITE(*,*) 'Physical dimensions of grid have changed'
          WRITE(*,*) 'Disregarding old values; proceed with caution!'
       END IF
       
       IF (nb.ne.nb_in) THEN
          WRITE(*,*) 'Geometry has changed from initial condition.'
          WRITE(*,*) 'Disregarding old geometry and forces'
          ALLOCATE( xold(2*nb_in), vold(2*nb_in), fold(2*nb_in), codeold(nb_in))
          READ(100) omega,xold,vold,fold,codeold,rhs_old, q, q0p, q0r
          CLOSE(100)
       ELSE
          ALLOCATE( xold(2*nb_in), vold(2*nb_in), fold(2*nb_in), codeold(nb_in))
          READ(100) omega,xold,vold,fold,codeold,rhs_old, q, q0p, q0r
          IF (ANY(xold.ne.xb).or.ANY(codeold.ne.codeb)) THEN
             WRITE(*,*) 'Geometry has changed from initial condition.'
             WRITE(*,*) 'Disregarding old geometry and forces'
          ELSE
             vb = vold
             fb = fold
          END IF
       END IF
 
       READ(100) s
       READ(100) rot_angle, rox, roy
       CLOSE(100)

       q0 = q0p + q0r

    ELSE    ! no initial.var or adjoint backward
      q = 0.d0
      rhs_old = 0.d0
      q0 = q0p + q0r
      rot_angle = 0.d0
      uv = motion_grid(0)
      rox = uv(4)
      roy = uv(5)
      omega = 0.d0






      IF(FORWARD_SIM) THEN ! no initial.var and forward

        !  ------- OSEEN VORTEX ---------------------   
        !         DO i=2,m
        !         DO j=2,n
        !           xx = (REAL(i)-1.d0)*delta-offsety
        !           yy = (REAL(j)-1.d0)*delta-offsetx
        !           rc=0.05d0
        !           gamma0=0.001d0
        !           omega(i,j,1)=gamma0/pi/rc*(exp(-((xx+0.75d0)**2.d0+(yy-0.5d0)**2.d0)/rc**2.d0))
        !         ENDDO
        !         ENDDO
        ! --------------------------------------------
        
        ! Grid movement
        q0p = motion_potential(0)
        q0r = motion_rotation(0)
        uv = motion_grid( 0 )
        
        IF(LINEAR)THEN ! no initial.var and forward and linear
          ! Initialize base flow variables
          omega_base_old = 0.d0
          omega_base_new = 0.d0
          omega_base = 0.d0
          q_base_old = 0.d0
          q_base_new = 0.d0
          q_base = 0.d0
          q0_base_old = 0.d0
          q0_base_new = 0.d0
          q0_base = 0.d0
          rhs_old_base = 0.d0 
        ENDIF       
        
      ELSE !backward adjoint
        
        ! No adjoint 'grid movement' possible
        q0p = 0.d0 
        q0r = 0.d0
        uv =  0.d0
        
        ! Initialize base flow variables
        omega_base_old = 0.d0
        omega_base_new = 0.d0
        omega_base = 0.d0
        q_base_old = 0.d0
        q_base_new = 0.d0
        q_base = 0.d0
        q0_base_old = 0.d0
        q0_base_new = 0.d0
        q0_base = 0.d0
        rhs_old_base = 0.d0
      ENDIF

    END IF

  END SUBROUTINE initial_condition
  !======================================================================================!
  
    FUNCTION motion_potential( it ) RESULT(qref)

    ! motion_potential = - u_b

    REAL(KIND(0.D0)), DIMENSION(Nq,mgridlev) :: qref
    INTEGER                                  :: it, k, i, j
    REAL(KIND(0.D0)) :: fac
    REAL(KIND(0.D0)), DIMENSION(5) :: uv

    uv = motion_grid(it)
    
    IF(ADJOINT.and.FORWARD_SIM)THEN
      DO k=1,mgridlev  ! one for each grid 
       fac                  = delta*2.d0**(k-1)  ! cell face length on each grid
       qref(1:(m+1)*n,k)    = -(uv(1)+control(it,2*nbodyf+2*nact+1))
       qref((m+1)*n+1:Nq,k) = -(uv(2)+control(it,2*nbodyf+2*nact+2))
       qref(:,k) = fac*qref(:,k)
      ENDDO
    ELSE
      DO k=1,mgridlev  ! one for each grid 
       fac                  = delta*2.d0**(k-1)  ! cell face length on each grid
       qref(1:(m+1)*n,k)    = -uv(1)
       qref((m+1)*n+1:Nq,k) = -uv(2)
       qref(:,k) = fac*qref(:,k)
      ENDDO
    ENDIF
    
  END FUNCTION motion_potential

  !======================================================================================!

  FUNCTION motion_rotation( it ) RESULT(qref)

    ! motion_rotation = - cross(omegab,x-ro)

    REAL(KIND(0.D0)), DIMENSION(Nq,mgridlev) :: qref
    INTEGER                                  :: it, k, i, j
    REAL(KIND(0.D0)) :: fac , xx, yy, omegab
    REAL(KIND(0.D0)), DIMENSION(5) :: uv

    uv = motion_grid(it)
    omegab = uv(3)
    IF(ADJOINT.AND.FORWARD_SIM)THEN  
      DO k=1,mgridlev  ! one for each grid 
         fac                  = delta*2.d0**(k-1)  ! cell face length on each grid
         DO i=1,m+1
            DO j=1,n
               xx = (REAL(i)-1-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
               yy = (REAL(j)-0.5D0-n/2)*delta*(2**(k-1)) + n/2*delta -offsety
               qref(u(i,j),k) = qref(u(i,j),k) + (omegab+control(it,2*nbodyf+2*nact+3))*( yy - roy )
            END DO
         END DO
         DO i=1,m
            DO j=1,n+1
               xx = (REAL(i)-0.5d0-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
               yy = (REAL(j)-1-n/2)*delta*(2**(k-1)) + n/2*delta -offsety
               qref(v(i,j),k) = qref(v(i,j),k) - (omegab+control(it,2*nbodyf+2*nact+3))*( xx - rox )
            END DO
         END DO    
         qref(:,k) = fac*qref(:,k)
      ENDDO
    ELSE
      DO k=1,mgridlev  ! one for each grid 
         fac                  = delta*2.d0**(k-1)  ! cell face length on each grid
         DO i=1,m+1
            DO j=1,n
               xx = (REAL(i)-1-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
               yy = (REAL(j)-0.5D0-n/2)*delta*(2**(k-1)) + n/2*delta -offsety
               qref(u(i,j),k) = qref(u(i,j),k) + omegab*( yy - roy )
            END DO
         END DO
         DO i=1,m
            DO j=1,n+1
               xx = (REAL(i)-0.5d0-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
               yy = (REAL(j)-1-n/2)*delta*(2**(k-1)) + n/2*delta -offsety
               qref(v(i,j),k) = qref(v(i,j),k) - omegab*( xx - rox )
            END DO
         END DO    
         qref(:,k) = fac*qref(:,k)
      ENDDO    
    
    ENDIF
    
  END FUNCTION motion_rotation

  !======================================================================================!
  
  
END MODULE variables
