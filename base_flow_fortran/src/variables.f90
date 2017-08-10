MODULE variables

  ! c Kunihiko (Sam) Taira and Tim Colonius
  ! December 2005
  ! calling program for ibfs v1.0
  ! ibfs v2.1 modified by Hsieh-Chen Tsai
  ! January 2013
  ! - now solves flows with arbitrary body motion in body-fixed frame

  USE grid
  USE user
  IMPLICIT NONE

  ! in what follows, the last index refers to the grid level, first 1 (or 2) indices to point in space
  
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: omega, s     ! vorticity 
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: rhs_old   ! nonlinear terms at previous time-step (for AB2 integration)
  REAL(KIND(0.D0)), DIMENSION(:,:)  , ALLOCATABLE :: q, q0, q0p, q0r  
      ! fluxes and additional potential flow and solid body rotation be added

  ! variables for Cholesky (stationary body)
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: cholmat
  REAL(KIND(0.D0)), DIMENSION(:)  , ALLOCATABLE :: cholvec


CONTAINS
  
  !======================================================================================!

  SUBROUTINE setup_variables

    USE parameters
    USE myfft
    INTEGER :: k
    
    CALL setup_fft

    ALLOCATE( omega(2:m,2:n,mgridlev), rhs_old(2:m,2:n,mgridlev), s(2:m,2:n,mgridlev) )

    ALLOCATE( q(Nq,mgridlev), q0(Nq,mgridlev), q0p(Nq,mgridlev), q0r(Nq,mgridlev) )

    !IF (stationary) THEN
       ALLOCATE( cholmat(Nf,Nf), cholvec(Nf) )
    !END IF




    





    IF (istart==0) THEN
       CALL initial_condition
       CALL setup_reg(xb)
    ELSE
       CALL read_variables(istart)
       CALL setup_reg(xb)
       !IF (stationary) THEN
       !   CALL read_cholesky
       !END IF
    END IF

  END SUBROUTINE setup_variables

  !======================================================================================!
 
  SUBROUTINE destroy_variables

    USE parameters
    USE myfft

    CALL destroy_fft
    DEALLOCATE( omega, rhs_old, q, q0, q0p, q0r, s)
    IF (stationary) THEN
       DEALLOCATE( cholmat, cholvec )
    END IF

  END SUBROUTINE destroy_variables

  !======================================================================================!
 
  SUBROUTINE destroy_grid

    USE parameters

    DEALLOCATE( f, u, v )
    DEALLOCATE( xb,vb,fb, codeb )

  END SUBROUTINE destroy_grid

  !================================================================
 
  SUBROUTINE write_cholesky

    USE parameters

    OPEN(unit=100,file="output/ib.chd",form="unformatted",status="unknown")
    WRITE(100) cholmat, cholvec
    CLOSE(100)

  END SUBROUTINE write_cholesky
  
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
    OPEN(unit=100,file="output/ib"//charit//".var",form="unformatted",status="unknown")

    WRITE(100) m,n,mgridlev,nb
    WRITE(100) re,dt,len,offsetx,offsety
    WRITE(100) omega,xb,vb,fb,codeb,rhs_old, q, q0p, q0r
    WRITE(100) s
    WRITE(100) rot_angle, rox, roy
    write(100) xb0, u_ib, ud_ib, udd_ib
    write(100) f_rdst
    CLOSE(100)
 
  END SUBROUTINE write_variables
 
 !======================================================================================!
 
  SUBROUTINE read_variables(it)

    USE parameters
    CHARACTER(7) :: charit
    INTEGER :: it
    REAL(KIND(0.D0)) :: dt_old

    write(*,*) 'reading variables at it=',it
    WRITE(charit,"(I7.7)") it
    OPEN(unit=100,file="output/ib"//charit//".var",form="unformatted",status="unknown")

    READ(100) m,n,mgridlev,nb
    READ(100) re,dt_old,len,offsetx,offsety
    READ(100) omega,xb,vb,fb,codeb,rhs_old, q, q0p, q0r
    READ(100) s
    READ(100) rot_angle, rox, roy
    read(100) xb0, u_ib, ud_ib, udd_ib
    read(100) f_rdst
    CLOSE(100)




    q0 = q0p + q0r

  END SUBROUTINE read_variables
  
  !======================================================================================!
 
  SUBROUTINE initial_condition

    USE parameters
    USE user
    CHARACTER(7) :: charit
    REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: xold,vold,fold
    INTEGER, DIMENSION(:), ALLOCATABLE :: codeold
    INTEGER :: m_in,n_in,mgridlev_in,nb_in
    REAL(KIND(0.D0)) :: re_in,dt_in,len_in,offsetx_in,offsety_in
    REAL(KIND(0.D0)), DIMENSION(5) :: uv
    LOGICAL :: readic

    INQUIRE(file="input/initial.var",exist=readic)

    IF (readic) THEN
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

    ELSE

       q = 0.d0
       omega = 0.d0
       rhs_old = 0.d0
       q0p = motion_potential(0)
       q0r = motion_rotation(0)
       q0 = q0p + q0r
       rot_angle = 0.d0
       uv = motion_grid( 0 )
       rox = uv(4)
       roy = uv(5)

    END IF

  END SUBROUTINE initial_condition

  !=======================================================================================

  FUNCTION motion_potential( it ) RESULT(qref)

    ! motion_potential = - u_b

    REAL(KIND(0.D0)), DIMENSION(Nq,mgridlev) :: qref
    INTEGER                                  :: it, k, i, j
    REAL(KIND(0.D0)) :: fac
    REAL(KIND(0.D0)), DIMENSION(5) :: uv

    uv = motion_grid(it)
    
    qref = 0.d0
    DO k=1,mgridlev  ! one for each grid 
       fac                  = delta*2.d0**(k-1)  ! cell face length on each grid
       qref(1:(m+1)*n,k)    = -uv(1)
       qref((m+1)*n+1:Nq,k) = -uv(2)
       qref(:,k) = fac*qref(:,k)
    ENDDO
    
  END FUNCTION motion_potential

!=======================================================================================

  FUNCTION motion_rotation( it ) RESULT(qref)

    ! motion_rotation = - cross(omegab,x-ro)

    REAL(KIND(0.D0)), DIMENSION(Nq,mgridlev) :: qref
    INTEGER                                  :: it, k, i, j
    REAL(KIND(0.D0)) :: fac , xx, yy, omegab
    REAL(KIND(0.D0)), DIMENSION(5) :: uv

    uv = motion_grid(it)
    omegab = uv(3)
    
    qref = 0.d0
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
    
  END FUNCTION motion_rotation

!=======================================================================================

END MODULE variables
