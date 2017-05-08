 MODULE operators_pressure

  !******************************************************************!
  !*   Module to solve for the total pressure (p+.5*u^2), from a    *!
  !*   given velocity and vorticity field                           *!
  !*   The pressure is defined at the cell centers                  *!
  !*   New features                                                 *!
  !*     - Calculate pressure using Sine FFT                        *!
  !*   New updates 2013.01.17 by Hsieh-Chen Tsai                    *!
  !*     - Calculate pressure in the rotating frame                 *!
  !******************************************************************!

  USE parameters
  USE grid
  IMPLICIT NONE

  REAL(KIND(0.D0)), DIMENSION(:,:,:), allocatable :: pressure !global variable to store static pressure

CONTAINS

!************************************************************************************!
   subroutine setup_pressure()

      allocate( pressure(m,n,mgridlev) )
      pressure = 0.0d0

   end subroutine setup_pressure
!************************************************************************************!


!************************************************************************************!
   subroutine destroy_pressure()

      deallocate( pressure )        

   end subroutine destroy_pressure
!************************************************************************************!


!************************************************************************************!
 SUBROUTINE calculate_pressure( itime )

   !***************************************************************!
   !*  Multiscale method to solve D D^T pressure = pressure_rhs   *!
   !*  D D^T pressure = D N(q) - D E^T f - d/dt(mdot) + L'mdot    *!
   !***************************************************************!
   USE myfft
   USE variables
   USE operators
   USE user

   REAL(KIND(0.D0)), DIMENSION(1:m+1,1:n+1,2,mgridlev) :: omegaq ! (at vertices)
   REAL(KIND(0.D0)), DIMENSION(1:Nq,mgridlev)          :: nonlinearq      ! (at edges)
   REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1))        :: omegaq_ubc, omegaq_vbc ! omegaq at the boundaries
   REAL(KIND(0.D0)), DIMENSION(1:m,1:n)                :: pressure_rhs      ! (at centers)
   REAL(KIND(0.D0)), DIMENSION(2*(m)+2*(n))            :: pressure_bc   ! pressure at the boundaries
   REAL(KIND(0.D0)), DIMENSION(1:m,1:n,mgridlev)       :: pressure_rhs_temp
!   REAL(KIND(0.D0)), DIMENSION(1:m,1:n,mgridlev)       :: pressure_stag, pressure_dyn ! pressure_stag = pressure + pressure_dyn (at centers)
   REAL(KIND(0.D0))                                    :: del
   REAL(KIND(0.D0)), DIMENSION(5) :: uv
   REAL(KIND(0.D0)) :: xx, yy , omegab
   INTEGER                                             :: i,j,k,itime
   
   ! ==================================================================
   ! ================ Nonlinear term N(q)= q x omega ==================
   ! ==================================================================	
   ! omegaq(1:m+1, 1:n+1, 1:2, 1:mgridlev)

   k          = mgridlev  ! === largest domain === !
   omegaq_ubc = 0.d0
   omegaq_vbc = 0.d0

   omegaq(:,:,:,k) = nonlinear_corner( omega(:,:,k), q(:,k), q0(:,k), omegaq_ubc, omegaq_vbc )
   
   nonlinearq(:,k) = nonlinear_corner2edge( omegaq(:,:,:,k) )

   DO k=mgridlev-1,1,-1 ! === telescope into smaller domain === !
      CALL get_bc( omegaq(:,:,1,k+1), omegaq_ubc, 1.d0 )
      CALL get_bc( omegaq(:,:,2,k+1), omegaq_vbc, 1.d0 )
      omegaq(:,:,:,k) = nonlinear_corner( omega(:,:,k), q(:,k), q0(:,k), omegaq_ubc, omegaq_vbc )
      nonlinearq(:,k) = nonlinear_corner2edge( omegaq(:,:,:,k) )
   END DO

   DO k=1,mgridlev
      del = delta * REAL( 2**(k-1) )
      nonlinearq(:,k) = nonlinearq(:,k) / (del**2) ! since omega=omega/(delta**2)
   END DO	

   pressure_rhs_temp = 0.0d0
   pressure_rhs_temp(:,:,1) = 0.d0 - divergence( reg(fb) )/dt ! E^T(fb)
   pressure_rhs_temp(:,:,1) = 0.d0 + divergence( rhs_forcing(itime,q(:,1)) ) 

   DO k=2,mgridlev   ! coarsify onto bigger domain
      pressure_rhs_temp(:,:,k) = coarsify_pressure(pressure_rhs_temp(:,:,k), pressure_rhs_temp(:,:,k-1), 1.d0)
   ENDDO

   DO k=1,mgridlev   ! Add D(N(q)) to RHS at each level
      pressure_rhs_temp(:,:,k) = pressure_rhs_temp(:,:,k) + divergence( nonlinearq(:,k) )
   ENDDO

   ! =====================================================================
   ! ================ Solve for pressure =================================
   ! =====================================================================
   ! ==== DD^T (0.5*|u|^2 + P + 0.5*omegab^2*|x|^2) = RHS           ==== !
   ! == At boundary, P=0 --> BC = 0.5*|q0/del|^2 -0.5*omegab^2*|x|^2  == !

   k           = mgridlev
   del         = 0.125D0 / ( delta * 2.d0**(k-1) )**2
   pressure_bc = 0.0d0
   
   ! user-defined angular velocity of the grid  
   uv = motion_grid(itime)
   omegab = uv(3)

   ! TOP and BOTTOM
   DO i=1,m
      xx = (REAL(i)-0.5D0-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
      j = 1
      yy = (REAL(j)-0.5D0-n/2)*delta*(2**(k-1)) + n/2*delta - offsety
      pressure_bc(bottom_phi + i) = del * ( (q(u(i+1,j),k)+q(u(i,j),k)+q0(u(i+1,j),k)+q0(u(i,j),k))**2 + & 
                                      (q(v(i,j+1),k)+q(v(i,j),k)+q0(v(i,j+1),k)+q0(V(i,j),k))**2 ) - &
                                      0.5*omegab**2*( (xx-rox)**2+(yy-roy)**2 )
      j = n
      yy = (REAL(j)-0.5D0-n/2)*delta*(2**(k-1)) + n/2*delta - offsety
      pressure_bc(top_phi + i)    = del * ( (q(u(i+1,j),k)+q(u(i,j),k)+q0(u(i+1,j),k)+q0(u(i,j),k))**2 + & 
                                      (q(v(i,j+1),k)+q(v(i,j),k)+q0(v(i,j+1),k)+q0(V(i,j),k))**2 ) - &
                                      0.5*omegab**2*( (xx-rox)**2+(yy-roy)**2 )
   ENDDO
   ! LEFT and RIGHT
   DO j=1,n
      yy = (REAL(j)-0.5D0-n/2)*delta*(2**(k-1)) + n/2*delta - offsety
      i = 1
      xx = (REAL(i)-0.5D0-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
      pressure_bc(left_phi + j)  = del * ( (q(u(i+1,j),k)+q(u(i,j),k)+q0(u(i+1,j),k)+q0(u(i,j),k))**2 + & 
                                      (q(v(i,j+1),k)+q(v(i,j),k)+q0(v(i,j+1),k)+q0(V(i,j),k))**2 ) - &
                                      0.5*omegab**2*( (xx-rox)**2+(yy-roy)**2 )
      i = m
      xx = (REAL(i)-0.5D0-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
      pressure_bc(right_phi + j) = del * ( (q(u(i+1,j),k)+q(u(i,j),k)+q0(u(i+1,j),k)+q0(u(i,j),k))**2 + & 
                                      (q(v(i,j+1),k)+q(v(i,j),k)+q0(v(i,j+1),k)+q0(V(i,j),k))**2 ) - &
                                      0.5*omegab**2*( (xx-rox)**2+(yy-roy)**2 )
   ENDDO

   ! ======== solve for pressure at mgridlev (coarsest grid) ========= !
   !before here pressure_rhs_temp should be grid size independent
   pressure_rhs = pressure_rhs_temp(:,:,k)         ! pressure_rhs at k level
   CALL apply_bc_phi( pressure_rhs, pressure_bc)   ! apply pressure_bc    
   pressure(:,:,mgridlev) = ddti( pressure_rhs )   ! pressure=(D D^T)-1 rhs(:,:,mgridlev)

   ! ========= telescope in =========================== !	
   DO k=mgridlev-1,1,-1  ! === Telescope in to solve for pressure on SMALLER DOMAIN

      CALL get_bc_phi( pressure(:,:,k+1), pressure_bc, 1.d0) ! get pressure_bc from k+1 level
      pressure_rhs = pressure_rhs_temp(:,:,k)       ! pressure_rhs at k level
      CALL apply_bc_phi( pressure_rhs, pressure_bc) ! apply pressure_bc
      pressure(:,:,k) = ddti( pressure_rhs )        ! compute new pressure at k level
     
   END DO
   ! ================================================== !	

   !pressure_stag = pressure ! pressure_stag = pressure_static + pressure_dynamic [== 0.5 rho |u|^2 ]

   DO k=1,mgridlev
      del = 0.125D0 / ( delta * 2.d0**(k-1) )**2
      DO j=1,n
         DO i=1,m                      ! We solved for P* : P* = P + (1/2) |u|^2 
            !pressure_dyn(i,j,k) = del*( ( q(u(i+1,j),k) +q(u(i,j),k) +q0(u(i+1,j),k) +q0(u(i,j),k) )**2 + &
            !                            ( q(v(i,j+1),k) +q(v(i,j),k) +q0(v(i,j+1),k) +q0(v(i,j),k) )**2 )
            xx = (REAL(i)-0.5D0-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
            yy = (REAL(j)-0.5D0-n/2)*delta*(2**(k-1)) + n/2*delta - offsety
            pressure(i,j,k) = pressure(i,j,k) - &
                              del * ( (q(u(i+1,j),k)+q(u(i,j),k)+q0(u(i+1,j),k)+q0(u(i,j),k))**2 + & 
                                      (q(v(i,j+1),k)+q(v(i,j),k)+q0(v(i,j+1),k)+q0(V(i,j),k))**2 ) + &
                              0.5*omegab**2*( (xx-rox)**2+(yy-roy)**2 ) 
         END DO
      END DO
   END DO

   ! scale pressure to the pressure coefficient
   pressure = 2.d0*pressure

END SUBROUTINE calculate_pressure
!************************************************************************************!


!************************************************************************************!
FUNCTION nonlinear_corner( omega, q, q0, omegaq_ubc, omegaq_vbc ) RESULT(omegaq)

   !***************************************************************!
   !*   nonlinear terms in rotational form (q x omega)            *!
   !*   (q x omega) at the corner [1:m+1, 1:n+1, 1:2]             *!
   !***************************************************************!
   REAL(KIND(0.D0)), DIMENSION(2:m,2:n)         :: omega
   REAL(KIND(0.D0)), DIMENSION(Nq)              :: q, q0
   REAL(KIND(0.D0)), DIMENSION(1:m+1,1:n+1,2)   :: omegaq 
   REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: omegaq_ubc, omegaq_vbc
   INTEGER                                      :: i,j

   DO j=2,n
      DO i=2,m
         ! x-dir :  v*omega
         omegaq(i,j,1) =  0.5D0 * ( q(v(i,j))+q(v(i-1,j)) + q0(v(i,j))+q0(v(i-1,j)) ) * omega(i,j)
         ! y-dir : -u*omega
         omegaq(i,j,2) = -0.5D0 * ( q(u(i,j))+q(u(i,j-1)) + q0(u(i,j))+q0(u(i,j-1)) ) * omega(i,j)
      ENDDO
    ENDDO

    DO j=1,n+1
       i = 1   ! ==== left boundary==== !
       omegaq(i,j,1) = omegaq_ubc(left+j)
       omegaq(i,j,2) = omegaq_vbc(left+j)
       i = m+1  ! ==== right boundary ==== !
       omegaq(i,j,1) = omegaq_ubc(right+j)
       omegaq(i,j,2) = omegaq_vbc(right+j)
    ENDDO

    DO i=1,m+1
       j = 1    ! ==== bottom boundary ==== !
       omegaq(i,j,1) = omegaq_ubc(bottom+i)
       omegaq(i,j,2) = omegaq_vbc(bottom+i)
       j = n+1  ! ==== top boundary ==== !
       omegaq(i,j,1) = omegaq_ubc(top+i)
       omegaq(i,j,2) = omegaq_vbc(top+i)
    ENDDO

END FUNCTION nonlinear_corner
!************************************************************************************!


!************************************************************************************!
FUNCTION nonlinear_corner2edge( omegaq ) RESULT(nonlinearq)

    !***************************************************************!
    !*   nonlinear terms in rotational form (q x omega)            *!
    !*   taking average from the corner to the edge [Nq]           *!
    !***************************************************************!

    REAL(KIND(0.D0)), DIMENSION(1:m+1,1:n+1,2) :: omegaq 
    REAL(KIND(0.D0)), DIMENSION(Nq)            :: nonlinearq
    INTEGER                                    :: i,j

    DO j=1,n
       DO i=1,m+1
           nonlinearq(u(i,j)) = 0.5D0*( omegaq(i,j+1,1) + omegaq(i,j,1) )
       END DO
    END DO

    DO j=1,n+1
       DO i=1,m
           nonlinearq(v(i,j)) = 0.5D0*( omegaq(i+1,j,2) + omegaq(i,j,2) )
       END DO
    END DO

END FUNCTION nonlinear_corner2edge
!************************************************************************************!


!************************************************************************************!
FUNCTION divergence( x )

  !***************************************************************!
  !*   divergence(x), x at edge & divergence at center           *!
  !***************************************************************!

   REAL(KIND(0.D0)), DIMENSION(Nq)      :: x
   REAL(KIND(0.D0)), DIMENSION(1:m,1:n) :: divergence
   INTEGER                              :: i,j

   DO j=1,n
      DO i=1,m
         divergence(i,j) = x(u(i+1,j)) - x(u(i,j)) + x(v(i,j+1)) - x(v(i,j))
      ENDDO
   ENDDO

END FUNCTION divergence
!************************************************************************************!


!************************************************************************************!
SUBROUTINE write_pressure(pressure, it)

  USE parameters
  USE variables 
  CHARACTER(7) :: charit
  INTEGER :: it
  REAL(KIND(0.D0)), DIMENSION(1:m,1:n,mgridlev) :: pressure
  
  WRITE(charit,"(I7.7)") it
  OPEN(unit=100,file="output/pressure"//charit//".var",form="unformatted",status="unknown")
  WRITE(100) pressure
  CLOSE(100)

END SUBROUTINE write_pressure
!************************************************************************************!


!************************************************************************************!
SUBROUTINE read_pressure(pressure, it)

    USE parameters
    USE variables
    CHARACTER(7) :: charit
    INTEGER :: it
    REAL(KIND(0.D0)), DIMENSION(1:m,1:n,mgridlev) :: pressure

    WRITE(charit,"(I7.7)") it
    OPEN(unit=100,file="output/pressure"//charit//".var",form="unformatted",status="unknown")
    READ(100) pressure
    CLOSE(100)

  END SUBROUTINE read_pressure
!************************************************************************************!

!***************************************************************************************!
FUNCTION coarsify_pressure( crhs, rhs, fac ) RESULT( arhs )
    
    !***************************************************************!
    !*   given mass source on a smaller, fine mesh, (rhs) interp.  *!
    !*   values to the center region of a larger, coarser mesh     *!
    !*   (crhs).  The values outside the center region are         *!
    !*   not unmodified. Result is placed in arhs                  *!
    !***************************************************************!
    REAL(KIND(0.D0)), DIMENSION(:,:)                     :: crhs, rhs
    REAL(KIND(0.D0)), DIMENSION(SIZE(rhs,1),SIZE(rhs,2)) :: arhs
    REAL(KIND(0.D0))                                     :: fac
    INTEGER                                              :: i,j,indi,indj
    
    arhs = crhs
    DO j=-n/4+1, n/4
       indj = n/2+2*j
       DO i=-m/4+1, m/4
          indi = m/2+2*i
          arhs(m/2+i,n/2+j) = ( rhs(indi-1, indj-1) + &  ! Left Bottom
                                rhs(indi-1, indj  ) + &  ! Left Top
                                rhs(indi,   indj-1) + &  ! Right Bottom
                                rhs(indi,   indj  )   )  ! Right Top
       ENDDO
    ENDDO

    arhs = fac * arhs

END FUNCTION coarsify_pressure
!***************************************************************************************!


!***************************************************************************************!
SUBROUTINE apply_bc_phi( r, rbc)
    
    !*****************************************************************!
    !*   given phi right outside of domain, phi_bc, (from larger,    *!
    !*   coraser mesh), add values to correct laplacian(D*D^T)of phi *!
    !*   , mass_rhs   on the (smaller, finer) domain, r.             *!
    !*****************************************************************!
    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:)   :: rbc
    INTEGER                          :: i,j

    ! add bc's from coarser grid
    ! TOP and BOTTOM
    DO i=1,m
       r(i,1) = r(i,1) - rbc(bottom_phi + i)
       r(i,n) = r(i,n) - rbc(   top_phi + i)
    END DO

    ! LEFT and RIGHT
    DO j=1,n
       r(1,j) = r(1,j) - rbc( left_phi + j)
       r(m,j) = r(m,j) - rbc(right_phi + j)
    END DO
    
  END SUBROUTINE apply_bc_phi
!************************************************************************************!


!***************************************************************************************!
SUBROUTINE get_bc_phi( r, rbc, fac)

    !******************************************************************!
    !*   given potential,phi on a larger, coarser mesh,               *!
    !*   interpolate it's values to the edge of a smaller, finer mesh *!
    !******************************************************************!

    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:)   :: rbc
    REAL(KIND(0.D0))                 :: fac
    INTEGER                          :: i,j

    ! get interpolated boundary conditions on finer grid
    ! TOP and BOTTOM
    DO i=0,m-2,2   
       rbc(bottom_phi + i+1) = .25d0*( r(m/4+i/2,n/4) + &
                               2.d0* r(m/4+1+i/2,n/4) + r(m/4+1+i/2,n/4+1) )
       rbc(top_phi    + i+1) = .25d0*( r(m/4+i/2,n/4+n/2+1) + &
                               2.d0* r(m/4+1+i/2,n/4+n/2+1) + r(m/4+1+i/2,n/4+n/2) )
    END DO

    DO i=1,m-1,2
       rbc(bottom_phi +i+1)  = .25d0*( r(m/4+2+(i-1)/2,n/4) + &
                                 2.d0* r(m/4+1+(i-1)/2,n/4) + r(m/4+1+(i-1)/2,n/4+1) )
       rbc(top_phi    +i+1)  = .25d0*( r(m/4+2+(i-1)/2,n/4+n/2+1) + &
                                 2.d0* r(m/4+1+(i-1)/2,n/4+n/2+1) + r(m/4+1+(i-1)/2,n/4+n/2) )
    END DO

    ! LEFT and RIGHT
    DO j=0,n-2,2
       rbc(left_phi  + j+1)  = .25d0*( r(m/4,n/4+j/2) + &
                                 2.d0* r(m/4,n/4+1+j/2) + r(m/4+1,n/4+1+j/2) )
       rbc(right_phi + j+1)  = .25d0*( r(m/4+m/2+1,n/4+j/2) + &
                                 2.d0* r(m/4+m/2+1,n/4+1+j/2) + r(m/4+m/2,n/4+1+j/2) )
    END DO
    DO j=1,n-1,2
       rbc(left_phi  + j+1)  = .25d0*( r(m/4,n/4+(j-1)/2) + &
                                 2.d0* r(m/4,n/4+1+(j-1)/2) + r(m/4+1,n/4+1+(j-1)/2) )
       rbc(right_phi + j+1)  = .25d0*( r(m/4+m/2+1,n/4+2+(j-1)/2) + &
                                 2.d0* r(m/4+m/2+1,n/4+1+(j-1)/2) + r(m/4+m/2,n/4+1+(j-1)/2) )
    END DO

    rbc = rbc * fac

  END SUBROUTINE get_bc_phi
!************************************************************************************!


END MODULE operators_pressure
