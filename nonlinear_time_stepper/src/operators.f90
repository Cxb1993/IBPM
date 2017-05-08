     MODULE operators

  !******************************************************************!
  !*                               IBFS 2.0                         *!
  !*                 c. Tim Colonius & Kunihiko (Sam) Taira         *!
  !*                 California Institute of Technology             *!
  !*                             March 1, 2007                      *!
  !******************************************************************!
  !*                                                                *!
  !*   New features                                                 *!
  !*     - completely revamped solution method                      *!
  !*     - discrete streamfunction/vorticity formulation            *!
  !*     - super grid solution of Lapalce eqn. for better           *!
  !*       far field boundary conditions                            *!
  !*     - Fast sine transform to invert Laplacians                 *!
  !*     - Cholesky factorization for solution of stationary        *!
  !*       bodies                                                   *!
  !*     - 2D for now                                               *!
  !*                                                                *!
  !******************************************************************!
  ! ibfs v2.1 modified by Hsieh-Chen Tsai
  ! January 2013 
  ! - now solves flows with arbitrary body motion in body-fixed frame

  USE parameters
  USE grid
  IMPLICIT NONE


CONTAINS

SUBROUTINE preprocess

USE variables

!***************************************************************!
!*   Cholesky factorization of EC(C^t C)^-1 C^t E^t            *!
!*      performed once only for stationary bodies              *!
!***************************************************************!

REAL(KIND(0.D0)), DIMENSION(Nf) :: z   ! temporary force vector
INTEGER :: i, ii
LOGICAL :: readchol
integer :: info, neqns, lda, j, lwork, jj
real(kind(0.d0)), dimension(:), allocatable :: work
integer, dimension(3*nb) :: ipiv_bg

!    INQUIRE(file='output/ib.chd',exist=readchol)  ! already done?
!    IF (readchol) THEN
!       CALL read_cholesky
!    ELSE

WRITE(*,*) 'precomputing body matrix for stationary geometry...please be patient'
DO i=1,nf   ! build matrix one column at a time
!WRITE(*,"(A,I5)",ADVANCE="NO") '...',i
z = 0.d0
z(i) = 1.d0
cholmat(1:Nf,i) = a_times( z )
END DO

do i = 1, nf
    do ii = 1, nf

        write(617,*) cholmat(i,ii)
    end do
end do





!WRITE(*,*) 'performing cholesky decomp'
CALL choldc
CALL write_cholesky  ! save results


!    END IF

END SUBROUTINE preprocess

  !*****************************************************************!

  SUBROUTINE advance(itime)

    !***************************************************************!
    !*    Advance to next time level                               *!
    !***************************************************************!

    USE myfft
    USE variables
    USE user

    INTEGER :: itime, k, kk

    ! intermediate results
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: rhs, rhsbc
    REAL(KIND(0.D0)), DIMENSION(Nf) :: accel, rhsf
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: omega_bc
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1),mgridlev) :: lastbc
    REAL(KIND(0.D0)), DIMENSION(nq) :: nl_temp
    real(kind(0.d0)), dimension(1:3*nb) :: Fint, F_bf, u0, ud0, &
            r_ub, r_cb, F_bg, d_xb, udd0
    real(kind(0.d0)), dimension(nf) :: F_sm
    real(kind(0.d0)), dimension(nb) :: f_rd_x, f_rd_y
    real(kind(0.d0)) :: err_fsi, tol_fsi
    integer :: fsi_count
    integer :: info, neqns, lda, lwork
    real(kind(0.d0)), dimension(:), allocatable :: work
    integer, dimension(3*nb) :: ipiv_bg



    itime = itime + 1

    IF ( MOD(itime,10).eq.0 ) THEN
       WRITE(*,*) "...Advancing to itime =",itime
    END IF

    q0p = motion_potential(itime-1)
    q0r = motion_rotation(itime-1)
    q0 = q0p + q0r
    rot_angle = rot_angle + delta_angle(itime-1)

    IF (itime==1 .and. stationary) THEN
       CALL preprocess
    END IF

    ! Step 1: solve intermediate curl(momentum) eq. on each grid
 
    DO k=1,mgridlev-1
       CALL get_bc(omega(:,:,k+1), lastbc(:,k), 0.25d0)
    END DO
    lastbc(:,mgridlev) = 0.d0
    
    DO k = mgridlev,1,-1

       IF (k.lt.mgridlev) THEN
          CALL get_bc( omega(:,:,k+1), omega_bc, 0.25d0 )
       ELSE
          omega_bc = 0.d0
       END IF

       ! Computing the nonliner term (nonlinear term mult. by dt/2)
       nl_temp = nonlinear( omega(:,:,k), q(:,k), q0(:,k), lastbc(:,k) )


       ! add user-defined RHS forcing term to momentum eq.
       IF (k==1) THEN
           nl_temp = nl_temp + rhs_forcing( itime, q(:,k) )
        END IF

       !rhs = rot( nonlinear( omega(:,:,k), q(:,k), q0(:,k), lastbc(:,k) ) )
       rhs = rot( nl_temp )
    
       ! If this is the very first time step, we need to use explicit Euler
       IF (itime==1) THEN
          rhs_old(:,:,k) = rhs(:,:)
       END IF

       rhsbc = 0.d0

       CALL apply_bc( rhsbc, lastbc(:,k), vfac(k) )
       CALL apply_bc( rhsbc, omega_bc,    vfac(k) )

       omega(:,:,k) = dst( lam1i(:,:,k) * &
          ( dst( con1(k)*rhs(:,:)       + &
                 con2(k)*rhs_old(:,:,k) + &
                         rhsbc(:,:)   ) + & 
                lam1(:,:,k)*dst(omega(:,:,k)) ) )

       ! update for next time step
       rhs_old(:,:,k) = rhs(:,:)

    END DO

    ! we now have vorticity evolution on all domains.  For smallest domain
    ! we add the vorticity created by the body

    ! first see if there is boundary motion and update regularization kernel if nec.
    IF (.not.stationary) THEN

        !Initialize for FSI loop:
        FSI_count = 0
        err_fsi = 10.d0
        tol_fsi = 1.d-5

        call get_M( Mmat )
        call var_update( u_ib, Fint )
        call get_K( Kmat )
        CALL setup_reg(xb)

        !Body force on solid:
            call bf_sol( f_rd_x, f_rd_y )
            F_bf = 0.d0
            call build_force( f_rd_x, f_rd_y, F_bf )

        !Store pos and vel at previous time step for iteration:
        u0 = u_ib
        ud0 = ud_ib


        !Need to get a consistent acceleration if at 1st time step:
        if (itime .eq. 1) then

            print *, "computing consistent initial body acceleration..."

            info = 0
            lwork = (3*nb)** 2
            allocate( work( lwork ) )
            neqns = 3*nb
            lda = 3*nb
            ipiv_bg = 0

            call dgetrf( neqns, lda, Mmat, lda, ipiv_bg, info)
            call dgetri( neqns, Mmat, lda, ipiv_bg, work, lwork, info)
            deallocate( work )

            udd_ib = matmul( Mmat, -Fint + F_bf + QWx( fb ) )

            call get_M( Mmat )

            print *, "initial acceleration computed"
        end if

        udd0 = udd_ib

        do while ( (err_fsi .ge. tol_fsi) .and. (fsi_count .le. 1000) )



                !-- step 1 : compute guess for surface stress:
                !---invert solid matrix:
                sol_mat = Kmat + 4.d0/dt**2.d0 * Mmat

                info = 0
                lwork = (3*nb)** 2
                allocate( work( lwork ) )
                neqns = 3*nb
                lda = 3*nb
                ipiv_bg = 0

                call dgetrf( neqns, lda, sol_mat, lda, ipiv_bg, info)
                call dgetri( neqns, sol_mat, lda, ipiv_bg, work, lwork, info)

                deallocate( work )

                !-----

                !---rhs contribution from fluid eqns:
                r_cb = 2.d0/dt * (u_ib - u0) - ud0
                !---


                !---rhs contribution from fluid eqns:

                !explicit terms from position treatment
                r_ub = matmul( Mmat, udd0 + 4.d0/dt * ud0 + 4.d0/dt**2.d0 * &
                    (u0 - u_ib) ) - Fint + F_bf

                r_ub = matmul( sol_mat, r_ub )
                !---

                !-- Put rhs together and truncate to be size of surface stress vector:
                F_bg = -1.d0 * delta * ( 2.d0/dt * r_ub + r_cb )

                F_sm = v_trunc( F_bg )

                if (mgridlev .ge. 2) then
                    CALL vort2flux( q, omega, s, mgridlev)
                else
                    CALL vort2flux( q, omega, s, mgridlev)
                end if

                rhsf = regT(q(:,1) ) + regT(q0(:,1)) + F_sm
                !---


                !--solve linear system for surface stress
                    !CALL cjgr( fb, rhsf )
                    call bicgstab( fb, rhsf )
                    !Get physical surface stress:
                    f_rdst = Redistribute( fb )

                !--
            !---

            !---step 2: compute increment in body position (dx)

                d_xb = r_ub + matmul( sol_mat, QWx( fb ) )

                !Check error:
                err_fsi = maxval(abs(d_xb))/maxval(abs(u_ib))

!                print *, "FSI err = ", err_fsi
            !---

            !---step 3 : update body positions and velocities and truncate to
            !-------     be size of fluid stresses

                u_ib = u_ib + d_xb
                ud_ib = -ud0 + 2.d0/dt * (u_ib - u0)
                udd_ib = 4.d0/dt**2.d0 * (u_ib - u0) - 4.d0/dt * ud0 - udd0

                xb = xb0 + v_trunc( u_ib )
                vb = vb0 + v_trunc( ud_ib )
                ab = ab0 + v_trunc( udd_ib )

            !-----


            !----final step: update terms for next iteration if necessary:

                if (err_fsi .ge. tol_fsi) then
                    call var_update( u_ib, Fint )
                    call get_K( Kmat )
                    call setup_reg(xb)
                end if
            !----



        end do  !end of FSI loop



        !-----

    ELSE


        CALL actuators(itime,xb,xb,vb)
        if (mgridlev .ge. 2) then
            CALL vort2flux( q, omega, s, mgridlev)
        else
            CALL vort2flux( q, omega, s, mgridlev)
        end if
        rhsf = regT(q(:,1)) + regT(q0(:,1)) + delta * vb
        fb = cholsl(rhsf)

        !Redistribute the force
        f_rdst = redistribute(fb)


    END IF

    ! complete omega on grid 1
    omega(:,:,1) = omega(:,:,1) - ainv( rot( reg( fb ) ) )


    ! coarsify final omega and correct velocities on all grids
    CALL vort2flux( q, omega, s, mgridlev )
    CALL write_slip( itime, q, delta * vb - regT(q0(:,1)) )
    CALL write_udtdx( itime, q, q0 )

  END SUBROUTINE advance

  !*****************************************************************!

!**********************************************!
function QWx( x )

real(kind(0.d0)), dimension(nf) :: x, frdst, fx, fy
real(kind(0.d0)), dimension(3*nb) :: fbg, QWx


!Redistribute the force
frdst = redistribute( x )

fx = frdst(1 : nb)
fy = frdst(nb + 1 : nf)

fbg = 0.d0
call build_force( fx, fy, fbg )

QWx = fbg * 0.5d0 / dt

end function
!**********************************************

!***********************************************
function v_trunc( v_sol )

integer :: kk

real(kind(0.d0)), dimension(nf) :: v_trunc
real(kind(0.d0)), dimension(3*nb) :: v_sol

do kk = 1, nb

    v_trunc(kk) = v_sol(3*(kk-1) + 1)
    v_trunc(nb + kk) = v_sol(3*(kk-1) + 2)

end do

end function

!*****************************************************************!
function redistribute(f_vector) result(f_redistributed)

!*****************************************************************!
!*Redistributes the force in a more physically pleasing way
!*****************************************************************!

real(kind(0.d0)), dimension(nf), intent(in) :: f_vector
real(kind(0.d0)), dimension(nf) :: f_redistributed, one_vect
real(kind(0.d0)), dimension(nq) :: wght, frc_reg, f_inter
integer :: i

!Define a vector of ones to get the weights:
one_vect = 1.d0

!Get appropriate weights for redistribution
wght = reg(one_vect)

!now redistribute the force: f_rdst = E*M^-1*E^T*f, where M^-1 is
!a diagonal matrix containing the inverse of the nonzero weights wght
frc_reg = reg(f_vector)

!initialize an intermediate force vector called f_inter:
f_inter = 0.d0
do i = 1, nq

    if ( wght(i) > 1.d-10 ) then  !Only take the reciprocal of weights if they
                                  !are above a certain tolerance

        f_inter(i) = 1.d0/wght(i) * frc_reg(i)

    end if

end do


!Now that we have obtained the appropriate weights, redistribute this back onto the IB:

f_redistributed = regT(f_inter)

end function redistribute
!*****************************************************************!

  SUBROUTINE vort2flux( vel, omega, sn, nlev )

    !***************************************************************!
    !*    Multiscale method to solve C^T C s = omega               *!
    !*    and return the velocity, C s.                            *!
    !*    Results are returned in vel on each of the first nlev    *!
    !*    grids.                                                   *!
    !*    Warning: the vorticity field on all but the finest mesh  *!
    !*     is modified by the routine in the following way:        *!
    !*     the value in the center of the domain is interpolated   *!
    !*     from the next finer mesh (the value near the edge is    *!
    !*     not changed.                                            *!
    !***************************************************************!

    USE myfft
    INTEGER :: nlev
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n,mgridlev) :: omega, sn
    REAL(KIND(0.D0)), DIMENSION(Nq,mgridlev) :: vel

    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: s, vort ! streamfun & vort
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: sbc  ! streamfun

    INTEGER :: k
    ! find coarsifieid omegas

!    write(67,*) 1, SUM( omega(:,:,1))

    if (nlev .ge. 2) then
        DO k=2,nlev
           omega(:,:,k) = coarsify(omega(:,:,k),omega(:,:,k-1))
    !       write(67,*) k, SUM( omega(:,:,k))
        END DO

    end if
!    write(67,*)

    ! invert laplacian on largest grid
    ! zero bc's on s assumed for largest grid
    sbc = 0.d0
    s = 0.d0
    ! compute s
    vort = omega(:,:,nlev)
    s = ctci( vort )
    sn(:,:,nlev) = s
    ! compute vel
    vel(:,nlev) = curl(  s, sbc )

    ! now telescope in

    DO k=nlev-1,1,-1

       CALL get_bc( s, sbc, 1.d0)       ! get bc's from pervious s
       vort = omega(:,:,k)
       CALL apply_bc( vort, sbc, 1.d0)  ! apply them 
       s = ctci( vort )                 ! compute new s
       sn(:,:,k) = s
       IF (nlev.ge.k) THEN              ! compute vel if desired
          vel(:,k) = curl( s, sbc )
       END IF
      
    END DO

  END SUBROUTINE vort2flux

  !*****************************************************************!

  FUNCTION curl( x, xbc )

    !***************************************************************!
    !*   returns curl(x) given x and the values of s'fun on bdy    *!
    !***************************************************************!

    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: x
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: xbc
    REAL(KIND(0.D0)), DIMENSION(Nq) :: curl
    INTEGER :: i,j

    DO j=2,n-1
       DO i=2,m 
          curl(u(i,j)) = x(i,j+1) - x(i,j)
       ENDDO
    ENDDO
    
    DO i=2,m   
       j=1
       curl(u(i,j)) = x(i,j+1) - xbc( bottom + i )
       j=n
       curl(u(i,j)) = xbc( top + i ) - x(i,j)
    ENDDO
    
    DO j=1,n
       i = 1
       curl(u(i,j)) = xbc( left + j + 1 ) -xbc( left + j )
       i = m+1
       curl(u(i,j)) = xbc( right + j + 1 ) -xbc( right + j )
    ENDDO

    DO j=2,n
       i = 1
       curl(v(i,j)) = - x(i+1,j) + xbc( left + j )
       DO i=2,m-1
          curl(v(i,j)) = x(i,j) - x(i+1,j)
       ENDDO
       i = m
       curl(v(i,j)) = - xbc( right + j) +  x(i,j)
    ENDDO
    
    DO i=1,m
       j = 1
       curl(v(i,j)) = xbc( bottom + i ) -xbc( bottom + i + 1 )
       j = n+1
       curl(v(i,j)) = xbc( top + i ) -xbc( top + i + 1 )
    ENDDO

  END FUNCTION curl

  !*****************************************************************!

 FUNCTION rot( x )

   !***************************************************************!
   !*   Transpose of curl                                         *!
   !***************************************************************!

   REAL(KIND(0.D0)), DIMENSION(Nq) :: x
   REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: rot
   INTEGER :: i,j

   DO j=2,n
      DO i=2,m
         rot(i,j) = x(v(i,j)) - x(v(i-1,j)) - x(u(i,j)) + x(u(i,j-1))
      ENDDO
   ENDDO

 END FUNCTION rot

  !*****************************************************************!

 FUNCTION coarsify( crhs, rhs ) RESULT( arhs )
    
   !***************************************************************!
   !*   given vorticity on a smaller, fine mesh, (rhs) interp.    *!
   !*   values to the center region of a larger, coarser mesh     *!
   !*   (crhs).  The values outside the center region are         *!
   !*   not unmodified. Result is placed in arhs                  *!
   !***************************************************************!
   REAL(KIND(0.D0)), DIMENSION(:,:)                     :: crhs, rhs
   REAL(KIND(0.D0)), DIMENSION(SIZE(rhs,1),SIZE(rhs,2)) :: arhs
   INTEGER                                              :: i,j,indi,indj
    
   arhs = crhs
   DO j=-n/4+1,n/4-1
      indj = n/2+2*j
      DO i=-m/4+1,m/4-1
         indi = m/2+2*i
         arhs(m/2+i,n/2+j) = rhs(indi  ,indj)   + &
                     0.5d0*( rhs(indi+1,indj)   + rhs(indi  ,indj+1)   + &
                             rhs(indi-1,indj)   + rhs(indi  ,indj-1) ) + &
                    0.25d0*( rhs(indi+1,indj+1) + rhs(indi+1,indj-1)   + &
                             rhs(indi-1,indj-1) + rhs(indi-1,indj+1) )
       ENDDO
    ENDDO

 END FUNCTION coarsify

  !*****************************************************************!

  SUBROUTINE get_bc( r, rbc, fac)
    
    !***************************************************************!
    !*   given vorticity on a larger, coarser mesh, interpolate    *!
    !*   it's values to the edge of a smaller, finer mesh          *!
    !***************************************************************!

    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:) :: rbc
    REAL(KIND(0.D0)) :: fac
    INTEGER :: i,j
    
    ! get interpolated boundary conditions on finer grid

    DO i=0,m,2
       rbc(bottom + i+1) = r(m/4+i/2,n/4)
       rbc(top + i+1) = r(m/4+i/2,3*n/4)
    END DO
    DO i=1,m-1,2
       rbc(bottom +i+1)  = 0.5d0*( r(m/4+(i+1)/2,n/4) + r(m/4-1+(i+1)/2,n/4) )
       rbc(top + i+1) = 0.5d0*( r(m/4+(i+1)/2,3*n/4) + r(m/4-1+(i+1)/2,3*n/4) )
    END DO

    DO j=0,n,2
       rbc(left + j+1) = r(m/4, n/4+j/2)
       rbc(right + j+1) = r(3*m/4, n/4+j/2)
    END DO
    DO j=1,n-1,2
       rbc(left + j+1) = 0.5d0*( r(m/4, n/4+(j+1)/2) + r(m/4, n/4-1+(j+1)/2) )
       rbc(right + j+1) = 0.5d0*( r(3*m/4, n/4+(j+1)/2) + r(3*m/4, n/4-1+(j+1)/2) )
    END DO

    rbc = rbc*fac

  END SUBROUTINE get_bc

  !*****************************************************************!

  SUBROUTINE apply_bc( r, rbc, fac)
    
    !***************************************************************!
    !*   given vorticity at edges of domain, rbc, (from larger,    *!
    !*   coraser mesh), add values to correct laplacian of         *!
    !*   vorticity  on the (smaller, finer) domain, r.             *!
    !***************************************************************!
    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:) :: rbc
    REAL(KIND(0.D0)) :: fac
    INTEGER :: i,j
    
    ! add bc's from coarser grid

    DO i=1,m-1
       r(i,1) = r(i,1) + fac* rbc( bottom + i+1 )
       r(i,n-1) = r(i,n-1) + fac* rbc( top + i+1 )
    END DO
    DO j=1,n-1
       r(1,j) = r(1,j) + fac* rbc( left + j+1 )
       r(m-1,j) = r(m-1,j) + fac* rbc( right + j+1 )
    END DO
    
  END SUBROUTINE apply_bc

  !*****************************************************************!

  SUBROUTINE cjgr( x, b ) 

    !***************************************************************!
    !*   Conjugate gradient solver for A x = b                     *!
    !*   Uses the routine a_times(x) below                         *!
    !*   no preconditioner (yet)
    !***************************************************************!

    IMPLICIT NONE

    REAL(KIND(0.d0)), DIMENSION(nf), INTENT(IN)  :: b
    REAL(KIND(0.d0)), DIMENSION(nf), INTENT(OUT) :: x
    REAL(KIND(0.d0)), DIMENSION(SIZE(b)) :: r,d,q,s
    REAL(KIND(0.d0)) :: delta_new,delta_old,eps, alpha, beta
    INTEGER :: iter, kk

    iter = 0


    r = b - a_times(x)

    r = r - b_times(x)


    delta_new = SQRT( DOT_PRODUCT(r,r) ) ! for check
    IF (delta_new < cgtol*cgtol ) THEN
       !      WRITE(99,*) '  cg: no iter required '
    END IF

    d = r !m_inverse_times(r)
    delta_new = DOT_PRODUCT(r,d) ! same as above for check
    eps = cgtol*cgtol
    !   write(99,*) "  cg: initial residual", delta_new

    DO WHILE ((iter.lt.cg_max_iter).and.(delta_new.gt.eps))

       iter=iter+1
!       write(*,*) 'cg...iter, delta_new',iter,delta_new
       q = a_times( d )

        q = q + b_times( d )

       alpha = delta_new / DOT_PRODUCT( d, q )
       x = x + alpha * d

!        print *, "maxval new x = ", maxval(abs(x))

       IF (MOD(iter,50).eq.0) THEN
          r = b - a_times( x )

          r = r - b_times( x )


       ELSE
          r = r - alpha * q
       END IF

       s = r !m_inverse_times(r)
       delta_old = delta_new
       delta_new = DOT_PRODUCT(r,s)

       beta = delta_new/delta_old
       d = s + beta * d

    END DO



!    print *, "maxval new x = ", maxval(abs(x))


    IF (iter.eq.cg_max_iter) THEN
       WRITE(*,*)  "......WARNING, cg used max iterations"
       WRITE(*,*)    "......itmax = ",iter,", res = ",delta_new
    END IF
    WRITE(*,*) "......CG used ",iter," iterations"

  END SUBROUTINE cjgr

  !*****************************************************************!


!*****************************************************************!

!*****************************************************************!
subroutine bicgstab( x, b)

integer :: j, iter
real(kind(0.d0)), dimension(nf), intent(in) :: b
real(kind(0.d0)), dimension(nf), intent(inout) :: x
real(kind(0.d0)), dimension(nf) :: r, rhat, p, nu, h, sv, tv
real(kind(0.d0)) :: rho_o, rho_n, alpha, om, eps, err, bta

!initialize


err = 1.d0
eps = cgtol * cgtol
iter = 0

r = b - a_times( x) - b_times( x)

rhat = r

rho_o = 1.d0
alpha = 1.d0
om = 1.d0

nu = 0.d0
p = 0.d0

do while ((iter .le. cg_max_iter) .and. (err .ge. eps))

    rho_n = dot_product(rhat, r)

    bta = (rho_n/rho_o) * (alpha/om)

    !upate rho:
    rho_o = rho_n

    p = r + bta * (p - om * nu)

    nu = a_times(p) + b_times( p )

    alpha = rho_n/dot_product(rhat, nu)

    h = x + alpha * p

    sv = r - alpha * nu

    tv = a_times(sv ) + b_times( sv)

    om = dot_product( tv, sv)/ dot_product(tv, tv)

    x = h + om * sv

    r = sv - om * tv

    err = dot_product( r, r)

!    err = dot_product( b - a_times( x) - b_times( x ), b - a_times(x) - b_times(x ) )


    iter = iter + 1

end do

IF (iter.eq.cg_max_iter) THEN
    WRITE(*,*)  "......WARNING, bicg used max iterations"
    WRITE(*,*)    "......itmax = ",iter,", res = ", err
!Else
!    WRITE(*,*) "......biCG used ",iter," iterations"
end if

end subroutine
!*****************************************************************!


  FUNCTION a_times(x)

    !***************************************************************!
    !*   Cholesky factorization of EC(C^t C)^-1 C^t E^t            *!
    !*      performed once only for stationary bodies              *!
    !***************************************************************!
    USE myfft
    REAL(KIND(0.D0)), DIMENSION(Nf) :: x, a_times
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n,mgridlev) :: vort, s
    REAL(KIND(0.D0)), DIMENSION(Nq, mgridlev ) :: vel

    ! zero all levels of temporary vorticity
    vort = 0.d0
    ! regularize forces for level 1 vorticity
    vort(:,:,1) = ainv( rot( reg(x) ) )
    ! coarsify vorticity and find curl(laplace^-1)

    if (mgridlev .ge. 2) then
    	CALL vort2flux( vel, vort, s, mgridlev )
    else
        call vort2flux( vel, vort, s, mgridlev )
    end if 


   ! regularize

    a_times = regT( vel(:,1) )

  END FUNCTION a_times

  !*****************************************************************!

 !*****************************************************************!
function b_times(x)

 !Perform 2*delta/ dt * Khatinv * Q * W * x

REAL(KIND(0.d0)), DIMENSION(nf), INTENT(in) :: x
real(kind(0.d0)), dimension(nf) :: v_tp
real(kind(0.d0)), dimension(nb) :: v_x, v_y
real(kind(0.d0)), dimension(3*nb) :: v_bg
real(kind(0.d0)), dimension(nf) :: b_times
integer :: kk

v_tp = redistribute( x )

v_x = v_tp( 1 : nb )
v_y = v_tp( nb + 1 : nf )

v_bg = 0.d0
call build_force( v_x, v_y, v_bg)

v_bg = delta * matmul( sol_mat, v_bg ) / dt**2.d0

v_tp = v_trunc( v_bg )

b_times = v_tp

end function


 !*****************************************************************!


  FUNCTION nonlinear( omega, q, q0, bc ) RESULT(fq)

    !***************************************************************!
    !*   nonlinear terms in rotational form                        *!
    !***************************************************************!
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: omega
    REAL(KIND(0.D0)), DIMENSION(Nq) :: q, q0, fq
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: bc
    REAL(KIND(0.D0)), DIMENSION(1:m+1,2:n) :: uavg
    REAL(KIND(0.D0)), DIMENSION(2:m,1:n+1) :: vavg
    INTEGER :: i,j

    DO j=2,n
       DO i=1,m+1
          uavg(i,j) = 0.5D0*( q(u(i,j))+q(u(i,j-1)) + q0(u(i,j))+q0(u(i,j-1)))
       END DO
    END DO
    DO j=1,n+1
       DO i=2,m
          vavg(i,j) = 0.5D0*( q(v(i,j))+q(v(i-1,j)) + q0(v(i,j))+q0(v(i-1,j)))
       END DO
    END DO

    DO j=2,n-1
       DO i=2,m
           fq(u(i,j)) = 0.5D0*( vavg(i,j+1)*omega(i,j+1) + vavg(i,j)*omega(i,j))
       ENDDO
    ENDDO
    
    DO i=2,m
       j = 1
       fq(u(i,j)) = 0.5D0*( vavg(i,j+1)*omega(i,j+1) + vavg(i,j)*bc(bottom+i) )
       j = n
       fq(u(i,j)) = 0.5D0*( vavg(i,j+1)*bc(top+i) + vavg(i,j)*omega(i,j) )
    END DO
    ! note...we don't need result for i=1 or 1=m+1 since it isn't needed by rot

    DO j=2,n
       i = 1
       fq(v(i,j))    = -0.5D0*( uavg(i+1,j)*omega(i+1,j) + uavg(i,j)*bc(left+j) )
       DO i=2,m-1
          fq(v(i,j)) = -0.5D0*( uavg(i+1,j)*omega(i+1,j) + uavg(i,j)*omega(i,j) )
       END DO
       i = m
       fq(v(i,j))    = -0.5D0*( uavg(i+1,j)*bc(right+j)  + uavg(i,j)*omega(i,j) )
    END DO
    ! note...we don't need result for j=1 or j=n+1 since it isn't needed by rot

 END FUNCTION nonlinear

!*****************************************************************!


!*****************************************************************!
FUNCTION reg( h0 ) RESULT( h )

!***************************************************************!
!*   regularization of body force                              *!
!***************************************************************!
USE grid
IMPLICIT NONE
REAL(KIND(0.D0)), DIMENSION(Nf), INTENT(IN) :: h0
REAL(KIND(0.D0)), DIMENSION(Nq)             :: h
INTEGER ::i, j, k, l, p, next

! initalize regularization field
h = 0.D0

DO k=1,nb

    i=indexx(k)
    j=indexx(k+nb)
    next=0
    DO l=-support2,support2
        DO p=-support2,support2

            next=next+1
            h(u(i+p,j+l)) = h(u(i+p,j+l)) &
            + weight(next,k)*h0(k)

            h(v(i+p,j+l)) = h(v(i+p,j+l)) &
            + weight(next,k+nb)*h0(k+nb)
        END DO
    END DO
END DO

END FUNCTION reg
!*****************************************************************!



!*****************************************************************!
FUNCTION regT( h ) RESULT( h0 )

!***************************************************************!
!*  interpolation to body point  (Transpose of reg)            *!
!***************************************************************!
USE grid
IMPLICIT NONE
REAL(KIND(0.D0)), DIMENSION(Nq), INTENT(IN) :: h
REAL(KIND(0.D0)), DIMENSION(Nf)             :: h0
INTEGER ::i, j, k, l, p, next

h0 = 0.D0

DO k=1,nb
    i=indexx(k)
    j=indexx(k+nb)
    next=0
    DO l=-support2,support2
        DO p=-support2,support2
            next=next+1
            h0(k) = h0(k) + weight(next,k)*h(u(i+p,j+l))
            h0(k+nb) = h0(k+nb) + weight(next,k+nb)*h(v(i+p,j+l))
        END DO
    END DO
END DO

END FUNCTION regT


!*****************************************************************!




!-------------------------------------------------------------!
subroutine get_M( Mmat )

implicit none

integer :: j, count
real(kind(0.d0)), dimension(3*(nel + 1), 3*(nel + 1)), intent(inout) :: Mmat
real(kind(0.d0)), dimension(6, 6) :: M_e
integer, dimension(6) :: index

integer :: jjj, jjjj

Mmat = 0.d0
M_e = 0.d0

do j = 1 , nel

    M_e(1,1) = 140.d0
    M_e(1,4) = 70.d0

    M_e(2,2) = 156.d0
    M_e(2,3) = 22.d0*h0_el(j)
    M_e(2,5) = 54.d0
    M_e(2,6) = -13.d0*h0_el(j)

    M_e(3,2) = 22.d0*h0_el(j)
    M_e(3,3) = 4.d0*h0_el(j)**2
    M_e(3,5) = 13.d0*h0_el(j)
    M_e(3,6) = -3.d0*h0_el(j)**2

    M_e(4,1) = 70.d0
    M_e(4,4) = 140.d0

    M_e(5,2) = 54.d0
    M_e(5,3) = 13.d0*h0_el(j)
    M_e(5,5) = 156.d0
    M_e(5,6) = -22.d0*h0_el(j)

    M_e(6,2) = -13.d0*h0_el(j)
    M_e(6,3) = -3.d0*h0_el(j)**2
    M_e(6,5) = -22.d0*h0_el(j)
    M_e(6,6) = 4.d0*h0_el(j)**2



    M_e = R_rho*R_th*h0_el(j)/420.d0 * M_e



    call get_ind( j, index )

    call M_assemble( Mmat, M_e, index )



end do




!Apply BCs:
if (STANDARD ) then
    count = 1
else
    count = 3*nb - 2
end if


do j = 1 , 3

    if (bc_type(j) .eq. 1) then !If Dirichlet condition

        Mmat(count,:) = 0.d0
        Mmat(:,count) = 0.d0
!        Mmat(j,j) = 1.d0
        count = count + 1
    end if

end do


end subroutine
!-------------------------------------------------------------!

!-------------------------------------------------------------
subroutine get_ind( iel, index )

implicit none

integer, intent( in ) :: iel
integer, dimension(6), intent( out ) :: index
integer :: nnel, ndof, edof, start, ii

nnel = 2 ! # of nodes per element
ndof = 3 ! # of DOF per node

edof = nnel * ndof !degrees of freedom for an element

start = (iel - 1) * (nnel - 1) * ndof

index = 0
do ii = 1 , edof

    index(ii) = start + ii

end do

end subroutine

!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine M_assemble( Mmat, M_e, index )

implicit none

real(kind(0.d0)), dimension(3*(nel + 1), 3*(nel + 1)) :: Mmat
real(kind(0.d0)), dimension(6, 6), intent(in) :: M_e
integer, dimension(6), intent( in ) :: index
integer :: j, jj, r, rr

do r = 1 , 6

    rr = index( r )

    do j = 1 , 6

        jj = index(j)

        Mmat( rr, jj ) = Mmat( rr, jj ) + M_e( r , j )

    end do

end do

end subroutine
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine var_update( udisp, Fint )

implicit none

integer :: j
real(kind(0.d0)), dimension(nel +1) :: x, y
real(kind(0.d0)), dimension(3*(nel + 1)), intent(in) :: udisp
!real(kind(0.d0)), dimension(nel) :: c, s, h_el
!real(kind(0.d0)), dimension(3*nel) :: qL
real(kind(0.d0)), dimension(3*(nel+1)), intent(out) :: Fint
real(kind(0.d0)) :: dx, dy
integer, dimension(6) :: index
integer, dimension(3) :: qind
real(kind(0.d0)) :: uL, theta1L, theta2L, beta1, beta2
real(kind(0.d0)), dimension(3,6) :: Bmat
real(kind(0.d0)), dimension(6) :: qint

!Update geometry based on original coordinates (x,y) in global frame and
!generalized displacement vector u:


x = xb0( 1 : nb )
y = xb0( nb + 1 : nf )

c_b = 0.d0
s_b = 0.d0
h_el = 0.d0
qL = 0.d0
Fint = 0.d0

do j = 1 , nel

    call get_ind(j, index)

    dx = x(j+1) + udisp(index(4)) - ( x(j) + udisp(index(1)) )
    dy = y(j+1) + udisp(index(5)) - ( y(j) + udisp(index(2)) )

    h_el(j) = sqrt( dx**2 + dy**2 )

    c_b(j) = dx/h_el(j)
    s_b(j) = dy/h_el(j)

    uL = ( h_el(j)**2 - h0_el(j)**2 )/( h_el(j) + h0_el(j) )

    qind = (/3*(j-1)+1 , 3*(j-1) + 2, 3*(j-1)+3 /)
    qL(qind(1)) = R_E*R_th*uL/(h0_el(j))

    beta1 = udisp(index(3)) + beta0( j )
    beta2 = udisp(index(6)) + beta0( j )

    theta1L = atan2( c_b(j)*sin(beta1) - s_b(j)*cos(beta1), &
    c_b(j)*cos(beta1) + s_b(j)*sin(beta1) )

    theta2L = atan2( c_b(j)*sin(beta2) - s_b(j)*cos(beta2), &
    c_b(j)*cos(beta2) + s_b(j)*sin(beta2) )

    qL(qind(2)) = 2.d0*R_E*R_sh/h0_el(j) * (2.d0*theta1L + theta2L)
    qL(qind(3)) = 2.d0*R_E*R_sh/h0_el(j) * (theta1L + 2.d0*theta2L)

    Bmat = 0.d0

    Bmat(1,1) = -c_b(j)
    Bmat(1,2) = -s_b(j)
    Bmat(1,4) = c_b(j)
    Bmat(1,5) = s_b(j)

    Bmat(2,1) = -s_b(j)/h_el(j)
    Bmat(2,2) = c_b(j)/h_el(j)
    Bmat(2,3) = 1.d0
    Bmat(2,4) = s_b(j)/h_el(j)
    Bmat(2,5) = -c_b(j)/h_el(j)

    Bmat(3,1) = -s_b(j)/h_el(j)
    Bmat(3,2) = c_b(j)/h_el(j)
    Bmat(3,4) = s_b(j)/h_el(j)
    Bmat(3,5) = -c_b(j)/h_el(j)
    Bmat(3,6) = 1.d0

    qint = matmul(transpose(Bmat), qL(qind(1):qind(3) ) ) !Internal forces in global frame

    call F_assemble( Fint, qint, index )

end do


call apply_bcs( Fint )

end subroutine
!-------------------------------------------------------------


!-------------------------------------------------------------
subroutine F_assemble( Fvect, F_e, index )

implicit none

integer :: j, jj
real(kind(0.d0)), dimension(3*(nel + 1)), intent(inout) :: Fvect
real(kind(0.d0)), dimension(6) :: F_e !element forces
integer, dimension(6) :: index

do j = 1 , 6

jj = index(j)

Fvect( jj ) = Fvect( jj ) + F_e( j );

end do

end subroutine
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine build_force( fx, fy, Fvect )

implicit none

real(kind(0.d0)), dimension(nel +1), intent(in) :: fx, fy
real(kind(0.d0)), dimension(3*(nel + 1)), intent(out) :: Fvect
integer :: j
real(kind(0.d0)), dimension(2) :: fxsmall, fysmall
real(kind(0.d0)), dimension(6) :: F_e
integer, dimension(6) :: index
real(kind(0.d0)), dimension(6) :: fsmall

!Take a vector of nodal x and y forces and turn them into the global force
!vector for the FEM formulation...

Fvect = 0.d0

do j = 1 , nel

    !Build element force vector:
    fxsmall = fx(j : j + 1) !Get x forces at that element
    fysmall = fy(j : j + 1) !Get y forces at that element


    fsmall = (/ fxsmall(1), fysmall(1), 0.d0, &
        fxsmall(2), fysmall(2), 0.d0 /)



    call get_Fe( fsmall, h0_el(j), F_e )

    !Assemble force vector:
    call get_ind( j, index )
    call F_assemble( Fvect, F_e, index )

end do

call apply_bcs( Fvect )

end subroutine
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine get_Fe( fs, hspace, F_e)

implicit none

real(kind(0.d0)), dimension(2) :: fxs, fys
real(kind(0.d0)), intent(in) :: hspace
real(kind(0.d0)), dimension(6) :: F_e
real(kind(0.d0)), dimension(6) :: fs

fxs(1) = fs(1)
fxs(2) = fs(4)
fys(1) = fs(2)
fys(2) = fs(5)

!Get element-wise rhs vector F_e (6x1 vector)
F_e = 0.d0

!section corresponding to x forces:
F_e(1) = hspace/3.d0*fxs(1) + hspace/6.d0*fxs(2)
F_e(2) = hspace/6.d0*fxs(1) + hspace/3.d0*fxs(2)

!section corresponding to y forces:
F_e(3) = 26.d0*hspace/70.d0*fys(1) + 9.d0*hspace/70.d0*fys(2);
F_e(4) = -11.d0*hspace**2/210.d0*fys(1) - hspace**2*13.d0/420.d0*fys(2);
F_e(5) = 9.d0*hspace/70.d0*fys(1) + 26.d0*hspace/70.d0*fys(2);
F_e(6) = 13.d0*hspace**2/420.d0*fys(1) + hspace**2*11/210.d0*fys(2);

!Rearrange vector to match structure of desired solution vector:
call vect_rearrange( F_e )

end subroutine
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine vect_rearrange( vect )

implicit none

real(kind(0.d0)), dimension(6), intent(inout) :: vect
real(kind(0.d0)), dimension(6) :: v_store

!Rearrange element matrices to correspond to the following arrangement of
!unknowns and eqns:
!           [u1, w1, theta1, u2, w2, theta2]
!(the matrices are originally arranged as [u1, u2, w1, theta1, w2, theta2]

v_store = vect;

vect(2) = v_store(3);
vect(3) = v_store(4);
vect(4) = v_store(2);

end subroutine
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine get_K( K_mat )

implicit none

integer :: j, jj, jjj, count
real(kind(0.d0)), dimension(3*(nel+1), 3*(nel+1)), intent(out) :: K_mat
integer, dimension(6) :: index
real(kind(0.d0)) :: r_c, Nf, Mf1, Mf2
real(kind(0.d0)), dimension(3,3) :: CL
real(kind(0.d0)), dimension(3,6) :: Bmat
real(kind(0.d0)), dimension(6,6) :: K1, K2, opmat1, opmat2, opmat3, K_e
integer, dimension(3) :: qind
real(kind(0.d0)), dimension(6) :: rvect, zvect


!Build global stiffness matrix from current beam configuration...

K_mat = 0.d0
do j = 1 , nel

    r_c = R_sh/R_th
    CL = 0.d0
    CL(1,1) = 1.d0
    CL(2,2) = 4.d0*r_c
    CL(2,3) = 2.d0*r_c
    CL(3,2) = 2.d0*r_c
    CL(3,3) = 4.d0*r_c
    CL = R_E*R_th/h0_el(j)* CL

    Bmat = 0.d0

    Bmat(1,1) = -c_b(j)
    Bmat(1,2) = -s_b(j)
    Bmat(1,4) = c_b(j)
    Bmat(1,5) = s_b(j)

    Bmat(2,1) = -s_b(j)/h_el(j)
    Bmat(2,2) = c_b(j)/h_el(j)
    Bmat(2,3) = 1.d0
    Bmat(2,4) = s_b(j)/h_el(j)
    Bmat(2,5) = -c_b(j)/h_el(j)

    Bmat(3,1) = -s_b(j)/h_el(j)
    Bmat(3,2) = c_b(j)/h_el(j)
    Bmat(3,4) = s_b(j)/h_el(j)
    Bmat(3,5) = -c_b(j)/h_el(j)
    Bmat(3,6) = 1.d0

    K1 = matmul(transpose(Bmat), matmul( CL , Bmat) )

    qind = (/ 3*( j - 1 ) + 1, 3*(j-1) + 2, 3*(j-1) + 3 /)
    Nf = qL( qind( 1 ) );
    Mf1 = qL( qind( 2 ) );
    Mf2 = qL( qind( 3 ) );

    zvect = (/s_b(j), -c_b(j), 0.d0, -s_b(j), c_b(j), 0.d0/)
    rvect = -1.d0*(/ c_b(j),  s_b(j), 0.d0, -c_b(j), s_b(j), 0.d0 /)


    !Build outer products necessary for building K2:
    opmat1 = 0.d0
    opmat2 = 0.d0
    opmat3 = 0.d0
    do jj  = 1, 6
        do jjj = 1, 6
            opmat1(jj,jjj) = zvect(jj)*zvect(jjj)
            opmat2(jj,jjj) = rvect(jj)*zvect(jjj)
            opmat3(jj,jjj) = zvect(jj)*rvect(jjj)
        end do
    end do

    K2 = Nf/h_el(j)*opmat1 + (Mf1 + Mf2)/( h_el(j)**2 )&
    *( opmat2 + opmat3)

    K_e = K1 + K2

    !Assemble the stiffness matrix:
    call get_ind( j, index )
    call K_assemble( K_mat, K_e, index )

end do

!BC's for K matrix:
if (STANDARD) then
    count = 1
else
    count = 3*nb - 2
end if
do j = 1 , 3

    if ( bc_type(j) .eq. 1) then !If Dirichlet condition

        K_mat(count,:) = 0.d0
        K_mat(:,count) = 0.d0
        K_mat(count, count) = 1.d0

        count = count + 1

    end if

end do

end subroutine
!-------------------------------------------------------------


!-------------------------------------------------------------
subroutine K_assemble( K_mat, K_e, index )

implicit none

integer, dimension(6), intent(in) :: index
real(kind(0.d0)), dimension(3*(nel+1),3*(nel+1)), intent(inout) :: K_mat
real(kind(0.d0)), dimension(6,6), intent(in) :: K_e
integer :: j, jj, r, rr

do r = 1 , 6

    rr = index( r )

    do j = 1 , 6

        jj = index(j)

        K_mat( rr, jj ) = K_mat( rr, jj ) + K_e( r , j )

    end do

end do


end subroutine
!-------------------------------------------------------------


!-------------------------------------------------------------
subroutine apply_bcs(  Fvect )

implicit none

integer :: j, count
!real(kind(0.d0)), dimension(3*(nel+1),3*(nel+1)), intent(inout) :: Kmat
real(kind(0.d0)), dimension(3*(nel+1)), intent(inout) :: Fvect

! bc_type -- 3x1 vector containing info on the type of BCs at node 1
!            (1 = Dirichlet, 2 = Neumman)
! bc_val -- 3x1 vector containing the values of the BCs (only set to
!           nonzero if the corresponding BC type is Dirichlet)

if (STANDARD) then
    count = 1
else
    count = 3*nb - 2
end if
do j = 1 , 3

    if ( bc_type(j) .eq. 1) then !If Dirichlet condition

        Fvect(count) = bc_val(j)

        count = count + 1
    end if

end do


end subroutine
!-------------------------------------------------------------


  SUBROUTINE write_slip( it, x, xb )

    !***************************************************************!
    !*   write actual slip at body points to monitor error         *!
    !***************************************************************!
     USE parameters
    INTEGER          :: it
    REAL(KIND(0.D0)), DIMENSION(Nq) :: x
    REAL(KIND(0.D0)), DIMENSION(Nf) :: xb
    REAL(KIND(0.D0)) :: slip

    ! computing slip

    slip = MAXVAL( ABS (regT(x) - xb) )

    IF (it==1) THEN
       OPEN(unit=106,file="output/slip.dat",form="formatted",status="replace")
    ELSE
       OPEN(unit=106,file="output/slip.dat",form="formatted",status="old",position="append")
    END IF
    WRITE(106,*) it, slip
    CLOSE(106)

  END SUBROUTINE write_slip

  !*****************************************************************!

  SUBROUTINE write_udtdx( it, x, y )

    !***************************************************************!
    !*   write maximum u*dt/dx                                     *!
    !***************************************************************!
     USE parameters
    INTEGER          :: it, i, j
    REAL(KIND(0.D0)), DIMENSION(Nq) :: x, y, z
    REAL(KIND(0.D0)), DIMENSION(1:m,1:n) :: va
    REAL(KIND(0.D0)) :: para

    ! computing u*dt/dx 
    Do i = 1,Nq
       z(i)=( x(i)+y(i) )*m/len
    End Do

    DO j=1,n
       DO i=1,m
          va(i,j) = 0.5D0*( ( z(u(i,j))+z(u(i+1,j)) )**2 + & 
               ( z(v(i,j))+z(v(i,j+1)) )**2 )**0.5
       END DO
    END DO
 
    para = MAXVAL(va)*dt*m/len

    IF (it==1) THEN
       OPEN(unit=106,file="output/cfl.dat",form="formatted",status="replace")
    ELSE
       OPEN(unit=106,file="output/cfl.dat",form="formatted",status="old",position="append")
    END IF

    !velocity is normalized by (a0*D)^0.5
    WRITE(106,*) it, para
    CLOSE(106)

  END SUBROUTINE write_udtdx

!*****************************************************************!

  SUBROUTINE write_force( it, frc )

    !***************************************************************!
    !*   write body forces to file                                 *!
    !***************************************************************!
  
    USE parameters
    USE user
    INTEGER          :: it, i, j, k
    REAL(KIND(0.D0)), DIMENSION(nf) :: frc
    REAL(KIND(0.D0)), DIMENSION(n_body+n_actuator):: forcex, forcey, forcex_lab, forcey_lab


    DO i=1,n_body
       forcex(i) = 2.d0*delta*SUM( getbody(frc, 1, i, bdy(i)%npts) )
       forcey(i) = 2.d0*delta*SUM( getbody(frc, 2, i, bdy(i)%npts) ) 
       forcex_lab(i) = forcex(i)*cos(rot_angle) - forcey(i)*sin(rot_angle)
       forcey_lab(i) = forcex(i)*sin(rot_angle) + forcey(i)*cos(rot_angle)
    END DO

    DO i=1,n_actuator
       forcex(n_body+i) = 2.d0*delta*SUM( getact(frc, 1, i, act(i)%npts) )
       forcey(n_body+i) = 2.d0*delta*SUM( getact(frc, 2, i, act(i)%npts) )
       forcex_lab(n_body+i) = forcex(n_body+i)*cos(rot_angle) - forcey(n_body+i)*sin(rot_angle)
       forcey_lab(n_body+i) = forcex(n_body+i)*sin(rot_angle) + forcey(n_body+i)*cos(rot_angle)
    END DO

    IF (it==1) THEN
       OPEN(unit=100,file="output/force.dat",form="formatted",status="replace")
    ELSE
       OPEN(unit=100,file="output/force.dat",form="formatted",status="old",position="append")
    END IF
   ! WRITE(100,*) it, (forcex(k), k=1,n_body+n_actuator), (forcey(k), k=1,n_body+n_actuator),&
   !                  (forcex_lab(k), k=1,n_body+n_actuator), (forcey_lab(k), k=1,n_body+n_actuator)
    WRITE(100,*) it, (forcex(k), k=1,n_body+n_actuator), (forcey(k), k=1,n_body+n_actuator)
    CLOSE(100)
 
  END SUBROUTINE write_force

  !*****************************************************************!

!*****************************************************************!

SUBROUTINE write_total_force( it, frc )

!***************************************************************!
!*   write body forces to file                                 *!
!***************************************************************!

use grid
USE parameters
INTEGER          :: it, j, k
REAL(KIND(0.D0)), DIMENSION(nf) :: frc
REAL(KIND(0.D0)), DIMENSION(nb):: forcex, forcey
real(kind(0.d0)), dimension(nb) :: yhat, xhat
CHARACTER(7) :: force_file_num

forcex = ( getbody(frc, 1, 1, bdy(1)%npts) )
forcey = ( getbody(frc, 2, 1, bdy(1)%npts) )

xhat  = getbody(xb, 1, 1, bdy(1)%npts)
yhat = getbody(xb, 2, 1, bdy(1)%npts)


WRITE(force_file_num,"(I7.7)") it

OPEN(unit=101,file="output/total_force_"//force_file_num//".dat",form="formatted",status="replace")

DO k = 1, nb
WRITE(101,*) xhat(k), yhat(k), vb(k), vb(k + nb), forcex(k), forcey(k)
END DO

CLOSE(101)

END SUBROUTINE write_total_force

!*****************************************************************!


!*****************************************************************!
subroutine write_tip_disp(it, disp)



USE parameters
USE user
INTEGER          :: it
REAL(KIND(0.D0)), dimension(2) :: disp


IF (it==1) THEN
OPEN(unit=10001,file="output/disp.dat",form="formatted",status="replace")
ELSE
OPEN(unit=10001,file="output/disp.dat",form="formatted",status="old",position="append")
END IF
! WRITE(100,*) it, (forcex(k), k=1,n_body+n_actuator), (forcey(k), k=1,n_body+n_actuator),&
!                  (forcex_lab(k), k=1,n_body+n_actuator), (forcey_lab(k), k=1,n_body+n_actuator)
WRITE(10001,*) it, disp(1), disp(2)
CLOSE(10001)


end subroutine


!*****************************************************************!



  SUBROUTINE choldc

    !***************************************************************!
    !*   Cholesky factorization of A                               *!
    !***************************************************************!
    USE variables
    REAL(KIND(0.D0)) :: sum

    INTEGER :: i,j,k

    DO i=1,Nf
     ! WRITE(*,"(A,I5)",ADVANCE="NO") '...',i
       DO j=i,Nf
          sum=cholmat(i,j)
          DO k=i-1,1,-1
             sum=sum-cholmat(i,k)*cholmat(j,k)
          END DO
          IF(i.EQ.j)THEN
            IF(sum.LE.0.) STOP 'choldc failed'
            cholvec(i)=SQRT(sum)
          ELSE
            cholmat(j,i)=sum/cholvec(i)
          ENDIF
       END DO
    END DO

  END SUBROUTINE choldc

  !*****************************************************************!

  FUNCTION cholsl(b) RESULT(x)

    !***************************************************************!
    !*   Solve A x = b given it's Cholesky decomposition           *!
    !***************************************************************!
    USE variables
    REAL(KIND(0.D0)), DIMENSION(:) :: b
    REAL(KIND(0.D0)), DIMENSION(SIZE(b)) :: x

    INTEGER :: i,k
    REAL(KIND(0.D0)) ::  sum

    DO i=1,Nf
       sum=b(i)
       DO  k=i-1,1,-1
          sum=sum-cholmat(i,k)*x(k)
       END DO
       x(i)=sum/cholvec(i)
    END DO
    DO i=Nf,1,-1
       sum=x(i)
       DO k=i+1,Nf
          sum=sum-cholmat(k,i)*x(k)
       END DO
       x(i)=sum/cholvec(i)
    END DO

  END FUNCTION cholsl

 !*****************************************************************!

!  SUBROUTINE accel_bodies(it,vlocal,flocal,ab)
!
!    USE grid
!    USE user
!
!    INTEGER :: it,ibdy,istrt,iend
!    REAL(KIND(0.D0)), DIMENSION(nf) :: x, vlocal,flocal, ab
!    REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: ax,ay
!
!    ab = 0.D0 ! default is no acceleration
!
!    ! this routine returns the accelration for each body on the domain
!
!    DO ibdy=1,n_body
!       IF (bdy(ibdy)%moving) THEN
!          istrt = bdy(ibdy)%ipos
!          iend = bdy(ibdy)%ipos+bdy(ibdy)%npts - 1
!          ALLOCATE( ax(iend-istrt+1), ay(iend-istrt+1) )
!          CALL body_accel( ibdy, it, xb(f(istrt:iend,1)), xb(f(istrt:iend,2)), &
!                                     vlocal(f(istrt:iend,1)), vlocal(f(istrt:iend,2)), &
 !                                    flocal(f(istrt:iend,1)), flocal(f(istrt:iend,2)), &
!                                     ax, ay)
!          ab(f(istrt:iend,1)) = ax
!          ab(f(istrt:iend,2)) = ay
!          DEALLOCATE( ax,ay )
!       END IF
!    END DO
!
!  END SUBROUTINE accel_bodies

 !*****************************************************************!

  SUBROUTINE actuators(it,xold,xlocal,vlocal)

    USE user
    INTEGER :: it, iact,istrt,iend,islv,i
    REAL(KIND(0.D0)), DIMENSION(nf) :: xold, xlocal,vlocal, atemp
    REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: vx,vy
    REAL(KIND(0.D0)) :: xmove, ymove

    DO iact=1,n_actuator
       
       istrt = act(iact)%ipos
       iend = act(iact)%ipos+act(iact)%npts - 1
       ALLOCATE( vx(iend-istrt+1), vy(iend-istrt+1) )
       CALL actuator_vel( iact, it, xlocal(f(istrt:iend,1)), xlocal(f(istrt:iend,2)), &
                               vx,vy)
       vlocal(f(istrt:iend,1)) = vx
       vlocal(f(istrt:iend,2)) = vy
       DEALLOCATE( vx,vy )

       ! that was the easy part...now fix up location if actuator is slaved to a body
       ! we move the actuator to a new position found by moving the actuator by the
       ! same relative amount as the (point on the) body it is slaved to.
       IF (act(iact)%slaved) THEN
          islv = bdy(act(iact)%slavebody)%ipos + act(iact)%slavept - 1
          xmove = xlocal(f(islv,1)) - xold(f(islv,1))
          ymove = xlocal(f(islv,2)) - xold(f(islv,2))
          istrt = act(iact)%ipos
          iend = act(iact)%ipos+act(iact)%npts - 1
          DO i=istrt,iend
             xlocal(f(i,1)) = xold(f(i,1)) + xmove
             xlocal(f(i,2)) = xold(f(i,2)) + ymove
          END DO
       END IF
  
    END DO   

  END SUBROUTINE actuators
 !*****************************************************************!

  FUNCTION rhs_forcing( it, q ) RESULT(dq)

    USE user
 
    INTEGER :: it,i,j
    REAL(KIND(0.D0)), DIMENSION(Nq) :: q, dq
    REAL(KIND(0.D0)), DIMENSION(m,n) :: dqx,dqy
    REAL(KIND(0.D0)) :: xx,yy,uu,vv

    dqx = 0.d0
    dqy = 0.d0
    dq = 0.d0

    DO j=2,n-1
       DO i=2,m-1
          xx = delta*(REAL(i)-0.5d0) - offsetx
          yy = delta*(REAL(j)-0.5d0) - offsety
          uu = 0.5d0 * ( q(u(i,j)) + q(u(i+1,j)) ) / delta
          vv = 0.5d0 * ( q(v(i,j)) + q(v(i,j+1)) ) / delta
          dqx(i,j) = delta* bodyforcex(it,xx,yy,uu,vv)
          dqy(i,j) = delta* bodyforcey(it,xx,yy,uu,vv)
       END DO
    END DO

    DO j=2,n-1
       DO i=2,m
          dq(u(i,j)) = 0.5d0 * ( dqx(i-1,j)+dqx(i,j) )
       END DO
    END DO
    DO i=2,m-1
       DO j=2,n
          dq(v(i,j)) = 0.5d0 * ( dqy(i,j-1)+dqy(i,j) )
       END DO
    END DO

  END FUNCTION rhs_forcing

 !*****************************************************************!

FUNCTION delta_angle( it ) RESULT(ang)

    USE user
 
    INTEGER :: it,i,j
    REAL(KIND(0.D0)) :: ang
    REAL(KIND(0.D0)) :: k1, k2
    REAL(KIND(0.D0)), DIMENSION(5) :: uv
    
    ! calculate_angle calculates the rotating angle of the grid when the motion of the gris is specified

    uv = motion_grid(it)    
    k1 = dt*uv(3)
    uv = motion_grid(it+1)  
    k2 = dt*uv(3)
    ang = 0.5d0*(k1+k2)

  END FUNCTION delta_angle

 !*****************************************************************!
END MODULE operators

