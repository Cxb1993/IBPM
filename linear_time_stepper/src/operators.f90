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
  ! - now solves flows with arbitrary bodyforcey motion in body-fixed frame
  !
  !
  ! ibfs v.2.2 modified by Thibault Flinois
  ! February 2013
  ! Now includes forward linearized simulation, based on Won Tae Joe's implementation
  ! Now includes adjoint equations optimization, also based on Won Tae Joe's implementation
  !
  

  USE parameters
  USE grid
  USE optimization
  IMPLICIT NONE

CONTAINS

  SUBROUTINE preprocess

    USE variables

    !***************************************************************!
    !*   Cholesky factorization of EC(C^t C)^-1 C^t E^t            *!
    !*      performed once only for stationary bodies              *!
    !***************************************************************!

    REAL(KIND(0.D0)), DIMENSION(Nf) :: z   ! temporary immersed body force vector
    INTEGER :: i
    LOGICAL :: readchol
    real(kind(0.d0)), dimension(nf,nf) :: EHtild
    real(kind(0.d0)), dimension(nf,3*nb) :: eyetil
    real(kind(0.d0)), dimension(3*nb,nf) :: Btil
    integer :: info, neqns, nrhs, lda, ldb, j, lwork, jj
    real(kind(0.d0)), dimension(:), allocatable :: work
    integer, dimension(3*nb + nf) :: ipiv_bg
    integer, dimension(nf) :: ipiv_fl
    real(kind(0.d0)), dimension(nb) :: z1, z2
    real(kind(0.d0)), dimension(3*nb) :: fvect
    real(kind(0.d0)), dimension( nf,  nf) :: waste1, waste2
    LOGICAL :: readbmat, readJmat

       WRITE(*,*) 'precomputing body matrix for stationary geometry...please be patient'


    INQUIRE(file='output/bmat.chd',exist=readbmat)  ! already done?
    IF (readbmat) THEN
       CALL read_bmat
    ELSE


!    print *, "starting to build bigmat blocks "

    !----------------
    !build matrices for fluid part
        DO i=1,nf   ! build matrix one column at a time
          z = 0.d0
          z(i) = 1.d0
          cholmat(1:Nf,i) = a_times( z )

          z = 0.d0
          z(i) = 1.d0
          EHtild(1:nf, i) = redistribute( z )

        END DO


    !-------------------


    !--------------------
    !Build remaining matrices


        call build_solid_mats( Sol_mat, M_sol, K_sol )



        WRITE(*,*) 'precomputing and storing Khat matrix inverse'
        !-------------------
        ! get matrix inverse

        !Have to solve a linear system
        info = 0

        lwork = (3*nb )** 2
        allocate( work( lwork ) )

        neqns = 3*nb
        lda = 3*nb
        ipiv_bg = 0

        call dgetrf( neqns, lda, sol_mat, lda, ipiv_bg, info)
        call dgetri(neqns, sol_mat, lda, ipiv_bg, work, lwork, info)

        !-------------------

        WRITE(*,*) 'done with Khat matrix inverse'


        eyetil = 0.d0
        do j = 1, 3*nb

            fvect = 0.d0
            fvect(j) = 1.d0

            call sol2fl_ind( fvect, z )

            eyetil( 1 : nf, j) = z !* delta * (2.d0/dt)

        end do


        !Build matrix that projects f onto FEM basis...
        Btil = 0.d0
        do j = 1 , nf
            z = 0.d0
            z(j) = 1.d0

            z1 = z(1 : nb)
            z2 = z(nb + 1 : nf)

            fvect = 0.d0

            call build_force( z1, z2, fvect )

            Btil(1:3*nb, j) = fvect
        end do

        Btil = matmul( Btil, EHtild )

        !--------------------

        WRITE(*,*) 'precomputing and storing main matrix inverse'
        !-------------------
        ! get matrix inverse


        cholmat = cholmat + delta/(dt**2.d0) * matmul( &
            eyetil , matmul(sol_mat, Btil ) )

        !Have to solve a linear system
        info = 0

        lwork = (nf )** 2
        deallocate( work )
        allocate( work( lwork ) )



        neqns = nf
        lda = nf
        ipiv_fl = 0

        call dgetrf( neqns, lda, cholmat, lda, ipiv_fl, info)
        call dgetri(neqns, cholmat, lda, ipiv_fl, work, lwork, info)

        deallocate( work )

        !-------------------

        WRITE(*,*) 'done with main matrix inverse'

        CALL write_bmat  ! save results

    end if

    INQUIRE(file='output/Jmat.chd',exist=readJmat)  ! already done?
    IF (readJmat) THEN
        CALL read_Jmat

    else

        do j = 1, nf

            z = 0.d0
            z(j ) = 1.d0

            JHmat(:, j) = JHtimes( z )
            JWmat(:, j) = JWtimes( z )
            JEmat(:, j) = JEtimes( z )

        end do

        call write_Jmat


    end if




  END SUBROUTINE preprocess

  !*****************************************************************!

  SUBROUTINE advance(itime)

    !***************************************************************!
    !*    Advance to next time level                               *!
    !***************************************************************!

    USE myfft
    USE variables
    USE user

    INTEGER :: itime,k
    ! intermediate results
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: rhs, rhsbc
    REAL(KIND(0.D0)), DIMENSION(Nf) :: accel, vbhalf, xbhalf, xbnew, rhsf, fbhalf, motion
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: omega_bc
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1),mgridlev) :: lastbc
    REAL(KIND(0.D0)), DIMENSION(nq) :: nl_temp
    
    itime = itime + 1
    IF ( MOD(itime,10).eq.0 .or. itime.eq.istart+1) THEN
       WRITE(*,*) "...Advancing to itime =",itime
    END IF

    q0p = motion_potential(itime-1)
    q0r = motion_rotation(itime-1)
    q0 = q0p + q0r
    rot_angle = rot_angle + delta_angle(itime)
    
    IF (itime==1.and.stationary) THEN
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
       CALL accel_bodies(itime,vb,fb/dt,accel)
       vbhalf = vb + dt*accel
       xbhalf = xb + dt*vb
       CALL actuators(itime,xb,xbhalf,vbhalf)
       CALL setup_reg(xbhalf)
       motion = delta*vbhalf - regT(q0(:,1))
       CALL vort2flux( q(:,1), omega, s, 1 )
       rhsf = regT(q(:,1) ) - motion
       fbhalf = fb ! initial guess for first iteration in CJGR
       CALL cjgr(fbhalf,rhsf)
       fbhalf = 0.5d0*(fbhalf+fb)
       CALL accel_bodies(itime,vbhalf,fbhalf/dt,accel)
       xbnew = xb + 0.5d0*dt*(vb+vbhalf)
       vb = vb + dt*accel
       CALL actuators(itime,xb,xbnew,vb)
       xb = xbnew
       CALL setup_reg(xb)
       motion = delta*vb - regT(q0(:,1))
       rhsf = regT(q(:,1) ) - motion
       CALL cjgr(fb,rhsf)
    ELSE
       CALL actuators(itime,xb,xb,vb)
       motion = delta*vb - regT(q0(:,1))
       CALL vort2flux( q(:,1), omega, s, 1)
       rhsf = regT(q(:,1))-motion
       fb = cholsl(rhsf)
    END IF

    !Redistribute the force to make it more physically pleasing
    f_rdst = redistribute(fb)


    ! complete omega on grid 1
    omega(:,:,1) = omega(:,:,1) - ainv( rot( reg( fb ) ) )
    
    ! coarsify final omega and correct velocities on all grids
    CALL vort2flux( q, omega, s, mgridlev )

    CALL write_slip( itime, q, motion )
    CALL write_udtdx( itime, q, q0 )



  END SUBROUTINE advance
  !*****************************************************************!
  SUBROUTINE advance_lin(itime)

    !***************************************************************!
    !*    Advance to next time level for linearized code           *!
    !***************************************************************!

    USE myfft
    USE variables

    INTEGER :: itime, k, j

    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: rhs, rhsbc
    REAL(KIND(0.D0)), DIMENSION(Nf) :: accel, vbhalf, xbhalf, xbnew, rhsf, fbhalf, motion
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: omega_bc
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1),mgridlev) :: lastbc, lastbc_base
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1),mgridlev) :: sbc
    REAL(KIND(0.D0)), DIMENSION(nq) :: nl_temp
    integer :: kk
    real(kind(0.d0)), dimension(nf + 3*nb) :: bigvect, rhsbig
    real(kind(0.d0)), dimension(3*nb) :: u_old, ud_old, udd_old, rchi, rzeta
    real(kind(0.d0)), dimension(nf) :: rchi_sm, rzeta_sm
    real(kind(0.d0)), dimension(3*nb, 3*nb) :: Mtmp
    integer :: info, neqns, lda, lwork
    real(kind(0.d0)), dimension(:), allocatable :: work
    integer, dimension(3*nb) :: ipiv_bg



    ! Read in base flow variables for the first time step
    ! Note: Perturbation actual variables (=perturbations to base state) initialized to 0 in initial condition
    IF(itime==istart) THEN
        CALL update_base(itime)


        call setup_reg( xb_base )

        CALL preprocess

        !Overwrite q...
        open(unit=12341, file="output/q_init.var",convert='BIG_ENDIAN',&
            form="unformatted", access='STREAM', status="old")
        do j = 1, nq
            read(12341) q(j,:)
        end do
        close(12341)

        !Update omega to be consistent with input q:
        do k   = 1, mgridlev
            omega(:, :, k  ) = rot( q(:, k  ) )
        end do


        ! coarsify final omega and correct velocities on all grids
        CALL vort2flux( q, omega, s, mgridlev )


        !Overwrite uib...
        open(unit=12342, file="output/uib_init.var",convert='BIG_ENDIAN',&
        form="unformatted", access='STREAM', status="old")
        do j = 1, 3*nb
            read(12342) u_ib( j )
        end do
        close(12342)


        !Overwrite udib...
        open(unit=12343, file="output/udib_init.var",convert='BIG_ENDIAN',&
        form="unformatted", access='STREAM', status="old")
        do j = 1, 3*nb
            read(12343) ud_ib(j)
        end do
        close(12343)


!        !Overwrite fb...
!        open(unit=12342, file="output/fb_init.var",convert='BIG_ENDIAN',&
!        form="unformatted", access='STREAM', status="old")
!        do j = 1, 2*nb
!            read(12342) fb( j )
!        end do
!        close(12342)


        !Overwrite udib...
        open(unit=12343, file="output/uddib_init.var",convert='BIG_ENDIAN',&
        form="unformatted", access='STREAM', status="old")
        do j = 1, 3*nb
            read(12343) udd_ib(j)
        end do
        close(12343)


        !Extract perturbation pos, vel, and accel
        do kk = 1, nb

            xbp(kk) = u_ib(3*(kk-1) + 1)
            xbp(nb + kk) = u_ib(3*(kk-1) + 2)

            vb(kk) = ud_ib(3*(kk-1) + 1)
            vb(nb + kk) = ud_ib(3*(kk-1) + 2)

            ab(kk) = udd_ib(3*(kk-1) + 1)
            ab(nb + kk) = udd_ib(3*(kk-1) + 2)

        end do

    ENDIF





    ! Increment time
    itime = itime + 1
    IF ( MOD(itime,10).eq.0 .or. itime.eq.istart+1) THEN
       WRITE(*,*) "...Advancing to itime =",itime
    END IF
    
    ! Update q0 imposed perturbation
    q0p = motion_potential(itime-1)
    q0r = motion_rotation(itime-1)
    q0 = q0p + q0r
    rot_angle = rot_angle + delta_angle(itime)


    ! Step 1: solve intermediate curl(momentum) eq. on each grid

    
    ! Update old BC for AB2 for base variables
     DO k=1,mgridlev-1
       CALL get_bc(omega_base(:,:,k+1), lastbc_base(:,k), 0.25d0)
     END DO
     lastbc_base(:,mgridlev) = 0.d0
    
    ! Update old BC for AB2 for perturbation variables
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
       nl_temp = nonlinear(omega_base(:,:,k), q(:,k),      q0,           lastbc_base(:,k)) + &    ! omega_base x q  +  
                 nonlinear(omega(:,:,k),      q_base(:,k), q0_base(:,k), lastbc(:,k))             ! omega x q_base 
          
       ! add user-defined body force to momentum eq.
       IF (k==1) THEN
           nl_temp = nl_temp + rhs_forcing( itime, q(:,k) )                        
        END IF
       
       ! Take curl of convective term
       rhs = rot( nl_temp )

        !Add RHS contribution from base surface stress:
        if (k .eq. 1) then

            rhs = rhs - rot( matmul( JHmat, xbp ) ) / (dt/delta**2.d0)
            !The dt/delta**2 is to have consistent scaling with the time-stepper

        end if



    
       ! If this is the very first time step, we need to use explicit Euler
       ! Otherwise rhs_old computed at previous time step
       IF (itime==istart+1) THEN
          rhs_old(:,:,k) = rhs(:,:)
       END IF

       rhsbc = 0.d0  ! BC to be added on as they are the boundary terms, so far considered to be zero

       CALL apply_bc( rhsbc, lastbc(:,k), vfac(k) ) ! add bc from previous time step to rhsbc
       CALL apply_bc( rhsbc, omega_bc,    vfac(k) ) ! add bc from current time step to rhsbc
       

        omega(:,:,k) = dst( lam1i(:,:,k) * &
        ( dst( con1(k)*rhs(:,:)       + &
        con2(k)*rhs_old(:,:,k) + &
        rhsbc(:,:)             ) + &
        lam1(:,:,k)*dst( omega(:,:,k)) ) )

       ! update for next time step
       rhs_old(:,:,k) = rhs(:,:)


    END DO
       IF(.NOT.STEADY) CALL update_base(itime)  ! Update value of non linear simulation variables
    ! we now have vorticity evolution on all domains.  For smallest domain
    ! we add the vorticity created by the body


    ! first see if there is boundary motion and update regularization kernel if nec
    IF (.not.stationary) THEN
       WRITE (*,*) '************************  WARNING   **************************************'
       WRITE (*,*) 'MOVING BODY RELATIVE TO GRID: FEATURE IS NOT YET SUPPORTED FOR LINEAR SIMULATION'
       WRITE (*,*) '**************************************************************************'

       CALL accel_bodies(itime,vb,fb/dt,accel)
       vbhalf = vb + dt*accel
       xbhalf = xb + dt*vb
       CALL actuators(itime,xb,xbhalf,vbhalf)
       CALL setup_reg(xbhalf)
       motion = delta*vbhalf - regT(q0(:,1))
       CALL vort2flux( q(:,1), omega, s, 1 )
       rhsf = regT(q(:,1) ) - motion
       fbhalf = fb ! initial guess for first iteration in CJGR
       CALL cjgr(fbhalf,rhsf)
       fbhalf = 0.5d0*(fbhalf+fb)
       CALL accel_bodies(itime,vbhalf,fbhalf/dt,accel)
       xbnew = xb + 0.5d0*dt*(vb+vbhalf)
       vb = vb + dt*accel
       CALL actuators(itime,xb,xbnew,vb)
       xb = xbnew
       CALL setup_reg(xb)
       motion = delta*vb - regT(q0(:,1))
       rhsf = regT(q(:,1) ) - motion
       CALL cjgr(fb,rhsf)
    ELSE


        !Store pos, vel, and accel from previous time step
        u_old = u_ib
        ud_old = ud_ib


        !--RHS for flow equations

            CALL actuators(itime,xb,xb,vb)
            CALL vort2flux( q(:,1), omega, s, 1)

            rhsf = regT(q(:,1)) + regT(q0(:,1))


            !contribution from solid:
            rchi = ud_old + 2.d0/dt * u_old

            rchi_sm = 0.d0
            call sol2fl_ind( rchi, rchi_sm )


            rchi_sm = delta * (rchi_sm ) + matmul( JEmat, xbp)


            rzeta = matmul( M_sol, udd_ib ) + &
                matmul( M_sol, 4.d0/dt * ud_ib ) + &
                matmul( M_sol, 4.d0/(dt**2.d0) * u_ib )

            rzeta = rzeta + matmul( JWmat, xbp )

            rzeta = matmul( sol_mat, rzeta )

            rzeta_sm = 0.d0
            call sol2fl_ind( rzeta, rzeta_sm)


            rhsf = rhsf + rchi_sm - 2.d0*delta/ dt * rzeta_sm





        !-----


        !-----Initialize body acceleration if needed:
!            if (itime .eq. 1 ) then
!
!                Mtmp = 0.d0
!                call get_M0( Mtmp )
!
!                print *, "computing consistent initial body acceleration..."
!
!                info = 0
!                lwork = (3*nb)** 2
!                allocate( work( lwork ) )
!                neqns = 3*nb
!                lda = 3*nb
!                ipiv_bg = 0
!
!                call dgetrf( neqns, lda, Mtmp, lda, ipiv_bg, info)
!                call dgetri( neqns, Mtmp, lda, ipiv_bg, work, lwork, info)
!                deallocate( work )
!
!                udd_ib = matmul( Mtmp, -matmul(K_sol, u_ib) + QWx( fb ) + matmul( JWmat, xbp ) )
!
!                print *, "initial acceleration computed"
!
!            end if


            udd_old = udd_ib

        !-----


        !--- soln to linear system

            fb = matmul( cholmat, rhsf )
            !Get redistributed force:
            f_rdst = redistribute( fb )


            !Get body position
            u_ib = rzeta + matmul( sol_mat, QWx( fb ) )


            !Extract vel and accel from pos
            udd_ib = 4.d0/(dt**2)*(u_ib - u_old) - 4.d0/dt * ud_old - udd_old
            ud_ib = -ud_old + 2.d0/dt * (u_ib - u_old )

        !---

        !Extract perturbation pos, vel, and accel
        do kk = 1, nb

            xbp(kk) = u_ib(3*(kk-1) + 1)
            xbp(nb + kk) = u_ib(3*(kk-1) + 2)

            vb(kk) = vb0(kk) + ud_ib(3*(kk-1) + 1)
            vb(nb + kk) = vb0(nb +kk) + ud_ib(3*(kk-1) + 2)

            ab(kk) =  ab0(kk) + udd_ib(3*(kk-1) + 1)
            ab(nb + kk) =  ab0(nb +kk) + udd_ib(3*(kk-1) + 2)

        end do


    END IF


   ! complete omega on grid 1
       omega(:,:,1) = omega(:,:,1) - ainv( rot( reg( fb )) )

    ! coarsify final omega and correct velocities on all grids
    CALL vort2flux( q, omega, s, mgridlev )
    
    CALL write_slip( itime, q, vb )
    CALL write_udtdx(itime, q, q0 )

    if (standard ) then
        call write_tip_disp(itime, u_ib(3*nb - 2 : 3*nb - 1))
    else
        call write_tip_disp(itime, u_ib(1 : 2))

    end if

  END SUBROUTINE advance_lin


    !*****************************************************************!
  SUBROUTINE advance_adj(itime)

    !***************************************************************!
    !*    Advance to next time level for adjoint code              *!
    !***************************************************************!

    USE myfft
    USE variables

    INTEGER :: itime, i,j,k

    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: rhs, rhsbc
    REAL(KIND(0.D0)), DIMENSION(Nf) :: accel, vbhalf, xbhalf, xbnew, rhsf, fbhalf, motion
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: omega_bc
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1),mgridlev) :: lastbc, lastbc_base
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n,mgridlev) :: rhs2, temp
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1),mgridlev) :: sbc, tempbc
    REAL(KIND(0.D0)), DIMENSION(nq) :: nl_temp
        
    ! Increment time BACKWARDS
    itime = itime - 1
    IF ( MOD(itime,10).eq.0 .or. itime.eq.istop) THEN
       WRITE(*,*) "...Advancing BACKWARDS to itime =",itime
    END IF
    
    ! Read in base flow variables for the first time step
    IF(itime==istop) THEN
      CALL vort2flux( q, omega, s, mgridlev ) ! Get q given initial condition
      CALL update_base(itime) ! get the first base flow field
    ENDIF

    
    ! The adjoint does not have a superimposed grid motion
    q0p = 0.d0
    q0r = 0.d0
    q0 = 0.d0
    rot_angle = 0.d0
    
    IF (itime==istop.and.stationary) THEN
       CALL preprocess
    END IF
    
    
    ! Step 1: solve intermediate curl(momentum) eq. on each grid

    
    ! Update old BC for AB2 for base variables
     DO k=1,mgridlev-1
       CALL get_bc(omega_base(:,:,k+1), lastbc_base(:,k), 0.25d0)
     END DO
     lastbc_base(:,mgridlev) = 0.d0
    
    ! Update old BC for AB2 for adjoint variables
    DO k=1,mgridlev-1
       CALL get_bc(omega(:,:,k+1), lastbc(:,k), 0.25d0)
    END DO
    lastbc(:,mgridlev) = 0.d0

     CALL vort2flux(q,omega,s,mgridlev)

     DO k=1,mgridlev
       temp(:,:,k) = qcrossq( q(:,k), q_base(:,k), q0_base(:,k) )
       tempbc(:,k) = qcrossq_bc( q(:,k), q_base(:,k), q0_base(:,k) )
     END DO

     CALL stream2vort(temp,rhs2)
     CALL apply_bc(rhs2(:,:,mgridlev),tempbc(:,mgridlev),-1.d0)

    DO k = mgridlev,1,-1

       IF (k.lt.mgridlev) THEN
          CALL get_bc( omega(:,:,k+1), omega_bc, 0.25d0 )
       ELSE
          omega_bc = 0.d0
       END IF

       ! Computing the nonliner term (nonlinear term mult. by dt/2)
       rhs = rot( -nonlinear(omega_base(:,:,k), q(:,k), q0(:,k), lastbc_base(:,k) ) ) &
             + rhs2(:,:,k)

                        
       ! add adjoint rhs forcing term to momentum eq.
       DO i=mmin+1,mmax
         DO j=nmin+1,nmax
           rhs(i,j) = rhs(i,j) + dt*CGAM*omega_base(i,j,k)
         ENDDO
       ENDDO
     
       ! If this is the very first time step, we need to use explicit Euler
       ! Otherwise rhs_old computed at previous time step
       IF (itime==istop) THEN
          rhs_old(:,:,k) = rhs(:,:)
       END IF

       rhsbc = 0.d0  ! BC to be added on as they are the boundary terms, so far considered to be zero

       CALL apply_bc( rhsbc, lastbc(:,k), vfac(k) ) ! add bc from previous time step to rhsbc
       CALL apply_bc( rhsbc, omega_bc,    vfac(k) ) ! add bc from current time step to rhsbc
       
       omega(:,:,k) = dst( lam1i(:,:,k) * &
            ( dst( con1(k)*rhs(:,:)       + &
            con2(k)*rhs_old(:,:,k) + &
            rhsbc(:,:)             ) + &
            lam1(:,:,k)*dst( omega(:,:,k)) ) )

       ! update for next time step
       rhs_old(:,:,k) = rhs(:,:)
    END DO
      CALL update_base(itime)  ! Update value of non linear simulation variables
      
    ! we now have vorticity evolution on all domains.  For smallest domain
    ! we add the vorticity created by the body


    ! first see if there is boundary motion and update regularization kernel if nec.
       
    IF (.not.stationary) THEN
       WRITE (*,*) '************************  WARNING   **************************************'
       WRITE (*,*) 'MOVING BODY RELATIVE TO GRID: FEATURE IS NOT YET SUPPORTED FOR ADJOINT OPT'
       WRITE (*,*) '**************************************************************************'
       CALL accel_bodies(itime,vb,fb/dt,accel)
       vbhalf = vb + dt*accel
       xbhalf = xb + dt*vb
       CALL actuators(itime,xb,xbhalf,vbhalf)
       CALL setup_reg(xbhalf)
       motion = delta*vbhalf - regT(q0(:,1))
       CALL vort2flux( q(:,1), omega, s, 1 )
       rhsf = regT(q(:,1) ) - motion
       fbhalf = fb ! initial guess for first iteration in CJGR
       CALL cjgr(fbhalf,rhsf)
       fbhalf = 0.5d0*(fbhalf+fb)
       CALL accel_bodies(itime,vbhalf,fbhalf/dt,accel)
       xbnew = xb + 0.5d0*dt*(vb+vbhalf)
       vb = vb + dt*accel
       CALL actuators(itime,xb,xbnew,vb)
       xb = xbnew
       CALL setup_reg(xb)
       motion = delta*vb - regT(q0(:,1))
       rhsf = regT(q(:,1) ) - motion
       CALL cjgr(fb,rhsf)
    ELSE
    
       CALL vort2flux( q(:,1), omega, s, 1)
       CALL calculate_rhs_adj           ! adjoint "slip velocity" at the body surface and actuators
       rhsf = regT(q(:,1)) + rhs_adj 
       fb = cholsl(rhsf)
    END IF

    !Redistribute the force to make it more physically pleasing
    f_rdst = redistribute(fb)

   ! complete omega on grid 1
       omega(:,:,1) = omega(:,:,1) - ainv( rot( reg( fb )) )
    ! coarsify final omega and correct velocities on all grids
    CALL vort2flux( q, omega, s, mgridlev )
    
    CALL write_slip( itime, q, motion )
    CALL write_udtdx(itime, q, q0 )
    
  END SUBROUTINE advance_adj


!*****************************************************************!


!**********************************************!
function JHtimes( x )

real(kind(0.d0)) :: eps
real(kind(0.d0)), dimension( nf), intent( in ) :: x
!real(kind(0.d0)), dimension( 2:m, 2:n) :: JHtimes
real(kind(0.d0)), dimension( nq ) :: reg_der, JHtimes

reg_der = reg( fb_base )
eps = 1.d-6

call setup_reg( xb_base + eps * x )

reg_der = ( reg(fb_base) - reg_der )/eps

!JHtimes = - rot( reg_der )

JHtimes = reg_der

!Reset regularization/ interp operators to be about base state
call setup_reg( xb_base )

end function
!**********************************************!


!**********************************************!
function JWtimes( x )

real(kind(0.d0)) :: eps
real(kind(0.d0)), dimension( nf ), intent( in) :: x
real(kind(0.d0)), dimension( 3*nb ) :: JWtimes

eps = 1.d-6
JWtimes = QWx( fb_base )

call setup_reg( xb_base + eps * x )

JWtimes = (QWx(fb_base) - JWtimes)/eps


!Reset regularization/ interp operators to be about base state
call setup_reg( xb_base )

end function
!**********************************************!

!**********************************************!
function JEtimes( x )

real(kind(0.d0)) :: eps
real(kind(0.d0)), dimension( nf ), intent( in) :: x
real(kind(0.d0)), dimension( nf ) :: JEtimes

eps = 1.d-6
JEtimes = ( regT(q_base(:,1) )  + regT(q0_base(:,1)) )

call setup_reg( xb_base + eps * x )

JEtimes = ( regT(q_base(:,1) )  + regT(q0_base(:,1)) - JEtimes)/eps

!Reset regularization/ interp operators to be about base state
call setup_reg( xb_base )

end function
!**********************************************!

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
    REAL(KIND(0.D0)), DIMENSION(Nq,nlev) :: vel

    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: s, vort ! streamfun & vort
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: sbc  ! streamfun

    INTEGER :: k
    ! find coarsifieid omegas

!    write(67,*) 1, SUM( omega(:,:,1))

    DO k=2,mgridlev
       omega(:,:,k) = coarsify(omega(:,:,k),omega(:,:,k-1))
!       write(67,*) k, SUM( omega(:,:,k))
    END DO
!    write(67,*)

    ! invert laplacian on largest grid
    ! zero bc's on s assumed for largest grid
    sbc = 0.d0
    ! compute s
    vort = omega(:,:,mgridlev)
    s = ctci( vort )
    sn(:,:,mgridlev) = s
    ! compute vel if desired
    IF (nlev.ge.mgridlev) THEN
       vel(:,mgridlev) = curl(  s, sbc )
    END IF

    ! now telescope in

    DO k=mgridlev-1,1,-1

       CALL get_bc( s, sbc, 1.d0)       ! get bc's from previous s
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
   
  SUBROUTINE stream2vort( stream, omg )

    !***************************************************************!
    !*  Multiscale method to compue S.Lambda.S(s): C^TC operator   *!
    !*  (Applied to the "q x qbase" part of the nonlinear term        *!
    !***************************************************************!

    USE myfft
    INTEGER :: nlev
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n,mgridlev) :: omg, stream

    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: vort ! streamfun & vort
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: omegabc  ! streamfun

    INTEGER :: k

    ! compute laplacian on largest grid
    ! zero bc's on s assumed for largest grid
    omegabc = 0.d0
    ! compute vort
    vort = dst(lam*(1.d0/normalize)*dst( stream(:,:,mgridlev) ))
    omg(:,:,mgridlev) = vort

    ! now telescope in
    DO k=mgridlev-1,1,-1
       CALL get_bc( stream(:,:,k+1), omegabc, 1.d0)  ! get bc's from previous s
       vort = dst(lam*(1.d0/normalize)*dst( stream(:,:,k) ))
       CALL apply_bc( vort, omegabc, -0.25d0)  ! apply them
       omg(:,:,k) = vort
    END DO

  END SUBROUTINE stream2vort

  ! *****************************************************************************************

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
    !*   its values to the edge of a smaller, finer mesh           *!
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
    !*   no preconditioner (yet)                                   *!
    !***************************************************************!

    IMPLICIT NONE

    REAL(KIND(0.d0)), DIMENSION(:), INTENT(IN)  :: b
    REAL(KIND(0.d0)), DIMENSION(:), INTENT(OUT) :: x
    REAL(KIND(0.d0)), DIMENSION(SIZE(b)) :: r,d,q,s
    REAL(KIND(0.d0)) :: delta_new,delta_old,eps, alpha, beta
    INTEGER :: iter

    iter = 0

    r = b - a_times(x) 

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
       alpha = delta_new / DOT_PRODUCT( d, q )
       x = x + alpha * d

       IF (MOD(iter,50).eq.0) THEN
          r = b - a_times( x ) 
       ELSE
          r = r - alpha * q
       END IF

       s = r !m_inverse_times(r)
       delta_old = delta_new
       delta_new = DOT_PRODUCT(r,s)

       beta = delta_new/delta_old
       d = s + beta * d

    END DO

    IF (iter.eq.cg_max_iter) THEN
       WRITE(*,*)  "......WARNING, cg used max iterations"
       WRITE(*,*)    "......itmax = ",iter,", res = ",delta_new
    END IF
    WRITE(*,*) "......CG used ",iter," iterations"

  END SUBROUTINE cjgr

  !*****************************************************************!

  FUNCTION a_times(x) 
    
    !***************************************************************!
    !*   Cholesky factorization of EC(C^t C)^-1 C^t E^t            *!
    !*      performed once only for stationary bodies              *!
    !***************************************************************!
    USE myfft
    REAL(KIND(0.D0)), DIMENSION(Nf) :: x, a_times
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n,mgridlev) :: vort, s
    REAL(KIND(0.D0)), DIMENSION(Nq) :: vel

    ! zero all levels of temporary vorticity
    vort = 0.d0
    ! regularize immersed body forces for level 1 vorticity
    vort(:,:,1) = ainv( rot( reg(x) ) )
    ! coarsify vorticity and find curl(laplace^-1)
    CALL vort2flux( vel, vort, s, 1 )
    ! regularize
    a_times = regT( vel )

  END FUNCTION a_times

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

!-------------------------------------------------------------!
subroutine build_solid_mats( Smat, Mmat, Kmat )

real(kind(0.d0)), dimension(3*(nel + 1), 3*(nel + 1)), intent(out) :: Smat, Kmat, Mmat
real(kind(0.d0)), dimension(nb) :: x, y
real(kind(0.d0)), dimension(nel) :: h0_el, h_el, c, s, beta0
real(kind(0.d0)), dimension(3*nel) :: qL
real(kind(0.d0)), dimension(3*(nel+ 1)) :: Fint
integer, dimension(3*(nel+1)) :: ipiv
integer :: info, neqns, nrhs, lda, ldb, j, lwork, jj
real(kind(0.d0)) :: dts
real(kind(0.d0)), dimension(:), allocatable :: work


Mmat = 0.d0
Kmat = 0.d0


dts = dt

x = xb0_base( 1 : nb)
y = xb0_base( nb + 1 : nf)


call var_init( x, y , beta0 )


call get_M( Mmat )


!Initial configuration of beam elements
call var_update( x, y, u_ib_base, beta0, c, s, qL, Fint )


!Calculate the stiffness matrix for the current load increment and
!current guess for nodal positions:
call get_K_lin( c, s, Kmat )


smat = Kmat + 4.d0/(dts**2)*Mmat


end subroutine

!-------------------------------------------------------------!
subroutine get_K_lin( c, s, K_mat )

real(kind(0.d0)), dimension( nel ), intent( in ) :: c, s
real(kind(0.d0)), dimension(3*nb, 3*nb), intent( out ) :: K_mat
real(kind(0.d0)), dimension(6,6) :: Rmat, K_e
integer :: j, count
integer, dimension(6) :: index

K_mat = 0.d0

do j = 1, nel

    K_e = 0.d0
    K_e(1,1) = R_E*R_th/h_el(j)
    K_e(1,4) = -R_E*R_th/h_el(j)

    K_e(2,2) = 12.d0*R_E*R_sh/(h_el(j)**3.d0)
    K_e(2,3) = 6.d0*R_E*R_sh/(h_el(j)**2.d0)
    K_e(2,5) = -12.d0*R_E*R_sh/(h_el(j)**3.d0)
    K_e(2,6) = 6.d0*R_E*R_sh/(h_el(j)**2.d0)

    K_e(3,2) = K_e(2,3)
    K_e(3,3) = 4.d0*R_e * R_sh/h_el(j)
    K_e(3,5) = -6.d0*R_e * R_sh/(h_el(j)**2.d0)
    K_e(3,6) = 2.d0*R_e * R_sh/h_el(j)

    K_e(4,1) = K_e(1,4)
    K_e(4,4) = R_E*R_th/h_el(j)

    K_e(5,2) = K_e(2,5)
    K_e(5,3) = K_e(3,5)
    K_e(5,5) = 12.d0*R_e * R_sh/(h_el(j)**3.d0)
    K_e(5,6) = -6.d0*R_e * R_sh/(h_el(j)**2.d0)

    K_e(6,2) = K_e(2,6)
    K_e(6,3) = K_e(3,6)
    K_e(6,5) = K_e(5,6)
    K_e(6,6) = 4.d0*R_e * R_sh/h_el(j)


    Rmat = 0.d0
    Rmat(1,1) = c(j)
    Rmat(1,2) = s(j)
    Rmat(2,1) = -s(j)
    Rmat(2,2) = c(j)
    Rmat(3,3) = 1.d0

    Rmat(4,4) = c(j)
    Rmat(4,5) = s(j)
    Rmat(5,4) = -s(j)
    Rmat(5,5) = c(j)
    Rmat(6,6) = 1.d0

    K_e = matmul( K_e, Rmat)
    K_e = matmul( transpose(Rmat), K_e )

    !Assemble the stiffness matrix:
    call get_ind( j, index )
    call K_assemble( K_mat, K_e, index )

end do



!Apply BCs:
count = 0
if (standard) then
    do j = 1, 3

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            K_mat(j,:) = 0.d0
            K_mat(:,j) = 0.d0
            K_mat(j,j) = 1.d0

        end if


    end do



    else
        do j = 3 * nb - 2, 3 * nb

            count = count + 1

            if (bc_type(count) .eq. 1) then !If Dirichlet condition

                K_mat(j,:) = 0.d0
                K_mat(:,j) = 0.d0
                K_mat(j,j) = 1.d0
            end if


        end do
end if




end subroutine


!-------------------------------------------------------------!





!-------------------------------------------------------------!
subroutine sol2fl_ind( ubg, usm )

real(kind(0.d0)), dimension( 3*nb ) :: ubg
real(kind(0.d0)), dimension( nf ) :: usm
integer :: kk


do kk = 1, nb

    usm(kk) = ubg(3*(kk-1) + 1)
    usm(nb + kk) =  ubg(3*(kk-1) + 2)

end do

end subroutine
!-------------------------------------------------------------!


!-------------------------------------------------------------!
subroutine get_M( Mmat )

implicit none

integer :: j, count
real(kind(0.d0)), dimension(3*(nel + 1), 3*(nel + 1)), intent(inout) :: Mmat
real(kind(0.d0)), dimension(6, 6) :: M_e, Rmat
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

    !print *, "M_e = ", M_e


    M_e = R_rho*R_th*h0_el(j)/420.d0 * M_e


    call get_ind( j, index )

    call M_assemble( Mmat, M_e, index )



end do


count = 0

!Apply BCs:
if (standard) then
    do j = 1 , 3

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            Mmat(j,:) = 0.d0
            Mmat(:,j) = 0.d0
        end if


    end do



else
    do j = 3 * nb - 2, 3 * nb

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            Mmat(j,:) = 0.d0
            Mmat(:,j) = 0.d0

        end if


    end do
end if

end subroutine
!-------------------------------------------------------------!

!-------------------------------------------------------------!
subroutine get_M0( Mmat )

implicit none

integer :: j, count
real(kind(0.d0)), dimension(3*(nel + 1), 3*(nel + 1)), intent(inout) :: Mmat
real(kind(0.d0)), dimension(6, 6) :: M_e, Rmat
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

    !print *, "M_e = ", M_e


    M_e = R_rho*R_th*h0_el(j)/420.d0 * M_e


    call get_ind( j, index )

    call M_assemble( Mmat, M_e, index )



end do


count = 0

!Apply BCs:
if (standard) then
    do j = 1, 3

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            Mmat(j,:) = 0.d0
            Mmat(:,j) = 0.d0
            Mmat(j,j) = 1.d0

        end if


    end do



else
    do j = 3 * nb - 2, 3 * nb

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            Mmat(j,:) = 0.d0
            Mmat(:,j) = 0.d0
            Mmat(j,j) = 1.d0
        end if


    end do
end if

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
subroutine var_init(x, y, beta0)

implicit none

integer :: j
real(kind(0.d0)), dimension(nel +1), intent(in) :: x, y
real(kind(0.d0)), dimension(nel), intent(out) :: beta0
real(kind(0.d0)) :: dx, dy

!Update geometry based on original coordinates (x,y) in global frame and
!generalized displacement vector u:

beta0 = 0.d0
h0_el = 0.d0

do j = 1 , nel

    dx = x(j+1) - x(j)
    dy = y(j+1) - y(j)

    h0_el(j) = sqrt( dx**2 + dy**2 )

    beta0( j ) = atan2( dy, dx )

end do

end subroutine
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine var_update( x, y, udisp, beta0, c, s, qL, Fint )

implicit none

integer :: j
real(kind(0.d0)), dimension(nel +1), intent(in) :: x, y
real(kind(0.d0)), dimension(3*(nel + 1)), intent(in) :: udisp
real(kind(0.d0)), dimension(nel), intent(in) :: beta0
real(kind(0.d0)), dimension(nel), intent(out) :: c, s
real(kind(0.d0)), dimension(3*nel), intent(out) :: qL
real(kind(0.d0)), dimension(3*(nel+1)), intent(out) :: Fint
real(kind(0.d0)) :: dx, dy
integer, dimension(6) :: index
integer, dimension(3) :: qind
real(kind(0.d0)) :: uL, theta1L, theta2L, beta1, beta2
real(kind(0.d0)), dimension(3,6) :: Bmat
real(kind(0.d0)), dimension(6) :: qint

!Update geometry based on original coordinates (x,y) in global frame and
!generalized displacement vector u:

c = 0.d0
s = 0.d0
h_el = 0.d0
qL = 0.d0
Fint = 0.d0

do j = 1 , nel

    call get_ind(j, index)

    dx = x(j+1) + udisp(index(4)) - ( x(j) + udisp(index(1)) )
    dy = y(j+1) + udisp(index(5)) - ( y(j) + udisp(index(2)) )

    h_el(j) = sqrt( dx**2 + dy**2 )

    c(j) = dx/h_el(j)
    s(j) = dy/h_el(j)

    uL = ( h_el(j)**2 - h0_el(j)**2 )/( h_el(j) + h0_el(j) )

    qind = (/3*(j-1)+1 , 3*(j-1) + 2, 3*(j-1)+3 /)
    qL(qind(1)) = R_E*R_th*uL/(h0_el(j))

    beta1 = udisp(index(3)) + beta0( j )
    beta2 = udisp(index(6)) + beta0( j )

    theta1L = atan2( c(j)*sin(beta1) - s(j)*cos(beta1), &
    c(j)*cos(beta1) + s(j)*sin(beta1) )

    theta2L = atan2( c(j)*sin(beta2) - s(j)*cos(beta2), &
    c(j)*cos(beta2) + s(j)*sin(beta2) )

    qL(qind(2)) = 2.d0*R_E*R_sh/h0_el(j) * (2.d0*theta1L + theta2L)
    qL(qind(3)) = 2.d0*R_E*R_sh/h0_el(j) * (theta1L + 2.d0*theta2L)

    Bmat = 0.d0

    Bmat(1,1) = -c(j)
    Bmat(1,2) = -s(j)
    Bmat(1,4) = c(j)
    Bmat(1,5) = s(j)

    Bmat(2,1) = -s(j)/h_el(j)
    Bmat(2,2) = c(j)/h_el(j)
    Bmat(2,3) = 1.d0
    Bmat(2,4) = s(j)/h_el(j)
    Bmat(2,5) = -c(j)/h_el(j)

    Bmat(3,1) = -s(j)/h_el(j)
    Bmat(3,2) = c(j)/h_el(j)
    Bmat(3,4) = s(j)/h_el(j)
    Bmat(3,5) = -c(j)/h_el(j)
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
real(kind(0.d0)) :: alpha1, alpha2, c1, c2, s1, s2, phi1, phi2
real(kind(0.d0)), dimension(6,6) :: Rmat
real(kind(0.d0)), dimension(6) :: fsmall

!Take a vector of nodal x and y forces and turn them into the global force
!vector for the FEM formulation...

Fvect = 0.d0

do j = 1 , nel

    !Build element force vector:
    fxsmall = fx(j : j + 1) !Get x forces at that element
    fysmall = fy(j : j + 1) !Get y forces at that element


    !    alpha1 = atan2( fysmall(1), fxsmall(1))
    !    alpha2 = atan2( fysmall(2), fxsmall(2))
    !
    !    phi1 = alpha1 - bta(j)
    !    phi2 = alpha2 - bta(j)
    !
    !    c1 = cos( phi1 )
    !    c2 = cos( phi2 )
    !
    !    s1 = sin( phi1 )
    !    s2 = sin( phi2 )


    fsmall = (/ fxsmall(1), fysmall(1), 0.d0, fxsmall(2), fysmall(2), 0.d0 /)


    call get_Fe( fsmall, h_el(j), F_e )

    !if (j .eq. 1) then
    !    F_e(1:3) = 2.d0 * F_e(1:3)
    !elseif (j .eq. nel) then
    !    F_e(4:6) = 2.d0 * F_e(4:6)
    !end if

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
subroutine get_K( c, s, qL, Kmat )

implicit none

integer :: j, jj, jjj, count
real(kind(0.d0)), dimension(nel), intent(in) :: c, s
real(kind(0.d0)), dimension(3*nel), intent(in) :: qL
real(kind(0.d0)), dimension(3*(nel+1), 3*(nel+1)), intent(out) :: Kmat
integer, dimension(6) :: index
real(kind(0.d0)) :: r_c, Nf, Mf1, Mf2
real(kind(0.d0)), dimension(3,3) :: CL
real(kind(0.d0)), dimension(3,6) :: Bmat
real(kind(0.d0)), dimension(6,6) :: K1, K2, opmat1, opmat2, opmat3, K_e
integer, dimension(3) :: qind
real(kind(0.d0)), dimension(6) :: rvect, zvect


!Build global stiffness matrix from current beam configuration...

Kmat = 0.d0
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

    Bmat(1,1) = -c(j)
    Bmat(1,2) = -s(j)
    Bmat(1,4) = c(j)
    Bmat(1,5) = s(j)

    Bmat(2,1) = -s(j)/h_el(j)
    Bmat(2,2) = c(j)/h_el(j)
    Bmat(2,3) = 1.d0
    Bmat(2,4) = s(j)/h_el(j)
    Bmat(2,5) = -c(j)/h_el(j)

    Bmat(3,1) = -s(j)/h_el(j)
    Bmat(3,2) = c(j)/h_el(j)
    Bmat(3,4) = s(j)/h_el(j)
    Bmat(3,5) = -c(j)/h_el(j)
    Bmat(3,6) = 1.d0

    K1 = matmul(transpose(Bmat), matmul( CL , Bmat) )

    qind = (/ 3*( j - 1 ) + 1, 3*(j-1) + 2, 3*(j-1) + 3 /)
    Nf = qL( qind( 1 ) );
    Mf1 = qL( qind( 2 ) );
    Mf2 = qL( qind( 3 ) );

    zvect = (/s(j), -c(j), 0.d0, -s(j), c(j), 0.d0/)
    rvect = -1.d0*(/ c(j),  s(j), 0.d0, -c(j), s(j), 0.d0 /)


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
    call K_assemble( Kmat, K_e, index );

end do



!Apply BCs:
count = 0
if (standard) then
    do j = 1, 3

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            Kmat(j,:) = 0.d0
            Kmat(:,j) = 0.d0
            Kmat(j,j) = 1.d0

        end if


    end do



else
    do j = 3 * nb - 2, 3 * nb

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            Kmat(j,:) = 0.d0
            Kmat(:,j) = 0.d0
            Kmat(j,j) = 1.d0
        end if


    end do
end if


end subroutine
!-------------------------------------------------------------


!-------------------------------------------------------------
subroutine K_assemble( Kmat, K_e, index )

implicit none

integer, dimension(6), intent(in) :: index
real(kind(0.d0)), dimension(3*(nel+1),3*(nel+1)), intent(inout) :: Kmat
real(kind(0.d0)), dimension(6,6), intent(in) :: K_e
integer :: j, jj, r, rr

do r = 1 , 6

    rr = index( r )

    do j = 1 , 6

        jj = index(j)

        Kmat( rr, jj ) = Kmat( rr, jj ) + K_e( r , j )

    end do

end do


end subroutine
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine apply_bcs( Fvect )

implicit none

integer :: j, count
real(kind(0.d0)), dimension(3*(nel+1)), intent(inout) :: Fvect

! bc_type -- 3x1 vector containing info on the type of BCs at node 1
!            (1 = Dirichlet, 2 = Neumman)
! bc_val -- 3x1 vector containing the values of the BCs (only set to
!           nonzero if the corresponding BC type is Dirichlet)

count = 0

!Apply BCs:
if (standard) then
    do j = 1, 3

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            Fvect(j) = bc_val(count)

        end if


    end do



else
    do j = 3 * nb - 2, 3 * nb

        count = count + 1

        if (bc_type(count) .eq. 1) then !If Dirichlet condition

            Fvect(j) = bc_val(count)

        end if


    end do
end if



end subroutine
!-------------------------------------------------------------


   ! *****************************************************************************************

  FUNCTION qcrossq( q, q_base, q0_base ) RESULT(fomega)

    !***************************************************************!
    !*   Computes the "q x qbase" part of nonlinear term           *!
    !***************************************************************!
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: fomega
    REAL(KIND(0.D0)), DIMENSION(Nq) :: q, q_base, q0_base
    REAL(KIND(0.D0)), DIMENSION(1:m+1,2:n) :: uavg
    REAL(KIND(0.D0)), DIMENSION(2:m,1:n+1) :: vavg
    INTEGER :: i,j

    fomega = 0.d0

    DO j=2,n
       DO i=1,m+1
          uavg(i,j) = 0.5D0*( q_base(u(i,j))+q_base(u(i,j-1)) + q0_base(u(i,j))+q0_base(u(i,j-1)))
       END DO
    END DO
    DO j=1,n+1
       DO i=2,m
          vavg(i,j) = 0.5D0*( q_base(v(i,j))+q_base(v(i-1,j)) + q0_base(v(i,j))+q0_base(v(i-1,j)))
       END DO
    END DO

    DO i=2,m
       DO j=2,n-1
          fomega(i,j+1) = fomega(i,j+1) + 0.5D0*( vavg(i,j+1)* q(u(i,j)) )
          fomega(i,j) = fomega(i,j) + 0.5D0*( vavg(i,j)* q(u(i,j)) )
       END DO
       j = 1
          fomega(i,j+1) = fomega(i,j+1) + 0.5D0*( vavg(i,j+1)*q(u(i,j)) )
       j = n
          fomega(i,j) = fomega(i,j) + 0.5D0*(  vavg(i,j)* q(u(i,j)) )
    END DO

    DO j=2,n
       DO i=2,m-1
           fomega(i+1,j) = fomega(i+1,j) - 0.5D0*( uavg(i+1,j)*q(v(i,j)) )
           fomega(i,j) = fomega(i,j) - 0.5D0*( uavg(i,j)*q(v(i,j)) )
       END DO
       i = 1
           fomega(i+1,j) = fomega(i+1,j) - 0.5D0*( uavg(i+1,j)*q(v(i,j)))
       i = m
           fomega(i,j) = fomega(i,j) - 0.5D0*( uavg(i,j)*q(v(i,j)) )
    END DO

  END FUNCTION qcrossq

  ! *****************************************************************************************

  FUNCTION qcrossq_bc( q, q_base, q0_base ) RESULT(fomegabc)

    !***************************************************************!
    !*   Computes the boundary conditions for the "q x qbase"      *!
    !*   part of the  nonlinear term                               *!
    !***************************************************************!
    REAL(KIND(0.D0)), DIMENSION(Nq) :: q, q_base, q0_base
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: fomegabc
    REAL(KIND(0.D0)), DIMENSION(1:m+1,2:n) :: uavg
    REAL(KIND(0.D0)), DIMENSION(2:m,1:n+1) :: vavg
    INTEGER :: i,j

    fomegabc = 0.D0

    DO j=2,n
       DO i=1,m+1
          uavg(i,j) = 0.5D0*( q_base(u(i,j))+q_base(u(i,j-1)) + q0_base(u(i,j))+q0_base(u(i,j-1)))
       END DO
    END DO
    DO j=1,n+1
       DO i=2,m
          vavg(i,j) = 0.5D0*( q_base(v(i,j))+q_base(v(i-1,j)) + q0_base(v(i,j))+q0_base(v(i-1,j)))
       END DO
    END DO

    DO i=2,m
       j = 1
                 
           fomegabc(bottom+i) = 0.5D0*( vavg(i,j)*q(u(i,j)) )
           !fomegabc(bottom+i) = fomegabc(bottom+i) - 0.5D0*( 0.5D0*( q(v(i,j)) + q(v(i-1,j)) )*q_base(u(i,j)) )
       j = n
           fomegabc(top+i) = 0.5D0*( vavg(i,j+1)*q(u(i,j)) )
           !fomegabc(top+i) = fomegabc(top+i) - 0.5D0*( 0.5D0*( q(v(i,j+1)) + q(v(i-1,j+1)) )*q_base(u(i,j)) )
    END DO
    ! note...we don't need result for i=1 or 1=m+1 since it isn't needed by rot

    DO j=2,n
       i = 1
           fomegabc(left+j) = -0.5D0*( uavg(i,j)*q(v(i,j)) )
           !fomegabc(left+j) = fomegabc(left+j) + 0.5D0*( 0.5D0*( q(u(i,j)) + q(u(i,j-1)) )*q_base(v(i,j)) )
       i = m
           fomegabc(right+j) = -0.5D0*( uavg(i+1,j)*q(v(i,j)) )
           !fomegabc(right+j) = fomegabc(right+j) + 0.5D0*( 0.5D0*( q(u(i+1,j)) + q(u(i+1,j-1)) )*q_base(v(i,j)) )
    END DO
    ! note...we don't need result for j=1 or j=n+1 since it isn't needed by rot

  END FUNCTION qcrossq_bc


  !*****************************************************************!
 

  FUNCTION reg( h0 ) RESULT( h )

    !***************************************************************!
    !*   regularization of immersed body force                     *!
    !***************************************************************!
    REAL(KIND(0.D0)), DIMENSION(Nf), INTENT(IN) :: h0
    REAL(KIND(0.D0)), DIMENSION(Nq)             :: h
    INTEGER :: i,k 

    h = 0.D0
    DO i=1,Nq
       DO k=1,Nsmear
          h(i) = h(i) + smear(i,k)*h0(ismear(i,k))
       END DO
    END DO

  END FUNCTION reg
 !*****************************************************************!

  FUNCTION regT( h ) RESULT( h0 ) 

    !***************************************************************!
    !*  interpolation to body point  (Transpose of reg)            *!
    !***************************************************************!
    REAL(KIND(0.D0)), DIMENSION(Nq), INTENT(IN) :: h
    REAL(KIND(0.D0)), DIMENSION(Nf)             :: h0
    INTEGER :: i,k 

    h0 = 0.D0

    DO i=1,Nq
       DO k=1,Nsmear
          h0(ismear(i,k)) = h0(ismear(i,k)) + smear(i,k)*h(i)
       END DO
    END DO
    
  END FUNCTION regT

  !*****************************************************************!

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

    IF(FORWARD_SIM)THEN

      IF(LINEAR)THEN !linear simulation
        IF (it==istart+1) THEN
           OPEN(unit=106,file="output/slip_lin.dat",form="formatted",status="replace")
        ELSE
           OPEN(unit=106,file="output/slip_lin.dat",form="formatted",status="old",position="append")
        END IF
      ELSEIF(ADJOINT)THEN !forward simulation preceding the adjoint simulation
        IF (it==istart+1) THEN
           OPEN(unit=106,file="output/slip.dat",form="formatted",status="replace")
        ELSE
           OPEN(unit=106,file="output/slip.dat",form="formatted",status="old",position="append")
        END IF
      ELSE !normal forward simulation
        IF (it==1) THEN
           OPEN(unit=106,file="output/slip.dat",form="formatted",status="replace")
        ELSE
           OPEN(unit=106,file="output/slip.dat",form="formatted",status="old",position="append")
        END IF
      ENDIF  

    ELSE ! backward simulation (adjoint)
    
      IF(it==istop)THEN
        OPEN(unit=106,file="output/slip_adj.dat",form="formatted",status="replace")
      ELSE
        OPEN(unit=106,file="output/slip_adj.dat",form="formatted",status="old",position="append")
      END IF
 
    ENDIF
    
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
    
    IF(FORWARD_SIM)THEN

      IF(LINEAR)THEN !linear simulation
        IF (it==istart+1) THEN
           OPEN(unit=106,file="output/cfl_lin.dat",form="formatted",status="replace")
        ELSE
           OPEN(unit=106,file="output/cfl_lin.dat",form="formatted",status="old",position="append")
        END IF
      ELSEIF(ADJOINT)THEN !forward simulation preceding the adjoint simulation
        IF (it==istart+1) THEN
           OPEN(unit=106,file="output/cfl.dat",form="formatted",status="replace")
        ELSE
           OPEN(unit=106,file="output/cfl.dat",form="formatted",status="old",position="append")
        END IF      
      ELSE !normal forward simulation
        IF (it==1) THEN
           OPEN(unit=106,file="output/cfl.dat",form="formatted",status="replace")
        ELSE
           OPEN(unit=106,file="output/cfl.dat",form="formatted",status="old",position="append")
        END IF
      ENDIF  

    ELSE ! backward simulation (adjoint)
    
      IF(it==istop)THEN
        OPEN(unit=106,file="output/cfl_adj.dat",form="formatted",status="replace")
      ELSE
        OPEN(unit=106,file="output/cfl_adj.dat",form="formatted",status="old",position="append")
      END IF
 
    ENDIF

    !velocity is normalized by (a0*D)^0.5
    WRITE(106,*) it*dt, para
    CLOSE(106)

  END SUBROUTINE write_udtdx


!*****************************************************************!

SUBROUTINE write_total_force( it, frc )

!***************************************************************!
!*   write body forces to file                                 *!
!***************************************************************!

use grid
USE parameters
INTEGER          :: it, j, k
REAL(KIND(0.D0)), DIMENSION(nf) :: frc
REAL(KIND(0.D0)), DIMENSION(nb):: forcex, forcey, xhat, yhat
CHARACTER(7) :: force_file_num

forcex = ( getbody(frc, 1, 1, bdy(1)%npts) )
forcey = ( getbody(frc, 2, 1, bdy(1)%npts) )

xhat  = getbody(xb, 1, 1, bdy(1)%npts)
yhat = getbody(xb, 2, 1, bdy(1)%npts)


WRITE(force_file_num,"(I7.7)") it

OPEN(unit=101,file="output/total_force_lin_"//force_file_num//".dat",form="formatted",status="replace")

DO k = 1, nb
WRITE(101,*) xhat(k), yhat(k), forcex(k), forcey(k)
END DO

CLOSE(101)


OPEN(unit=102,file="output/total_disp_lin_"//force_file_num//".dat",form="formatted",status="replace")

DO k = 1, 3*nb
WRITE(102,*) u_ib( k ), ud_ib(k), udd_ib( k )
END DO

CLOSE(102)

END SUBROUTINE write_total_force

!*****************************************************************!
!*****************************************************************!
subroutine write_tip_disp(it, disp)



USE parameters
USE user
INTEGER          :: it
REAL(KIND(0.D0)), dimension(2) :: disp


IF (it==1) THEN
OPEN(unit=10001,file="output/disp_lin.dat",form="formatted",status="replace")
ELSE
OPEN(unit=10001,file="output/disp_lin.dat",form="formatted",status="old",position="append")
END IF
! WRITE(100,*) it, (forcex(k), k=1,n_body+n_actuator), (forcey(k), k=1,n_body+n_actuator),&
!                  (forcex_lab(k), k=1,n_body+n_actuator), (forcey_lab(k), k=1,n_body+n_actuator)
WRITE(10001,*) it, disp(1), disp(2)
CLOSE(10001)


end subroutine


!*****************************************************************!

!*****************************************************************!

  SUBROUTINE write_force( it, frc )

    !*****************************************************************!
    !*   write immersed body forces to file for forward simulations  *!
    !*****************************************************************!
  
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
    
    IF(LINEAR)THEN ! linearized simulation
      IF (it==1) THEN
         OPEN(unit=100,file="output/force_lin.dat",form="formatted",status="replace")
      ELSE
         OPEN(unit=100,file="output/force_lin.dat",form="formatted",status="old",position="append")
      END IF
    ELSE !normal simulation
      IF (it==1) THEN
         OPEN(unit=100,file="output/force.dat",form="formatted",status="replace")
      ELSE
         OPEN(unit=100,file="output/force.dat",form="formatted",status="old",position="append")
      END IF
    ENDIF  
    WRITE(100,*) it, (forcex(k), k=1,n_body+n_actuator), (forcey(k), k=1,n_body+n_actuator),&
                     (forcex_lab(k), k=1,n_body+n_actuator), (forcey_lab(k), k=1,n_body+n_actuator)
    CLOSE(100)
 
  END SUBROUTINE write_force

!*****************************************************************!

  SUBROUTINE write_force_adj( it, frc ,iopt)

    !******************************************************************!
    !*   write immersed body forces to file for adjoint optimization  *!
    !******************************************************************!
  
    USE parameters
    USE user
    INTEGER          :: it, i, j, k, iopt
    REAL(KIND(0.D0)), DIMENSION(nf) :: frc
    REAL(KIND(0.D0)), DIMENSION(n_body+n_actuator):: forcex, forcey, forcex_lab, forcey_lab
    CHARACTER(3) :: char_iopt

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
    
    WRITE(char_iopt,"(I3.3)") iopt
    IF (it==istart+1) THEN
      OPEN(unit=100,file="output/force"//char_iopt//".dat",form="formatted",status="replace")
    ELSE
      OPEN(unit=100,file="output/force"//char_iopt//".dat",status="old",position="append")
    END IF
   
    WRITE(100,*) it, (forcex(k), k=1,n_body+n_actuator), (forcey(k), k=1,n_body+n_actuator),&
                     (forcex_lab(k), k=1,n_body+n_actuator), (forcey_lab(k), k=1,n_body+n_actuator)
    CLOSE(100)
 
  END SUBROUTINE write_force_adj

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

  SUBROUTINE accel_bodies(it,vlocal,flocal,abloc)

    USE grid
    USE user

    INTEGER :: it,ibdy,istrt,iend
    REAL(KIND(0.D0)), DIMENSION(nf) :: x, vlocal,flocal, abloc
    REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: ax,ay

    abloc = 0.D0 ! default is no acceleration

    ! this routine returns the acceleration for each body on the domain

    DO ibdy=1,n_body
       IF (bdy(ibdy)%moving) THEN
          istrt = bdy(ibdy)%ipos
          iend = bdy(ibdy)%ipos+bdy(ibdy)%npts - 1
          ALLOCATE( ax(iend-istrt+1), ay(iend-istrt+1) )
          CALL body_accel( ibdy, it, xb(f(istrt:iend,1)), xb(f(istrt:iend,2)), &
                                     vlocal(f(istrt:iend,1)), vlocal(f(istrt:iend,2)), &
                                     flocal(f(istrt:iend,1)), flocal(f(istrt:iend,2)), &
                                     ax, ay)
          abloc(f(istrt:iend,1)) = ax
          abloc(f(istrt:iend,2)) = ay
          DEALLOCATE( ax,ay )
       END IF
    END DO

  END SUBROUTINE accel_bodies

 !*****************************************************************!

  SUBROUTINE actuators(it,xold,xlocal,vlocal)
    USE variables
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
       IF(ADJOINT.AND.FORWARD_SIM)THEN
         vlocal(f(istrt:iend,1)) = vx + control(it,nbodyf*2+iact)
         vlocal(f(istrt:iend,2)) = vy + control(it,nbodyf*2+nact+iact)
       ELSE
         vlocal(f(istrt:iend,1)) = vx
         vlocal(f(istrt:iend,2)) = vy
       ENDIF
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
 
    INTEGER :: it,i,j,k
    REAL(KIND(0.D0)), DIMENSION(Nq) :: q, dq
    REAL(KIND(0.D0)), DIMENSION(m,n) :: dqx,dqy
    REAL(KIND(0.D0)) :: xx,yy,uu,vv

    dqx = 0.d0
    dqy = 0.d0
    dq = 0.d0

  IF(ADJOINT.and.FORWARD_SIM)THEN
    DO j=2,n-1
       DO i=2,m-1
          xx = delta*(REAL(i)-0.5d0) - offsetx
          yy = delta*(REAL(j)-0.5d0) - offsety
          uu = 0.5d0 * ( q(u(i,j)) + q(u(i+1,j)) ) / delta
          vv = 0.5d0 * ( q(v(i,j)) + q(v(i,j+1)) ) / delta
          DO k=1,nbodyf
            dqx(i,j) = delta* ( bodyforcex(it,xx,yy,uu,vv) + control(it,k) * gaussianforce(xx,yy,gx(k),gy(k)) )
            dqy(i,j) = delta* ( bodyforcey(it,xx,yy,uu,vv) + control(it,nbodyf+k) * gaussianforce(xx,yy,gx(k),gy(k)) )
          ENDDO
       END DO
    END DO
  ELSE
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
  ENDIF
    

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

  ! *****************************************************************************************  
  SUBROUTINE get_motion_deriv(itime)

    !************************************************************************!
    !*    Computes the derivative of the Hamiltonian with respect           *!
    !*    to the grid translation and rotation                              *!
    !************************************************************************!   
  REAL(KIND(0.D0)),DIMENSION(nq) :: ux,uy,urot,q0temp
  REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: bc_base
  REAL(KIND(0.D0)) :: xx, yy
  INTEGER :: itime,k,i,j

  ! Q * U terms
  DO i=nbodyf*2+nact*2+1,nbodyf*2+nact*2+3
    gradient(itime,i) = CU(i)*control(itime,i)
  ENDDO
  
  ! Create ux, uy and urot vectors 
  q0temp=0.d0
  ux=0.d0
  uy=0.d0
  urot=0.d0

    DO i=1,m+1
      DO j=1,n
        ux(u(i,j)) = 1
        yy = (REAL(j)-0.5d0)*delta-offsety
        urot(u(i,j)) = -(yy-roy)        
      ENDDO
    ENDDO
    
    DO i=1,m
      DO j=1,n+1 
        uy(v(i,j)) = 1
        xx = (REAL(i)-0.5d0) - offsetx
        urot(v(i,j)) = (xx-rox)       
      ENDDO
    ENDDO
  ! Derivatives with respect to qc
  gradient(itime,nbodyf*2+nact*2+1) = gradient(itime,nbodyf*2+nact*2+1) + &
                                             delta*qTq(ux,reg(fb),1,m,1,n)
  gradient(itime,nbodyf*2+nact*2+2) = gradient(itime,nbodyf*2+nact*2+2) + &
                                             delta*qTq(uy,reg(fb),1,m,1,n)
  gradient(itime,nbodyf*2+nact*2+3) = gradient(itime,nbodyf*2+nact*2+3) + &
                                             delta*qTq(urot,reg(fb),1,m,1,n)
  
  !    print*, gradient(itime,nbodyf*2+nact*2+1), gradient(itime,nbodyf*2+nact*2+2), &
  !  gradient(itime,nbodyf*2+nact*2+3)
  
  ! Advection term:
  bc_base = 0.d0
  IF(mgridlev.gt.1)THEN
    CALL get_bc(omega_base(:,:,2),bc_base,0.25d0)
  ENDIF
  
  gradient(itime,nbodyf*2+nact*2+1) = gradient(itime,nbodyf*2+nact*2+1) - &
                                             1/delta* qTq(   nonlinear(omega_base(:,:,1),ux,q0temp,bc_base),&
                                                             q(:,1),1,m,1,n)
  gradient(itime,nbodyf*2+nact*2+2) = gradient(itime,nbodyf*2+nact*2+2) - &
                                             1/delta* qTq(   nonlinear(omega_base(:,:,1),uy,q0temp,bc_base),&
                                                             q(:,1),1,m,1,n)
  gradient(itime,nbodyf*2+nact*2+3) = gradient(itime,nbodyf*2+nact*2+3) - &
                                             1/delta* qTq(   nonlinear(omega_base(:,:,1),urot,q0temp,bc_base),&
                                                             q(:,1),1,m,1,n)  
   !print*, gradient(itime,nbodyf*2+nact*2+1), gradient(itime,nbodyf*2+nact*2+2), &
   !gradient(itime,nbodyf*2+nact*2+3)
  END SUBROUTINE get_motion_deriv

  ! *****************************************************************************************
   SUBROUTINE calculate_gradient(itime)
   
    !************************************************************************!
    !*    Computes the Gradient of the Hamiltonian with respect to the      *!  
    !*    control variables                                                 *!
    !************************************************************************!

  INTEGER :: itime,i,j,k, nbpts
  REAL(KIND(0.D0)) :: xx, yy
  
  ! X and Y body force
  DO k=1,nbodyf
    gradient(itime,k) = CU(k)*control(itime,k)
    gradient(itime,k+nbodyf) = CU(k+nbodyf)*control(itime,k+nbodyf)
    DO j=1,n
       DO i=1,m+1
          xx = delta*(REAL(i)-1d0) - offsetx
          yy = delta*(REAL(j)-0.5d0) - offsety
          gradient(itime,k) = gradient(itime,k) + delta*gaussianforce(xx,yy,gx(k),gy(k))*q(u(i,j),1)
       END DO
    END DO      
  
    DO j=1,n+1
       DO i=1,m
          xx = delta*(REAL(i)-0.5d0) - offsetx
          yy = delta*(REAL(j)-1.d0) - offsety
          gradient(itime,k+nbodyf) = gradient(itime,k+nbodyf) + delta*gaussianforce(xx,yy,gx(k),gy(k))*q(v(i,j),1)
       END DO
    END DO      
  ENDDO
    
  !Actuators
  
  ! Find how many immersed body points there are in total
  nbpts = 0    
  DO i=1,n_body 
    nbpts=nbpts + bdy(i)%npts       
  END DO

  DO i=nbpts+1,nb       ! i = index of actuator in fb
    j=i-nbpts           ! j = used for index of actuator in CU, gradient and control
    ! X actuator velocity
    gradient(itime,nbodyf*2+j) = CU(nbodyf*2+j)*control(itime,nbodyf*2+j)+fb(f(i,1))
    ! Y actuator velocity
    gradient(itime,nbodyf*2+nact+j) = CU(nbodyf*2+nact+j)*control(itime,nbodyf*2+nact+j)+fb(f(i,2))
  ENDDO
  
  !grid translation and rotation
  CALL get_motion_deriv(itime)  
    
  END SUBROUTINE calculate_gradient

  ! ************************************************************************************************************************************************************

        
  SUBROUTINE eval_cost(distance,totalcost)

  !******************************************************************!
  !*     - Given generalized "distance", r                          *!
  !*       Calculate "totalcost" from forward simulation            *!
  !******************************************************************!
    INTEGER :: it, i
    REAL(KIND(0.d0)):: distance, totalcost
    
    CALL update_control(distance)

    it = istart           ! always start from istart  
    FORWARD_SIM=.TRUE.
    CALL setup_variables
    DO WHILE (it < istop )
      IF(.NOT.LINEAR) THEN
        CALL advance(it)
        CALL integ_cost(it)
      ELSE
        CALL advance_lin(it)
        CALL integ_cost(it)
      ENDIF  
    END DO

    CALL calculate_totalcost(totalcost)
    
    IF(totalcost+1.d0.eq.totalcost)THEN
      WRITE(*,*) '**** WARNING: Totalcost is a NaN: Setting value to 1.e20 ****'
      totalcost = 1.e20
    ENDIF
    
    CALL destroy_variables
    
  END SUBROUTINE eval_cost
  
  ! *****************************************************************************************  
  
  SUBROUTINE minimization_1d(distance_opt)

  !******************************************************************!
  !*     - Find optimal distance r, "distance_opt"                  *!
  !*       which gives minimum cost functional J, "totalcost_opt    *!
  !******************************************************************!
  
  USE variables
  
    REAL(KIND(0.d0)) :: distance_opt, totalcost_opt, tol
    REAL(KIND(0.d0)) :: ax,  bx,  cx  
    REAL(KIND(0.d0)) :: afun, bfun, cfun 
   
    
    control_old = control
    
    ax=2.0d0
    bx=2.2d0
    cx=0.d0

    afun=0.d0
    bfun=0.d0
    cfun=0.d0
    CALL mnbrak(ax,bx,cx,afun,bfun,cfun)

    distance_opt = bx
    totalcost_opt = bfun


    tol=0.1d0
    CALL brent(ax, bx, cx, tol, distance_opt, totalcost_opt)

    WRITE(*,*) 'DISTANCE'
    WRITE(*,*) distance_opt
    WRITE(*,*) '*************************'

  END SUBROUTINE minimization_1d

  ! *****************************************************************************************  
  
  SUBROUTINE mnbrak(ax, bx, cx, afun, bfun, cfun)

  !******************************************************************!
  !*     - Numerical Recipes, pp. 393
  !*     - Given distinct initial points ax and bx                  *!
  !*       Search in the downhill direction and returns new points  *!
  !*       ax, bx, cx that bracket a minimum of the function        *!
  !*       fa, fb, fc are function values at three points           *!
  !******************************************************************!
    REAL(KIND(0.d0)) :: ax, bx, cx
    REAL(KIND(0.d0)) :: afun, bfun, cfun 
    REAL(KIND(0.d0)) :: GOLD, GLIMIT, TINY  
    REAL(KIND(0.d0)) :: dum, uufun, qq, rr, uu, uulim
   
   
        WRITE(*,*) ' ---------  Bracketing Distance  -------------------------'
    GOLD = 1.618034
    GLIMIT = 10.d0
    TINY = 1.e-20

    CALL eval_cost(ax, afun)
    CALL eval_cost(bx, bfun)
    
    IF (bfun .GT. afun) then   ! Switch roles of a and b if necessary
        dum = ax               ! s.t. we cann go downhill from a to b
         ax = bx
         bx = dum
        dum = bfun
       bfun = afun
       afun = dum 
    END IF
    cx = bx + GOLD * (bx-ax)

    CALL eval_cost(cx, cfun)

    DO WHILE (bfun .GE. cfun)
 
       rr = (bx-ax)*(bfun-cfun)
       qq = (bx-cx)*(bfun-afun)
       uu = bx-((bx-cx)*qq-(bx-ax)*rr) / (2.d0*SIGN(MAX(abs(qq-rr),TINY),qq-rr))
       uulim = bx+GLIMIT*(cx-bx)

       IF ((bx-uu)*(uu-cx) .GT. 0.d0) THEN

          CALL eval_cost(uu, uufun)
          IF (uufun .LT. cfun) THEN
             ax   = bx
             afun = bfun
             bx   = uu
             bfun = uufun
    GO TO 1
          ELSE IF (uufun .GT. bfun) THEN
             cx   = uu
             cfun = uufun
    GO TO 1
          END IF
          uu = cx + GOLD * (cx-bx)
          CALL eval_cost(uu, uufun)

       ELSE IF ((cx-uu)*(uu-uulim) .GT. 0.d0) THEN

          CALL eval_cost(uu, uufun)
          IF (uufun .LT. cfun) THEN
             bx = cx
             cx = uu
             uu = cx + GOLD*(cx-bx)
             bfun = cfun
             cfun = uufun
             CALL eval_cost(uu, uufun)
          END IF

       ELSE IF ((uu-uulim)*(uulim-cx) .GE. 0.d0) THEN

          uu = uulim
          CALL eval_cost(uu, uufun)

       ELSE
          
          uu = cx + GOLD*(cx-bx)
          CALL eval_cost(uu, uufun)

       END IF
       ax = bx
       bx = cx
       cx = uu
       afun = bfun
       bfun = cfun
       cfun = uufun

    END DO
    1  RETURN
  END SUBROUTINE mnbrak

  ! *****************************************************************************************  
  
  SUBROUTINE brent(aax, bbx, ccx, tol, xmin, fxmin)

  !******************************************************************!
  !*     - Brent's Method to find 1D minimum                        *!
  !*     - Numerical Recipes, pp. 397                               *!
  !*     - Given bracketing triplet of obscissas ax,bx,cx           *!
  !*       ax<bx<cx, f(bx)<f(ax), f(bx)<f(cx)                       *!
  !*       Isolate minimum within tol                               *!
  !*       xmin = abscissa of min, fxmin=f(xmin)                    *!
  !******************************************************************!
    REAL(KIND(0.d0)) :: aax, bbx, ccx, tol, xmin, fxmin, CGOLD, ZEPS
    REAL(KIND(0.d0)) :: aa, bb, dd, ee, eetemp, fuu, fvv, fww, fxx  
    REAL(KIND(0.d0)) :: pp, qq, rr, tol1, tol2, uu,  vv,  ww,  xx, xxm  
    INTEGER :: itmax, iter
    
    WRITE(*,*) ' ---------  Performing Minimization  ------------------'
    
    itmax = 100
    CGOLD = 0.3819660d0
    ZEPS = 1.0e-10

    aa = min(aax, ccx)
    bb = max(aax, ccx)
    vv = bbx
    ww = vv
    xx = vv
    ee = 0.d0
    CALL eval_cost(xx, fxx)
    fvv = fxx
    fww = fxx

    DO iter = 1,itmax
       
       xxm = 0.5d0 * (aa+bb)
       tol1 = tol * abs(xx) + ZEPS
       tol2 = 2.d0 * tol1
       IF (abs(xx-xxm) .LE. (tol2-0.5d0*(bb-aa))) THEN
           GO TO 3
       END IF
       IF (abs(ee) .GT. tol1) THEN
          rr=(xx-ww)*(fxx-fvv)
          qq=(xx-vv)*(fxx-fww)
          pp=(xx-vv)*qq - (xx-ww)*rr
          qq=2.d0*(qq-rr)
          IF (qq .GT. 0.d0) THEN
             pp=-pp
          END IF
          qq=abs(qq)
          eetemp=ee
          ee=dd
          IF ((abs(pp) .GE. abs(0.5d0*qq*eetemp)) .OR. (pp .LE. qq*(aa-xx)) .OR. (pp .GE. qq*(bb-xx))) THEN
             GO TO 1
          END IF
          dd = pp/qq
          uu = xx+dd
          IF ((uu-aa .LT. tol2) .OR. (bb-uu .LT. tol2)) THEN
             dd = sign(tol1, xxm-xx)
          END IF
          GO TO 2
       END IF
    1  IF (xx .GE. xxm) THEN
          ee=aa-xx
       ELSE
          ee=bb-xx
       END IF
       dd=CGOLD*ee
    2  IF (abs(dd) .GE. tol1) THEN
          uu=xx+dd
       ELSE
          uu=xx+sign(tol1,dd)
       END IF
       CALL eval_cost(uu, fuu)
       IF (fuu .LE. fxx) THEN
          IF (uu .GE. xx) THEN
             aa=xx
          ELSE
             bb=xx
          END IF
          vv=ww
          fvv=fww
          ww=xx
          fww=fxx
          xx=uu
          fxx=fuu
       ELSE
          IF (uu .LT. xx) THEN
             aa = uu
          ELSE
             bb = uu
          END IF
          IF ((fuu .LE. fww) .OR. (ww .EQ. xx)) THEN
              vv = ww
             fvv = fww
              ww = uu
             fww = fuu
          ELSE IF ((fuu .LE. fvv) .OR. (vv .EQ. xx) .OR. (vv .EQ. ww)) THEN
             vv=uu
             fvv=fuu
          END IF
       END IF

    END DO 
    STOP 'brent exceed maximum iterations!!!!!!'
 3  xmin = xx
    fxmin = fxx


  END SUBROUTINE brent


  ! *****************************************************************************************
  
   
END MODULE operators

