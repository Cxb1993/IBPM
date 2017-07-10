PROGRAM main

  ! c Kunihiko (Sam) Taira and Tim Colonius
  ! December 2005
  ! calling program for ibfs v1.0
  
  ! Immersed Boundary Fractional Step method
  ! Current features/limitations
  !    - 2D, Incompressible Navier-Stokes
  !    - 2nd-order-accurate except for IB
  !    - 1st-order-accurate IB
  !    - Requires nearly uniform grid near immersed boundary
  !    - serial only

  ! Update: ibfs v1.2
  !    - actuation added
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
  USE grid
  USE variables
  USE operators
  USE operators_pressure
  USE optimization

  IMPLICIT NONE
  INTEGER       :: it,iopt,i,j,k
  CHARACTER(10) :: date, time
  REAL(KIND(0.D0)) :: distance ! Distance for the control variables update
  REAL(KIND(0.D0)) :: totalcost
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: init_adj_circ
  
  CALL input
  CALL setup

    IF (compute_pressure) CALL setup_pressure
  



    IF(ADJOINT) THEN
     
      WRITE(*,*) ' ---------------------------------------------------------' 
      WRITE(*,*) ' ---------  Running Adjoint Optimization  ----------------'
      WRITE(*,*) ' ---------------------------------------------------------'

      DO iopt=iopt_start,iopt_stop
        CALL setup_adj_variables
        ALLOCATE(init_adj_circ(2:m,2:n,mgridlev))
        WRITE(*,*) '----------------------------------------------------------'
        WRITE(*,*) 'Optimization iteration:',iopt
        WRITE(*,*) '----------------------------------------------------------'
! ---------------- Get optimization parameters ---------------------------------
        WRITE(*,*) ' ---------  Loading control parameters  ------------------'
        CALL read_control(iopt)    
        IF(iopt.eq.iopt_start) CALL write_control(iopt)     
! ---------------- Forward Simulation ------------------------------------------
        WRITE(*,*) ' ---------  Running Forward Simulation  ------------------'
        FORWARD_SIM=.TRUE.
        it = istart
        CALL setup_variables
        DO WHILE (it < istop )
          IF(.NOT.LINEAR) THEN
            CALL advance(it)
            CALL integ_cost(it)
          ELSE
            WRITE(*,*) '**** WARNING : NONLINEAR SIMULATION NOT YET SUPPORTED FOR ADJOINT EQUATIONS ***'
            CALL advance_lin(it)
            CALL integ_cost(it)
          ENDIF
          ! save a restart
          IF ((MOD(it,isave).eq.0).or.(it==istart+1).or.(it==istop)) THEN
            CALL write_variables(it)
            call write_total_force(it, f_rdst / dt)
            IF (compute_pressure) THEN 
              CALL calculate_pressure( it )
              CALL write_pressure(pressure,it)
            END IF
          END IF
          CALL write_force_adj( it, fb/Dt,iopt )
        END DO
      
! ---------------- Evaluate cost & initial condition for adjoint simulation ----
        WRITE(*,*) ' ---- Evaluating initial condition for adj simulation ----'       

        init_adj_circ = 0.d0
        DO i=mmin+1,mmax
          DO j=nmin+1,nmax
            init_adj_circ(i,j,1) = R*omega(i,j,1)
          ENDDO
        ENDDO
      
        WRITE(*,*) ' ---- Evaluating total cost from results -----------------'
        CALL calculate_totalcost(totalcost)
        CALL write_cost(iopt,totalcost)
      
      
! ---------------- Adjoint Simulation ------------------------------------------
        WRITE(*,*) ' ---------  Running Adjoint Simulation  ------------------' 
        FORWARD_SIM=.FALSE.
        CALL check_adj_files
        CALL destroy_variables ! clear all variables from forward simulation
        CALL setup_variables ! for adjoint simulation     
        it = istop+1
        omega = init_adj_circ ! set initial condition
        CALL vort2flux( q(:,1), omega, s, 1 )
        CALL write_variables(it-1)
        DO WHILE(it.gt.istart+1) 
          CALL advance_adj(it)
          CALL calculate_gradient(it)
           ! save a restart
           IF ((MOD(it,isave).eq.0).or.(it==istart+1)) then
                CALL write_variables(it)
                call write_total_force(it, f_rdst / dt)
            end if
        ENDDO
        
        CALL write_gradient(iopt)
        CALL destroy_variables
! ---------------- Update Control ----------------------------------------------
        IF(iopt_start.ne.iopt_stop)THEN ! Only update control if several optimization iterations are required
          WRITE(*,*) ' ---------  Updating control  ----------------------------'

          ! Get distance for each control variable
          CALL  minimization_1d(distance)
        
          ! Update control for each control variable
          CALL update_control(distance)
          
          ! Write new control to file
          CALL write_control(iopt+1)
        ENDIF
        
! ---------------- Clean up ----------------------------------------------------           
        CALL destroy_adj_variables
        DEALLOCATE(init_adj_circ)
      ENDDO
      
! ---------------- ONLY A FORWARD SIMULATION REQUIRED --------------------------
    ELSE

    FORWARD_SIM=.TRUE.

      CALL setup_variables
      IF(LINEAR)THEN
       IF(STEADY)THEN
         WRITE(*,*) ' ---------------------------------------------------------' 
         WRITE(*,*) ' ---------  Running Linearized STEADY Simulation  --------'
         WRITE(*,*) ' ---------------------------------------------------------' 
       ELSE
         WRITE(*,*) ' ---------------------------------------------------------' 
         WRITE(*,*) ' ---------  Running Linearized UNSTEADY Simulation  ------'
         WRITE(*,*) ' ---------------------------------------------------------' 
       ENDIF
       CALL check_lin_files
     ELSE
         WRITE(*,*) ' ---------------------------------------------------------' 
         WRITE(*,*) ' ---------  Running Non-Linear Simulation  ---------------'
         WRITE(*,*) ' ---------------------------------------------------------'
     ENDIF
      
      it=istart
      DO WHILE (it < istop )
         IF(.NOT.LINEAR) THEN
           CALL advance(it)
         ELSE
           CALL advance_lin(it)
         ENDIF
         ! save a restart
         IF ((MOD(it,isave).eq.0).or.(it==istop)) THEN
            CALL write_variables(it)
            call write_total_force(it, f_rdst / dt)
            IF (compute_pressure) THEN 
               CALL calculate_pressure( it )
               CALL write_pressure(pressure,it)
            END IF
         END IF

        CALL write_force( it, fb/Dt )

      END DO
      CALL destroy_variables
    ENDIF
    
    

    IF (compute_pressure) CALL destroy_pressure
    CALL destroy_grid
    
      WRITE(*,*) ' ---------------------------------------------------------' 
      WRITE(*,*) ' ---------  RUN SUCCESSFULLY COMPLETED  ------------------'
      WRITE(*,*) ' ---------------------------------------------------------'

END PROGRAM main
