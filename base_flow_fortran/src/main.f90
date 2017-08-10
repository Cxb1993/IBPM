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

  USE parameters
  USE grid
  USE variables
  USE operators
  USE operators_pressure

  IMPLICIT NONE
  INTEGER       :: it
  CHARACTER(10) :: date, time

  CALL input
  CALL setup
  CALL setup_variables
  IF (compute_pressure) CALL setup_pressure
  
  it = istart
 
  DO WHILE (it < istop )

     CALL advance(it)
     ! save a restart

     IF ((MOD(it,isave).eq.0).or.(it==istop)) THEN
        CALL write_variables(it)
        CALL write_total_force(it, f_rdst/dt)
        IF (compute_pressure) THEN
           CALL calculate_pressure( it )
           CALL write_pressure(pressure,it)
        END IF
     END IF
    
    CALL write_force( it, f_rdst/Dt )
    if (standard) then
        call write_tip_disp(it, u_ib(3*nb - 2 : 3*nb - 1))
    else
        call write_tip_disp(it, u_ib(1 : 2) )
    end if

  END DO

  CALL destroy_variables
  IF (compute_pressure) CALL destroy_pressure
  CALL destroy_grid

END PROGRAM main
