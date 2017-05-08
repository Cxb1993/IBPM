MODULE user
  
  USE parameters
  USE grid
  IMPLICIT NONE

  ! module for user-defined routines to move bodies relative to the grid, add actuators, or move bodies with the grid.
  ! 
  ! see specific examples in the test cases

CONTAINS

  !--------------------------------------------------------------------------------

  SUBROUTINE body_accel( bdyno, it, x, y, vx, vy, fx, fy, ax, ay )

    INTEGER, INTENT(IN) :: bdyno, it
    REAL(KIND(0.D0)), DIMENSION(:) :: x,y,vx, vy, fx,fy, ax, ay

    ! DO NOT MODIFY x,y,vx,vy,fx,fy

print *, "bdyno", bdyno

    ax = 0.D0 ! do not remove; overwrite values below inside if block
    ay = 0.D0 ! do not remove; overwrite values below inside if block

    ! for each body (that you want to move relative to the grid), specify the acceleration, potentially as a function of the forces on the body.

    ! inputs
    !  bdyno : the body number (the number x in the filename body.00x.inp
    !  it : the time step (multiply by dt to get physical time)
    !  x,y  : the positions of each point on the body
    !  fx,fy : the forces (per unit of grid spacing) on the body
    !  
    ! outputs
    !  ax, ay : the acceleration components

!    IF (bdyno==1) THEN
!
!
!        ay = -0.5D0*sin(dt*real(it)) !-omega^2*R*cos(omega*t)
!        ax = -0.5d0*cos(dt*real(it)) !-omega^2*R*sin(omega*t)
!
!    END IF


  END SUBROUTINE body_accel

 !*****************************************************************!

  SUBROUTINE actuator_vel( actno, it, x, y, vx, vy )

    INTEGER, INTENT(IN) :: actno, it
    REAL(KIND(0.D0)), DIMENSION(:) :: x,y,vx,vy
    
    ! DO NOT MODIFY x,y

    vx = 0.D0 ! do not remove; overwrite values below inside if block
    vy = 0.D0 ! do not remove; overwrite values below inside if block

    ! for each actuator, specify the imposed velocity as a function of time

    ! inputs
    !  actno : the actuator number (the number x in the filename actuator.00x.inp
    !  it : the time step (multiply by dt to get physical time)
    !  x,y  : the positions of each point defining the actuator
    !  
    ! outputs
    !  vx, vy : the velocity components
    
    IF (actno==1) THEN

       vx = 0.D0 ! Steady actuation with speed 0 in x dirn.
       vy = 0.D0 

    END IF


  END SUBROUTINE actuator_vel

 !*****************************************************************!

  FUNCTION motion_grid( it )

    REAL(KIND(0.D0)), DIMENSION(5) :: motion_grid
    REAL(KIND(0.D0)) :: rx, ry
    INTEGER  :: it

    ! motion_grid specifies the (possibly time-varying) u component, v component of the velocity, angular velocity, and the coordinates of the center of rotation when the body is moving with the grid.
   
    rx = 0.D0
    ry = 0.D0

    ! default is no flow
    motion_grid(1) = -1.D0   ! x-component of velocity in the body-fixed frame
    motion_grid(2) = 0.D0   ! y-component of velocity in the body-fixed frame
    motion_grid(3) = 0.D0   ! angular velocity
    motion_grid(4) = rx    ! x-coordinate of the center of rotation
    motion_grid(5) = ry    ! y-coordinate of the center of rotation

  END FUNCTION motion_grid

 !*****************************************************************!

  FUNCTION bodyforcex(it, x, y, u, v ) 

   REAL(KIND(0.D0)) :: x,y,u,v, bodyforcex
   INTEGER  :: it
  
   ! specify a body force (per unit mass) on RHS of u-momentum equation,
   ! as a function of location (x,y) and possibly velocity (u,v) and possibly time

   ! the specified function should decay to zero near the edges of the (mgrid=1)
   ! inner domain
   !
   ! units are in terms of velocity (result should be consistent
   ! with U^2/L units

   bodyforcex = 0.d0 ! do not remove; overwrite values below inside if block

 END FUNCTION bodyforcex

 !*****************************************************************!

 FUNCTION bodyforcey(it, x, y, u, v ) 

   REAL(KIND(0.D0)) :: x,y,u,v, bodyforcey, xc, yc
   INTEGER  :: it
  
   ! specify a body force (per unit mass) on RHS of v-momentum equation,
   ! as a function of location (x,y) and possibly velocity (u,v) and possibly time

   ! the specified function should decay to zero near the edges of the (mgrid=1)
   ! inner domain
   !
   ! units are in terms of velocity (result should be consistent
   ! with U^2/L units

   bodyforcey = 0.d0 ! do not remove; overwrite values below inside if block


!Add small body force for a small unit of time
xc = 0.5D0 ! center of the bodyforce
yc = -0.2d0

IF ((real(it)*dt) .ge. 115.0d0 .and. real(it)*dt .le. 116.0d0) then

    IF ( ( (x-xc)**2.d0+(y-yc)**2.d0 ) .le. (10.D0*delta)**2.d0 ) THEN
        bodyforcey = 5.d-4
    end if
end IF

 END FUNCTION bodyforcey

 !*****************************************************************!

!******************
subroutine bf_sol( fx, fy )

real(kind(0.d0)), dimension(nb) :: fx, fy

fx = 0.d0
fy = 0.d0


end subroutine

END MODULE user

