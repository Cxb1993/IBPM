MODULE optimization
  
  ! c Kunihiko (Sam) Taira and Tim Colonius
  ! December 2005
  ! calling program for ibfs v1.0

  ! Updated: March 2006 
  !   ibfs v1.2 - now includes actuation
  !
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
  USE variables
  USE grid

  
  IMPLICIT NONE 

CONTAINS

 ! *****************************************************************************************
   SUBROUTINE integ_cost(it)

    !************************************************************************!
    !*    Increments the interpolation of the cost functional based on      *!  
    !*    current value of state variables.                                 *!
    !*    The dt factor will be added later                                 *!
    !*    int = int + f(itime)                                              *!
    !************************************************************************!


    REAL(KIND(0.D0)) :: cost_gam, cost_bx, cost_by, cost_ax, cost_ay
    INTEGER :: it,i, nbpts
    
    IF(it.EQ.istart+1)THEN
      tempcost = 0.d0
    ENDIF
    cost_gam=0.d0
    cost_bx=0.d0
    cost_by=0.d0
    cost_ax=0.d0
    cost_ay=0.d0
    nbpts=0

    
    ! Find how many body points there are in total
    DO i=1,n_body 
      nbpts=nbpts + bdy(i)%npts       
    END DO
  
    ! Get all incremental costs 
     
    ! state cost
    cost_gam = cost_gam + CGAM*qTq(q(:,1),q(:,1),mmin,mmax,nmin,nmax)  

    ! x and y body immersed body force cost
    DO i=1,nbpts
      cost_bx = cost_bx + CBX*fb(f(i,1))**2.d0
      cost_by = cost_by + CBY*fb(f(i,2))**2.d0
    ENDDO
    
    ! x and y actuator force cost
    DO i=nbpts+1,nb
      cost_ax = cost_ax + CAX*fb(f(i,1))**2.d0
      cost_ay = cost_ay + CAY*fb(f(i,2))**2.d0
    ENDDO
    
    ! Add up all the contributions
    tempcost = tempcost + cost_gam + cost_bx + cost_by + cost_ax + cost_ay

    END SUBROUTINE integ_cost
 
  ! *****************************************************************************************
  SUBROUTINE calculate_totalcost(totalcost)
  
    !************************************************************************!
    !*    Adds up all the contributions from the cost funtional to get      *!
    !*    its final value for the last forward simulation that was run      *!    
    !************************************************************************!
  
    REAL(KIND(0.D0)) :: totalcost, controlcost
    INTEGER :: i,j
     
    ! Add cost from the forward simulation integration
    totalcost = tempcost
    WRITE(*,*)  '*** COST SUMMARY ************'
    WRITE(*,*)  'Cost from integration'
    WRITE(*,*)  tempcost/2*dt
    ! Add control Costs
    controlcost=0.d0
    DO j=1,2*nbodyf+2*nact+3
      DO i=istart+1,istop  
        controlcost= controlcost +CU(j)*control(i,j)**2.d0
        totalcost = totalcost + CU(j)*control(i,j)**2.d0
      ENDDO
    ENDDO
    
    WRITE(*,*)  'Control Cost'
    WRITE(*,*)  controlcost /2*dt
    
    ! Add dt factor to complete integration
    totalcost = totalcost * dt

    ! Add final state cost
    totalcost = totalcost + R*qTq(q(:,1),q(:,1),mmin,mmax,nmin,nmax)
    
    WRITE(*,*)  'Final State Cost'
    WRITE(*,*)  0.5*R*qTq(q(:,1),q(:,1),mmin,mmax,nmin,nmax)
    
    ! Divide total by 2
    totalcost = totalcost*0.5d0
    
    WRITE(*,*)  'Total Cost'
    WRITE(*,*)  totalcost
    WRITE(*,*)  '*****************************'
    
  END SUBROUTINE calculate_totalcost
 
  ! *****************************************************************************************
   SUBROUTINE calculate_rhs_adj
   
    !************************************************************************!
    !*    Computes the RHS "slip" velocity terms required in the            *!
    !*    adjoint equations at the body and actuator points                 *!    
    !************************************************************************!
 
  
    INTEGER :: i, nbpts
    
    nbpts = 0
    
    ! Find how many body points there are in total
    DO i=1,n_body 
      nbpts=nbpts + bdy(i)%npts       
    END DO

    ! zero rhs_adj
    rhs_adj =0.d0
    
    ! x and y body point contribution
    DO i=1,nbpts
      rhs_adj(i) = CBX*fb_base(f(i,1))
      rhs_adj(nb+i) = CBY*fb_base(f(i,2))
    ENDDO

    ! x and y actuator contribution
    DO i=nbpts+1,nb
      rhs_adj(i) = CAX*fb_base(f(i,1))
      rhs_adj(nb+i) = CAY*fb_base(f(i,2))
    ENDDO
   END SUBROUTINE calculate_rhs_adj  
    

  ! *****************************************************************************************  
  FUNCTION qTq(q1,q2,xmin,xmax,ymin,ymax)

    !************************************************************************!
    !*    Computes the matrix multiplication (dot product) (q1)^T(q2)       *!
    !*    where q1 and q2 are two fluxes defined over mgridlev grids        *!
    !************************************************************************!  

  REAL(KIND(0.D0)), DIMENSION(nq) :: q1,q2
  REAL(KIND(0.D0)) :: qTq
  INTEGER i,j,k,xmin,xmax,ymin,ymax

  qTq = 0.d0
  ! x fluxes
  DO i=xmin,xmax+1
    DO j=ymin,ymax
      qTq = qTq + q1(u(i,j))*q2(u(i,j))
    ENDDO
  ENDDO
  
  ! y fluxes
  DO i=xmin,xmax
    DO j=ymin,ymax+1
      qTq = qTq + q1(v(i,j))*q2(v(i,j))
    ENDDO
  ENDDO
  
! -------------------------------------------------------------------------!
! The following commented lines are to be used if one wants to sum the state
! cost over all the grids (would need updating to be used)
! -------------------------------------------------------------------------!
!  ! Dot product of entire region
!  DO k=1,mgridlev
!    DO i=1,nq
!      qTq = qTq + q1(i,k)*q2(i,k)/2.d0**(k-1)
!    ENDDO
!  ENDDO
!  
!  ! Now remove overlapping region  
!  DO k=2,mgridlev
!    DO i=1,m+1
!      DO j=1,n
!      
!        IF(u(i,j).ge.m/4+1.and.u(i,j).le.3*m/4)THEN
!          qTq = qTq - q1(u(i,j),k)*q2(u(i,j),k)/2.d0**(k-1)
!        ENDIF
!      ENDDO
!    ENDDO
!        
!    DO i=1,m
!      DO j=1,n+1
!        IF(v(i,j).ge.n/4+1.and.v(i,j).le.3*n/4)THEN
!          qTq = qTq - q1(v(i,j),k)*q2(v(i,j),k)/2.d0**(k-1)
!        ENDIF     
!      ENDDO
!    ENDDO
!    
!  ENDDO
  
  END FUNCTION qTq

  ! *****************************************************************************************
   FUNCTION gaussianforce(x, y, xc, yc)
   
    !************************************************************************!
    !*    creates a 2D Gaussian distribution of amplitude 1, centred at     *!
    !*    xc, yc, and with width defined by parameter a. It is used as a    *!
    !*    waveform for the control body force which multiplies it by an     *!
    !*    amplitude determined by the array control                         *!
    !************************************************************************!  

   REAL(KIND(0.D0)) :: x,y,u,v, gaussianforce
   REAL(KIND(0.D0)) :: a, xc, yc

  
   gaussianforce = 0.d0
   a = 0.15D0/2.64665D0
   ! a is a constant that determines the width of the gaussian and should be chosen to make sure that the Gaussian is smooth enough
   
   IF ( ( (x-xc)**2.d0+(y-yc)**2.d0 ) .le. (10.D0*delta)**2.d0 ) THEN  ! Only > 0 in a radius < 10 cells around centre.
      gaussianforce = EXP(-a*((x-xc)/delta)**2.D0) ! x-dir function
      gaussianforce = gaussianforce*EXP(-a*((y-yc)/delta)**2.D0) ! y-dir function
   END IF
   
 END FUNCTION gaussianforce
 
   ! ***************************************************************************************** 
   
  SUBROUTINE update_control(distance)
   
    !************************************************************************!
    !*    Updates the values of the control for the next iteration          *!
    !************************************************************************!
    INTEGER :: i,j
    REAL(KIND(0.D0)) :: distance
    
    WRITE(*,*) 'Updating the following control(s):'
    DO j=1,2*nbodyf+2*nact+3
      IF(CU(j).ne.0.d0)THEN
      WRITE(*,*) j
        DO i=istart+1,istop
         control(i,j) = control_old(i,j) - distance * gradient(i,j)
        END DO
      ENDIF
    ENDDO
  END SUBROUTINE update_control  
    

  ! *****************************************************************************************  
  

  
END MODULE optimization
