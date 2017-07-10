MODULE parameters

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
  
  
  IMPLICIT NONE
  
  ! parameters
  INTEGER :: istart                ! initial time index
  INTEGER :: istop                ! last time index
  INTEGER :: isave                 ! save a restart every isave steps
  INTEGER :: m                    ! cells in x
  INTEGER :: n                   ! cells in y
  REAL(KIND(0.0D0)) :: dt     ! time step
  REAL(KIND(0.0D0)) :: Re      ! Reynolds number
  REAL(KIND(0.0D0)) :: cgtol   ! tol. for cg convergence (poission eq)
  INTEGER :: cg_max_iter           ! max. iterations for any cg iteration
  INTEGER :: n_actuator                 ! number of actuators 
  INTEGER :: n_body          ! number of moving bodies
                                        ! (will be counted in grid.f90) 
  REAL(KIND(0.D0)) :: len   ! length scale for grid
  REAL(KIND(0.D0)) :: offsetx   ! offset for grid in x
  REAL(KIND(0.D0)) :: offsety    ! offset for grid in y
  INTEGER :: mgridlev
  REAL(KIND(0.D0)) :: pi

  REAL(KIND(0.D0)) :: rot_angle  ! rotating angle of the grid
  REAL(KIND(0.D0)) :: rox  ! x-coord. of center of rotation rotation
  REAL(KIND(0.D0)) :: roy  ! y-coord. of center of rotation rotation

  LOGICAL :: stationary       ! all stationary bodies w.r.t the grid
  LOGICAL :: compute_pressure  ! whether to output pressure

    real(kind(0.d0)) :: R_rho
    real(kind(0.d0)) :: R_E
    real(kind(0.d0)) :: R_sh
    real(kind(0.d0)) :: R_th

    logical :: standard
    logical :: pinned

    integer, dimension(3) :: bc_type
    real(kind(0.d0)), dimension(3) :: bc_val

  LOGICAL :: linear           ! linearized simulation?
  LOGICAL :: steady           ! if linear simulation: steady?
  INTEGER :: ibase            ! base flow field is read every ibase iterations
  
  LOGICAL :: adjoint          ! adjoint optimization?
  INTEGER :: iopt_start       ! iteration to start with
  INTEGER :: iopt_stop        ! max iteration number
  INTEGER :: nact             ! number of actuators in cost function
  INTEGER :: nbodyf           ! number of body forces in cost function
  REAL(KIND(0.D0)),DIMENSION(:),ALLOCATABLE :: CINIT,CU        ! Control default values and weights
  REAL(KIND(0.D0)) :: CGAM               ! Weight parameter for state cost
  REAL(KIND(0.D0)) :: CBX                ! Weight parameter for x immersed body force cost
  REAL(KIND(0.D0)) :: CBY                ! Weight parameter for y immersed body force cost
  REAL(KIND(0.D0)) :: CAX                ! Weight parameter for x actuator force cost
  REAL(KIND(0.D0)) :: CAY                ! Weight parameter for y actuator force cost
  REAL(KIND(0.D0)) :: R                  ! Weight parameter for final state cost
  REAL(KIND(0.D0)),DIMENSION(:),ALLOCATABLE :: XCTRL,YCTRL,GX,GY              ! Coordinates of the control region for the state
  INTEGER :: mmin,mmax,nmin,nmax
  
CONTAINS
  
  SUBROUTINE input

    LOGICAL :: readinput,readopt
    INTEGER :: i,j,nacttemp,nbodyftemp
    REAL(KIND(0.D0)) :: dx,xx,yy,xxm1,yym1
    REAL(KIND(0.D0)),DIMENSION(2) :: XCREG,YCREG
    REAL(KIND(0.D0)),DIMENSION(3) :: CINITU,CUU
    REAL(KIND(0.D0)),DIMENSION(100) :: CINITGX,CINITGY,CINITAX,CINITAY    ! Default steady values for the controls      
    REAL(KIND(0.D0)),DIMENSION(100) :: CUGX,CUGY,CUAX,CUAY       
    REAL(KIND(0.D0)),DIMENSION(100) :: XG,YG             
   
    NAMELIST /read_parameters/  istart,istop,isave,m,n,dt,Re,cgtol,              &
                                cg_max_iter,len,offsetx,offsety,                &
                                mgridlev, compute_pressure, R_rho, &
                                R_E, R_sh, R_th, standard, pinned,   &
                                linear, steady, ibase, adjoint, nact, nbodyf
                                
    NAMELIST /read_optimization/ iopt_start,iopt_stop,                          &                         
                                 CINITGX,CINITGY,CINITAX,CINITAY,CINITU,        &
                                 CGAM,R,XCREG,YCREG,                            &
                                 CBX,CBY,CAX,CAY,                               &
                                 CUGX,CUGY,XG,YG,CUAX,CUAY,CUU                                                   
    pi  = 4.0d0*atan(1.0d0)
                                 
    ! read input
    INQUIRE(file='input/ib.inp',exist=readinput)

    IF (readinput) THEN
       OPEN(unit=3,file='input/ib.inp',form='formatted',status='old')
       READ(unit=3,nml=read_parameters)
       CLOSE(3)
       WRITE(*,*) 'read input file'
      ELSE
       STOP 'cannot find input file'
    END IF


    if (pinned) then
        bc_type = (/ 1, 1, 0 /)
    ELSE
        bc_type = (/ 1, 1, 1 /)
    end if

        bc_val = (/ 0.d0, 0.d0, 0.d0 /)


    IF(ADJOINT)THEN
      ALLOCATE(XCTRL(2),YCTRL(2),GX(nbodyf),GY(nbodyf))
      ALLOCATE(CU(2*nbodyf+2*nact+3),CINIT(2*nbodyf+2*nact+3))

      
      ! read optimization input
      INQUIRE(file='input/opt.inp',exist=readopt)
      IF (readopt) THEN
         OPEN(unit=3,file='input/opt.inp',form='formatted',status='old')
         READ(unit=3,nml=read_optimization)
         CLOSE(3)    
         XCTRL=XCREG
         YCTRL=YCREG
         DO i=1,nbodyf            
           GX(i)=XG(i)           
           GY(i)=YG(i)             
         ENDDO        
         ! make the CU and CINIT vectors from the values read
         DO i=1,nbodyf
           CU(i)=CUGX(i)
           CU(nbodyf+i)=CUGY(i)
           CINIT(i)=CINITGX(i)
           CINIT(nbodyf+i)=CINITGY(i)
         ENDDO
         
         DO i=1,nact
           CU(2*nbodyf+i)=CUAX(i)
           CU(2*nbodyf+nact+i)=CUAY(i)
           CINIT(2*nbodyf+i)=CINITAX(i)
           CINIT(2*nbodyf+nact+i)=CINITAY(i)
         ENDDO 
         
         CU(2*nact+2*nbodyf+1)= CUU(1) !X Translation
         CU(2*nact+2*nbodyf+2)= CUU(2) !Y Translation
         CU(2*nact+2*nbodyf+3)= CUU(3) !Rotation
         
         CINIT(2*nact+2*nbodyf+1)=CINITU(1)
         CINIT(2*nact+2*nbodyf+2)=CINITU(2)
         CINIT(2*nact+2*nbodyf+3)=CINITU(3)
         
         WRITE(*,*) '-------- Cost function information: -----------------------'
         WRITE(*,*) 'CGAM, R =',CGAM,R
         WRITE(*,*) 'XCTRL   =',xctrl
         WRITE(*,*) 'YCTRL   =',yctrl
         WRITE(*,*) 'CBX,CBY =',CBX,CBY
         WRITE(*,*) 'CAX,CAY =',CAX,CAY
         WRITE(*,*) 'CUGX    =',(CU(i),i=1,nbodyf)
         WRITE(*,*) 'CUGY    =',(CU(i+nbodyf),i=1,nbodyf)
         WRITE(*,*) 'GX      =',GX
         WRITE(*,*) 'GY      =',GY
         WRITE(*,*) 'CUAX    =',(CU(2*nbodyf+i),i=1,nact)
         WRITE(*,*) 'CUAY    =',(CU(2*nbodyf+nact+i),i=1,nact)                                                              
         WRITE(*,*) 'CUU     =',(CU(2*nbodyf+2*nact+i),i=1,3) 
         WRITE(*,*) '-------- Default Control Values: --------------------------'
         WRITE(*,*) 'CINIT   =',CINIT
         WRITE(*,*) '-------- Lumped Control Weight Vector: --------------------'
         WRITE(*,*) 'CU      =',CU
         
         ! Find the indices corresponding to the control region
         dx=len/real(m)
         mmin=0
         mmax=0
         nmin=0
         nmax=0
         
         
         WRITE(*,*) '-------- Control Region: ----------------------------------'
         ! X limits 
         IF(xctrl(1).gt.xctrl(2).or. &
            xctrl(1).lt.-offsetx.or. &
            xctrl(1).gt.m*dx-offsetx.or. &
            xctrl(2).lt.-offsetx.or. &
            xctrl(2).gt.m*dx-offsetx)THEN        
            WRITE(*,*) 'Using Default x control region: full domain'
            mmin=1
            mmax=m
          ELSE
            DO i=2,m
              xx = (REAL(i)-1.d0)*dx - offsetx
              xxm1 = (REAL(i)-2.d0)*dx - offsetx    
              IF(xx.ge.xctrl(1).and.xxm1.le.xctrl(1))THEN
                mmin=i-1
              ENDIF
              IF(xx.ge.xctrl(2).and.xxm1.le.xctrl(2))THEN
                mmax=i
              ENDIF
            ENDDO
          ENDIF

          WRITE(*,*) 'X Range =', (REAL(mmin)-1.d0)*dx - offsetx, (REAL(mmax)-1.d0)*dx - offsetx
    
          ! Y limits 
          IF(yctrl(1).gt.yctrl(2).or. &
             yctrl(1).lt.-offsety.or. &
             yctrl(1).gt.n*dx-offsety.or. &
             yctrl(2).lt.-offsety.or. &
             yctrl(2).gt.n*dx-offsety)THEN        
             WRITE(*,*) 'Using Default y control region: full domain'
             nmin=1
             nmax=n
          ELSE
            DO i=2,n
             yy = (REAL(i)-1.d0)*dx - offsety
             yym1 = (REAL(i)-2.d0)*dx - offsety
             IF(yy.ge.yctrl(1).and.yym1.le.yctrl(1))THEN
               nmin=i-1
             ENDIF
             IF(yy.ge.yctrl(2).and.yym1.le.yctrl(2))THEN
               nmax=i
             ENDIF
           ENDDO
         ENDIF
         WRITE(*,*) 'Y Range =', (REAL(nmin)-1.d0)*dx - offsety, (REAL(nmax)-1.d0)*dx - offsety 
         
         WRITE(*,*) 'read optimization input file'
      ELSE
         STOP 'cannot find optimization input file'
      END IF !endif readopt
         WRITE(*,*) '-----------------------------------------------------------'
    ENDIF !endif adjoint
    

  END SUBROUTINE input

END MODULE parameters


