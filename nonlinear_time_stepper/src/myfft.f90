MODULE myfft

   IMPLICIT NONE

   INCLUDE 'fftw3.f'           ! needed for defining the plan

   REAL(KIND(0.D0)), DIMENSION(:),     ALLOCATABLE :: viscfac,vfac,con1,con2
   REAL(KIND(0.D0)), DIMENSION(:,:),   ALLOCATABLE :: in, laminv
   REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: lam1, lam1i, lam1inv
   REAL(KIND(0.D0)), DIMENSION(:,:),   ALLOCATABLE :: in_ddti, laminv_ddti 
!   REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: intfac1,intfac2,intfac3

    !!!!NCP stuff
    complex(kind(0.d0)), dimension(:), allocatable :: ncp_out
    real(kind(0.d0)), dimension(:), allocatable :: ncp_in
    integer*8 :: ncp_plan
    integer :: nf_div_2

    !!!


   INTEGER*8 :: forward, inverse, forward_ddti
    INTEGER :: mm, nn
   REAL*8  :: normalize

CONTAINS

  SUBROUTINE setup_fft

    USE grid
    USE parameters
    
    INTEGER                :: i,j,k
    REAL*8                 :: del2, del22, normalize_ddti
    real*8, dimension(m,n) :: lam_ddti
    real*8, dimension(m-1,n-1) :: lam
    
    mm = m
    nn = n
    ALLOCATE( viscfac(mgridlev), con1(mgridlev), con2(mgridlev) )
    ALLOCATE( vfac(mgridlev) )
    ALLOCATE( in(mm-1,nn-1), laminv(mm-1,nn-1) )
    !              lam(mm-1,nn-1) )
    ALLOCATE( lam1(mm-1,nn-1,mgridlev), &
              lam1i(mm-1,nn-1,mgridlev), &
              lam1inv(mm-1,nn-1,mgridlev) )
!    ALLOCATE( intfac1(mm-1,nn-1,mgridlev), intfac2(mm-1,nn-1,mgridlev), intfac3(mm-1,nn-1,mgridlev)  )	
    ALLOCATE( in_ddti(mm,nn), laminv_ddti(mm,nn)  )          

    CALL dfftw_plan_r2r_2d(forward, mm-1,nn-1,in,in, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)
    CALL dfftw_plan_r2r_2d(forward_ddti, mm,nn,in_ddti,in_ddti, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)
    
    normalize_ddti = 4.d0*REAL( (mm+1)*(nn+1) )                         
    normalize      = 4.d0*REAL(     mm*nn     )    

    ! eigenvalues for inverse of C^T C
    del2 = delta*delta
    DO k=1,mgridlev
       del22      =  del2*4.d0**(k-1)
       viscfac(k) =  dt/Re/del22
       vfac(k)    =  0.5d0*dt/Re/del22/normalize
       con1(k)    =  1.5d0*dt/del22/normalize
       con2(k)    = -0.5d0*dt/del22/normalize
    ENDDO

    DO j=1,nn-1
       DO i=1,mm-1
          ! turn lam into local variable, because it is only used here
          lam(i,j) = -2.d0*( COS( pi*REAL(i)/REAL(mm) ) + &
                             COS( pi*REAL(j)/REAL(nn) ) - 2.d0 )
          laminv(i,j) = 1.d0/lam(i,j)/normalize
          ! integrating factors for viscous terms
          DO k=1,mgridlev
             del22 = del2* 4.d0**(k-1)
             lam1(i,j,k)    =      (1.d0 - 0.5d0*dt*lam(i,j)/del22/Re)/normalize ! original
             lam1i(i,j,k)   = 1.d0/(1.d0 + 0.5d0*dt*lam(i,j)/del22/Re) ! original

! it appears that this HUGE array is NEVER used anywhere, so why compute/store it????
!             lam1inv(i,j,k) = 1.d0/(1.d0 - 0.5d0*dt*lam(i,j)/del22/Re)

!computed, but never USED !!????
!             intfac2(i,j,k) =  1.5*dt* EXP(      -dt*lam(i,j) / del22 / Re )/del22/normalize
!             intfac3(i,j,k) = -0.5*dt* EXP( -2.d0*dt*lam(i,j) / del22 / Re )/del22/normalize
!             intfac1(i,j,k) =          EXP(      -dt*lam(i,j) / del22 / Re )      /normalize
          END DO
      END DO
    END DO
    DO j=1,nn                                                          
       DO i=1,mm 
          lam_ddti(i,j)    = 2.d0*( COS( pi*REAL(i)/REAL(mm+1) ) + &          
                                    COS( pi*REAL(j)/REAL(nn+1) ) - 2.d0 )  
          laminv_ddti(i,j) = 1.d0/lam_ddti(i,j)/normalize_ddti       
       END DO                                                          
    END DO


!!!! NCP stuff

    nf_div_2 = floor(nf/2.d0) + 1

    !print *, "nf/2, fl(nf/2), fl(nf/2) + 1, nf_div_2 ", nf/2.d0, floor(nf/2.d0), floor(nf/2.d0) + 1, nf_div_2
    allocate(ncp_in(nf),ncp_out(nf_div_2))

    call dfftw_plan_dft_r2c_1d(ncp_plan, nf, ncp_in, ncp_out, FFTW_EXHAUSTIVE)


!!!!

  END SUBROUTINE setup_fft
! ********************************************************************************** !
  SUBROUTINE destroy_fft

    CALL dfftw_destroy_plan(forward)
    DEALLOCATE( in, laminv, viscfac,vfac,con1,con2 )
!    deallocate( intfac1,intfac2,intfac3 )
    DEALLOCATE( lam1, lam1i, lam1inv )

    CALL dfftw_destroy_plan(forward_ddti)
    DEALLOCATE( in_ddti, laminv_ddti )

    call dfftw_destroy_plan(ncp_plan)
    deallocate( ncp_in, ncp_out )

  END SUBROUTINE destroy_fft
! ********************************************************************************** !
  FUNCTION dst( psi)

    REAL(KIND(0.D0)), DIMENSION(:,:) :: psi
    REAL(KIND(0.D0)), DIMENSION(2:mm,2:nn) :: dst

    ! discrete sine transform
    ! careful...two calls of dst need to be divided by "normalize"
    ! to return original vector

    in  = psi
    CALL dfftw_execute_r2r(forward,in,in)
    dst = in

  END FUNCTION dst
! ********************************************************************************** !
  FUNCTION ctci( omega ) 

    REAL(KIND(0.D0)), DIMENSION(:,:) :: omega
    REAL(KIND(0.D0)), DIMENSION(2:mm,2:nn) :: ctci
        
    in =  omega
    CALL dfftw_execute_r2r(forward,in,in)
    in = laminv * in
    CALL dfftw_execute_r2r(forward,in,in)
    ctci =  in
    
  END  FUNCTION ctci
 ! ********************************************************************************** ! 
  FUNCTION ainv( omega ) 

    REAL(KIND(0.D0)), DIMENSION(:,:)       :: omega
    REAL(KIND(0.D0)), DIMENSION(2:mm,2:nn) :: ainv
    
    in = omega
    CALL dfftw_execute_r2r(forward,in,in)
    in = lam1i(:,:,1) * in / normalize
    CALL dfftw_execute_r2r(forward,in,in)
    ainv = in

  END  FUNCTION ainv
! ********************************************************************************** !
  FUNCTION ddti( phi )                  

    REAL(KIND(0.D0)), DIMENSION(:,:) :: phi
    REAL(KIND(0.D0)), DIMENSION(1:mm,1:nn) :: ddti
        
    in_ddti = phi                    
    CALL dfftw_execute_r2r(forward_ddti,in_ddti,in_ddti)    
    in_ddti = laminv_ddti * in_ddti     
    CALL dfftw_execute_r2r(forward_ddti,in_ddti,in_ddti)    
    ddti = in_ddti                     
    
  END  FUNCTION ddti                   
! ********************************************************************************** !

! ********************************************************************************** !
FUNCTION dft_1d( in_var )

REAL(KIND(0.D0)), DIMENSION(:) :: in_var
complex(KIND(0.D0)), DIMENSION(1 : nf_div_2) :: dft_1d

! discrete sine transform
! careful...two calls of dst need to be divided by "normalize"
! to return original vector

ncp_in  = in_var
ncp_out = 0.d0
CALL dfftw_execute_dft_r2c(ncp_plan,ncp_in,ncp_out)

dft_1d = ncp_out

END FUNCTION dft_1d
! ********************************************************************************** !

END MODULE myfft
