MODULE GLOBVAR
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N_ind=10000   ! Number of individuals
  INTEGER, PARAMETER :: N_int=100    ! Number of integration points
  INTEGER, PARAMETER :: N_mix=2      ! Number of mixture components for factors
  INTEGER, PARAMETER :: N_x=2        ! Number of X's
  INTEGER, PARAMETER :: N_z=1 

  ! Mixture
  REAL(8) :: m_fac(N_mix),v_fac(N_mix),p_fac(N_mix)

  !! NOTICE THAT I AM ASSUMING THE X's ARE THE SAME IN ALL TESTS
  REAL(8) :: uni(N_ind,N_int),nor(N_ind,N_int)
  REAL(8) :: x(N_ind,N_x),Z(N_ind,N_z),M(N_ind),Y(N_ind)
  REAL(8) :: bm(N_x),b1(N_x),b0(N_x),gx(N_x),gz(N_z),vm,v1,v0
  REAL(8) :: am,a1,a0,av
  INTEGER :: d(N_ind)

  ! Variables for MPI
  INTEGER :: myid,numprocs
  
END MODULE GLOBVAR
