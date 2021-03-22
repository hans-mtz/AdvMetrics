MODULE GLOBVAR
  use MATRIX

  IMPLICIT NONE
  ! This modlue declares the variables that are going to be used throughout the program
  ! higher level languages do not require you to declare your variables ahead of time

  INTEGER, PARAMETER :: N_individuals=50000, &	! number of people in the dataset
                        N_covariates=3, &       ! number of covariates in the regression (X's)
                        s=100,n_boot=100!, &                 ! number of simulations
                        ! b_h(N_covariates)!=OLS
  REAL(8) ::  x(N_individuals,N_covariates),&
              y(N_individuals),&!theta_hat(N_covariates),&
              vcov(N_covariates,N_covariates),&
              D(N_individuals), nor(N_individuals,s),&
              uni(N_individuals,n_boot),b_h(N_covariates), nor_v(N_individuals)
  real(8) :: k=0.1d0
  INTEGER :: myid,numprocs,smooth,ind
  ! logical:: hal, smooth!, sigma
  ! b_h=OLS(x,D)

END MODULE GLOBVAR
