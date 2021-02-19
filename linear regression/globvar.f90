MODULE GLOBVAR
  
  IMPLICIT NONE
  ! This modlue declares the variables that are going to be used throughout the program
  ! higher level languages do not require you to declare your variables ahead of time
  
  INTEGER, PARAMETER :: N_individuals=50253, &	! number of people in the dataset
       N_covariates=2                           ! number of covariates in the regression (X's)
  REAL(8) :: x(N_individuals,N_covariates),y(N_individuals),beta(N_covariates),variance
  INTEGER :: myid,numprocs
  
END MODULE GLOBVAR
