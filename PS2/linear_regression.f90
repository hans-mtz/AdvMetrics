PROGRAM MAIN
  USE GLOBVAR
  USE INIT
  USE LIKELIHOOD
  USE MINIMIZATION
  USE ANNEALING, ONLY : SIMULATED_ANNEALING
  use io
  IMPLICIT NONE
  REAL(8) :: theta(N_covariates+1), vars(3)
  INTEGER :: j,dummy_argument

  ! Reading in the data
  CALL READDATA()

  ! Read in, initial value for parameters
  OPEN(1,file='point.out')
  DO j=1,N_covariates+1
     READ(1,fmt=*) theta(j)
  END DO
  CLOSE(1)

  CALL SIMULATED_ANNEALING(1000.0d0,0.85d0,12,theta,loglikelihood,100,0)
  CALL BFGS(theta,1.0d-6,1.0d-6,loglikelihood,0,.TRUE.)

  ! After it is done lets transform everything back
  beta=theta(1:N_covariates)
  ! Since variance needs to be positive, instead of maximizing with respect to variance
  ! lets maximize with respect to a transformation. Since variance=EXP(LOG(variance)) if we call
  ! alpha=LOG(variance), and maximize with respect to alpha we are ok.
  variance=EXP(theta(N_covariates+1))

  ! and write the output
  OPEN(2,file='output.txt')
  WRITE(2,fmt='(A,40F16.8)') 'beta',beta
  WRITE(2,fmt='(A,F16.8)') 'variance',variance
  CLOSE(2)

  !vars=(/beta(1),beta(2),variance/)
  call write_file('out.txt',theta)

END PROGRAM MAIN
