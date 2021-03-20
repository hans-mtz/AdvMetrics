! include 'likelihood.f90'
! include 'io.f90'

PROGRAM MAIN
  USE GLOBVAR
  USE INIT
  USE MATRIX
  ! USE SIMPLEX
  use wmod
  ! USE LIKELIHOOD
  ! USE MINIMIZATION
  ! USE ANNEALING, ONLY : SIMULATED_ANNEALING
  use io
  IMPLICIT NONE
  REAL(8) :: theta(N_covariates), theta_i(N_covariates),theta_s(N_covariates),&
             theta_is(N_covariates),theta_ss(N_covariates), results(5,3)
                !,std_e(N_covariates),&
               !t_hat(N_covariates),se_hat(N_covariates),t_bs_u(N_covariates),se_bs_u(N_covariates)
  ! INTEGER :: i
  integer, PARAMETER :: maxiter=300
  REAL(8), PARAMETER :: tol=1.0d-6 !theta_bs(bs,N_covariates), std_e_bs(bs,N_covariates),
  ! real(8) :: dist!,b_h(N_covariates),v !, nor(N_individuals,s)

  ! boost=.FALSE.
  ! bs=100
  ! Reading in the data
  CALL READDATA()

  call initiate()

  b_h=OLS(x,d)
  ! print*, "OLS_0 = ",b_h

  ! vcov=varcovar(x,d)
  ! print*,"Final OLS VCOV = "
  ! print*,vcov
  ! call write_file("OLS_vcov.txt",vcov)

  ! print*,
  ! print*,"Normal draws = ", sum(nor, DIM=1)
  ! print*,"Uniform draws = ", sum(uni, DIM=1)
  ! print*,"Starting varcovar"
  !
  call vcov_bs(x,d,vcov)

  ! call read_file("sigma_inv.txt",vcov)

  ! call READDATA()

  ! print*,"x = ", sum(x,dim=1)
  ! Read in, initial value for parameters

  ! ! *****
  call read_file("point.out",theta)

  print*, "theta from file = ",theta

  theta_i=rec_msle(theta,tol,maxiter,0,1)
  print*, "MLSE theta (indicator)= ",theta_i

  ! print*, "Is theta the same after MSLE indicator?", theta

  ! b_h=OLS(x,D)
  theta_s=rec_msle(theta,tol,maxiter,0,2)
  print*, "MLSE theta (smooth)= ",theta_s

  theta_is=rec_msle(theta,tol,maxiter,0,3)
  print*, "MLSE theta (indicator) with weights= ",theta_is

  ! print*, "Is theta the same after MSLE indicator?", theta

  ! b_h=OLS(x,D)
  theta_ss=rec_msle(theta,tol,maxiter,0,4)
  print*, "MLSE theta (smooth) with weights= ",theta_ss
  !
  ! print*, "OLS = ",b_h
  !
  ! ! theta_i=theta
  ! ! theta_s=theta
  ! ! theta_is=theta
  ! ! theta_ss=theta
  ! !
  ! ! smooth=.FALSE.
  ! ! print*,"theta_i",theta_i
  ! ! call Nelder_Meade(theta_i,tol,objf,1,maxiter)
  ! ! print*, "MLSE theta (indicator) NM= ",theta_i
  ! !
  ! ! ! CALL READDATA()
  ! !
  ! ! smooth=.TRUE.
  ! ! print*,"theta_s",theta_i
  ! ! call BFGS(theta_s,tol,tol,objf,1,.TRUE.)
  ! ! print*, "MLSE theta (smooth) BFGS= ",theta_s
  !
  ! ! CALL READDATA()
  ! !
  ! ! smooth=.FALSE.
  ! ! call Nelder_Meade(theta_is,tol,objfvcov,1,maxiter)
  ! ! print*, "MLSE theta (indicator) NM W = ",theta_is
  ! !
  ! ! CALL READDATA()
  ! !
  ! ! smooth=.TRUE.
  ! ! call BFGS(theta_ss,tol,tol,objfvcov,1,.TRUE.)
  ! ! print*, "MLSE theta (smooth) BFGS W= ",theta_ss
  !
  !
  results=reshape([b_h,theta_i,theta_s,theta_is,theta_ss],[5,3],order=[2,1])
  call write_file("results.txt",results)
  ! !******


END PROGRAM MAIN
