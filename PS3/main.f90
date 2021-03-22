! include 'likelihood.f90'
! include 'io.f90'
!## My program -----------------
PROGRAM MAIN
  USE GLOBVAR
  USE INIT
  USE MATRIX
  USE SIMPLEX
  use wmod
  ! USE LIKELIHOOD
  USE MINIMIZATION
  ! USE ANNEALING, ONLY : SIMULATED_ANNEALING
  use io
  IMPLICIT NONE
  REAL(8) :: theta(N_covariates), theta_i(N_covariates),theta_s(N_covariates),&
             theta_is(N_covariates),theta_ss(N_covariates), results(4,3),&
             theta_bfgs(N_covariates)
             !beta(N_covariates),b_h2(N_covariates)!,&
               !t_hat(N_covariates),se_hat(N_covariates),t_bs_u(N_covariates),se_bs_u(N_covariates)
  ! INTEGER :: i
  integer, PARAMETER :: maxiter=300
  REAL(8), PARAMETER :: tol=1.0d-9 !theta_bs(bs,N_covariates), std_e_bs(bs,N_covariates),
  ! real(8) :: dist!,b_h(N_covariates),v !, nor(N_individuals,s)

  ! boost=.FALSE.
  ! bs=100
  ! Reading in the data
  CALL READDATA()

  call initiate()

  ! print*,uni(1:10,3:5)
  ! print*,nor(1:10,3:5)
  ! print*,nor_v(1:10)

  b_h=OLS(x,d)
  print*, "OLS_0 = ",b_h
  ! print*, "x",x(5:8,:)

  ! vcov=varcovar(x,d)
  ! call bs_m(x,d,varbeta,vcov)
  call bs_m(x,d,varcovar,vcov)
  ! vcov=varbeta(x,d)
  ! print*,"VCOV = "
  ! print*,vcov
  ! call write_file("OLS_vcov.txt",vcov)

  ! print*,
  ! print*,"Normal draws = ", sum(nor, DIM=1)
  ! print*,"Uniform draws = ", sum(uni, DIM=1)
  ! print*,"Starting varcovar"
  !
  ! call vcov_bs(x,d,vcov)

  ! call read_file("sigma_inv.txt",vcov)

  ! call READDATA()

  ! print*,"x = ", sum(x,dim=1)
  ! Read in, initial value for parameters
  ! smooth=.FALSE.
  ! ! *****
  !
  ! call read_file("point.out",theta)
  ! call read_file("pointBFGS.out",theta_s)
  ! smooth=0
  ! theta_i=b_star(theta,smooth)
  ! print*,theta_i
  ! ! call READDATA()
  ! smooth=1
  ! theta_s=b_star(theta,smooth)
  ! print*,theta_s
  ! print*, "theta from file = ",theta
  ! call READDATA()
  ! !
  ! ! b_h2=OLS(x,d)
  ! print*, "OBJ = ",objf(theta)
  ! print*, "OBJ2 = ",objf(theta_i)
  ! print*, "OBJ3 = ",objf(theta_s)
  ! print*, "x",x(5:8,:)


  ! !********


  !
  ! ! call b_star(theta_i,smooth)
  !
  ! ! print*, "Same theta after b_star?", theta_i
  !
  ! ! b_h=OLS(x,D)
  ! smooth=.TRUE.
  ! theta_s=rec_msle(theta,tol,maxiter,0,.TRUE.,2)
  ! print*, "MLSE theta (smooth)= ",theta_s
  !
  ! smooth=.FALSE.
  ! theta_is=rec_msle(theta,tol,maxiter,0,.FALSE.,3)
  ! print*, "MLSE theta (indicator) with weights= ",theta_is
  !
  ! ! print*, "Is theta the same after MSLE indicator?", theta
  !
  ! ! b_h=OLS(x,D)
  ! smooth=.TRUE.
  ! theta_ss=rec_msle(theta,tol,maxiter,0,.TRUE.,4)
  ! print*, "MLSE theta (smooth) with weights= ",theta_ss
  ! !
  !*********
  ! print*, "OLS = ",b_h
  call read_file("point.out",theta)
  call read_file("pointBFGS.out",theta_bfgs)
  print*,"theta_NM",theta
  print*,"theta_BFGS",theta_bfgs
  theta_i=theta
  theta_s=theta_bfgs
  theta_is=theta
  theta_ss=theta_bfgs

  smooth=0
  ind=1
  ! print*,"theta_i",theta_i
  call Nelder_Meade(theta_i,tol,objf,1,maxiter)
  print*, "MLSE theta (indicator) NM (sampling u=1)= ",theta_i
  ! beta=b_star(theta_i,smooth)
  ! print*,beta

  ! theta_s=theta_i
  smooth=1
  ind=1
  ! print*,"theta_s",theta_s
  call BFGS(theta_s,tol,tol,objf,1,.TRUE.)
  print*, "MLSE theta (smooth) BFGS (sampling u=1)= ",theta_s
  ! beta=b_star(theta_i,smooth)
  ! print*,beta

  smooth=0
  ind=1
  call Nelder_Meade(theta_is,tol,objfvcov,1,maxiter)
  print*, "MLSE theta (indicator) NM W  (sampling u=1)= ",theta_is
  ! beta=b_star(theta_i,smooth)
  ! print*,beta

  smooth=1
  ind=1
  call BFGS(theta_ss,tol,tol,objfvcov,1,.TRUE.)
  print*, "MLSE theta (smooth) BFGS W (sampling u=1)= ",theta_ss
  ! beta=b_star(theta_i,smooth)
  ! print*,beta

  results=reshape([theta_i,theta_s,theta_is,theta_ss],[4,3],order=[2,1])
  ! results=reshape([b_h,theta_i,theta_s],[5,3],order=[2,1])
  call write_file("results.txt",results)
  ! !******
  !
  ! call read_file("point.out",theta)
  ! call read_file("pointBFGS.out",theta_s)
  theta_i=theta
  theta_s=theta_bfgs
  theta_is=theta
  theta_ss=theta_bfgs

  smooth=0
  ind=0
  ! print*,"theta_i",theta_i
  call Nelder_Meade(theta_i,tol,objf,1,maxiter)
  print*, "MLSE theta (indicator) NM (sampling u=100)= ",theta_i
  ! beta=b_star(theta_i,smooth)
  ! print*,beta

  ! theta_s=theta_i
  smooth=1
  ind=0
  ! print*,"theta_s",theta_s
  call BFGS(theta_s,tol,tol,objf,1,.TRUE.)
  print*, "MLSE theta (smooth) BFGS (sampling u=100)= ",theta_s
  ! beta=b_star(theta_i,smooth)
  ! print*,beta

  smooth=0
  ind=0
  call Nelder_Meade(theta_is,tol,objfvcov,1,maxiter)
  print*, "MLSE theta (indicator) NM W  (sampling u=100)= ",theta_is
  ! beta=b_star(theta_i,smooth)
  ! print*,beta

  smooth=1
  ind=0
  call BFGS(theta_ss,tol,tol,objfvcov,1,.TRUE.)
  print*, "MLSE theta (smooth) BFGS W (sampling u=100)= ",theta_ss
  ! beta=b_star(theta_i,smooth)
  ! print*,beta

  results=reshape([theta_i,theta_s,theta_is,theta_ss],[4,3],order=[2,1])
  ! results=reshape([b_h,theta_i,theta_s],[5,3],order=[2,1])
  call write_file("results_s100.txt",results)
  ! !******


END PROGRAM MAIN
