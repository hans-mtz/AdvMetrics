! include 'likelihood.f90'
! include 'io.f90'

PROGRAM MAIN
  USE GLOBVAR
  USE INIT
  USE LIKELIHOOD
  USE MINIMIZATION
  USE ANNEALING, ONLY : SIMULATED_ANNEALING
  use io
  IMPLICIT NONE
  REAL(8) :: theta(N_covariates), vars(N_covariates),std_e(N_covariates),&
               t_hat(N_covariates),se_hat(N_covariates),t_bs_u(N_covariates),se_bs_u(N_covariates)
  INTEGER :: j,dummy_argument,i,v(9),m(3,3),v1(9),m2(3,3),sim
  REAL(8) :: theta_bs(bs,N_covariates), std_e_bs(bs,N_covariates)

  ! boost=.FALSE.
  ! bs=100
  ! Reading in the data
  CALL READDATA()

  ! Read in, initial value for parameters
  OPEN(1,file='point.out')
  DO j=1,N_covariates
     READ(1,fmt=*) theta(j)
  END DO
  CLOSE(1)

  ! v=[((i+j,i=1,3),j=0,2)]
  ! FORALL(i=1:3,j=5:9:2) m(i,)=i+j
  ! ! m=RESHAPE(v,[3,3],[2,1])
  !
  ! print*,m(1,:)
  ! print*,m(2,:)
  ! print*,m(3,:)

  ! m=m+m
  ! print*,m(1,:)
  ! print*,m(2,:)
  ! print*,m(3,:)
  ! print*, (("x[",trim(i),",:]= ",x(i,:),"y(i)=", y(i)), i=1,5)
  ! print*, "Theta = ",theta
  ! CALL SIMULATED_ANNEALING(1000.0d0,0.85d0,12,theta,loglikelihood,100,0)
  ! beta_hat=probit(x,theta,loglikelihood)
  CALL BFGS(theta,1.0d-6,1.0d-6,loglikelihood,0,.TRUE.)



  !
  ! ! After it is done lets transform everything back
  ! beta=theta(1:N_covariates)
  ! ! Since variance needs to be positive, instead of maximizing with respect to variance
  ! ! lets maximize with respect to a transformation. Since variance=EXP(LOG(variance)) if we call
  ! ! alpha=LOG(variance), and maximize with respect to alpha we are ok.
  ! variance=EXP(theta(N_covariates+1))

  ! variance=var_probit(x,theta)
  std_e=se(x,theta)
  !
  ! ! and write the output
  ! OPEN(2,file='output.txt')
  ! WRITE(2,fmt='(A,40F16.8)') 'beta',beta
  ! WRITE(2,fmt='(A,F16.8)') 'variance',variance
  ! CLOSE(2)

  print*, "Theta = ",theta
  print*, "S.E. = ", std_e

  call write_file("theta_hat.txt", [(theta(i),i=1,3),(std_e(j),j=1,3)])

  t_hat=theta
  se_hat=std_e
  !Bootstrap

  !resample
  !hal=.FALSE. !if .TRUE. halton sequences are used to sample in boostrap
            ! else uniform distribution is used.
   print*,"What kind of bootstrap? .TRUE. for Halton, .FALSE. for Uniform."
   read*,hal

  do j=1,bs
     theta=[0.0d0,0.0d0,0.0d0]
     sim=j
     call shuffle(sim)

     CALL BFGS(theta,1.0d-6,1.0d-6,loglikelihood,1,.TRUE.)

     forall(i=1:3 )theta_bs(j,i)=theta(i)

     ! theta_bs(j,:)=theta
     std_e=se(x,theta)
     forall(i=1:3 )std_e_bs(j,i)=std_e(i)

  end do

  !Bias corrected estimator θ⋆=2θ_hat-θ_bar⋆
  print*, "Theta (BS)= ",(2.0d0*t_hat(i)-sum(theta_bs(:,i))/DBLE(bs),i=1,3)
  print*, "S.E. (BS)= ", (2.0d0*se_hat(i)-sum(std_e_bs(:,i))/DBLE(bs),i=1,3)

  t_bs_u=[(2.0d0*t_hat(i)-sum(theta_bs(:,i))/DBLE(bs),i=1,3)]
  se_bs_u=[(2.0d0*se_hat(i)-sum(std_e_bs(:,i))/DBLE(bs),i=1,3)]

  call write_file("theta_bs_u.txt", [(t_bs_u(i),i=1,3),(se_bs_u(j),j=1,3)])
  !
  ! !vars=(/beta(1),beta(2),variance/)
  call write_file('theta_bs.txt',theta_bs)
  call write_file('std_bs.txt',std_e_bs)

END PROGRAM MAIN
