program MAIN
    USE GLOBVAR
    USE INIT
    USE MINIMIZATION
    USE MATRIX
    USE RANDOM
    USE HW3
    USE PROBABILITY, ONLY : CDF_NORMAL,PDF_NORMAL
    use SIMPLEX
    implicit none

    REAL(8) :: theta(N_covariates), A(N_covariates,N_covariates), temp(N_covariates,N_covariates), xb
    REAL(8) :: test,theta_a(N_covariates), B_hat_BS(BS,N_covariates)
    REAL(8) :: inde, N_individuals_R, theta_BS(N_covariates), BS_data(BS,N_covariates),SE(N_covariates)
    REAL(8) :: BS_se(N_covariates),theta_BS_avg(N_covariates),theta_b(N_covariates),theta_c(N_covariates),theta_d(N_covariates)
    REAL(8), PARAMETER :: tol=1.0d-8
    integer, PARAMETER :: max_NM=300
    INTEGER :: j,i,k

    ! Reading in the data
    CALL READDATA()

    ! Read in, initial value for parameters
    OPEN(1,file='point1.out')
    DO j=1,N_covariates
       READ(1,fmt=*) theta(j)
    END DO
    CLOSE(1)
    theta(1)=1.0d0
    theta(2)=1.0d0
    theta(3)=-1.0d0
    !WRITE(*,*) y(5),x(5,1),x(5,2),x(5,3),D(5)
    DD = REAL(D)
    B_hat=OLS(x,DD)
    WRITE(*,*) B_hat

    hah=0
    !a)
    mode=0
    call set_seed(1,2,3,4)
    DO i=1,N_individuals
        epis(i)=Sample_Normal(0.0d0,1.0d0)
    END DO
    I_hat=0.0d0
    ! WRITE(*,*) epis(1:10)
    !test= GMM(theta)
    theta_a(1)=theta(1)
    theta_a(2)=theta(2)
    theta_a(3)=theta(3)
    WRITE(*,*) "B_hat:", B_hat
    call Nelder_Meade(theta_a,tol,GMM,1,max_NM)
    !WRITE(*,*) I_hat(49990:50000)


    !b)
    mode=1
    I_hat=0.0d0
    WRITE(*,*) "B_hat:", B_hat
    theta_b(1)=theta(1)
    theta_b(2)=theta(2)
    theta_b(3)=theta(3)
    CALL BFGS(theta_b,1.0d-6,1.0d-6,GMM,0,.TRUE.)
    WRITE(*,*) "        BFGS Theta: ", theta_b
    WRITE(*,*) "Nelder Meade Theta: ", theta_a


    !c)
    temp_x=x
    temp_y=y
    temp_D=D
    N_individuals_R = REAL(N_individuals)
    BS_data = 0.0d0
    call set_seed(1,2,3,4)
    do i=1,BS

        DO j=1,N_individuals
            inde = Sample_Uniform(1.0d0, N_individuals_R)
            !WRITE(*,*) inde
            x(j,:)=temp_x(int(inde),:)
            D(j)=temp_D(int(inde))

        END DO
        DD = REAL(D)
        B_hat_BS(i,:) = OLS(x,DD)

    END DO
    theta_BS_avg=0.0d0
    !WRITE(*,*) "BS Beta ", B_hat_BS(99,:)
    !WRITE(*,*) "BS Beta ", B_hat_BS(5,:)
    theta_BS_avg = SUM (B_hat_BS, DIM=1)/BS
    B_hat_BS(:,1)=B_hat_BS(:,1)-theta_BS_avg(1)
    B_hat_BS(:,2)=B_hat_BS(:,2)-theta_BS_avg(2)
    B_hat_BS(:,3)=B_hat_BS(:,3)-theta_BS_avg(3)

    BS_data_se = matmul(TRANSPOSE(B_hat_BS),B_hat_BS)
    BS_data_se = BS_data_se/(BS-1.0d0)
    !WRITE(*,*) "BS_data_se 1", BS_data_se(1,:)
    !WRITE(*,*) "BS_data_se 2", BS_data_se(2,:)
    !WRITE(*,*) "BS_data_se 3", BS_data_se(3,:)
    hah=1
    mode=0
    test= GMM(theta)

    theta_c(1)=theta(1)
    theta_c(2)=theta(2)
    theta_c(3)=theta(3)
    call Nelder_Meade(theta_c,tol,GMM,1,max_NM)

    mode=1
    theta_d(1)=theta(1)
    theta_d(2)=theta(2)
    theta_d(3)=theta(3)
    CALL BFGS(theta_d,1.0d-6,1.0d-6,GMM,0,.TRUE.)

    WRITE(*,*) "                 NM Theta: ", theta_a
    WRITE(*,*) "               BFGS Theta: ", theta_b
    WRITE(*,*) "  NM Theta Optimal Weight: ", theta_c
    WRITE(*,*) "BFGS Theta Optimal Weight: ", theta_d
    WRITE(*,*) "                    SE_BS: ",sqrt(BS_data_se(1,1)),sqrt(BS_data_se(2,2)),sqrt(BS_data_se(3,3))











END PROGRAM MAIN
