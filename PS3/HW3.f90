MODULE HW3
    USE GLOBVAR
    USE INIT
    USE MINIMIZATION
    USE LIKELIHOOD
    USE MATRIX
    USE RANDOM
    USE PROBABILITY, ONLY : CDF_NORMAL,PDF_NORMAL
    USE SIMPLEX
    implicit none

    CONTAINS

  REAL(8) FUNCTION GMM(theta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: aux, tempp(N_covariates,N_covariates)
    INTEGER :: i
    REAL(8) :: B_star(N_covariates),D_star(N_individuals), diff(N_covariates)

    ! Let's recover the parameters from the vector
    theta_hat=theta(1:N_covariates)

    GMM =0.0d0
    D_star=0.0d0
    I_hat = theta_hat(1)+theta_hat(2)*x(:,2)+theta_hat(3)*x(:,3)+epis
    IF(mode==0) THEN
        DO i = 1, N_individuals
            IF(I_hat(i)>=0) THEN
                D_star(i)=1.0d0
            ELSE
                D_star(i)=0.0d0
            END IF
        END DO
    ELSE IF (mode==1) THEN
        !WRITE(*,*) "hahahahahahahahahhaahhahahahah"
        DO i = 1, N_individuals
            D_star(i)=exp(I_hat(i)/0.1)/(1+exp(I_hat(i)/0.1))
        END DO
    END IF

    B_star=OLS(x,D_star)
    diff = B_star-B_hat
    IF(hah==0) THEN
        GMM = diff(1)*diff(1)+diff(2)*diff(2)+diff(3)*diff(3)
        !WRITE(*,*) "GMM:", GMM
        !WRITE(*,*) "B_star:", B_star
    ELSE IF (hah==1) THEN
        tempp = 0
        tempp(:,1) = diff
        tempp = matmul(matmul(transpose(tempp),Matrix_Inverse(BS_data_se)),tempp)
        WRITE(*,*) "tempp 1:", tempp(1,:)
        WRITE(*,*) "tempp 2:", tempp(2,:)
        WRITE(*,*) "tempp 3:", tempp(3,:)
        GMM = tempp(1,1)

    END IF



  END FUNCTION GMM


END MODULE HW3
