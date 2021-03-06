MODULE LIKELIHOOD
  USE GLOBVAR
  USE LIKELIHOOD_EVALUATION, ONLY : LIKE_EVAL
  IMPLICIT NONE
  
CONTAINS 
  
  REAL(8) FUNCTION LOGLIKELIHOOD(theta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: aux
    INTEGER :: i
    
    ! Let's recover the parameters from the vector
    beta=theta(1:n_covariates)
    ! Since variance needs to be positive, instead of maximizing with respect to variance
    ! lets maximize with respect to a transformation. Since variance=EXP(LOG(variance)) if we call 
    ! alpha=LOG(variance), and maximize with respect to alpha we are ok. 
    variance=EXP(theta(n_covariates+1))
    
    LOGLIKELIHOOD =0.0d0
    DO i = 1, N_individuals
       aux = LIKE_EVAL(i)
       IF (aux<1.0d-300) THEN
          LOGLIKELIHOOD = LOGLIKELIHOOD + 690.77552789821368151024216786026954650879d0 ! -Log of 1.0d-300
       ELSE
          LOGLIKELIHOOD = LOGLIKELIHOOD - LOG(aux) ! - because we are minimizing
       END IF
    END DO    
  END FUNCTION LOGLIKELIHOOD
  
END MODULE LIKELIHOOD

