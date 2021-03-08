! include 'likelihood_evaluation.f90'

MODULE LIKELIHOOD
  USE GLOBVAR
  USE INIT
  USE LIKELIHOOD_EVALUATION, ONLY : LIKE_EVAL
  USE PROBABILITY
  USE MATRIX
  USE LAPACK95
  USE MINIMIZATION
  IMPLICIT NONE

  PUBLIC :: LOGLIKELIHOOD, var_probit,se!, probit

  PRIVATE

CONTAINS

  REAL(8) FUNCTION LOGLIKELIHOOD(theta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: aux
    INTEGER :: i

    ! Let's recover the parameters from the vector
    beta=theta!(1:N_covariates)
    ! Since variance needs to be positive, instead of maximizing with respect to variance
    ! lets maximize with respect to a transformation. Since variance=EXP(LOG(variance)) if we call
    ! alpha=LOG(variance), and maximize with respect to alpha we are ok.
    ! variance=EXP(theta(N_covariates))

    LOGLIKELIHOOD =0.0d0
    DO i = 1, N_individuals
      aux = LIKE_EVAL(i)
      if ( D(i)==1 ) then

        IF ((1.0d0-aux)<1.0d-300) THEN
           LOGLIKELIHOOD = LOGLIKELIHOOD - log(1.0d0-1.0d-300)! 690.77552789821368151024216786026954650879d0 ! -Log of 1.0d-300
        ELSE
           LOGLIKELIHOOD = LOGLIKELIHOOD - LOG(1.0d0 - aux) ! - because we are minimizing
        END IF
      else
         ! aux = LIKE_EVAL(i)
         IF (aux<1.0d-300) THEN
            LOGLIKELIHOOD = LOGLIKELIHOOD + 690.77552789821368151024216786026954650879d0 ! -Log of 1.0d-300
         ELSE
            LOGLIKELIHOOD = LOGLIKELIHOOD - LOG(aux) ! - because we are minimizing
         END IF
       end if
    END DO
  END FUNCTION LOGLIKELIHOOD

  function var_probit(x,beta)
    ! Function that obtains the probit variance
    ! B=E[xx'φ(x'β)^2/Φ(x'β)Φ(-x'β)]
    ! takes
    real(8), intent(in) :: x(:,:),beta(:)
    ! integer :: n,m
    INTEGER :: i,j,k
    real(8) :: var_probit(size(x(1,:)),size(x(1,:))),xv(size(x(1,:)))
    real(8) :: xb,xx(size(x(1,:)),size(x(1,:))),phi,phi_s,den,xb_m
    real(8), PARAMETER :: mu=0.0d0, var=1.0d0
    integer :: n, IPIV(size(x(1,:)))
    n=size(x(:,1))
    ! m=size(x(1,:))
    print*,2
    var_probit=0.0d0
    do i = 1,size(x(:,1))
      xb = -dot_product(x(i,:),beta)
      ! xb_m=-xb
      ! x_p=RESHAPE(x(i,:),[3,1],order=[2,1])
      ! ! xx=MATMUL(x_p,x(i,:))
      ! xv=[((x(i,j)*x(i,k),j=1,3),k=1,3)]
      ! xx=reshape(xv,[3,3],order=[2,1])
      FORALL(j=1:3,k=1:3) xx(j,k)=x(i,j)*x(i,k)

      phi=PDF_NORMAL(xb,mu,var)
      phi_s=phi*phi
      den=CDF_NORMAL(xb,mu,var)*(1.0d0-CDF_NORMAL(xb,mu,var))
      var_probit=var_probit+xx*(phi_s/den)
    end do

    ! print*,"before dividing"
    ! print*,var_probit
    var_probit=var_probit/DBLE(n)
    ! print*,"before inverting"
    ! print*, var_probit
    ! var_probit=Matrix_Inverse_symmetric(var_probit)
    ! call SYTRF(var_probit,"U",IPIV)
    ! call SYTRI(var_probit,IPIV,"U")
    !
    ! call GETRF(var_probit,IPIV)
    ! call GETRI(var_probit,IPIV)

    call POTRF(var_probit)
    call POTRI(var_probit)

    ! var_probit=Matrix_Inverse(var_probit)

  end function var_probit

  function se(x,beta)
    real(8), INTENT(IN) :: x(:,:),beta(:)
    real(8) :: vcov(size(x(1,:)),size(x(1,:))),se(size(x(1,:)))
    integer :: i,n,l
    vcov=var_probit(x,beta)
    ! print*,"varcovar= "
    ! print*,vcov
    n=size(x(:,1))
    l=size(vcov(1,:))
    se=[ (SQRT(abs(vcov(i,i))/DBLE(n)),i=1,l) ]

  end function se

  ! function probit(x,theta,fun)
  !   interface
  !     function fun(b) result(f_res)
  !       real(8), INTENT(IN) :: b(:)
  !       real(8) :: f_res
  !     end function fun
  !   end interface
  !   real(8), INTENT(IN) :: x(:,:),theta(N_covariates)
  !   REAL, ALLOCATABLE, DIMENSION(:) :: probit
  !   real :: b(N_covariates)
  !   b=theta
  !   CALL BFGS(b,1.0d-6,1.0d-6,fun,0,.TRUE.)
  !
  !   probit=theta
  ! end function probit

END MODULE LIKELIHOOD
