! include 'lapack.f90'
module wmod
  use globvar
  ! use init
  use matrix
  use simplex
  use LAPACK95 ,only : POTRF,POTRI,sytrf,sytri,GETRF,GETRI
  use io
  use MINIMIZATION
  ! use lapack

  implicit none

  public :: vcov_bs,rec_msle, varcovar
  private

contains

  subroutine b_star(theta,smooth)
    !This function samples from U, generates I, and D
    !estimates β⋆ from sampled data and updates theta.
      implicit none
    ! use GLOBVAR
    ! use init
    ! real(8), INTENT(IN) :: x(N_individuals,N_covariates)
    real(8), INTENT(INOUT) :: theta(N_covariates)
    real(8) :: xb(N_individuals), big_i(N_individuals,s), &
               d_s(N_individuals,s),&
               b_s(s,N_covariates)

    integer :: i,j
    logical, INTENT(IN) :: smooth
    d_s=0.0d0
    ! generate fake data set
    xb=MATMUL(x,theta)
    ! print*,"Sum xb = ",sum(xb)/DBLE(N_individuals)
    forall (i=1:N_individuals, j=1:s) big_i(i,j) = xb(i) + nor(i,j)

    ! If smooth (global) is TRUE then it uses D=exp(xβ)/1+exp(xβ)
    !instead of the indicator function.
    if (smooth) then
      forall (i=1:N_individuals,j=1:s) d_s(i,j) = exp(big_i(i,j)/k)/(1.0d0+exp(big_i(i,j)/k))
    else
      forall (i=1:N_individuals,j=1:s,big_i(i,j).GT.0.0d0) d_s(i,j) = 1.0d0
      ! do i=1,N_individuals
      !   do j=1,s
      !     if (big_i(i,j).GT.0.0d0) then
      !       d_s(i,j)=1.0d0
      !     else
      !       d_s(i,j)=0.0d0
      !     endif
      !   enddo
      ! enddo
    end if
    ! end forall

    ! Run OLS for each simiulation
    ! forall(k=1:s)
    do j=1,s
      b_s(j,:)=OLS(x,d_s(:,j))
    end do
    ! get b_star
    theta=sum(b_s, dim=1)/DBLE(s)

  end subroutine b_star
  !
  real(8) function objf(theta)
    !This function takes in current theta values and
    !returns the value of the objective function
    implicit none
    real(8), INTENT(IN) :: theta(:)
    ! real(8) :: b_h(N_covariates), D(N_)
    real(8) :: diff(N_covariates), idmat(N_covariates,N_covariates),&
               p(N_covariates), beta(N_covariates)
    integer :: i,j
    ! Sample from U, get I, get D, and estimate β ^⋆
    beta=theta
    ! print*,"beta",beta
    ! call b_star(beta)
    ! print*,"updated beta ",beta
    !
    ! Get vector of differences
    diff=beta-b_h
    ! Generate identity matrix
    idmat=0.0d0
    forall(i=1:size(idmat(1,:)), j=1:size(idmat(1,:)), i==j) idmat(i,j)=1.0d0
    ! Estimate objective function
    p=matmul(idmat,diff)
    objf= DOT_PRODUCT(p,diff)

  end function objf

  real(8) function objfvcov(theta)
    !Same as objf but using optimal weight matrix vcov (global)
    implicit none
    real(8), INTENT(IN) :: theta(:)
    ! real(8) :: b_h(N_covariates), D(N_)
    real(8) :: diff(N_covariates), idmat(N_covariates,N_covariates),&
               p(N_covariates), beta(N_covariates)
    integer :: i,j

    beta=theta
    ! call b_star(beta)
    ! b_h=OLS(x,D)
    ! print*,"Target beta = ",b_h
    diff=beta-b_h
    ! print*,"beta - beta hat = ",diff
    ! dt=transpose(d)
    ! idmat=0.0d0
    ! forall(i=1:size(idmat(1,:)), j=1:size(idmat(1,:)), i==j) idmat(i,j)=1.0d0
    ! print*,idmat
    p=matmul(vcov,diff)

    objfvcov= DOT_PRODUCT(p,diff)

  end function objfvcov

  recursive function rec_msle(theta,tol,maxiter,iter,model) result(reslt)
    implicit none
    real(8), INTENT(IN) :: tol
    real(8), INTENT(IN) :: theta(:)
    INTEGER, INTENT(IN) :: maxiter, iter
    ! logical, intent(IN) :: smooth,sigma
    real(8) :: reslt(N_covariates)
    INTEGER :: i
    real(8) :: t(N_covariates), dist,beta(N_covariates)
    integer, intent(in) :: model
    ! print*, "theta in = ",theta
    ! t=theta
    beta=theta
    i=iter + 1

    call b_star(beta,smooth)

    ! if (sigma) then
    !   call Nelder_Meade(beta,tol,objf,1,maxiter)
    ! else
    !   call Nelder_Meade(beta,tol,objfvcov,1,maxiter)
    ! endif

    select case (model)
    case (1)
        print*,"NMI"
        call Nelder_Meade(beta,tol,objf,1,maxiter)
      case (2)
        print*,"BFGSI"
        call BFGS(beta,tol,tol,objf,1,.TRUE.)
      case (3)
        print*,"NMW"
        call Nelder_Meade(beta,tol,objfvcov,1,maxiter)
      case (4)
        print*,"BFGSW"
        call BFGS(beta,tol,tol,objfvcov,1,.TRUE.)
    end select

    ! print*,"updated theta = ",beta
    dist=ABS(MAXVAL(beta/b_h- 1.0d0))

    if (dist<tol) then
      print*, "Tolerance reached in ",iter,&
              " iteraions. Distance = ",dist
      reslt=beta
    elseif (i>maxiter) then
      print*, "Max iter reached. Distance = ",dist
      reslt=beta
    else
      ! reslt=theta
      if (modulo(i,10)==0) then
        print*,"Iteration ",i
        print*,"Distance ",dist
      endif
      ! print*,"Theta out = ",beta

      reslt=rec_msle(beta,tol,maxiter,i,model)
    endif
  end function rec_msle


  subroutine vcov_bs(x,d,vcov)
    !This subroutine estimates the bootstrap variacnce covariance matrix
    !takes in matrix of covariates x and vector of outcome variable d
    !and returns the bootstrapped varcovar matrix
    real(8) :: x_bs(N_individuals,N_covariates), d_bs(N_individuals),&
               t1(N_individuals,N_covariates), t2(N_individuals),&
               varcovar_bs(N_covariates,N_covariates),&
               olsv(N_covariates,N_covariates),&
               xxpu(N_covariates,N_individuals)!,AB(N_covariates,N_covariates)
   real(8), intent(in) :: x(:,:),d(:)
   real(8), intent(out) :: vcov(:,:)
   integer :: i!,j,k
   t1=x
   t2=d

   olsv=varcovar(t1,t2) !getting OLS varcovar matrix

   vcov=0.0d0
   do i=1,s
     !Bootstrapping data
     x_bs=t1(int(uni(:,i)),:)
     d_bs=t2(int(uni(:,i)))

     ! print*,"x=",sum(x,dim=1)
     ! print*,"x_bs=", sum(x_bs,dim=1)

     !Get varcovar for current draw
     varcovar_bs=varcovar(x_bs,d_bs)

     ! print*,"covariance matrix done for s=",i
     vcov=vcov+varcovar_bs
    enddo

    vcov=vcov/DBLE(s) !Getting the bootstrap estimator
    vcov=2.0d0*olsv-vcov !Getting the unbiased bootstrap estimator

    print*,"Final VCOV = ",vcov
    call write_file("sigma.txt",vcov)

    call POTRF(vcov)
    call POTRI(vcov)

    print*,"Final VCOV inverted = ",vcov
    call write_file("sigma_inv.txt",vcov)

  end subroutine vcov_bs

  function varcovar(x,d) result(varcovar_bs)
    !This function estimates the variance covariance matrix of OLS assuming
    !Homosckedasticity. Takes in a matrixx of covariates x
    !and outcomes variable vector d and returns a symmetric matrix

    real(8) :: x_bs(N_individuals,N_covariates), d_bs(N_individuals),&
               xx_bs(N_covariates,N_covariates),&!b_bs(N_covariates),&
               xb_bs(N_individuals),residuals(N_individuals),&
               x_bs_t(N_covariates,N_individuals),u_bar
   real(8), intent(in) :: x(:,:),d(:)
   real(8) :: varcovar_bs(N_covariates,N_covariates)
   integer :: i,ipiv(size(xx_bs))
   x_bs=x
   d_bs=d

    ! Getting residuals
    ! xux=0.0d0
    xx_bs=0.0d0
    ! print*,"OLS from global"
    ! print*,b_h
    do i=1,N_individuals
      xb_bs(i)=DOT_PRODUCT(x_bs(i,:),b_h)
      residuals(i)=d_bs(i)-xb_bs(i)
      residuals(i)=residuals(i)*residuals(i)
    enddo

    u_bar=sum(residuals)/dble(N_individuals)
    ! print*, "residuals = ",u_bar


    ! !Get Variance Covariance Matrix
    ! Getting x'x
    x_bs_t=transpose(x_bs)
    xx_bs=matmul(x_bs_t,x_bs)

    ! Inverting x'x. Getting (x'x)^(-1)

    call POTRF(xx_bs)
    call POTRI(xx_bs)
    !
    ! call SYTRF(xx_bs,'U',ipiv)
    ! call SYTRI(xx_bs,ipiv)
    !
    ! call GETRF(xx_bs,ipiv)
    ! call GETRI(xx_bs,ipiv)

    ! print*,"A inverted = "
    ! print*,xx_bs

    !Getting variance matrix
    varcovar_bs=xx_bs*u_bar


    ! print*,"Final OLS VCOV = "
    ! print*,varcovar_bs
    ! call write_file("OLS_vcov.txt",varcovar_bs)

  end function varcovar

end module wmod
