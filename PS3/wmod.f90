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

  public :: b_star,varcovar,objf,objfvcov,bs_m,varbeta!,rec_msle
  private

contains

  function b_starm(theta,smooth) result(theta_hat)
    !This function samples from U, generates I, and D
    !estimates β⋆ from sampled data and updates theta.
    implicit none
    real(8), INTENT(IN) :: theta(N_covariates)
    real(8) :: xb(N_individuals), big_i(N_individuals,s), &
               d_s(N_individuals,s),&
               b_s(s,N_covariates), theta_hat(N_covariates),t1(N_individuals,N_covariates),&
               t2(N_individuals),t3(N_covariates)

    integer :: i,j,smooth

    ! generate fake data set
    t1=x
    t2=d
    t3=theta
    xb=MATMUL(t1,t3)

    do j=1,s
      big_i(:,j) = xb + nor(:,j)
    enddo

    if (smooth==1) then
      ! print*,"Smoothing D"
      do i=1,N_individuals
        do j=1,s
          d_s(i,j) = exp(big_i(i,j)/k)/(1.0d0+exp(big_i(i,j)/k))
        enddo
      enddo
    else if (smooth==0) then
      ! print*,"Indicator D"
      do i=1,N_individuals
        do j=1,s
          if (big_i(i,j)>=0.0d0) then
            d_s(i,j)=1.0d0
          else
            d_s(i,j)=0.0d0
          endif
        enddo
      enddo
    end if

    ! Run OLS for each simiulation

    do j=1,s
      b_s(j,:)=OLS(x,d_s(:,j))
    end do

    theta_hat=sum(b_s, dim=1)/DBLE(s)

  end function b_starm

  function b_star(theta,smooth) result(theta_hat)
    !This function samples from U, generates I, and D
    !estimates β⋆ from sampled data and updates theta.
    implicit none

    real(8), INTENT(IN) :: theta(N_covariates)
    real(8) :: xb(N_individuals), big_i(N_individuals), &
               d_s(N_individuals),&
               theta_hat(N_covariates), t1(N_individuals,N_covariates),&
               t2(N_covariates)
    integer :: i,j,smooth
    ! logical, INTENT(IN) :: smooth

    ! generate fake data set
    t1=x
    t2=theta
    xb=MATMUL(t1,t2)
    ! print*,"Sum xb = ",sum(xb)/DBLE(N_individuals)
    big_i = xb + nor(:,1)

    ! beta=0.0d0
    if (smooth==1) then
      d_s = exp(big_i/k)/(1.0d0+exp(big_i/k))

      ! do i=1,N_individuals
      !   d_s(i) = exp(big_i(i)/k)/(1.0d0+exp(big_i(i)/k))
      ! enddo

    else if (smooth==0) then
      ! print*,"Indicator D"
      do i=1,N_individuals
        if (big_i(i)>=0.0d0) then
          d_s(i)=1.0d0
        else
          d_s(i)=0.0d0
        endif
      enddo
    end if

    ! Run OLS
    theta_hat=OLS(x,d_s)

  end function b_star
  !
  real(8) function objf(theta)
    !This function takes in current theta values and
    !returns the value of the objective function
    implicit none
    real(8), INTENT(IN) :: theta(:)
    real(8) :: diff(N_covariates),&! idmat(N_covariates,N_covariates),&
               theta_hat(N_covariates),b_s(N_covariates)

    theta_hat=theta
    ! print*,"beta",beta
    ! Sample from U, get I, get D, and estimate β ^⋆
    if (ind==1) then
      ! print*,"sampling u=1"
      b_s=b_star(theta_hat,smooth)
    else if (ind==0) then
      ! print*,"sampling u=100"
      b_s=b_starm(theta_hat,smooth)
    end if

    !
    ! Get vector of differences
    ! print*,b_s
    ! print*,b_h

    diff=b_s-b_h

    objf=dot_product(diff,diff)


  end function objf

  real(8) function objfvcov(theta)
    !Same as objf but using optimal weight matrix vcov (global)
    implicit none
    real(8), INTENT(IN) :: theta(:)
    real(8) :: diff(N_covariates),&! idmat(N_covariates,N_covariates),&
               b_s(N_covariates),theta_hat(N_covariates)
    integer :: i,j

    theta_hat=theta

    if (ind==1) then
      b_s=b_star(theta_hat,smooth)
    else if (ind==0) then
      b_s=b_starm(theta_hat,smooth)
    end if
      ! b_h=OLS(x,D)
    ! print*,"Target beta = ",b_h
    diff=b_s-b_h
    ! print*,"beta - beta hat = ",diff

    objfvcov= DOT_PRODUCT(matmul(diff,vcov),diff)

  end function objfvcov

  ! recursive function rec_msle(theta,tol,maxiter,iter,smooth,model) result(reslt)
  !   implicit none
  !   real(8), INTENT(IN) :: tol
  !   real(8), INTENT(IN) :: theta(:)
  !   INTEGER, INTENT(IN) :: maxiter, iter
  !   logical, intent(IN) :: smooth!,sigma
  !   real(8) :: reslt(N_covariates)
  !   INTEGER :: i
  !   real(8) :: t(N_covariates), dist,beta(N_covariates),beta_s(N_covariates)
  !   integer, intent(in) :: model
  !   print*, "theta in = ",theta
  !   ! t=theta
  !   beta=theta
  !   i=iter + 1
  !
  !   ! call b_star(beta,smooth)
  !
  !   ! if (sigma) then
  !   !   call Nelder_Meade(beta,tol,objf,1,maxiter)
  !   ! else
  !   !   call Nelder_Meade(beta,tol,objfvcov,1,maxiter)
  !   ! endif
  !
  !   select case (model)
  !   case (1)
  !       ! print*,"NMI"
  !       call Nelder_Meade(beta,1.0d-8,objf,1,maxiter)
  !     case (2)
  !       ! print*,"BFGSI"
  !       call BFGS(beta,1.0d-8,1.0d-8,objf,1,.TRUE.)
  !     case (3)
  !       ! print*,"NMW"
  !       call Nelder_Meade(beta,1.0d-8,objfvcov,1,maxiter)
  !     case (4)
  !       ! print*,"BFGSW"
  !       call BFGS(beta,1.0d-8,1.0d-8,objfvcov,1,.TRUE.)
  !   end select
  !
  !   print*,"updated theta = ",beta
  !   print*,"OLS = ",b_h
  !
  !   beta_s=b_star(beta,smooth)
  !   print*,"Beta star = ",beta_s
  !   dist=MAXVAL(ABS(beta_s-b_h))
  !   print*,"Distance = ",dist
  !
  !   if (dist<tol) then
  !     print*, "Tolerance reached in ",iter,&
  !             " iteraions. Distance = ",dist
  !     reslt=beta
  !   elseif (i>maxiter) then
  !     print*, "Max iter reached. Distance = ",dist
  !     reslt=beta
  !   else
  !     ! reslt=theta
  !     if (modulo(i,10)==0) then
  !       select case (model)
  !       case (1)
  !           print*,"NMI"
  !           ! call Nelder_Meade(beta,1.0d-8,objf,1,maxiter)
  !         case (2)
  !           print*,"BFGSI"
  !           ! call BFGS(beta,1.0d-8,1.0d-8,objf,1,.TRUE.)
  !         case (3)
  !           print*,"NMW"
  !           ! call Nelder_Meade(beta,1.0d-8,objfvcov,1,maxiter)
  !         case (4)
  !           print*,"BFGSW"
  !           ! call BFGS(beta,1.0d-8,1.0d-8,objfvcov,1,.TRUE.)
  !       end select
  !       print*,"Iteration ",i
  !       print*,"Distance ",dist
  !     endif
  !
  !     ! print*,"Theta out = ",beta
  !
  !     reslt=rec_msle(beta,tol,maxiter,i,smooth,model)
  !   endif
  ! end function rec_msle


  subroutine bs_m(x,d,f,varm)
    !This subroutine estimates the bootstrap variacnce covariance matrix
    !takes in matrix of covariates x and vector of outcome variable d
    !and returns the bootstrapped varcovar matrix
    real(8) :: x_bs(N_individuals,N_covariates), d_bs(N_individuals),&
               t1(N_individuals,N_covariates), t2(N_individuals),&
               varcovar_bs(N_covariates,N_covariates),&
               olsv(N_covariates,N_covariates)!,&
               ! varm(N_covariates,N_covariates)!,AB(N_covariates,N_covariates)
   real(8), intent(in) :: x(:,:),d(:)
   real(8), intent(out) :: varm(N_covariates,N_covariates)!,vcov(:,:)
   integer :: i!,j,k

   interface
     function f(g,h) result(f_result)
       real(8), intent(in) :: g(:,:), h(:)
       real(8) :: f_result(size(g,2),size(g,2))
     end function
   end interface

   t1=x
   t2=d

   olsv=f(t1,t2) !getting OLS varcovar matrix

   varm=0.0d0
   do i=1,n_boot
     !Bootstrapping data
     x_bs=t1(int(uni(:,i)),:)
     d_bs=t2(int(uni(:,i)))

     ! print*,"x=",sum(x,dim=1)
     ! print*,"x_bs=", sum(x_bs,dim=1)

     !Get varcovar for current draw
     varcovar_bs=f(x_bs,d_bs)

     ! print*,"covariance matrix done for s=",i
     varm=varm+varcovar_bs
    enddo

    varm=varm/DBLE(n_boot) !Getting the bootstrap estimator
    varm=2.0d0*olsv-varm !Getting the unbiased bootstrap estimator

    print*,"Final VCOV = "
    print*,varm
    call write_file("sigma.txt",varm)

    call POTRF(varm)
    call POTRI(varm)

    print*,"Final VCOV inverted = "
    print*,varm
    call write_file("sigma_inv.txt",varm)

  end subroutine bs_m

  function varcovar(xin,din) result(varcovar_bs)
    !This function estimates the variance covariance matrix of OLS assuming
    !Homosckedasticity. Takes in a matrixx of covariates x
    !and outcomes variable vector d and returns a symmetric matrix

    real(8) :: x_bs(N_individuals,N_covariates), d_bs(N_individuals),&
               xx_bs(N_covariates,N_covariates),b_bs(N_covariates),&
               xb_bs(N_individuals),residuals(N_individuals),&
               x_bs_t(N_covariates,N_individuals),u_bar
   real(8), intent(in) :: xin(:,:),din(:)
   real(8) :: varcovar_bs(size(xin,2),size(xin,2))
   integer :: i,ipiv(size(xx_bs))
   x_bs=xin
   d_bs=din

    ! Getting residuals
    ! xux=0.0d0
    ! xx_bs=0.0d0
    ! print*,"OLS from global"
    ! print*,b_h
    b_bs=OLS(x,d)
    do i=1,N_individuals
      xb_bs(i)=DOT_PRODUCT(x_bs(i,:),b_bs)
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

  function varbeta(xin,din) result(var_b)

    real(8),intent(in) :: xin(:,:), din(:)
    real(8) :: var_b(size(xin,2),size(xin,2))
    real(8) :: t1(size(xin,1),size(xin,2)),t2(size(din)),x_i(N_individuals,N_covariates)
    real(8) :: d_i(N_individuals),b_i(n_boot,N_covariates),mbi(N_covariates)
    integer :: i

    t1=xin
    t2=din

    do i=1,n_boot
      x_i=t1(uni(:,i),:)
      d_i=t2(uni(:,i))
      b_i(i,:)=OLS(x_i,d_i)
    enddo

    mbi=sum(b_i,1)/dble(n_boot)
    ! print*,"Average beta hat=",mbi

    do i=1,N_covariates
      b_i(:,i)=b_i(:,i)-mbi(i)
    enddo
    ! print*,sum(b_i,1)/dble(n_boot)

    var_b=matmul(transpose(b_i),b_i)
    var_b=var_b/dble(n_boot-N_covariates)

    print*,"Final VCOV = "
    print*,var_b
    call write_file("sigma.txt",var_b)

    call POTRF(var_b)
    call POTRI(var_b)

    print*,"Final VCOV inverted = "
    print*,var_b
    call write_file("sigma_inv.txt",var_b)

  end function varbeta

end module wmod
