
include "lapack.f90"

MODULE GLOBVAR
  USE INTERPOL
  ! USE omp_lib
  IMPLICIT NONE
  ! This modlue declares the variables that are going to be used throughout the program
  ! higher level languages do not require you to declare your variables ahead of time

  INTEGER, PARAMETER :: N_grid=20, N_sim=1008
  INTEGER, PARAMETER :: N_quad=20
  REAL(8), PARAMETER :: beta=0.97d0, &
       r=0.1d0, &
       imin=0.1d0, &
       sigma2_l=0.5d0

  REAL(8) :: y(1:2,1:N_sim),constemp,valuetemp!,test_data(1:5000,1:7)
  REAL(8) :: x_quad(1:N_quad),w_quad(1:N_quad),v_new(1:2,1:N_sim),c_new(1:2,1:N_sim)
  real(8) :: assets(1:3,1:N_sim), v2_1temp, c2_1temp, a3_1temp, res_1, res_0, &
              v2_0temp, c2_0temp, a3_0temp, v1_1temp, a2_1temp, c1_1temp, &
              v1_0temp, a2_0temp, c1_0temp, res1_1, res1_0
  integer :: lab(1:2, 1:N_sim), t_lab, lt, verb=1
  integer :: num_threads=20,num_procs, max_thr, act_lev
  real(8) :: rho_0(1:15), rho_h(1:15), l_v(1:5000), f_data(1:N_sim,1:7)!, eps1, eps2

  ! real(8) :: alpha=1.1d0, &
  !            mu=0.0d0, &
  !            gamma=2.0d0, &
  !            delta=0.2d0, &
  !            sigma2_e=1.0d0, &
  !            lambda=5.0d0
  real(8) :: alpha,&!=1.1d0, &
            mu,&!=0.0d0, &
            gamma,&!=2.0d0, &
            delta,&!=0.2d0, &
            sigma2_e,&!=1.0d0, &
            lambda!=5.0d0

  integer :: ierr, myid, numprocs
  TYPE(PLINEARSTRUC) :: inter_value

CONTAINS

  real(8) function u(c, l)
    implicit none
    real(8), INTENT(IN) :: c
    real(8) :: arg1
    INTEGER, INTENT(IN) :: l

    if (alpha==1.0d0) then
      arg1=-1.0d300
    elseif (c<1.0d-16) then
      arg1=-1.0d300
    else
      arg1=(c**(1.0d0-alpha)-1.0d0)/(1.0d0-alpha)
    endif
    u=arg1-lambda*DBLE(l)
    ! print*, "lambda",lambda
  end function u

  real(8) function w(a)
    implicit none
    real(8), INTENT(IN) :: a
    real(8) :: arg1
    !INTEGER, INTENT(IN) :: lab
    if (gamma==1.0d0) then
      arg1=-1.0d300
    elseif (a<1.0d-16) then
      arg1=-1.0d300
    else
      arg1=((a*(1.0d0+r))**(1.0d0-gamma)-1.0d0)/(1.0d0-gamma)
    endif
    w=arg1
  end function w

  real(8) elemental function y_fun(l,lp,eps)! result(yt)
    real(8), intent(in) :: eps
    integer, intent(in) :: l, lp
    real(8) :: yt

    if (lp==0) then
      y_fun=0.0d0
    elseif (lp==1) then
      y_fun=exp(mu+delta*dble(l)+eps)
    endif

  end function y_fun

  real(8) function rcs(a, l, lprime, eps)
    implicit none
    real(8), INTENT(in) :: a, eps
    integer, INTENT(IN) :: l, lprime
    REAL(8) :: yt, inc

    ! lny=mu+delta*DBLE(l)+eps
    ! yt=exp(lny)
    ! if (lprime==0) then
    !   inc=imin
    ! else
    !   inc=yt+MAX(imin-yt,0.0d0)
    ! endif
    yt=y_fun(l,lprime,eps)
    inc=yt+MAX(imin-yt,0.0d0)
    rcs=inc+(1.0d0+r)*a
  end function rcs

  ! We need to know what y is inside this function
  ! I have a guess interpolated of the value function somewhere
  ! and this function needs access to it
  ! Somehow this routine needs to know the quadrature points
  REAL(8) FUNCTION OBJ_V2_1(theta)
    ! USE PROBABILITY
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: cons,arg1,kprime,arg2
    INTEGER :: j
    cons=EXP(theta(1))
    cons=res_1/(1.0d0+cons)
    ! Consumption lives between 0 and y

    arg1=u(cons,1)
    ! print*, arg1
    kprime=res_1-cons
    ! integ=0.0d0
    arg2=beta*w(kprime)
    ! print*, arg2, cons, kprime
    ! lt=t_lab
    c2_1temp=cons
    a3_1temp=kprime
    OBJ_V2_1=-arg1-arg2
    v2_1temp=-OBJ_V2_1
  END FUNCTION OBJ_V2_1

  REAL(8) FUNCTION OBJ_V2_0(theta)
    ! USE PROBABILITY
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: cons,arg1,kprime,arg2
    INTEGER :: j
    cons=EXP(theta(1))
    cons=res_0/(1.0d0+cons)
    ! Consumption lives between 0 and y

    arg1=u(cons,0)
    ! print*, arg1
    kprime=res_0-cons
    ! integ=0.0d0
    arg2=beta*w(kprime)
    ! print*, arg2, cons, kprime
    ! lt=t_lab
    c2_0temp=cons
    a3_0temp=kprime
    OBJ_V2_0=-arg1-arg2
    v2_0temp=-OBJ_V2_0
  END FUNCTION OBJ_V2_0

  real(8) function V2(a2,l1,eps,i)
    use SIMPLEX1
    ! use omp_lib
    implicit none
    real(8), INTENT(IN) :: a2, eps
    integer, INTENT(IN) :: l1,i
    ! integer :: t_lab
    real(8) :: v2_1, v2_0, theta_2(1), a3_0, a3_1,&
                c2_0, c2_1, theta_1(1)

    res_1=rcs(a2,l1,1,eps)
    res_0=rcs(a2,l1,0,eps)
    ! print*, res


    theta_1(1)=0.0d0
    CALL NELDER_MEADE1(theta_1,1.0d-8,OBJ_V2_1,verb,1000)
    CALL NELDER_MEADE1(theta_1,1.0d-8,OBJ_V2_1,verb,1000)
    v2_1=v2_1temp
    a3_1=a3_1temp
    c2_1=c2_1temp

    theta_2(1)=0.0d0
    CALL NELDER_MEADE1(theta_2,1.0d-8,OBJ_V2_0,verb,1000)
    CALL NELDER_MEADE1(theta_2,1.0d-8,OBJ_V2_0,verb,1000)
    v2_0=v2_0temp
    a3_0=a3_0temp
    c2_0=c2_0temp

    if (v2_1>v2_0) then
      lab(2,i)=1
      assets(3,i)=a3_1
      c_new(2,i)=c2_1
      ! v_new(3,i)=v2_1
      V2=v2_1
      ! print*, "L=1"
    else
      lab(2,i)=0
      assets(3,i)=a3_0
      c_new(2,i)=c2_0
      ! v_new(2,i)=v2_0
      V2=v2_0
      ! print*, "L=0"
    endif

  end function V2

  real(8) function V2c(a2,l1,eps)
    use SIMPLEX1
    ! use omp_lib
    implicit none
    real(8), INTENT(IN) :: a2, eps
    integer, INTENT(IN) :: l1
    ! integer :: t_lab
    real(8) :: v2_1, v2_0, theta_2(1), a3_0, a3_1,&
                c2_0, c2_1, theta_1(1)

    res_1=rcs(a2,l1,1,eps)
    res_0=rcs(a2,l1,0,eps)
    ! print*, res

    theta_1(1)=0.0d0
    CALL NELDER_MEADE1(theta_1,1.0d-8,OBJ_V2_1,verb,1000)
    CALL NELDER_MEADE1(theta_1,1.0d-8,OBJ_V2_1,verb,1000)
    ! v2_1=v2_1temp
    ! a3_1=a3_1temp
    ! c2_1=c2_1temp

    theta_2(1)=0.0d0
    CALL NELDER_MEADE1(theta_2,1.0d-8,OBJ_V2_0,verb,1000)
    CALL NELDER_MEADE1(theta_2,1.0d-8,OBJ_V2_0,verb,1000)
    ! v2_0=v2_0temp
    ! a3_0=a3_0temp
    ! c2_0=c2_0temp

    if (v2_1temp>v2_0temp) then
      V2c=v2_1temp
    else
      V2c=v2_0temp
    endif

  end function V2c

  subroutine V2s(a2,l1,eps,i)
    use SIMPLEX1
    ! use omp_lib
    implicit none
    real(8), INTENT(IN) :: a2, eps
    integer, INTENT(IN) :: i,l1
    ! integer :: t_lab
    real(8) :: v2_1, v2_0, theta_2(1), a3_0, a3_1,&
                c2_0, c2_1, theta_1(1)
    !
    res_1=rcs(a2,l1,1,eps)
    res_0=rcs(a2,l1,0,eps)
    ! print*, res


    theta_1(1)=0.0d0
    CALL NELDER_MEADE1(theta_1,1.0d-8,OBJ_V2_1,verb,1000)
    CALL NELDER_MEADE1(theta_1,1.0d-8,OBJ_V2_1,verb,1000)
    ! v2_1=v2_1temp
    ! a3_1=a3_1temp
    ! c2_1=c2_1temp

    theta_2(1)=0.0d0
    CALL NELDER_MEADE1(theta_2,1.0d-8,OBJ_V2_0,verb,1000)
    CALL NELDER_MEADE1(theta_2,1.0d-8,OBJ_V2_0,verb,1000)
    ! v2_0=v2_0temp
    ! a3_0=a3_0temp
    ! c2_0=c2_0temp


    if (v2_1temp>v2_0temp) then
      lab(2,i)=1
      assets(3,i)=a3_1temp
      c_new(2,i)=c2_1temp
      v_new(2,i)=v2_1temp
      ! V2=v2_1
      ! print*, "L=1"
    else
      lab(2,i)=0
      assets(3,i)=a3_0temp
      c_new(2,i)=c2_0temp
      v_new(2,i)=v2_0temp
      ! V2=v2_0
      ! print*, "L=0"
    endif

  end subroutine V2s

  REAL(8) FUNCTION OBJ_V1_1(th)
    USE PROBABILITY
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: th(:)
    REAL(8) :: cons,arg1,kprime,arg2, eps, integ
    INTEGER :: j

    cons=EXP(th(1))
    cons=res1_1/(1.0d0+cons)
    ! Consumption lives between 0 and y

    arg1=u(cons,1)

    kprime=res1_1-cons
    integ=0.0d0

    DO j=1,N_quad
       eps=x_quad(j)
       integ=integ + V2c(kprime,1,eps)*PDF_NORMAL(x_quad(j),mu,sigma2_e)*w_quad(j)
    END DO

    arg2=beta*integ

    c1_1temp=cons
    a2_1temp=kprime
    OBJ_V1_1=-arg1-arg2
    v1_1temp=-OBJ_V1_1

  END FUNCTION OBJ_V1_1

  REAL(8) FUNCTION OBJ_V1_0(th)
    USE PROBABILITY
    ! use mpi
    ! use omp_lib
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: th(:)
    REAL(8) :: cons,arg1,kprime,arg2, eps, integ
    INTEGER :: j

    cons=EXP(th(1))
    cons=res1_0/(1.0d0+cons)
    ! Consumption lives between 0 and y

    arg1=u(cons,0)
    ! print*, arg1
    kprime=res1_0-cons

    integ=0.0d0

    DO j=1,N_quad
       eps=x_quad(j)
       integ=integ + V2c(kprime,0,eps)*PDF_NORMAL(x_quad(j),mu,sigma2_e)*w_quad(j)
    END DO

    arg2=beta*integ
    ! print*, arg2, cons, kprime
    ! lt=t_lab
    c1_0temp=cons
    a2_0temp=kprime
    OBJ_V1_0=-arg1-arg2
    v1_0temp=-OBJ_V1_0
  END FUNCTION OBJ_V1_0

  SUBROUTINE T(h,eps)
    USE SIMPLEX1
    use INTEGRATION
    ! use omp_lib
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: h
    real(8), INTENT(IN) :: eps
    ! integer, INTENT(IN) :: l1,i
    ! integer :: t_lab
    real(8) :: v1_1, v1_0, theta_1(1), a2_0, a2_1,&
                c1_0, c1_1, a1, theta_2(1)

    ! sigma2_e=sigma2_e
    res1_1=rcs(assets(1,h),0,1,eps)
    res1_0=rcs(assets(1,h),0,0,eps)

    theta_1(1)=0.0d0
    CALL NELDER_MEADE1(theta_1,1.0d-8,OBJ_V1_1,verb,1000)
    CALL NELDER_MEADE1(theta_1,1.0d-8,OBJ_V1_1,verb,1000)
    ! v1_1=v1_1temp
    ! a2_1=a2_1temp
    ! c1_1=c1_1temp

    theta_2(1)=0.0d0
    CALL NELDER_MEADE1(theta_2,1.0d-8,OBJ_V1_0,verb,1000)
    CALL NELDER_MEADE1(theta_2,1.0d-8,OBJ_V1_0,verb,1000)
    ! v1_0=v1_0temp
    ! a2_0=a2_0temp
    ! c1_0=c1_0temp

    if (v1_1temp>v1_0temp) then
      lab(1,h)=1
      assets(2,h)=a2_1temp
      c_new(1,h)=c1_1temp
      v_new(1,h)=v1_1temp
    else
      lab(1,h)=0
      assets(2,h)=a2_0temp
      c_new(1,h)=c1_0temp
      v_new(1,h)=v1_0temp
    endif

  END SUBROUTINE T

END MODULE GLOBVAR
