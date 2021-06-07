MODULE GLOBVAR
  USE INTERPOL
  IMPLICIT NONE
  ! This modlue declares the variables that are going to be used throughout the program
  ! higher level languages do not require you to declare your variables ahead of time
  
  INTEGER, PARAMETER :: N_grid=20
  INTEGER, PARAMETER :: N_quad=16
  REAL(8), PARAMETER :: alpha=0.4d0, &
       beta=0.96d0, &
       mu=0.0d0, &
       sigma=0.1d0
  
  REAL(8) :: y_grid(N_grid),v_old(N_grid),y,constemp,valuetemp,var,v_old_cons(N_grid)
  REAL(8) :: x_quad(N_quad),w_quad(N_quad),v_new(N_grid),c_new(N_grid)
  TYPE(PLINEARSTRUC) :: inter_value

CONTAINS

  SUBROUTINE T(h)
    USE SIMPLEX
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: h
    REAL(8) :: theta(1)
  
    y=inter_value%v(1)%x(h)
    theta(1)=0.0d0
    CALL NELDER_MEADE(theta,1.0d-8,OBJECTIVE,0,1000)
    CALL NELDER_MEADE(theta,1.0d-8,OBJECTIVE,0,1000)
    v_new(h)=valuetemp
    c_new(h)=constemp
  END SUBROUTINE T
  
  ! We need to know what y is inside this function
  ! I have a guess interpolated of the value function somewhere
  ! and this function needs access to it
  ! Somehow this routine needs to know the quadrature points
  REAL(8) FUNCTION OBJECTIVE(theta)
    USE PROBABILITY
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: cons,arg1,kprime,integ,yprime,arg2
    INTEGER :: j
    cons=EXP(theta(1))
    cons=y/(1.0d0+cons)
    ! Consumption lives between 0 and y
    IF (cons<1.0d-16) THEN
       arg1=-1.0d300
    ELSE
       arg1=LOG(cons)
    END IF
    kprime=y-cons
    integ=0.0d0
    DO j=1,N_quad
       yprime=(kprime**alpha)*EXP(x_quad(j))
       integ=integ + PLINEAR_INTER(inter_value,yprime,v_old)*PDF_NORMAL(x_quad(j),mu,var)*w_quad(j)
    END DO
    arg2=beta*integ
    constemp=cons
    OBJECTIVE=-arg1-arg2
    valuetemp=-OBJECTIVE
  END FUNCTION OBJECTIVE

  REAL(8) FUNCTION OBJECTIVE_CONS(cons)
    USE PROBABILITY
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: cons
    REAL(8) :: arg1,kprime,integ,yprime,arg2
    INTEGER :: j
    ! Consumption lives between 0 and y
    IF (cons<1.0d-16) THEN
       arg1=-1.0d300
    ELSE
       arg1=LOG(cons)
    END IF
    kprime=y-cons
    integ=0.0d0
    DO j=1,N_quad
       yprime=(kprime**alpha)*EXP(x_quad(j))
       integ=integ + PLINEAR_INTER(inter_value,yprime,v_old_cons)*PDF_NORMAL(x_quad(j),mu,var)*w_quad(j)
    END DO
    arg2=beta*integ
    OBJECTIVE_CONS=+arg1+arg2
  END FUNCTION OBJECTIVE_CONS

    END MODULE GLOBVAR
