MODULE NORMAL
  IMPLICIT NONE
CONTAINS
  FUNCTION FUN_PDF_NORMAL(x,mu,var)
    IMPLICIT NONE
    !This says that x,mu and var are inputs to the function
    !that are NOT ALLOWED to change. They cannot suffer any
    !alterations inside the function
    REAL(8), INTENT(IN) :: x,mu,var
    REAL(8) :: FUN_PDF_NORMAL ! This will be the output of the function
    ! Declaring SQRT(2*pi) as a parameter. I do this so that, again,
    ! I cannot change it by accident in the function. Not only that,
    ! this way I do not do a SQRT operation everytime I
    ! call this function
    REAL(8), PARAMETER :: root2pi=2.506628275D0
    REAL(8) :: arg,den
    arg = (x-mu)
    arg=arg**2
    arg=arg/var ! Create ((x-mu)^2)/var
    den = SQRT(var)*root2pi
    FUN_PDF_NORMAL = exp(-0.5*arg)/den
  END FUNCTION FUN_PDF_NORMAL

  SUBROUTINE SUB_PDF_NORMAL(x,mu,var,out)
    IMPLICIT NONE
    !This says that x,mu and var are inputs to the function
    !that are NOT ALLOWED to change. They cannot suffer any
    !alterations inside the function
    REAL(8), INTENT(IN) :: x,mu,var
    REAL(8), INTENT(OUT) :: out ! This will the output
    ! Declaring SQRT(2*pi) as a parameter. I do this so that, again,
    ! I cannot change it by accident in the function. Not only that,
    ! this way I do not do a SQRT operation every time I
    ! call this function
    REAL(8), PARAMETER :: root2pi=2.506628275D0
    REAL(8) :: arg,den
    arg = (x-mu)
    arg=arg**2
    arg=arg/var ! Create ((x-mu)^2)/var
    den = SQRT(var)*root2pi
    out = exp(-0.5*arg)/den
  END SUBROUTINE SUB_PDF_NORMAL
END MODULE NORMAL
