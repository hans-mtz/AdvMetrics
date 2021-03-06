PROGRAM MAIN
  USE RANDOM
  USE INTEGRATION
  ! USE sobol, only: r8i8_uniform_01
  IMPLICIT NONE

  INTEGER, PARAMETER :: n1=6,n2=50
  REAL(8), PARAMETER :: a=-1.0d0,b=2.0d0
  REAL(8) :: int1,int2,int3,int4,int5,int6,int7,int12
  INTEGER :: j
  integer(8) :: i
  REAL(8) :: quadx(n1),quadw(n1),c,m,t

  !Suppose we want to integrate x^2 between -1 and 1
  !First Gauss Legendre ignoring that my code already does the change of variables
  CALL GAUSS_LEGENDRE(quadx,quadw,-1.0d0,1.0d0)
  m=0.5d0*(b-a)
  c=0.5d0*(b+a)
  quadw=quadw*m
  int1=0
  DO j=1,n1
      t=c+m*quadx(j)
      int1=int1+(t*t*quadw(j))
  END DO

  CALL GAUSS_LEGENDRE(quadx,quadw,-1.0d0,2.0d0)
  ! m=0.5d0*(b-a)
  ! c=0.5d0*(b+a)
  ! quadw=quadw*m
  int12=0
  DO j=1,n1
      ! t=c+m*quadx(j)
      t=quadx(j)
      int12=int12+(t*t*quadw(j))
  END DO

  ! Now Gauss-Chebyshev
  CALL GAUSS_CHEBYSHEV(quadx,quadw)
  m=0.5d0*(b-a)
  c=0.5d0*(b+a)
  quadw=quadw*m
  int2=0.0d0
  DO j=1,n1
      t=c+m*quadx(j)
      int2=int2 + SQRT((1.0d0-quadx(j)*quadx(j)))*t*t*quadw(j)
  END DO

  ! Now Monte Carlo
  ! First with 10 points
  int3=0.0d0
  DO j=1,n1
      m=Sample_Uniform(-1.0d0,2.0d0)
      int3=int3+3.0d0*m*m
  END DO
  int3=int3/DBLE(n1)

  ! Now with 50 points
  int4=0.0d0
  DO j=1,n2
      m=Sample_Uniform(-1.0d0,2.0d0)
      int4=int4+3.0d0*m*m
  END DO
  int4=int4/DBLE(n2)

  ! And now with a Halton Sequence instead of a uniform
  int5=0.0d0
  DO j=100,n1+100
      CALL HALTON(j,m)
      m=m*3.0d0-1.0d0
      int5=int5+3.0d0*m*m
  END DO
  int5=int5/DBLE(n1)

  ! And now with a Halton Sequence instead of a uniform
  int6=0.0d0
  DO j=1000,n2+1000
      CALL HALTON(j,m)
      m=m*3.0d0-1.0d0
      int6=int6+3.0d0*m*m
  END DO
  int6=int6/DBLE(n2)

  WRITE(*,'(7F10.5)') 2.33d0,int1,int12!,int2,int3,int4,int5,int6

  ! And now with a Sobol Sequence instead of a uniform
  ! int7=0.0d0
  ! DO i=1,n2
  !     m=r8i8_uniform_01(i)
  !     m=m*3.0d0-1.0d0
  !     int7=int7+3.0d0*m*m
  ! END DO
  ! int7=int7/DBLE(n2)
  ! print*, "Now with Sobol sequence this is the integral",int7

END PROGRAM MAIN
