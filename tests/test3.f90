PROGRAM TEST3
  USE NORMAL ! So that it recognizes the functions I have in the module
  IMPLICIT NONE
  REAL(8) :: a,c,d,e(3)
  ! I do not know the size of b, but I know it will be one dimensional
  REAL(8), POINTER :: b(:)
  INTEGER :: j
  ! Now I allocate b
  ALLOCATE(b(3)) ! So now b is a vector of size 3
  b=(/1.0d0,1.10d0,2.0d0/)
  e=b
  a=2.0d0
  ! Two ways of using what I have in normal
  DO j=1,3
     c=FUN_PDF_NORMAL(0.0d0,0.0d0,1.0d0)*a
     CALL SUB_PDF_NORMAL(0.0d0,0.0d0,1.0d0,d)
     c=a*d
  END DO
  DEALLOCATE(b)
  ALLOCATE(b(2)) ! So now b is size 2
WRITE(*,*) c
END PROGRAM TEST3
