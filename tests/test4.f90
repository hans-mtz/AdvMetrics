PROGRAM TEST4
  IMPLICIT NONE
  REAL(8) :: a,b,c,d,e,f,g
  CHARACTER(len=100) :: t
  !Let's create some numbers
  a=1.0d0
  b=2.0d0
  c=3.0d0
  OPEN(1,file="output.txt")
  WRITE(1,fmt='(A)') "Writing a,b and b,c in the next line"
  WRITE(1,fmt='(2F16.8)') a,b
  WRITE(1,fmt='(2F16.8)') b,c
  CLOSE(1)
  ! Now let's read from it
  OPEN(1,file="output.txt")
  READ(1,fmt=*) t
  READ(1,fmt=*) d,e
  READ(1,fmt=*) f,g
  CLOSE(1)
  ! Finally write them to the screen so we are sure they worked
  WRITE(*,*) "a,d",a,d
  WRITE(*,*) "b,e",b,e
  WRITE(*,*) "b,f",b,f
  WRITE(*,*) "c,g",c,g
END PROGRAM TEST4
