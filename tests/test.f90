!My perverted version of Hello World
PROGRAM MAIN ! All programs start with the program declaration
IMPLICIT NONE ! Always use it, bad things happen if you do not
CHARACTER(len=100) :: name ! Declarations of vari- ables by type
!name = 'Hans'
WRITE(*,*) 'Please give me your name' ! Ask for name READ(*,*) name
READ*, name
WRITE(*,*) trim(name),' Stinks!, Jajaja'
END PROGRAM MAIN
