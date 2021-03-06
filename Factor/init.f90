MODULE INIT
USE GLOBVAR
!USE MPI
USE RANDOM
IMPLICIT NONE

CONTAINS
  
  SUBROUTINE INITIALIZE(id)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id
    INTEGER, PARAMETER :: N_vars=6     ! Number of columns in dataset
    REAL(8) :: dataset(N_ind,N_vars)
    INTEGER :: l_y,l_x(N_x),l_z,l_d,l_m
    INTEGER :: i,j,ierr,h,t
    
    ! Location of the Y's
    l_y=5
    ! Location of X's
    l_x=(/1,2/)
    l_z=3
    l_m=4
    l_d=6

    CALL SET_SEED(176765183,581076745,294342786,765109823)
    ! Now read in the dataset and fill in the data matrices
    IF (id==0) OPEN(1,file="data.raw")
    DO i=1,N_ind
       IF (id==0) READ(1,fmt=*) dataset(i,:)
!       CALL MPI_BCAST(dataset(i,:),N_vars,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
       x(i,:)=dataset(i,l_x)
       z(i,:)=dataset(i,l_z)
       y(i)=dataset(i,l_y)
       m(i)=dataset(i,l_m)
       d(i)=dataset(i,l_d)
       DO j=1,N_int
          uni(i,j)=Sample_Uniform(0.0d0,1.0d0)
          nor(i,j)=Sample_Normal(0.0d0,1.0d0)
       END DO
    END DO

    IF (id==0) THEN 
       CLOSE(1)
       WRITE(*,*) "dataset read"
       WRITE(*,*) 'done init'
    END IF
    
  END SUBROUTINE INITIALIZE
  
END MODULE INIT
