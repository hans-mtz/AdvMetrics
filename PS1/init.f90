MODULE INIT
  IMPLICIT NONE
CONTAINS
  
  SUBROUTINE READDATA()
    USE GLOBVAR
    IMPLICIT NONE	
    INTEGER, PARAMETER :: N_variables=4		! Number of variables in the dataset	
    REAL(8) :: dataset(N_individuals,N_variables)	! Array where we will read in the data
    INTEGER :: i,&
         l_covariates(N_covariates), &          ! space to assign which columns of the dataset contain the X's
         l_y					! space to assign which column contains Y
    
    OPEN(1,file='data.raw')	                ! opening the file dataset.txt where the dataset is saved in 
    DO i = 1, N_individuals		        ! ascii format and asigning it to handle 1
       READ(1,fmt=*)	dataset(i,:)	        ! reading from handle 1 into the array I created
    END DO					! That is, loading the dataset into memory
    CLOSE(1)					! Closing the file since now I have its contents in the dataset array
    WRITE(*,*) "dataset read"
    
    ! Now let's assign the data  
    l_covariates=(/2,3/)	                ! Telling it which columns contain the X's
    l_y=1	                                ! first column of the dataset contains the Y's
    x=dataset(:,l_covariates)	                ! assigning X
    y=dataset(:,l_Y)	                        ! assigning Y
    ! Done reading my data, now moving to the main program
    
  END SUBROUTINE READDATA
  
END MODULE INIT
