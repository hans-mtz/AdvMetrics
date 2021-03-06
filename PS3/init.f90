MODULE INIT
   USE GLOBVAR
   USE RANDOM
   use matrix
   ! use LAPACK95
   IMPLICIT NONE
   public :: READDATA, initiate
   private

CONTAINS

   SUBROUTINE READDATA()
       ! USE GLOBVAR
       ! IMPLICIT NONE
       INTEGER, PARAMETER :: N_variables=4		! Number of variables in the dataset
       REAL(8) :: dataset(N_individuals,N_variables), c(N_individuals,1)	! Array where we will read in the data
       INTEGER :: i,j,&
            l_covariates(N_covariates - 1), &          ! space to assign which columns of the dataset contain the X's
            l_y					! space to assign which column contains Y


       OPEN(1,file='data_ps2.out')	                ! opening the file dataset.txt where the dataset is saved in
       DO i = 1, N_individuals		        ! ascii format and asigning it to handle 1
          READ(1,fmt=*)	dataset(i,:)	        ! reading from handle 1 into the array I created
       END DO					! That is, loading the dataset into memory
       CLOSE(1)					! Closing the file since now I have its contents in the dataset array
       WRITE(*,*) "dataset read"

       ! Now let's assign the data
       l_covariates=(/2,3/)
       c=1	                ! Telling it which columns contain the X's
       ! c=reshape(c,(/N_individuals,1/))
       l_y=1	                                ! first column of the dataset contains the Y's
       x(:,l_covariates)=dataset(:,l_covariates)
       x(:,1)=c(:,1)	                ! assigning X
       y=dataset(:,l_Y)	                        ! assigning Y
       D=dataset(:,4)
       ! Done reading my data, now moving to the main program
       ! call initiate()

   END SUBROUTINE READDATA

   subroutine initiate()
      integer :: i,j

      CALL SET_SEED(3,4,2,1)

      ! Sampling from U_I ~ N(0,1)
      do i=1,N_individuals
        do j=1,s
           ! uni(i,j)=Sample_Uniform(1.0d0, DBLE(N_individuals))
           nor(i,j)=Sample_Normal(0.0d0,1.0d0)
        enddo
      enddo

      do i=1,N_individuals
         nor_v(i)=Sample_Normal(0.0d0,1.0d0)
      enddo

      do i=1,N_individuals
        do j=1,n_boot
           uni(i,j)=Sample_Uniform(1.0d0, DBLE(N_individuals))
           ! nor(i,j)=Sample_Normal(0.0d0,1.0d0)
        enddo
      enddo

   end subroutine initiate

  subroutine shuffle(sim, hal)
     ! use GLOBVAR
     integer, INTENT(IN) :: sim
     integer :: i,j,k
     REAL(8) :: up,s, temp(N_individuals,N_covariates), tempd(N_individuals)
     logical, intent(in) :: hal
     up= REAL(N_individuals)

     call READDATA()
     ! print*, size(x)
     ! print*, (sum(x(:,i))/N_individuals, i=1,3)
     temp=x
     ! x=0.0d0
     tempd=D
     ! xbar=[(sum(x(:,i)),i=1,N_covariates)]
    ! if (boost) then
      do i=1,N_individuals
         if (hal) then
            call HALTON(i+ 100+sim*33,s)
            s=1.0d0 +(50000.0d0-1.0d0)*s
         else
            call set_seed(sim*333,666,3505,i*111)
            s=Sample_Uniform(1.0d0,up)
         end if
         j=int(s)
         forall(k=1:3) x(i,k)=temp(j,k)!-xbar
         D(i)=tempd(j)
      end do
    ! endif
    ! print*,size(x)
    ! print*, size(x)
    ! print*, (sum(x(:,i))/N_individuals, i=1,3)
  end subroutine shuffle

END MODULE INIT
