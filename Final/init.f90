MODULE INIT
   ! USE GLOBVAR
   ! USE RANDOM
   ! use matrix
   ! use NRUTIL
   ! use LAPACK95, only : ORMQR, GELS, GEQRF
   ! use BLAS95, only : TRSM
   IMPLICIT NONE
   public :: READDATA, initiate, simulate, get_rho, logb
   private

CONTAINS

   SUBROUTINE READDATA(ass,tdf)
       ! USE GLOBVAR
       IMPLICIT NONE
       INTEGER, PARAMETER :: N_variables=7, N_obs=5000		! Number of variables in the dataset
       REAL(8) :: dataset(1:N_obs,1:N_variables)	! Array where we will read in the data
       real(8), intent(inout) :: ass(:,:), tdf(:,:)!, rho(15)
       ! real(8) :: y1(N_obs), lny1(N_obs), y2(N_obs)
       ! real(8) :: x1(N_obs,4), x2(N_obs,4), a2(N_obs), a3(N_obs)
       ! real(8) :: cov1(3), cov2(3)
       ! real(8) :: r1, r2, r3, r4
       INTEGER :: i,j
       ! real(8), DIMENSION(:,:), ALLOCATABLE :: x3
       ! real(8), DIMENSION(:), ALLOCATABLE :: lny2


       OPEN(1,file='data.out')	                ! opening the file dataset.txt where the dataset is saved in
       DO i = 1, N_obs	        ! ascii format and asigning it to handle 1
          READ(1,fmt=*)	dataset(i,:)	        ! reading from handle 1 into the array I created
       END DO					! That is, loading the dataset into memory
       CLOSE(1)					! Closing the file since now I have its contents in the dataset array
       WRITE(*,*) "dataset read"

       ass(1,:)=dataset(1:size(ass,2), 1)

       tdf=dataset

       ! call get_rho(rho, dataset)
       ! print*, "rho_0 ready"
   END SUBROUTINE READDATA

   subroutine initiate(v,sig2)
      use random
      integer :: i
      real(8), INTENT(INOUT) :: v(:)
      real(8), intent(in) :: sig2

      CALL SET_SEED(3,4,2,1)

      do i=1,5000
         v(i)=exp(Sample_Normal(0.0d0,sqrt(sig2)))
      enddo

   end subroutine initiate

   subroutine simulate(eps_v,sig2)
      use random
      implicit none
      real(8), intent(in) :: sig2
      real(8), INTENT(INOUT) :: eps_v(:,:)
      integer :: i

      CALL SET_SEED(361,631,691,196)
      ! CALL SET_SEED()

      ! var=sigma*sigma
      do i=1,size(eps_v(:,1))
         eps_v(i,1)=Sample_Normal(0.0d0,sqrt(sig2))
         eps_v(i,2)=Sample_Normal(0.0d0,sqrt(sig2))
      enddo

   end subroutine simulate

   subroutine get_rho(rho,df)
      use matrix
      use nrutil
      ! use wmod1
      implicit none
      real(8), intent(in) :: df(:,:)
      real(8), intent(inout) :: rho(15)
      integer :: nd !=N_sim!=size(df(:,1))
      real(8) :: y1(1:size(df(:,1))), lny1(1:size(df(:,1)))
      real(8) :: x1(1:size(df(:,1)),1:4), x2(1:size(df(:,1)),1:4), a2(1:size(df(:,1))), a3(1:size(df(:,1)))
      real(8) :: cov1(3), cov2(3), nu1(size(df,1)), x1b(size(df,1))
      real(8) :: nu2(size(df,1)), x2b(size(df,1))
      real(8) :: x3b(size(df,1)), pff(2)
      INTEGER :: i,j,l2dim
      logical :: pick(size(df,1))
      real(8), DIMENSION(:,:), ALLOCATABLE :: x3
      real(8), DIMENSION(:), ALLOCATABLE :: lny2, y2, nu3
      ! real(8) :: x3(size(df, dim=1),1:2), lny2(size(df, dim=1))
      ! Now let's assign the data

      ! nd=N_sim
      nd=size(df,1)
      print*, "df size", nd
      cov1=(/3,1,2/)
      x1(:,1)=1.0d0	                ! Telling it which columns contain the X's
      x1(:,2:4)=df(:,cov1)
      a2=df(:,4)	          ! first column of the df contains the Y's
      rho(1:4)=LM(x1,a2)

      ! print "(a,4f16.6)", "SN OLS", rho(1:4)
      ! print "(a,4f16.6)", "my OLS", myOLS(x1,a2)

      ! rho(5)=sum((a2-MATMUL(transpose(x1),rho(1:4)))**2.0d0)/DBLE(nd)
      x1b=mmul(x1,rho(1:4))
      nu1=a2-x1b
      rho(5)=sum((nu1)**2)/DBLE(nd)-(sum(nu1)/dble(nd))**2
      cov2=(/6,4,5/)
      x2(:,1)=1.0d0
      x2(:,2:4)=df(:,cov2)
      a3=df(:,7)
      rho(6:9)=LM(x2,a3)
      ! print*,"6-9 OK"
      x2b=mmul(x2,rho(6:9))
      nu2=a3-x2b
      rho(10)=sum((nu2)**2)/DBLE(nd)-(sum(nu2)/dble(nd))**2
      ! print*,"10 OK"
      rho(11)=sum(log(df(:,3)), mask=(df(:,2)>0))/sum(df(:,2))

      rho(12)=sum((log(df(:,3))-rho(11))**2.0d0,mask=(df(:,2)>0))/sum(df(:,2))

      pick=(df(:,5)>0)
      l2dim=count(pick)
      ALLOCATE(x3(l2dim,1:2))
      ALLOCATE(lny2(1:l2dim))
      ALLOCATE(y2(1:l2dim))
      ALLOCATE(nu3(1:l2dim))
      print*, "L2==1", l2dim
      x3(1:l2dim,1)=1.0d0
      x3(:,2)=pack(df(:,2), mask=pick)
      y2=pack(df(:,6), mask=pick)
      lny2=log(y2)

      ! print "(3f16.6)",(x3(i,:),lny2(i),i=1,20)

      ! print*, "size x3 and lny2", size(x3,1), size(lny2)
      ! rho(13:14)=OLS(transpose(x3),lny2)
      rho(13:14)=LM(x3,lny2)
      ! pff=myOLS(x3,lny2)
      ! x3(:,1)=1.0d0
      ! x3(:,2)=df(:,2)
      ! lny2=log(df(:,6))
      ! rho(13:14)=OLS(x3,lny2,pick,l2dim)

      ! print*,"13-14 OK", x3(14:18,2)
      ! print "(3f16.6)",(x3(i,:),lny2(i),i=1,20)

      ! print "(3(2f16.6,/))",rho(13:14), myOLS(x3,lny2), pff

      ! x3b=[(DOT_PRODUCT(x3(i,:),rho(13:14)),i=1,l2dim)]

      x3b=mmul(x3,rho(13:14))
      ! print*,"matmul OK"!, (x3b(i),i=1,10)

      nu3=lny2-x3b
      ! print*, "nu3 OK"
      rho(15)=sum((nu3)**2)/dble(l2dim)-(sum(nu3)/dble(l2dim))**2
      ! print*, "15 OK"

      DEALLOCATE(x3,y2,lny2,nu3)
      print*, "rho ready"

   end subroutine get_rho

   function LM(xin,yin) result(beta)
      use nrutil
      use LAPACK95, only : POTRF, POTRI, getrf, getri
      ! use BLAS95, only : trsm
      implicit none
      real(8), intent(in) :: xin(:,:), yin(:)
      real(8), DIMENSION(size(xin, dim=2),size(xin, dim=2)) :: xx
      real(8) :: beta(size(xin,2)), xt(size(xin,2),size(xin,1)), xy(size(xin,dim=2))
      real(8) :: xout(size(xin,1),size(xin,2)),yout(size(xin,1))
      real(8) :: tau(size(xin,2))
      integer :: info, rank, ipiv(size(xin,2))
      ! ALLOCATE(xx(size(xin, dim=2),size(xin, dim=2)))
      xt=transpose(xin)
      xx=mmul(xt,xin)
      xy=mmul(yin,xin)

      ! call POTRF(xx)
      ! call POTRI(xx)

      call getrf(xx, ipiv)
      call getri(xx, ipiv)

      beta=mmul(xx,xy)
      ! xout=xin
      ! yout=yin
      !
      ! ! call gelsd(xout,yout,rank=rank,s=s,info=info)
      ! ! call gels(xout,yout,"N",info=info)
      !
      ! call geqrf(xout, tau, info)
      ! call ORMQR(xout, tau, yout, "L", "T")
      ! call TRSM(xout,yout,"L","U")

      ! print*, yout

      ! if (info .ne. 0) stop "Unsuccessful regression"
      !
      ! beta=yout(1:size(xin,2))

   end function LM

   elemental real(8) function logb(x,b)
     implicit none
     real(8), intent(in) :: x, b
     ! integer, intent(in) :: b
     logb = log(x) / log(b)
   end function

END MODULE INIT
