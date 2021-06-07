module wmod
  use GLOBVAR
  use INTERPOL
  use mod_csv
  use INIT
  ! use mpi
  ! use omp_lib
  implicit none
  public :: q3, q4, q5, obj_ind,arrange_data
contains

  subroutine q3(r_t2)
    ! implicit none
    real(8), INTENT(INOUT) :: r_t2(N_sim,9)
    INTEGER :: j
    ! integer :: num_threads=4

    !$ call omp_set_num_threads(num_threads)
    !$omp parallel do

    DO j=1,N_sim
         assets(2,j)=inter_value%v(1)%x(j)
         v_new(2,j)=V2(assets(2,j),0,0.0d0,j)
     END DO
     !$omp end parallel do
     r_t2(:,1)=0
     r_t2(:,4)=v_new(2,:)
     r_t2(:,3)=assets(2,:)
     r_t2(:,6)=assets(3,:)
     r_t2(:,8)=lab(2,:)

     !$OMP PARALLEL DO
     DO j=1,N_sim
         ! assets(2,j)=inter_value%v(1)%x(j)
         v_new(2,j)=V2(assets(2,j),1,0.0d0,j)
     END DO
     !$OMP end parallel do
     r_t2(:,2)=1
     r_t2(:,5)=v_new(2,:)
     ! r_t2(:,6)=reshape(assets(2,:),[N_grid,1])
     r_t2(:,7)=assets(3,:)
     r_t2(:,9)=lab(2,:)

     ! call write_file("q3.csv",r_t2)

     call csvwrite("res",r_t2)
     print*, "Done writing to file results of q3"
 end subroutine q3

  subroutine q4(r_t1)
      ! implicit none
      real(8), INTENT(INOUT) :: r_t1(N_sim, 4)
      integer :: j
      ! integer :: num_threads=4

      !$ call omp_set_num_threads(num_threads)

      !$OMP PARALLEL DO
      DO j=1,N_sim
         assets(1,j)=inter_value%v(1)%x(j)
         ! res1_1=rcs(assets(1,j),0,1,0.0d0)
         ! res1_0=rcs(assets(1,j),0,0,0.0d0)
         CALL T(j, 0.0d0)
      END DO
      ! DO j=1,N_grid
      !     print*, j,v_new(1,j), assets(1,j), assets(2,j), lab(1,j)
      ! END DO

      r_t1(:,1)=assets(1,:)
      r_t1(:,2)=assets(2,:)
      r_t1(:,3)=v_new(1,:)
      r_t1(:,4)=lab(1,:)

      call csvwrite("res_t1",r_t1)
      print*, "Done writing results of q4"
      ! print*, u(c_new(2,1), 1),u(c_new(2,1), 0), w(assets(2,1)), rcs(assets(2,1),lab(2,1), 0.0d0)
      ! print*, OBJ_V2_0(c_new(2,:)), OBJ_V2_1(c_new(2,:)), V2(assets(1,1),0,0.0d0,0)

  end subroutine q4

  subroutine q5(theta)
      use SIMPLEX1
      implicit none
      real(8), intent(inout) :: theta(:)
      ! For each individual
      !1) Take lambda and A1 from data
      !2) for a guess of sigma_e,gamma, alpha, mu, and delta
      ! simulate fake data:
        !2.1) for t=1, given A1 and lambda, and the realisation
        ! of eps1, choose L1 and A2
        !2.2) for t=2, given A2, L1, and the realisation of
        !eps2, choose L2 and A3.
      !3) From fake data, estimate rho_h(15)
      !4) evaluate objective function
      ! use GLOBVAR
      ! use INIT

      ! assets(1,:)=dataset(N_sim, 1)

      CALL NELDER_MEADE1(theta,1.0d-4,obj_ind,3,1000)
      CALL NELDER_MEADE1(theta,1.0d-4,obj_ind,3,1000)
      ! if (myid .eq. 0) THEN
      !   theta([1,2,3,5])=1.4**(theta([1,2,3,5]))
      !   ! theta([1,2,3,5])=logb(theta([1,2,3,5]),1.4d0)
      !   print*, "theta = ", theta!, "from proc ",myid
      !   print*, "obj fun value = ", valuetemp!, "from proc ",myid
      !   ! print*, t4, " seconds for q5", "from proc ",myid
      !   ! flush(1)
      !   call csvwrite("theta",theta)
      !   call csvwrite("fake_df_mpi",f_data)
      !   print*, "Done writing results of q5"
      ! endif
      ! call csvwrite("theta",theta)
      ! call csvwrite("fake_df",f_data)
      !
      ! print*, "Obj value = ", valuetemp
      ! print*, "theta is =", theta
  end subroutine q5

  real(8) function obj_ind(theta)
      use mpi
      use INTEGRATION
      use BLAS95
      implicit none
      real(8), intent(in) :: theta(:)
      real(8):: eps_v(1:N_sim,2)!, sigma_e
      real(8):: df(1:N_sim,7), eps1, eps2,a2, diff(15)
      integer :: i, FIRST, LAST, nbar, l1


      ! print *, 'Before MPI_WTIME 2: Process ', myid, ' of ', numprocs, ' is alive'

      if (myid .eq. 0) then
          sigma2_e=1.4**(theta(1))
          gamma=1.4**(theta(2))
          ! gamma=0.99d0/(1.0d0+gamma)
          alpha=1.4**(theta(3))
          ! alpha=0.99d0/(1.0d0+alpha)
          mu=theta(4)
          delta=1.4**(theta(5))
          call simulate(eps_v,sigma2_e)

          CALL GAUSS_LEGENDRE(x_quad,w_quad,-4.0d0*SQRT(sigma2_e),4.0d0*sqrt(sigma2_e))
          ! print*, "Sending eps_v1",eps_v([1,20,21,40,41,60,61,80,81,100,101,120],1)
          ! print*, "Sending eps_v2",eps_v([1,20,21,40,41,60,61,80,81,100,101,120],2)
          ! print*, "Sending: ", sigma2_e, gamma, alpha, mu, delta
      endif

      call MPI_BCAST(sigma2_e,1,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gamma,1,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(alpha,1,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(mu,1,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(delta,1,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(eps_v(1:N_sim,1),N_sim,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(eps_v(1:N_sim,2),N_sim,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(x_quad,N_quad,MPI_REAL8,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(w_quad,N_quad,MPI_REAL8,0,MPI_COMM_WORLD, ierr)

      if (ierr .ne. mpi_success) then
        print*, "Error: MPI_BCAST not succesful obj fun"
        stop
      endif
      ! print *, 'MPI_BCAST: Process ', myid, ' of ', numprocs, ' is alive'
      nbar=N_sim/numprocs
      first=(myid*nbar) + 1
      last=nbar*(myid + 1)
      ! print*, "Sending from ", myid, "first ",FIRST,"-last ",LAST

      ! print*, "l_v: Sending from ", myid, "[",l_v(FIRST),",",l_v(LAST),"]"
      ! print*, "eps_v1: Sending from ", myid, "[",eps_v(FIRST,1),",",eps_v(LAST,1),"]"
      ! print*, "eps_v2: Sending from ", myid, "[",eps_v(FIRST,2),",",eps_v(LAST,2),"]"
      ! print*, "receiving to ",myid," : " ,sigma2_e, gamma, alpha, mu, delta
      ! do i=1,N_sim!, numprocs
      DO i=FIRST,LAST
          lambda=l_v(i)
          eps1=eps_v(i,1)
          call T(i,eps1)
          eps2=eps_v(i,2)
          a2=assets(2,i)
          l1=lab(1,i)
          call V2s(a2,l1,eps2,i)
      end do

      ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

      call mpi_gather(assets(2,FIRST:LAST),nbar,MPI_REAL8,f_data(:,4),nbar,MPI_REAL8,0,MPI_COMM_WORLD, IERR)
      call mpi_gather(DBLE(LAB(1,FIRST:LAST)),nbar,MPI_REAL8,f_data(:,2),nbar,MPI_REAL8,0,MPI_COMM_WORLD, IERR)
      call mpi_gather(assets(3,FIRST:LAST),nbar,MPI_REAL8,f_data(:,7),nbar,MPI_REAL8,0,MPI_COMM_WORLD, IERR)
      call mpi_gather(DBLE(LAB(2,FIRST:LAST)),nbar,MPI_REAL8,f_data(:,5),nbar,MPI_REAL8,0,MPI_COMM_WORLD, IERR)

      if (ierr .ne. mpi_success) then
        print*, "Error: MPI_gather not succesful obj fun"
        stop
      endif
      ! print *, 'Gather: Process ', myid, ' of ', numprocs, ' is alive'

      ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
      ! if (ierr .ne. mpi_success) then
      !   print*, "Error: MPI_BARRIER not succesful main"
      !   stop
      ! endif
      if (myid .eq. 0) then
          ! print*, "Received",f_data([1,20,21,40,41,60,61,80,81,100,101,120],4)
          ! y(1,:)=exp(mu+eps_v(:,1))
          ! y(2,:)=exp(mu+delta*DBLE(lab(1,:))+eps_v(:,2))
          f_data(:,1)=assets(1,:)
          f_data(:,3)=y_fun(0,INT(f_data(:,2)),eps_v(:,1))
          f_data(:,6)=y_fun(INT(f_data(:,2)),INT(f_data(:,5)),eps_v(:,2))

          ! call arrange_data(df)

          ! f_data=df

          call get_rho(rho_h, f_data)

          diff=rho_h-rho_0

          obj_ind=nrm2(diff)
          valuetemp=obj_ind
      end if

      call MPI_BCAST(f_data,N_sim*7, MPI_REAL8,0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(valuetemp,1, MPI_REAL8,0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(obj_ind,1, MPI_REAL8,0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(rho_h,15, MPI_REAL8,0, MPI_COMM_WORLD, ierr)
      if (ierr .ne. mpi_success) then
        print*, "Error: MPI_BCAST not succesful main "
        stop
      endif
      ! print *, 'Bcast 2: Process ', myid, ' of ', numprocs, ' is alive'

  end function obj_ind

  subroutine arrange_data(df)
      implicit none
      real(8), intent(inout) :: df(:,:)

      df(:,[1,4,7])=transpose(assets)
      df(:,2)=DBLE(lab(1,:))
      df(:,5)=DBLE(lab(2,:))
      df(:,3)=y(1,:)
      df(:,6)=y(2,:)

  end subroutine arrange_data


end module wmod
