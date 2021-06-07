PROGRAM MAIN
  USE GLOBVAR
  use INIT
  ! USE INTEGRATION
  USE INTERPOL
  USE IO
  use mod_csv
  use wmod
  use mpi
  IMPLICIT NONE

  INTEGER :: j, h, rc, i
  REAL(8) :: limits(2),arg1, a
  ! REAL(8) :: c1,c2,c3,c4,v_true(N_grid),vold2(N_grid)
  real(8) :: r_t2(N_sim,9), r_t1(N_sim,4), theta(5)
  real*8 :: t1, t2, t3,t4, test_data(5000,7)

  !
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
  ! print *, 'Process ', myid, ' of ', numprocs, ' is alive'


  ! limits(1)=0.0d0
  ! limits(2)=5.0d0
  ! CALL PLINEAR_INIT(inter_value,N_grid,limits,1)

  ! call CPU_TIME(t1)
  !
  ! call q3(r_t2)
  !
  ! call CPU_TIME(t2)
  !
  ! call CPU_TIME(t3)
  !
  ! call q4(r_t1)
  !
  ! call CPU_TIME(t4)
  !
  ! print*, t2-t1, " seconds for q3"
  ! print*, t4-t3, " seconds for q4"

  if (myid .eq. 0) then
    call READDATA(assets, test_data)
    !
    ! print*, "printing rho_0 = ",rho_0
    ! call get_rho(rho_h,test_data)
    ! print*, "rho_0 from fun =  ", rho_h
    call initiate(l_v,sigma2_l)
    ! print*, "Received",l_v([1,20,21,40,41,60,61,80,81,100,101,120])

    ! print*, "printing lambda vector = ", l_v(1:100:10)

    call read_file("point_h.out",theta)
    theta([1,2,3,5])=logb(theta([1,2,3,5]),1.4d0)

    call get_rho(rho_0, test_data(1001:(1000+N_sim),:))
  endif

  ! print *, 'Process ', myid, ' of ', numprocs, ' is alive'

  call MPI_BCAST(assets(1,1:N_sim),N_sim, MPI_REAL8,0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(rho_0,15, MPI_REAL8,0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(theta,15, MPI_REAL8,0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(l_v,5000, MPI_REAL8,0, MPI_COMM_WORLD, ierr)

  if (ierr .ne. mpi_success) then
    print*, "Error: MPI_BCAST not succesful main "
    stop
  endif

  ! print *, 'BCAST: Process ', myid, ' of ', numprocs, ' is alive'

  ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  ! if (ierr .ne. mpi_success) then
  !   print*, "Error: MPI_BARRIER not succesful main"
  !   stop
  ! endif
  ! print *, 'Barrier: Process ', myid, ' of ', numprocs, ' is alive'

  T1=MPI_WTIME()
  ! print *, 'Before MPI_WTIME 1: Process ', myid, ' of ', numprocs, ' is alive'

  a=obj_ind(theta)
  ! call q5(theta)

  ! print *, 'Before MPI_WTIME 2: Process ', myid, ' of ', numprocs, ' is alive'

  t2= MPI_WTIME()
  t3=t2-t1

  ! print *, 'Before Loop: Process ', myid, ' of ', numprocs, ' is alive'

  if (myid .eq. 0) THEN
    theta([1,2,3,5])=1.4**(theta([1,2,3,5]))
    ! theta([1,2,3,5])=logb(theta([1,2,3,5]),1.4d0)

    print "(a,5f16.2)", "theta = ", theta!, "from proc ",myid
    print "(a,f16.6)", "obj fun value = ", valuetemp!, "from proc ",myid
    ! print "(a,f16.6)", "obj fun value a = ", a
    print "(i3, 2f16.6)", (i, rho_0(i), rho_h(i), i=1,15)
    ! print*, t4, " seconds for q5", "from proc ",myid
    ! flush(1)
    call csvwrite("theta",theta)
    call csvwrite("fake_df_mpi_6b",f_data)
  endif

  call MPI_REDUCE(t3,t4,1, MPI_REAL8, mpi_max, 0, MPI_COMM_WORLD, ierr)
  if (ierr .ne. mpi_success) then
    print*, "Error: MPI_REDUCE not succesful main"
    stop
  endif

  ! print *, 'REDUCE: Process ', myid, ' of ', numprocs, ' is alive'

  if (myid==0) print*,  t4, " seconds for q5"
  ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  ! if (ierr .ne. mpi_success) then
  !   print*, "Error: MPI_BARRIER not succesful main"
  !   stop
  ! endif
  call MPI_FINALIZE(rc)

END PROGRAM MAIN
