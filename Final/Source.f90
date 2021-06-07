PROGRAM MAIN
  USE GLOBVAR
  use INIT
  ! USE INTEGRATION
  USE INTERPOL
  USE IO
  use mod_csv
  use wmod1

  IMPLICIT NONE

  INTEGER :: j,h, rc,i
  REAL(8) :: limits(2),arg1, a
  ! REAL(8) :: c1,c2,c3,c4,v_true(N_grid),vold2(N_grid)
  real(8) :: r_t2(N_sim,9), r_t1(N_sim,4), theta(5)
  real*8 :: t1, t2, t3,t4, rho_test(15), test_data(5000,7)!, rho_0(15)

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

  call READDATA(assets,test_data)
  !
  ! print*, "printing rho_0 = ",rho_0
  ! call get_rho(rho_h,test_data)
  ! print*, "rho_0 from fun =  ", rho_h
  call initiate(l_v,sigma2_l)

  ! print*, "printing lambda vector = ", l_v(1:100:10)

  call read_file("point.out",theta)

  call get_rho(rho_0,test_data(1:N_sim,:))

  ! print*, "rho_hat = ", rho_test
  call CPU_TIME(t1)

  a=obj_ind(theta)
  ! call q5(theta)


  call CPU_TIME(t2)
  ! t2= MPI_WTIME()
  t3=t2-t1


  print*, "theta = ", theta!, "from proc ",myid
  print*, "obj fun value = ", valuetemp!, "from proc ",myid
  ! print*, "obj fun value a = ", a
  print*, t3, " seconds for q5"!, "from proc ",myid
  print "(i3,2f16.6)",(i, rho_0(i),rho_h(i), i=1,15)
  !"(15(3F4.6,/))"! flush(1)
  call csvwrite("fake_df",f_data)


END PROGRAM MAIN
