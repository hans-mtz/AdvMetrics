program main
  use mpi
  implicit none
  integer :: ierr, nprocs, procid

  call MPI_INIT(ierr)

  call MPI_COMM_RANK(MPI_COMM_WORLD, procid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  print*, "Hello from procs ", procid, " out of ", nprocs

  call MPI_FINALIZE(ierr)

end program main 
