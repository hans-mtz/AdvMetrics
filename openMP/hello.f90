program hello
  use omp_lib
  implicit none

    integer :: num_threads = 4
    integer :: thread_num = 0
    !$ call omp_set_num_threads(num_threads)
    print*, "Display Hello world!"
    print*, "Number of threads used = ", num_threads

    !intendations are for reference only!

    !$omp parallel
        !$omp critical
            !$ thread_num = omp_get_thread_num()
            print *, "Hello world from thread number ", thread_num
        !$omp end critical
    !$omp end parallel


end program hello
