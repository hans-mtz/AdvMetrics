module boot
  use GLOBVAR
  use loglikelihood

  implicit none

contains
  subroutine bootstrap_sub_v(subr,dt,ns,bs_f)
    real(8), INTENT(IN) :: f,dt(:,:)
    integer, INTENT(IN) :: ns
    real(8), DIMENSION(ns,:),ALLOCATABLE, INTENT(OUT) :: bs_f

    interface
      subroutine subr(theta,a,b,fun,ir,bool)


      end subroutine subr




  end subroutine bootstrap

end module boot
