module Q2BM
  USE RANDOM
  USE INTEGRATION
  use NORMAL
  use io
  IMPLICIT NONE

  public :: q2b

  PRIVATE

contains

    subroutine q2b()
        !
      INTEGER, PARAMETER :: n=4
      !REAL(8), PARAMETER ::
      REAL(8) :: miu(n),sigma(n),intx2(n),intx3(n),intx4(n),f,int_r(n)
      INTEGER :: i,j,k,h,sim(n),skip,ns,initial,final,l,m
      REAL(8) :: x1(n),x2(n),x3(n),x4(n),z1,z2,z3,z4,rho,drho,dx1

      !Suppose we want to integrate normal ϕ(x) between -∞ and a ∈ (μ-3σ,μ+3σ)
      !Make grid of values for a
      miu=(/0.0d0,0.0d0,0.0d0,1.0d0/)
      sigma=(/1.0d0,1.0d0,1.0d0,1.0d0/)
      sim=(/20,25,50,100/)
      skip=100

      ! 1) get random [0,1] number for x1, transform to [-∞, ∞] with ρ(y)=ln(y/(1-y))
      ! for that value of x1, loop x2
      ! 2) x2 loop, set x1=ub  as upper bound [-∞, ub], transform from random [0,1]
      ! use ρ(y)=up+(y-1)/y dρ/dy=y^(-2) wih 𝞥(ρ(y))dρ/dy with μ and σ
      ! 3) do same for x3 and x4
      ! 4) repeat for different x1
      ! And now with a Halton Sequence instead of a uniform

      simloop: do i=1,size(sim)
          ns=sim(i)
          initial=skip
          final=ns+skip
          intx2(i)=0.0d0
          intx3(i)=0.0d0
          intx4(i)=0.0d0

          x1loop: do j=initial,final

              CALL HALTON(j,z1)
              x1(j)=log(z1/(1.0d0-z1))
              dx1=1.0d0/(z1*(1.0d0-z1))

              x2loop: DO k=initial,final

                  CALL HALTON(k,z2)
                  rho=x1(j)+(z2-1.0d0)/z2
                  drho=z2*z2
                  f=FUN_PDF_NORMAL(rho,miu(2),sigma(2))
                  intx2(i)=intx2(i)+(f/drho)!*dx1
                  ! print*, "x1=",x1(j)," x2=",rho, " sim=",sim(i)
              END DO x2loop

              x3loop: DO l=initial,final

                  CALL HALTON(l,z3)
                  rho=x1(j)+(z3-1.0d0)/z3
                  drho=z3*z3
                  f=FUN_PDF_NORMAL(rho,miu(3),sigma(3))
                  intx3(i)=intx3(i)+(f/drho)!*dx1
                  ! print*, "x1=",x1(j)," x2=",rho, " sim=",sim(i)
              END DO x3loop

              x4loop: DO m=initial,final

                  CALL HALTON(m,z4)
                  rho=x1(j)+(z4-1.0d0)/z4
                  drho=z4*z4
                  f=FUN_PDF_NORMAL(rho,miu(4),sigma(4))
                  intx4(i)=intx4(i)+(f/drho)!*dx1
                  ! print*, "x1=",x1(j)," x2=",rho, " sim=",sim(i)
              END DO x4loop

          enddo x1loop
          intx2(i)=intx2(i)/(DBLE(ns)*DBLE(ns))
          intx3(i)=intx3(i)/(DBLE(ns)*DBLE(ns))
          intx4(i)=intx4(i)/(DBLE(ns)*DBLE(ns))
          ! intx2(i)=intx2(i)/(DBLE(ns)*DBLE(ns))
          ! intx2(i)=1.0d0-intx2(i)
      enddo simloop

      do i=1,size(sim)
          int_r(i)=(1.0d0-intx2(i))*(1.0d0-intx3(i))*(1.0d0-intx4(i))
      enddo


      print*, "MC [20,25,50,100] P(x2<x1)= ",(intx2(i),i=1,size(sim))
      print*, "MC [20,25,50,100] P(x3<x1)= ",(intx3(i),i=1,size(sim))
      print*, "MC [20,25,50,100] P(x4<x1)= ",(intx4(i),i=1,size(sim))
      print*, "Monte Carlo Integral result is ", (int_r(i),i=1,size(sim))

  end subroutine q2b

  !
  !         x3loop: do
  ! int5=int5/DBLE(n1)
  !
  ! ! And now with a Halton Sequence instead of a uniform
  ! int6=0.0d0
  ! DO j=1000,n2+1000
  !     CALL HALTON(j,m)
  !     m=m*3.0d0-1.0d0
  !     int6=int6+3.0d0*m*m
  ! END DO
  ! int6=int6/DBLE(n2)

  ! Problem 2
  ! Using Guass-Legendre since it behaved better

  !
  ! !Use Gauss Legendre
  ! ! !Simple change of variables integrating in the interval [μ-4σ,a]
  ! CALL GAUSS_LEGENDRE(x1,w1,-1.0d0,1.0d0)
  ! x1loop: do i=1,size(x1)
  !     u=x1(i)
  !    ! Loop for integral P(x2<x1)
  !    CALL GAUSS_LEGENDRE(x2,w2,-1.0d0,1.0d0)
  !
  !     l=miu(2)-4.0d0*sigma(2)
  !     m=0.5d0*(u-l)
  !     c=0.5d0*(u+l)
  !     w2=w2*m
  !     intx2=0.0d0
  !     x2loop: do j=1,size(x2)
  !         t=c+m*x2(j)
  !         f=FUN_PDF_NORMAL(t,miu(2),sigma(2))
  !         intx2=intx2+(f*w2(j))
  !     END DO x2loop
  !
  !
  !     !Loop for integral P(x3<x1)
  !     CALL GAUSS_LEGENDRE(x3,w3,-1.0d0,1.0d0)
  !     ! u2=x1(i)
  !     l2=miu(3)-4.0d0*sigma(3)
  !     m2=0.5d0*(u-l2)
  !     c2=0.5d0*(u+l2)
  !     w3=w3*m2
  !     intx3=0.0d0
  !     x3loop: do h=1,n
  !         t2=c2+m2*x3(h)
  !         f2=FUN_PDF_NORMAL(t2,miu(3),sigma(3))
  !         intx3=intx3+(f2*w3(h))
  !     END DO x3loop
  !
  !     !Loop for integral P(x4<x1)
  !     CALL GAUSS_LEGENDRE(x4,w4,-1.0d0,1.0d0)
  !     ! u3=x1(i)
  !     l3=miu(4)-4.0d0*sigma(4)
  !     m3=0.5d0*(u-l3)
  !     c3=0.5d0*(u+l3)
  !     w4=w4*m3
  !     intx4=0.0d0
  !     x4loop: do k=1,n
  !         t3=c3+m3*x4(k)
  !         f3=FUN_PDF_NORMAL(t3,miu(4),sigma(4))
  !         intx4=intx4+(f3*w4(k))
  !     END DO x4loop
  ! enddo x1loop
  !
  ! quad_int=(1.0d0-intx2)*(1.0d0-intx3)*(1.0d0-intx4)
  !
  ! print*, "Gauss-Legendre Quadrature integrals are", intx2,intx3,intx4
  ! print*, "Gauss-Legendre Quadrature integral is ", quad_int
  !


  ! print*,"Starting loop to make grid"
  ! do i=1,nm
  !     do j=1,ns
  !         up=miu(i)+3.0d0*sigma(j)
  !         lb=miu(i)-3.0d0*sigma(j)
  !         delta=(up-lb)/DBLE(na-1)
  !         do k=1,na
  !           !print*, "iteration i=",i," j=",j," k=",k
  !           !print*, "miu=",miu(i),"sigma=",sigma(j),"delta=",delta
  !           a(k,i,j)=lb+DBLE(k-1)*delta
  !           !print*, "a=",a(k,i,j)
  !         enddo
  !     enddo
  ! enddo
  ! print*,"grid done"
  ! do i=1,nm
  !     do j=1,ns
  !         print*, "a=", a(:,i,j)
  !     enddo
  ! enddo

  !i=1; j=1; k=3
  ! miuloop: do i=1,size(miu)
  !     sigmaloop: do j=1,size(sigma)
  !         aloop: do k=1,size(a(:,1,1))
  !             print*,"Starting loop for:"
  !             print*, "a=",a(k,i,j),"miu=",miu(i),"sigma=",sigma(j)
  !             !a) First Gauss Legendre ignoring that my code already does the change of variables
  !             !Simple change of variables integrating in the interval [μ-4σ,a]
  !             CALL GAUSS_LEGENDRE(quadx,quadw,-1.0d0,1.0d0)
  !
  !             u=a(k,i,j)
  !             l=miu(i)-4.0d0*sigma(j)
  !             m=0.5d0*(u-l)
  !             c=0.5d0*(u+l)
  !             quadw=quadw*m
  !             res(k,i,j,1)=0.0d0
  !             DO h=1,n
  !                 t=c+m*quadx(h)
  !                 f=FUN_PDF_NORMAL(t,miu(i),sigma(j))
  !                 res(k,i,j,1)=res(k,i,j,1)+(f*quadw(h))
  !             END DO
  !
  !             CALL GAUSS_CHEBYSHEV(quadx,quadw)
  !             u=a(k,i,j)
  !             res(k,i,j,2)=0.0d0
  !             DO h=1,n
  !                 drho= (1.0d0+quadx(h))
  !                 W=SQRT(1.0d0-quadx(h)*quadx(h))
  !                 v=(quadw(h)*W)/drho
  !                 t=LOG(1.0d0 + quadx(h))-LOG(2.0d0)+u
  !                 f=FUN_PDF_NORMAL(t,miu(i),sigma(j))
  !                 g=f/W
  !                 res(k,i,j,2)= res(k,i,j,2) + g*v
  !             END DO
  !
  !             !print*,"Gauss-Chebyshev integral result=",res(k,i,j,2)
  !
  !             ! GAUSS HERMITE ρ(y)=-exp(-y)+a ∂ρ/∂y=exp(-y)
  !
  !             CALL GAUSS_HERMITE(quadx,quadw)
  !             u=a(k,i,j)
  !             res(k,i,j,3)=0.0d0
  !             DO h=1,n
  !                 drho= EXP(-quadx(h))
  !                 W=EXP(-(quadx(h)*quadx(h)))
  !                 v=quadw(h)!/W
  !                 t=-EXP(- quadx(h))+u
  !                 f=FUN_PDF_NORMAL(t,miu(i),sigma(j))
  !                 g=W*f*drho
  !                 res(k,i,j,3)= res(k,i,j,3) + g*v
  !             END DO
  !
  !             !print*,"Gauss-Hermite integral result=",res(k,i,j,3)
  !         end do aloop
  !     enddo sigmaloop
  ! enddo miuloop

  !printing results to file
  ! call write_res("results.txt",res,miu,sigma,a)

  ! Problem 2
  ! Using Guass-Legendre since it behaved better


  ! ! Now Monte Carlo
  ! ! First with 10 points
  ! int3=0.0d0
  ! DO j=1,n1
  !     m=Sample_Uniform(-1.0d0,2.0d0)
  !     int3=int3+3.0d0*m*m
  ! END DO
  ! int3=int3/DBLE(n1)
  !
  ! ! Now with 50 points
  ! int4=0.0d0
  ! DO j=1,n2
  !     m=Sample_Uniform(-1.0d0,2.0d0)
  !     int4=int4+3.0d0*m*m
  ! END DO
  ! int4=int4/DBLE(n2)
  !
  ! ! And now with a Halton Sequence instead of a uniform
  ! int5=0.0d0
  ! DO j=1,n1
  !     CALL HALTON(j,m)
  !     m=m*3.0d0-1.0d0
  !     int5=int5+3.0d0*m*m
  ! END DO
  ! int5=int5/DBLE(n1)
  !
  ! ! And now with a Halton Sequence instead of a uniform
  ! int6=0.0d0
  ! DO j=1,n2
  !     CALL HALTON(j,m)
  !     m=m*3.0d0-1.0d0
  !     int6=int6+3.0d0*m*m
  ! END DO
  ! int6=int6/DBLE(n2)
  !
  ! WRITE(*,'(7F10.5)') 2.33d0,int1,int2,int3,int4,int5,int6


end module Q2BM
