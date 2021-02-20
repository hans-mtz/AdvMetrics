module Q2AM
  USE RANDOM
  USE INTEGRATION
  use NORMAL
  use io
  IMPLICIT NONE

  PUBLIC :: q2a

  PRIVATE

contains

  subroutine  q2a()

      INTEGER, PARAMETER :: n=12,nm=4,ns=4
      !REAL(8), PARAMETER ::
      REAL(8) :: miu(nm),sigma(ns),quad_int,dx1,dm,dm2,dm3,f1,dt1,dt2,dt3,dt4
      INTEGER :: i,j,k,h,h1,h2,h3
      REAL(8) :: l,l1,l2,l3,u,u1,u2,u3,m,m1,m2,m3,c,c1,c2,c3,t,t1,t2,t3,t4,f,f2,f3
      REAL(8) :: x1(n),x2(n),x3(n),x4(n),w1(n),w2(n),w3(n),w4(n),intx2,intx3,intx4

      miu=(/0.0d0,0.0d0,0.0d0,1.0d0/)
      sigma=(/1.0d0,1.0d0,1.0d0,1.0d0/)
      ! sim=(/20,25,50,100/)
      ! skip=1000

      ! Problem 2
      ! Using Guass-Legendre since it behaved better
      ! !Use Gauss Legendre
      ! 
      intx2=0.0d0
      intx3=0.0d0
      intx4=0.0d0

      ! l1=miu(1) - 4.0d0*sigma(1)
      ! u1=miu(1) + 4.0d0*sigma(1)
      CALL GAUSS_LEGENDRE(x1,w1,-1.0d0,1.0d0)
      x1loop: do i=1,size(x1)
          ! l1=miu(1)-4.0d0*sigma(1)
          ! u1=miu(1)+4.0d0*sigma(1)
          ! m1=0.5d0*(u1-l1)
          ! c1=0.5d0*(u1+l1)
          ! u=c1+m1*x1(i)
          ! w1=w1/m1
          ! u1=x1(i)
          t1=x1(i)/(x1(i)*x1(i)-1.0d0)
          ! t1=log((x1(i)+ 1.0d0)/(x1(i)- 1.0d0))
          ! dt1=-2.0d0/(x1(i)*x1(i)- 1.0d0)
          ! dt1=-(x1(i)*x1(i)+1.0d0)/((x1(i)*x1(i)-1.0d0)**2.0d0)
          ! f1=FUN_PDF_NORMAL(x1(i),miu(1),sigma(1))*w1(i)

         ! Loop for integral P(x2<x1)
         ! u=miu(2) + 4.0d0*sigma(2)
         ! l=miu(2) - 4.0d0*sigma(2)
         CALL GAUSS_LEGENDRE(x2,w2,-1.0d0,1.0d0)

          ! m=0.5d0*(u-l)
          ! ! c=0.5d0*(u+l)
          ! dm=(x1(i)-l)/(u-l)
          w2=w2*w1(i)
          ! intx2=0.0d0
          x2loop: do j=1,size(x2)
              ! t=c+m*x2(j)
              t2=t1+((x2(j)- 1.0d0)/(x2(j)+ 1.0d0))
              dt2=2.0d0/((x2(j)+ 1.0d0)**2.0d0)
              ! dx1=(x2(j)-l)/(u-l)
              f=FUN_PDF_NORMAL(t2,miu(2),sigma(2))
              intx2=intx2+f*w2(j)/dt2!*dt1
          END DO x2loop


          !Loop for integral P(x3<x1)
          ! u2=miu(3) + 4.0d0*sigma(3)
          ! l2=miu(3) - 4.0d0*sigma(3)
          CALL GAUSS_LEGENDRE(x3,w3,-1.0d0,1.0d0)

          ! m2=0.5d0*(u2-l2)
          ! c2=0.5d0*(u+l2)
          ! dm2=(x1(i) - l2)/(u2 - l2)
          w3=w3*w1(i)
          ! intx3=0.0d0
          x3loop: do h=1,n
              ! t2=c2+m2*x3(h)
              t3=t1+((x3(h)- 1.0d0)/(x3(h)+ 1.0d0))
              dt3=2.0d0/((x3(h)+ 1.0d0)**2.0d0)
              ! dx1=(x3(h)-l2)/(u2-l2)
              f2=FUN_PDF_NORMAL(t3,miu(3),sigma(3))
              intx3=intx3+f2*w3(h)/dt3!*dx1)
          END DO x3loop

          !Loop for integral P(x4<x1)
          ! u3=miu(4) + 4.0d0*sigma(4)
          ! l3=miu(4) - 4.0d0*sigma(4)
          CALL GAUSS_LEGENDRE(x4,w4,-1.0d0,1.0d0)

          ! m3=0.5d0*(u-l3)
          ! c3=0.5d0*(u+l3)
          ! dm3=(x1(i)-l3)/(u3-l3)
          w4=w4*w1(i)
          ! intx4=0.0d0
          x4loop: do k=1,n
              ! t3=c3+m3*x4(k)
              t4=t1+((x4(k)- 1.0d0)/(x4(k)+ 1.0d0))
              dt4=2.0d0/((x4(k)+ 1.0d0)**2.0d0)
              ! dx1=(x4(k)-l3)/(u3-l3)
              f3=FUN_PDF_NORMAL(t4,miu(4),sigma(4))
              intx4=intx4+f3*w4(k)/dt4!*dx1)
          END DO x4loop
      enddo x1loop

      quad_int=(1.0d0-intx2)*(1.0d0-intx3)*(1.0d0-intx4)

      print*, "Gauss-Legendre Quadrature integrals are", intx2,intx3,intx4
      print*, "Gauss-Legendre Quadrature integral is ", quad_int

end subroutine  q2a

    ! !Suppose we want to integrate normal ϕ(x) between -∞ and a ∈ (μ-3σ,μ+3σ)
    ! !Make grid of values for a
    ! miu=(/0.0d0,0.0d0,0.0d0,1.0d0/)
    ! sigma=(/1.0d0,1.0d0,1.0d0,1.0d0/)
    ! sim=(/20,25,50,100/)
    ! skip=1000
    !
    ! ! 1) get random [0,1] number for x1, transform to [-∞, ∞] with ρ(y)=ln(y/(1-y))
    ! ! for that value of x1, loop x2
    ! ! 2) x2 loop, set x1=ub  as upper bound [-∞, ub], transform from random [0,1]
    ! ! use ρ(y)=up+(y-1)/y dρ/dy=y^(-2) wih 𝞥(ρ(y))dρ/dy with μ and σ
    ! ! 3) do same for x3 and x4
    ! ! 4) repeat for different x1
    ! ! And now with a Halton Sequence instead of a uniform
    !
    ! simloop: do i=1,size(sim)
    !     ns=sim(i)
    !     initial=skip
    !     final=ns+skip
    !     intx2(i)=0.0d0
    !
    !     x1loop: do j=initial,final
    !
    !         CALL HALTON(j,z1)
    !         x1(j)=log(z1/(1.0d0-z1))
    !
    !         x2loop: DO k=initial,final
    !
    !             CALL HALTON(k,z2)
    !             rho=x1(j)+(z2-1.0d0)/z2
    !             drho=1.0d0/(z2*z2)
    !             f=FUN_PDF_NORMAL(rho,miu(2),sigma(2))
    !             intx2(i)=intx2(i)+drho*f/DBLE(ns)
    !             print*, "x1=",x1(j)," x2=",rho, " sim=",sim(i)
    !         END DO x2loop
    !
    !
    !     enddo x1loop
    !     intx2(i)=1.0d0-intx2(i)
    ! enddo simloop
    !
    ! print*, "Integral result is ", (intx2(i),i=1,4)
    !
    ! !
    ! !         x3loop: do
    ! ! int5=int5/DBLE(n1)
    ! !
    ! ! ! And now with a Halton Sequence instead of a uniform
    ! ! int6=0.0d0
    ! ! DO j=1000,n2+1000
    ! !     CALL HALTON(j,m)
    ! !     m=m*3.0d0-1.0d0
    ! !     int6=int6+3.0d0*m*m
    ! ! END DO
    ! ! int6=int6/DBLE(n2)

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

END module Q2AM
