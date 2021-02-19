PROGRAM MAIN
  USE RANDOM
  USE INTEGRATION
  use NORMAL
  use io
  use Q2AM
  use Q2BM
  IMPLICIT NONE

  INTEGER, PARAMETER :: n=12,na=5,nm=2,ns=2,nmet=3
  !REAL(8), PARAMETER ::
  REAL(8) :: miu(nm),sigma(ns),y(n),delta,a(na,nm,ns),res(na,nm,ns,nmet)
  INTEGER :: i,j,k,h,h1,h2,h3
  REAL(8) :: quadx(n),quadw(n),jacob,up,lb,l,u,m,c,t,f,drho,v,W,g


  !Suppose we want to integrate normal ϕ(x) between -∞ and a ∈ (μ-3σ,μ+3σ)
  !Make grid of values for a
  miu=(/-1.0d0,1.0d0/)
  sigma=(/0.5d0,5.0d0/)
  print*,"Starting loop to make grid"
  do i=1,nm
      do j=1,ns
          up=miu(i)+3.0d0*sigma(j)
          lb=miu(i)-3.0d0*sigma(j)
          delta=(up-lb)/DBLE(na-1)
          do k=1,na
            !print*, "iteration i=",i," j=",j," k=",k
            !print*, "miu=",miu(i),"sigma=",sigma(j),"delta=",delta
            a(k,i,j)=lb+DBLE(k-1)*delta
            !print*, "a=",a(k,i,j)
          enddo
      enddo
  enddo
  print*,"grid done"
  ! do i=1,nm
  !     do j=1,ns
  !         print*, "a=", a(:,i,j)
  !     enddo
  ! enddo

  !i=1; j=1; k=3
  miuloop: do i=1,size(miu)
      sigmaloop: do j=1,size(sigma)
          aloop: do k=1,size(a(:,1,1))
              ! print*,"Starting loop for:"
              ! print*, "a=",a(k,i,j),"miu=",miu(i),"sigma=",sigma(j)
              !a) First Gauss Legendre ignoring that my code already does the change of variables
              !Simple change of variables integrating in the interval [μ-4σ,a]
              CALL GAUSS_LEGENDRE(quadx,quadw,-1.0d0,1.0d0)

              u=a(k,i,j)
              l=miu(i)-4.0d0*sigma(j)
              m=0.5d0*(u-l)
              c=0.5d0*(u+l)
              quadw=quadw*m
              res(k,i,j,1)=0.0d0
              DO h=1,n
                  t=c+m*quadx(h)
                  f=FUN_PDF_NORMAL(t,miu(i),sigma(j))
                  res(k,i,j,1)=res(k,i,j,1)+(f*quadw(h))
              END DO

              CALL GAUSS_CHEBYSHEV(quadx,quadw)
              u=a(k,i,j)
              res(k,i,j,2)=0.0d0
              DO h1=1,n
                  drho= 1.0d0/(1.0d0+quadx(h1))
                  ! drho= (2.0d0*u)/((1.0d0+quadx(h))**2.0d0)
                  W=SQRT(1.0d0-quadx(h1)*quadx(h1))
                  v=quadw(h1)*drho*W/W
                  t=LOG(1.0d0 + quadx(h1))-LOG(2.0d0)+u
                  ! t=(2.0d0*u*quadx(h))/(1.0d0+quadx(h))
                  f=FUN_PDF_NORMAL(t,miu(i),sigma(j))
                  g=f*W
                  res(k,i,j,2)= res(k,i,j,2) + g*v
              END DO

              !print*,"Gauss-Chebyshev integral result=",res(k,i,j,2)

              ! GAUSS HERMITE ρ(y)=-exp(-y)+a ∂ρ/∂y=exp(-y)

              CALL GAUSS_HERMITE(quadx,quadw)
              u=a(k,i,j)
              res(k,i,j,3)=0.0d0
              DO h2=1,n
                  drho= EXP(quadx(h2)*quadx(h2)-quadx(h2))
                  W=EXP(quadx(h2)*quadx(h2))
                  v=quadw(h2)*drho/W
                  t=-EXP(- quadx(h2))+u
                  f=FUN_PDF_NORMAL(t,miu(i),sigma(j))
                  g=W*f
                  res(k,i,j,3)= res(k,i,j,3) + g*v
              END DO

              !print*,"Gauss-Hermite integral result=",res(k,i,j,3)
          end do aloop
      enddo sigmaloop
  enddo miuloop

  !printing results to file
  call write_res("results.txt",res,miu,sigma,a)

  call q2a()

  call q2b()

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

END PROGRAM MAIN
