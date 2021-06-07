PROGRAM MAIN
  USE GLOBVAR
  USE INTEGRATION
  USE INTERPOL
  IMPLICIT NONE

  INTEGER :: j,h
  REAL(8) :: limits(2),arg1
  REAL(8) :: c1,c2,c3,c4,v_true(N_grid),vold2(N_grid)

  var=sigma*sigma
  
  CALL GAUSS_LEGENDRE(x_quad,w_quad,-4.0d0*sigma,4.0d0*sigma)
  
  limits(1)=0.1d0
  limits(2)=10.0d0
  CALL PLINEAR_INIT(inter_value,N_grid,limits,1)


  c1=LOG(1.0d0-(alpha*beta))/(1.0d0-beta)
  c2=(mu+(alpha*LOG(alpha*beta)))/(1.0d0-alpha)
  c3=1.0d0/(1.0d0-beta)
  c4=1.0d0/(1.0d0-(alpha*beta))
  
  DO j=1,N_grid
     y=inter_value%v(1)%x(j)
     v_old(j)=c1 + c2*(c3-c4) + c4*LOG(y)
     v_true(j)=v_old(j)
  END DO
  DO j=1,N_grid
     CALL T(j)
  END DO
  DO j=1,N_grid
      WRITE(*,*) j,v_true(j),v_new(j)
  END DO

  
  
  DO j=1,N_grid
     v_old(j)=DBLE(j)
  END DO
  DO 
     DO j=1,N_grid
        CALL T(j)
     END DO
     ! POlicy function iteration
     v_old_cons=v_new
     DO h=1,20
         DO j=1,N_grid
              y=inter_value%v(1)%x(j)
             v_new(j)=OBJECTIVE_CONS(c_new(j))
         END DO
         v_old_cons=v_new
     END DO
     
     arg1=SUM((v_new-v_old)**2.0d0)
     IF (arg1<1.0d-6) EXIT
     v_old=v_new
  END DO
  WRITE(*,*)
  WRITE(*,*)
  DO j=1,N_grid
     WRITE(*,*) j,v_true(j),v_new(j),(1.0d0-(alpha*beta))*(inter_value%v(1)%x(j)),c_new(j)
  END DO

END PROGRAM MAIN
  
