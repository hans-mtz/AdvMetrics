PROGRAM MAIN
  USE NRUTIL
  USE RANDOM
  USE INTERPOL
  IMPLICIT NONE
  
  TYPE(PLINEARSTRUC) :: pg,pg1
  TYPE(CHEBSTRUC) :: cg,cg1
  TYPE(MLINEARSTRUC) :: mg,mg1
  REAL(8) :: limits(3,2),point_ml(3),point_pl(3),point_ch(3),cc
  INTEGER :: j,n,k
  REAL(8), POINTER :: fx_ml(:),fx_pl(:),fx_ch(:)
  
  limits(:,1)=0.0d0
  limits(:,2)=4.0d0
  
  CALL MLINEAR_INIT(mg,(/4,5,3/),limits)
  CALL PLINEAR_INIT(pg,(/4,5,3/),limits)
  CALL CHEB_INIT(cg,(/4,5,3/),limits,.FALSE.)
  
  ALLOCATE(fx_ml(mg%lt),fx_ch(cg%pn),fx_pl(pg%lt))
  
  DO n=1,mg%lt
     DO j=1,mg%nv
        point_ml(j)=mg%v(j)%x(mg%i(n)%i(j))
     END DO
     fx_ml(n)=FUN(point_ml)
     DO j=1,pg%nv
        point_pl(j)=pg%v(j)%x(pg%i(n)%i(j))
     END DO
     fx_pl(n)=FUN(point_pl)
     DO j=1,cg%nv
        point_ch(j)=cg%v(j)%x(cg%si(n)%i(j))
     END DO
     fx_ch(n)=FUN(point_ch)
  END DO
  CALL CHEB_COEF(cg,fx_ch)
  
  OPEN(1,file="out.out")
  WRITE(1,'(A)') "Let's first check they give the right answer at the grid points"
  DO n=1,mg%lt
     DO j=1,mg%nv
        point_ml(j)=mg%v(j)%x(mg%i(n)%i(j))
        point_pl(j)=pg%v(j)%x(pg%i(n)%i(j))
        point_ch(j)=cg%v(j)%x(cg%si(n)%i(j))
     END DO
     !cc=MLINEAR_INTER(mg,point_ml,fx_ml)
     WRITE(1,'(40F16.8)') point_ml,PLINEAR_INTER(pg,point_pl,fx_pl),MLINEAR_INTER(mg,point_ml,fx_ml),fx_ml(n),&
          point_ch,CHEB_INTER(cg,point_ch),fx_ch(n)
  END DO
  WRITE(1,'(A)') 
  WRITE(1,'(A)') "Now let's check how they do off the grid points"
  DO n=1,10
     DO j=1,mg%nv
        point_pl(j)=Sample_Uniform(limits(j,1),limits(j,2))
     END DO
     WRITE(1,'(40F16.8)') point_pl,PLINEAR_INTER(pg,point_pl,fx_pl),MLINEAR_INTER(mg,point_pl,fx_ml),&
          CHEB_INTER(cg,point_pl),FUN(point_pl)
  END DO
  
CONTAINS
  
  REAL(8) FUNCTION FUN(x)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x(:)
    INTEGER :: i
    FUN=0.0d0
    FUN = (x(1)**3.0d0 + x(2)**2.5d0 + x(3)**3.5d0)**0.33d0 
  END FUNCTION FUN
  
  
END PROGRAM MAIN

