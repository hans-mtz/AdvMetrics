MODULE LIKELIHOOD_EVALUATION
  USE GLOBVAR
  USE PROBABILITY, ONLY : CDF_NORMAL,PDF_NORMAL
  IMPLICIT NONE
  
CONTAINS
  
  REAL(8) FUNCTION LIKE_EVAL(i)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    INTEGER :: j,k,ex,h,t
    REAL(8) :: th,xm,xb,indx,cdf

    xm=DOT_PRODUCT(x(i,:),bm)
    IF (d(i)==0) THEN
       xb=DOT_PRODUCT(x(i,:),b0)
    ELSE
       xb=DOT_PRODUCT(x(i,:),b1)
    END IF
    indx=DOT_PRODUCT(x(i,:),gx)+DOT_PRODUCT(z(i,:),gz)
    LIKE_EVAL=0.0d0    
    DO h=1,N_int
       th=0.0d0
       cdf=0.0d0
       ex=0
       DO k=1,N_mix
          IF ((uni(i,h)>cdf).AND.(uni(i,h)<=cdf+p_fac(k))) THEN
             th=m_fac(k) + SQRT(v_fac(k))*nor(i,h)
             ex=1
          END IF
          IF (ex==1) EXIT
          cdf=cdf+p_fac(k)
       END DO
       LIKE_EVAL = LIKE_EVAL + F_FAC(th,xm,xb,indx,i)
    END DO
    LIKE_EVAL=LIKE_EVAL/DBLE(N_int)
  END FUNCTION LIKE_EVAL

  REAL(8) FUNCTION F_FAC(th,xm,xb,indx,i)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    REAL(8), INTENT(IN) :: th,xm,xb,indx
    INTEGER :: j
    F_FAC=PDF_NORMAL(m(i)-xm-th*am,0.0d0,vm)
    IF (d(i)==0) THEN
       F_FAC=F_FAC*PDF_NORMAL(y(i)-xb-th*a0,0.0d0,v0)
       F_FAC=F_FAC*CDF_NORMAL(-indx-th*av,0.0d0,1.0d0)
    ELSE
       F_FAC=F_FAC*PDF_NORMAL(y(i)-xb-th*a1,0.0d0,v1)
       F_FAC=F_FAC*(1.0d0-CDF_NORMAL(-indx-th*av,0.0d0,1.0d0))
    END IF
  END FUNCTION F_FAC
  
END MODULE LIKELIHOOD_EVALUATION
