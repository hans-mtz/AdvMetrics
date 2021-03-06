MODULE LIKELIHOOD
  USE GLOBVAR
  USE LIKELIHOOD_EVALUATION, ONLY : LIKE_EVAL
  USE PROBABILITY, ONLY : LN_PDF_GAMMA
!  USE MPI
  IMPLICIT NONE
  
CONTAINS
  
  REAL(8) FUNCTION LOGLIKELIHOOD(theta)
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: myout,out,aux,pen,den,weight
    INTEGER :: i,j,k,ierr,n,t
    
    n=TRANSFORM(theta)

    ! OK no let's go to the likelihood
    myout=0.0d0
    DO i = myid+1, N_ind, numprocs ! The loop, notice it starts at myid+1 and strides numprocs, that way I split the work amongst numprocs processors
       ! A bunch of crap that needs to be done only once for each individual (so I do not repeat the operation)
       aux = LIKE_EVAL(i)
       IF (aux<1.0d-300) THEN ! In case I get a very small number (1.0d-300 or less) I cannot take the log so I assign -log(1.0d-300)
          myout = myout + 690.77552789821368151024216786026954650879d0
       ELSE
          myout = myout - LOG(aux)
       END IF
    END DO
!    CALL MPI_REDUCE(myout,out,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!    CALL MPI_BCAST(out,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    
    ! Penalizing the likelihood so that the mixture parameters "behave"
!    pen=0.0d0
    !    p_dirich=2.0d0
!    DO k=1,N_mix
!       pen = pen+LN_PDF_GAMMA(1.0d0/v_fac(k),2.0d0,1.0d0)
       !          pen = pen+LN_PDF_NORMAL(m_fac(j,k),0.0d0,4.0d0)
!    END DO
    !       pen = pen+LN_PDF_DIRICHLET(p_fac(j,:),p_dirich)
    
    !    LOGLIKELIHOOD = out -pen
    LOGLIKELIHOOD=myout
    RETURN

  END FUNCTION LOGLIKELIHOOD

  INTEGER FUNCTION TRANSFORM(theta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    INTEGER :: n,j,k,h,lo,hi
    n=0
    DO k=1,N_mix-1
       n=n+1
       m_fac(k)=theta(n)
       n=n+1
       v_fac(k)=1.4d0**theta(n)
       n=n+1
       p_fac(k)=1.4d0**theta(n)
    END DO
    n=n+1
    v_fac(N_mix)=1.4d0**theta(n)
    p_fac(N_mix)=1.0d0
    p_fac=p_fac/SUM(p_fac)
    m_fac(N_mix)=-DOT_PRODUCT(p_fac(1:N_mix-1),m_fac(1:N_mix-1))/p_fac(N_mix) ! Because the mean is normalized to zero
    DO k=1,N_x
       n=n+1
       bm(k)=theta(n)
    END DO
    n=n+1
    vm=1.4d0**theta(n)
    am=1.0d0
    DO k=1,N_x
       n=n+1
       b0(k)=theta(n)
    END DO
    n=n+1
    v0=1.4d0**theta(n)
    n=n+1
    a0=theta(n)
    DO k=1,N_x
       n=n+1
       b1(k)=theta(n)
    END DO
    n=n+1
    v1=1.4d0**theta(n)
    n=n+1
    a1=theta(n)
    DO k=1,N_x
       n=n+1
       gx(k)=theta(n)
    END DO
    DO k=1,N_z
       n=n+1
       gz(k)=theta(n)
    END DO
    n=n+1
    av=theta(n)

    TRANSFORM=n
    RETURN
  END FUNCTION TRANSFORM
  
END MODULE LIKELIHOOD
