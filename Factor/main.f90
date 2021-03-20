PROGRAM MAIN
  USE GLOBVAR
  USE INIT
  USE LIKELIHOOD
  USE LIKELIHOOD_EVALUATION, ONLY : LIKE_EVAL
  USE MINIMIZATION
  USE SIMPLEX
!  USE MPI
  
  IMPLICIT NONE
  
  REAL(8), ALLOCATABLE :: theta(:)
  INTEGER :: j,k,ierr,n,lo,hi
  REAL(8) :: val

  myid=0
  numprocs=1
  
!  CALL MPI_INIT( ierr ) ! Always call it at the begining of the program
!  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) ! Always call it at the begining of the program, myid is going to be the id for the processor
!  CALL MPI_COMM_SIZE( MPI_COMM_WORLD,numprocs,ierr) ! Always call it at the begining of the program, numprocs is the number of processors that you gave when you ran it     

  ! INITIALIZE THE MODEL (i.e. READ IN THE DATA)
  CALL INITIALIZE(myid)
  
  ! Impose normalizations
  am=1.0d0

  ALLOCATE(theta(10000))
  theta=0.0d0
  n=TRANSFORM(theta)
  DEALLOCATE(theta)
  ALLOCATE(theta(n))
  WRITE(*,*) n, SIZE(theta)
  IF (myid==0) THEN
     OPEN(1,file="point.out") ! File with initial values for the parameters
     DO j=1,n
        READ(1,fmt=*) theta(j)
     END DO
     CLOSE(1)
  END IF
!  CALL MPI_BCAST(theta,n,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  n=TRANSFORM(theta)

  IF (myid==0) WRITE(*,*) "Calling Minimization Routine"
!  CALL NELDER_MEADE(theta,1.0d-10,LOGLIKELIHOOD,myid,5000)
  CALL BFGS(theta,1.0d-4,1.0d-14,loglikelihood,myid,.TRUE.)
  val= LOGLIKELIHOOD(theta)
  
  n=TRANSFORM(theta)
  
  ! Ok, now write them out
  IF (myid==0) THEN
     OPEN(1,file='output.out')
     WRITE(1,'(A,F42.16)') "Value",val
     WRITE(1,'(A,I10)') "Number of parameter",n 
     WRITE(1,*)
     WRITE(1,*)
     WRITE(1,'(A)') "measurement equation"
     WRITE(1,'(A,<N_x>F16.8,A,F16.8,A,F16.8)') 'bm',bm,' am',am,'  variance',vm
     WRITE(1,*)
     WRITE(1,'(A)') "Outcome d=0"
     WRITE(1,'(A,<N_x>F16.8,A,F16.8,A,F16.8)') 'b0',b0,' a0',a0,'  variance',v0
     WRITE(1,*)
     WRITE(1,'(A)') "Outcome d=1"
     WRITE(1,'(A,<N_x>F16.8,A,F16.8,A,F16.8)') 'b1',b1,' a1',a1,'  variance',v1
     WRITE(1,*)
     WRITE(1,'(A)') "Choice equation"
     WRITE(1,'(A,<N_x>F16.8,A,<N_z>F16.8,A,F16.8,A,F16.8)') 'gx',gx,' gz',gz,' loading',av,'  variance',1.0d0 
     WRITE(1,*) 
     WRITE(1,'(A,I4)') 'Factor'
     WRITE(1,'(A,<N_mix>F16.8,A,<N_mix>F16.8,A,<N_mix>F16.8)') 'Means',m_fac,'  Variances',v_fac,'  Weights',p_fac
     CLOSE(1)
  END IF
  
!  CALL MPI_FINALIZE(ierr)
  STOP
  
END PROGRAM MAIN
