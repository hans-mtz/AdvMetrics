FC = mpifort
MPI = mpifort
EX = mpiexec
FFLAGS = -O3 -mkl -m64 -lmkl_blas95_ilp64 -lmkl_blas95_lp64 \
	-lmkl_lapack95_ilp64 -lmkl_lapack95_lp64 -xHost -lmylib \
	#-qopenmp -par-num-threads=5 -parallel #-init=snan # -gen-dep=deps.txt
DFLAGS = -g3 -check all -init=snan -o0
SRC =	hello_mpi.f90

OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(@F) $<

myp: $(OBJ)
	$(FC) $(FFLAGS) -r16 -o $@ $(notdir $(OBJ))

debug:
	$(FC) $(DFLAGS) $(FFLAGS) -o $@ $(notdir $(OBJ))


clean:
	rm -f *.o *.mod myp
