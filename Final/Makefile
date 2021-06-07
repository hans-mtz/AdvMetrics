FC = mpifort
# FC = ifort
FFLAGS = -O3 -mkl -m64 -lmkl_blas95_lp64 \
	-lmkl_lapack95_lp64 -xHost -lmylib \
	-parallel #-qopenmp #-par-num-threads=5  #-init=snan # -gen-dep=deps.txt
DFLAGS = -g3 -check all -init=snan -o0
SRC =	simplex1.f90 \
	globvar.f90 \
	init.f90 \
	wmod.f90 \
	Source1.f90

OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(@F) $<

myp: $(OBJ)
	$(FC) $(FFLAGS) -r16 -o $@ $(notdir $(OBJ))

debug:
	$(FC) $(DFLAGS) $(FFLAGS) -o $@ $(notdir $(OBJ))

clean:
	rm -f *.o *.mod myp
