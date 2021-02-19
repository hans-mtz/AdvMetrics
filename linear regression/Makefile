FC = ifort
FFLAGS = -O3 -mkl -m64 -lmkl_blas95_ilp64 -lmkl_blas95_lp64 \
	-lmkl_lapack95_ilp64 -lmkl_lapack95_lp64 -xHost -init=snan #-gen-dep=deps.txt
DFLAGS = -g3 -check all -init=snan
prefix = /Users/hans/F_code_mkl
# prefix1 = /opt/intel/oneapi/mkl/2021.1.1/include
SRC = $(prefix)/nrutil.f90 \
	$(prefix)/integration.f90 \
	$(prefix)/random.f90 \
	$(prefix)/probability.f90 \
	$(prefix)/minimization.f90 \
	$(prefix)/anneal.f90 \
	globvar.f90 \
	init.f90 \
	likelihood_evaluation.f90 \
	likelihood.f90 \
	io.f90 \
	linear_regression.f90

OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(@F) $<

myprog: $(OBJ)
	$(FC) $(FFLAGS) -r16 -o $@ $(notdir $(OBJ))

debug:
	$(FC) $(DFLAGS) $(FFLAGS) -o $@ $(notdir $(OBJ))

clean:
	rm -f *.o *.mod *.out myprog
