FC = ifort
FFLAGS = -O3 -mkl -m64 -lmkl_blas95_ilp64 -lmkl_blas95_lp64 -lmkl_lapack95_ilp64 \
	-lmkl_lapack95_lp64 -xHost
DFLAGS = -g3 -check bounds
prefix = /Users/hans/F_code_mkl
# prefix1 = /opt/intel/oneapi/mkl/2021.1.1/include
SRC = $(prefix)/nrutil.f90 \
	$(prefix)/integration.f90 \
	$(prefix)/random.f90 \
	$(prefix)/probability.f90 \
	main.f90
	# $(prefix)/sobol.f90

OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(@F) $<

myprog: $(OBJ)
	$(FC) $(FFLAGS) -r16 -i8 -o $@ $(notdir $(OBJ))

test:
	$(FC) $(FFLAGS) -o $@ $(SRC)

debug:
	$(FC) $(DFLAGS) $(FFLAGS) -o $@ $(notdir $(OBJ))

clean:
	rm -f *.o *.mod *.out myprog
