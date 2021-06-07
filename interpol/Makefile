IFORT = ifort
MGFORT = mygfor 
GFORT = gfortran

SRCS = 	../../F_code/nrutil.f90 \
	../../F_code/interpol.f90 \
	../../F_code/random.f90 \
	main.f90

serial:
	$(IFORT) $(SRCS) -ip -ipo -O3 -unroll -o Multiroy.exe
debug:
	$(IFORT) $(SRCS) -g -o debug.exe 
debugmg:
	$(MGFORT) $(SRCS) -g -fcheck=bounds -o debug.exe
debugg:
	$(GFORT) $(SRCS) -g -fcheck=bounds -o debug.exe

clean:
	rm -rf *.mod *.o *.exe *.ipo *.il *~ \#* PI*
 


