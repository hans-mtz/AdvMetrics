likelihood.mod : \
  likelihood.f90

likelihood.o : \
  likelihood.f90 likelihood_evaluation.mod globvar.mod

