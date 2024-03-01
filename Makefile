FC=gfortran
EXE=cash

.PHONY: default
default:
	$(FC) -o $(EXE) main.f90

