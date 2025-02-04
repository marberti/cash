CC=gcc
FC=gfortran
EXE=cash

.PHONY: default
default:
	$(CC) -c mergeball_kernel.c
	$(FC) -c main.f90
	$(FC) -o $(EXE) mergeball_kernel.o main.o

.PHONY: clean
clean:
	rm -f *.o *.mod

