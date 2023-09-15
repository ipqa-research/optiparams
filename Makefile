#FC=gfortran
#FFLAG=-Wall -Wextra
FC=ifort -fast

deps:
	$(FC) -c src/asa057.f90
	$(FC) -c src/praxis.f
	$(FC) -c src/Pure.for
	$(FC) -c src/RKPR.for
	$(FC) -c src/SRK_PR.for
	mv *.o obj/

sist: deps
	$(FC) -o sist ./Sistemas/OptimCMRKP2011.for obj/* -qmkl

clean:
	rm obj/*
