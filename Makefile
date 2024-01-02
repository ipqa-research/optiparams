#FC=gfortran
#FFLAG=-Wall -Wextra
FC=ifort -g -extend-source -fast # -fopenmp

deps:
	mkdir -p obj
	$(FC) -c src/asa057.f90
	$(FC) -c src/sa.f90
	$(FC) -c src/praxis.f
	$(FC) -c src/Pure.for
	$(FC) -c src/RKPR.for
	$(FC) -c src/SRK_PR.for
	mv *.o obj/

sist: deps
	$(FC) -o sist ./Sistemas/OptimCMRKP2011.for obj/* -qmkl

series: deps
	$(FC) -o series ./Series/OptimQMRseries2017.f90 obj/* -qmkl

clean:
	rm obj/*
