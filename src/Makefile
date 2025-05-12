FC = ifort

all: findsqs

findsqs: findsqs.f90 functions.o structure.o symmetry.o supercell.o cluster.o configurations.o outfiles.o
	$(FC) -o findsqs findsqs.f90 *.o
	cp findsqs ../bin/

supercell.o: supercell.f90 functions.o
	$(FC) -c supercell.f90

functions.o: functions.f90
	$(FC) -c functions.f90

structure.o: structure.f90 functions.o
	$(FC) -c structure.f90

symmetry.o: symmetry.f90 functions.o
	$(FC) -c symmetry.f90

cluster.o: cluster.f90 functions.o
	$(FC) -c cluster.f90 

configurations.o: configurations.f90 functions.o
	$(FC) -c configurations.f90

outfiles.o: outfiles.f90 functions.o structure.o
	$(FC) -c outfiles.f90

clean:
	rm *.mod *.o findsqs
