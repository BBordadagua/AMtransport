.SUFFIXES: .f90

F90=gfortran
F90FLAGS=-fdefault-real-8 -O2 -g -fcheck=all -Wall

SRCDIR=src

SOURCE= parameters.f90 mathfunctions_mod.f90 data_handling_mod.f90 coefficients_mod.f90 modeamplitude_mod.f90 rotationprofile_mod.f90 mixedmodesflux.f90
OBJECTS=$(SOURCE:.f90=.o)

all:
	cd $(SRCDIR); $(MAKE) -f ../Makefile program.exec

program.exec: $(OBJECTS)
	$(F90) $(F90FLAGS) $(OBJECTS) -o ../$@

.f90.o:
	$(F90) $(F90FLAGS) -c $<

clean:
	rm -f program.exec output/*.txt src/*.o src/*.mod src/*~

