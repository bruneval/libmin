FC=gfortran # -g -fbacktrace -fcheck=all
CC=gcc     # -g
CPP=-cpp
FORTRANLIB=-lgfortran

.PHONY: clean tarball

all: test_c test_fortran libmin.a

test_fortran: test_fortran.o libmin.a
	$(FC) -o test_fortran test_fortran.o -L. -lmin 

test_c: test_c.o libmin.a libmin.h
	$(CC) -o test_c test_c.o  -L. -lmin -lm  $(FORTRANLIB)

test_fortran.o: test_fortran.f90 libmin.f03
	$(FC) -c $(CPP) test_fortran.f90

test_c.o: test_c.c libmin.h
	$(CC) -c test_c.c

libmin.a: libmin.o lbfgs.o mcsrch.o libmin.h
	$(AR) -crs libmin.a libmin.o mcsrch.o lbfgs.o

lbfgs.o: lbfgs.f90 libmin.h
	$(FC) -c $(CPP) lbfgs.f90

mcsrch.o: mcsrch.f90 libmin.h
	$(FC) -c $(CPP) mcsrch.f90

libmin.o: libmin.c libmin.h
	$(CC) -c libmin.c 

tarball:
	cd ..; tar czf libmin/libmin.tgz libmin/Makefile libmin/*.h libmin/*.f90 libmin/*.c libmin/*.f03
clean:
	$(RM) -f *.o *.a test_fortran test_c

