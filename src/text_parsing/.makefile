test_c: cwraps.f90 cfun_parse_text.c f_call_c.f90
	rm -f *.mod *.o *.x
	gfortran -fcheck=all -pedantic -c cwraps.f90
	gcc -Wextra -pedantic -c cfun_parse_text.c
	gfortran -fcheck=all -pedantic -c f_call_c.f90
	gfortran -fcheck=all -pedantic f_call_c.o cfun_parse_text.o -o toto.x

all: cwraps.f90 cfun_parse_text.c f_call_c.f90
	rm -f *.mod *.o *.x
	gfortran -fcheck=all -pedantic -c cwraps.f90
	gcc -Wextra -pedantic -c cfun_parse_text.c
	gfortran -fcheck=all -pedantic -c ../misc/fortran_kinds.f90
	gfortran -fcheck=all -pedantic -c ../periodic_box/P3_equations.f90
	gfortran -fcheck=all -pedantic -c driver_P3_equations.f90
	gfortran -fcheck=all -pedantic driver_P3_equations.o cfun_parse_text.o P3_equations.o fortran_kinds.o -o toto.x
