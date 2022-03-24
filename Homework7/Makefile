O = obj

# gfortran compiler
F90 = gfortran
FLAGS = -g -O3 -s

# Portland compiler
#F90 = pgf90
#FLAGS = -fast -Mbounds -Mneginfo -Mdclchk -Mstandard -Knoieee
#FLAGS = -fast -Mnobounds -Mneginfo -Mdclchk -Munroll=c:6 -Mstandard -Knoieee

# Absoft compiler
#F90 = f90
##FLAGS = -O2 -W132 -YEXT_NAMES=LCS  -s -B108 -YCFRL=1
#FLAGS = -O2 -s

all: diffusion

diffusion: constants.h \
       $O/diffusion.o \
       $O/gll_library.o \
       $O/lagrange_poly.o \
       $O/define_derivative_matrix.o
	${F90} $(FLAGS) -o xdiffusion \
       $O/diffusion.o \
       $O/gll_library.o \
       $O/lagrange_poly.o \
       $O/define_derivative_matrix.o

clean:
	rm -f $O/*.o *.o seismogram xwave xdiffusion snapshots/*

$O/diffusion.o: constants.h diffusion.f90
	${F90} $(FLAGS) -c -o $O/diffusion.o diffusion.f90

$O/define_derivative_matrix.o: constants.h define_derivative_matrix.f90
	${F90} $(FLAGS) -c -o $O/define_derivative_matrix.o define_derivative_matrix.f90

$O/gll_library.o: gll_library.f90
	${F90} $(FLAGS) -c -o $O/gll_library.o gll_library.f90

$O/lagrange_poly.o: lagrange_poly.f90
	${F90} $(FLAGS) -c -o $O/lagrange_poly.o lagrange_poly.f90

