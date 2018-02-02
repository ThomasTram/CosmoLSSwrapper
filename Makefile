
#FC = gfortran

FC = gcc
CC = gcc

FFLAGS = -cpp -O3 -ffast-math -ffree-line-length-none -fopenmp -fmax-errors=10 -fPIC
LNK =

OUTPUT_DIR = source

vpath %.f90 source
vpath %.o source


UTILS = MiscUtils.o StringUtils.o MpiUtils.o FileUtils.o ObjectLists.o Interpolation.o Romberg.o

WRAPPER = LogLikeCosmoLSS.o pyLogLikeCosmoLSS.o

OBJ = $(addsuffix $(OUTPUT_DIR), $(UTILS) $(WRAPPER))

main: libCosmoLSS.a pyLogLikeCosmoLSS.pyx pyLogLikeCosmoLSS.h
	touch pyLogLikeCosmoLSS.c; rm pyLogLikeCosmoLSS.c
	alias ld=$(CC);CC=$(CC) python setup.py build_ext --inplace

libCosmoLSS.a: $(UTILS) $(WRAPPER)
	cd $(OUTPUT_DIR); ar rvs ../libCosmoLSS.a $(UTILS) $(WRAPPER)

%.o: %.f90
	cd $(OUTPUT_DIR); $(FC) $(FFLAGS) -c ../$< -o $*.o

.PHONY: clean

clean:
	rm -f *.c  $(OUTPUT_DIR)/*.o  $(OUTPUT_DIR)/*.mod libCosmoLSS.a

