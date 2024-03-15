
FC ?= gcc
CC ?= gcc

FFLAGS = -cpp -O3 -ffree-line-length-none -fopenmp -fmax-errors=10 -fbounds-check -g -fbacktrace -ffpe-trap=invalid
#FFLAGS = -cpp -fPIC -fbounds-check -g -O0 -Wall -Wextra -fbacktrace -ffpe-trap=invalid,zero,overflow
LNK =

OUTPUT_DIR = source

vpath %.f90 source
vpath %.o source


UTILS = MiscUtils.o MpiUtils.o StringUtils.o FileUtils.o ObjectLists.o Interpolation.o Romberg.o

WRAPPER = LogLikeCosmoLSS.o pyLogLikeCosmoLSS.o

OBJ = $(addsuffix $(OUTPUT_DIR), $(UTILS) $(WRAPPER))

main: libCosmoLSS.a pyLogLikeCosmoLSS.pyx pyLogLikeCosmoLSS.h
	touch pyLogLikeCosmoLSS.c; rm pyLogLikeCosmoLSS.c
	alias ld=$(CC);CC=$(CC) python setup.py build_ext --inplace


LogLikeCosmoLSS.o: Interpolation.o

pyLogLikeCosmoLSS.o: LogLikeCosmoLSS.o

StringUtils.o: MiscUtils.o MpiUtils.o

FileUtils.o: StringUtils.o MiscUtils.o MpiUtils.o

ObjectLists.o: FileUtils.o MpiUtils.o

Interpolation.o: FileUtils.o MpiUtils.o MiscUtils.o StringUtils.o ObjectLists.o

libCosmoLSS.a: $(UTILS) $(WRAPPER)
	cd $(OUTPUT_DIR); ar rvs ../libCosmoLSS.a $(UTILS) $(WRAPPER)

%.o: %.f90
	cd $(OUTPUT_DIR); $(FC) $(FFLAGS) -c ../$< -o $*.o

.PHONY: clean

clean:
	rm -f *.c  $(OUTPUT_DIR)/*.o  $(OUTPUT_DIR)/*.mod libCosmoLSS.a

