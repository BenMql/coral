#!/bin/bash

export CORAL_ROOT=/home/ben/dev/coral/
export FFTW_ROOT=/usr/local/
export DECOMP2D_ROOT=/home/ben/software/2decomp_light/
export MKLROOT=/opt/intel/compilers_and_libraries_2019.4.243/linux/mkl
export MPIFLAGS=-std=f2008\ -fmax-errors=6\ -fcheck=all\ -pedantic\ -O0\ -g\ -fbacktrace\ -Wextra\ -Wall\ 
#export MPIFLAGS=-pedantic\ -std=f2003
#export MPIFLAGS=-O2\ -std=f2003\ -march=skylake
export MOD_DIR_FLAG=-J
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MKLROOT}/lib/intel64
#MPIFLAGS = -O3 -march=core-avx2
export GCC_COLORS='locus=01;32'

