#!/bin/bash

export CORAL_ROOT=/home/bmiquel/dev/coral/
export DECOMP2D_ROOT=/home/bmiquel/software/2decomp_light/ 
export FFTW_ROOT=/usr/local/
export FFTW_ROOT=${HOME}/software/fftw-3.3.8/
export MKLROOT=/opt/intel/compilers_and_libraries_2019.1.144/linux/mkl
export MPIFLAGS=-std=gnu\ -fmax-errors=6\ -fcheck=all\ -pedantic\ -g\ -fbacktrace\ 
export MPIFLAGS=-std=f2008\ -fmax-errors=6\ -fcheck=all\ -pedantic\ -g\ -fbacktrace\ -Wextra\ -Wall\ 
export MPIFLAGS=-std=gnu\ -fmax-errors=6\ -fcheck=all\ -pedantic\ -g\ -fbacktrace\ -Wextra\ -Wall\ 
#export MPIFLAGS=-std=gnu\ -fmax-errors=6\ -O2\ -march=native
export MOD_DIR_FLAG=-J
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MKLROOT}/lib/intel64
export GCC_COLORS='locus=01;32'

