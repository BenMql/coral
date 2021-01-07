#!/bin/bash
module purge
module load intel
module load impi
module load fftw

export CORAL_ROOT=/home/bemi5354/coral/
export FFTW_ROOT=${CURC_FFTW_ROOT}/
export MKLROOT=${CURC_INTEL_ROOT}/compilers_and_libraries_2017.4.196/linux/mkl
#export MPIFLAGS=-stand f95 
#export MPIFLAGS=-O3\ -xCORE-AVX2
export MPIFLAGS=-O0\ -check\ all\ -stand\ f03
export MOD_DIR_FLAG=-module\ 
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MKLROOT}/lib/intel64
#MPIFLAGS = -O3 -march=core-avx2

