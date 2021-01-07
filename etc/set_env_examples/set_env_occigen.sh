#!/bin/bash
module purge
module load intel/17.0
module load openmpi/intel/2.0.4
module load fftw3/3.3.5

export CORAL_ROOT=${HOMEDIR}/dev/pencils_coral_v4/
export FFTW_ROOT=${FFTW3_ROOT}/
export DECOMP2D_ROOT=${HOMEDIR}/software/2decomp_light/
export MPIFLAGS=-O0\ -march=core-avx2
export MOD_DIR_FLAG=-module\ 
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MKLROOT}/lib/intel64
export GCC_COLORS='locus=01;32'

