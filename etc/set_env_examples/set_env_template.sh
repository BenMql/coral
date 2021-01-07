#!/bin/bash

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Coral directory:
export CORAL_ROOT=${HOME}/dev/coral/
#===========================================================================


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2Decomp_fft library path:                                                    
export DECOMP2D_ROOT=${HOME}/software/2decomp_fft/
#===========================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FFTW directory (by default, this is where the FFTW installer put things):
# --- by default, the FFTW installer installs the library in /usr/local/
# --- if you have installed fftw with "./configure --prefix=/my/custom/directory/", 
#     then this prefix is the needed information   
# --- may be defined by an environment variable on clusters.                           
export FFTW_ROOT=/usr/local/
#===========================================================================



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MKL directory
# --- by default, in /opt/intel/... (with the right year/version numbers)
# --- may be defined by an environment variable on clusters.                           
export MKLROOT=/opt/intel/compilers_and_libraries_2019.1.144/linux/mkl
#===========================================================================




# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Flag for compilation. Pick one (comment/uncomment) below (or create your own)
# ....... pedantic compilation, gnu compiler:
#export MPIFLAGS=-std=gnu\ -fmax-errors=6\ -fcheck=all\ -pedantic\ -g
# ....... optimisations, gnu compiler:                 
export MPIFLAGS=-O2\ -march=native
#===========================================================================




# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Flag for indicating where to put the .mod files.
# Pick the right one for your compiler.
export MOD_DIR_FLAG=-J      # gnu   compiler
export MOD_DIR_FLAG=-module # intel compiler
#===========================================================================




# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Update your LD_LIBRARY_PATH with the path to MKL
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MKLROOT}/lib/intel64
#===========================================================================





# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Miscellanous                                            
# --- Use colors to display errors at compilation
export GCC_COLORS='locus=01;32'
#===========================================================================

