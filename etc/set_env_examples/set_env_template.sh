#!/bin/bash

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Coral directory:
export CORAL_ROOT=${HOME}/dev/coral_241025_a48/
#===========================================================================


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2Decomp_fft library path:                                                    
export DECOMP2D_ROOT=${HOME}/software/2decomp_emptyFourierFourier/
#===========================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FFTW directory (by default, this is where the FFTW installer put things):
# --- by default, the FFTW installer installs the library in /usr/local/
# --- if you have installed fftw with "./configure --prefix=/my/custom/directory/", 
#     then this prefix is the needed information   
# --- may be defined by an environment variable on clusters.                           
export FFTW_ROOT=${HOME}/software/fftw/fftw-3.3.10/
#===========================================================================



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MKL directory
# --- by default, in /opt/intel/... (with the right year/version numbers)
# --- may be defined by an environment variable on clusters.                           
export MKLROOT=${HOME}/software/intel/oneapi/mkl/2024.0/
#===========================================================================



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Flag for compilation. Pick one (comment/uncomment) below (or create your own)
# ....... pedantic compilation, gnu compiler:
export MPIFLAGS=-std=gnu\ -fmax-errors=6\ -fcheck=all\ -pedantic\ -g
# ....... optimisations, gnu compiler:                 
#export MPIFLAGS=-O2\ -march=native
#===========================================================================




# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Flag for indicating where to put the .mod files.
# Pick the right one for your compiler.
export MOD_DIR_FLAG=-J      # gnu   compiler
#export MOD_DIR_FLAG=-module # intel compiler
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



## ............. BELOW ARE AUTOMATED CHECKS.
## ............. USERS TYPICALLY DO NOT NEED TO TWEAK
## ............. THE CONTENT BELOW THIS LINE.


# ~~~~  
# ~~~~  search for MKL Lapack sources 
# ~~~~  


export LAPACK_SOURCES=${MKLROOT}/include/lapack.f90
echo "[..] Searching for Lapack sources in:"
echo "[..] $LAPACK_SOURCES"
if [ -f "$LAPACK_SOURCES" ]; then
	echo "[OK] Lapack sources have been found!"
else
	export LAPACK_SOURCES=${MKLROOT}/include/mkl_lapack.f90
	echo "[..] Searching for Lapack sources somewhere else:"
	echo "[..] $LAPACK_SOURCES"
	if [ -f "$LAPACK_SOURCES" ]; then
		echo "[OK] Lapack sources have been found."
	else
		echo "/!!\\ Lapack sources have not been found."
		echo "/!!\\ Compilation is likely to fail."
	fi
fi

# ~~~~  
# ~~~~  search for fftw sources 
echo "[..] Searching for FFTW sources and libs in"
echo "[..] $FFTW_ROOT"
# ~~~~  
if [ -f "$FFTW_ROOT/include/fftw3-mpi.f03" ]; then
	echo "[OK] fftw3-mpi.f03 has been found!"
else
	echo "/!!\\ fftw3-mpi.f03 has not been found in $FFTW_ROOT/lib."
	echo "/!!\\ Compilation is likely to fail."
fi
if [ -f "$FFTW_ROOT/lib/libfftw3.a" ]; then
	echo "[OK] libfftw3.a has been found!"
else
	echo "/!!\\ libfftw3.a has not been found in $FFTW_ROOT/lib."
	echo "/!!\\ Compilation is likely to fail."
fi
if [ -f "$FFTW_ROOT/lib/libfftw3_mpi.a" ]; then
	echo "[OK] libfftw3_mpi.a has been found!"
else
	echo "/!!\\ libfftw3_mpi.a has not been found in $FFTW_ROOT/lib."
	echo "/!!\\ Compilation is likely to fail."
fi



