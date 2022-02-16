# Compiler Options: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MPIFC = mpif90 -cpp
MPICC = mpicc 


# Directories: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DECOMP_SRC = $(HOME)/software/2decomp_fft/src/
MODDIR     = $(CORAL_ROOT)/.mod/
BUILD_DIR  = $(CORAL_ROOT)/build/
SPARSE_SRC = $(CORAL_ROOT)/src/sparse_tools/
TEXT_SRC   = $(CORAL_ROOT)/src/text_parsing/
CHEBY_SRC  = $(CORAL_ROOT)/src/cheby_tools/
TSTEP_SRC  = $(CORAL_ROOT)/src/timesteppers/
FFTW_SRC   = $(CORAL_ROOT)/src/fftw_tools/
MISC_SRC   = $(CORAL_ROOT)/src/misc/
OUT_SRC    = $(CORAL_ROOT)/src/output_pack/
PENCILS    = $(CORAL_ROOT)/src/layer_pencils/
TRIPLY     = $(CORAL_ROOT)/src/triply_periodic/
SLABS      = $(CORAL_ROOT)/src/slab_layer/
MPI_SRC    = $(CORAL_ROOT)/src/MPI_tools/
MKL_LIB    = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64 -I${MKLROOT}/include
FFTW_INC:=$(FFTW_ROOT)/include/
FFTW_LIB:=$(FFTW_ROOT)lib/
JMODDIR = $(MOD_DIR_FLAG)$(MODDIR)
MPI_FFTW_link = -lfftw3_mpi -lfftw3 -lm

$(info $(shell mkdir -p $(BUILD_DIR)))
$(info $(shell mkdir -p $(MODDIR)))

$(CHEBY_SRC)lapack_module.o: ${MKLROOT}/include/lapack.f90
	$(MPIFC) $(MPIFLAGS) -c $^ -o $@ $(JMODDIR)

$(FFTW_SRC)fftw3_wrap.o: $(FFTW_SRC)/fftw3_wrap.f90
	$(MPIFC) $(MPIFLAGS) -c $^ -o $@ $(JMODDIR) -I$(FFTW_INC)

$(FFTW_SRC)fftw3_mpi.o: $(FFTW_SRC)/fftw3_mpi.f90
	$(MPIFC) $(MPIFLAGS) -c $^ -o $@ $(JMODDIR) -I$(FFTW_INC)


GIT_VERSION:="$(shell git describe --abbrev=4 --dirty --always --tags)"
.PHONY: get_git_version
get_git_version:
	printf "Module include_git_version\nImplicit None\ncontains\nsubroutine display_git_version()\nprint *, \"$(GIT_VERSION)\"\nend subroutine\nend module" > $(MISC_SRC)include_git_version.f90

$(MISC_SRC)include_git_version.o: get_git_version
	$(MPIFC) $(MPIFLAGS) -c $(MISC_SRC)include_git_version.f90 -o $@ $(JMODDIR)

%.o: %.f90
	$(MPIFC) $(MPIFLAGS) -c $^ -o $@ $(JMODDIR) -I$(DECOMP2D_ROOT)/include

%.o: %.c  
	$(MPICC) -Wextra -pedantic -c $^ -o $@ 


clean:
	rm -f $(CORAL_ROOT)src/*/*.o
	rm -f $(MODDIR)*.mod

clean_o: 
	rm -f $(CORAL_ROOT)src/*/*.o

clean_mod: 
	rm -f $(MODDIR)*mod



                        ##################
                        ##################
# ~~~~~~~~~~~~~~~~~~~   #  LAYER PENCILS #  ~~~~~~~~~~~~~~~~~~~
                        ##################
                        ##################


Layer_pencil_Objects := $(MISC_SRC)chdir_mod.o    
Layer_pencil_Objects += $(MISC_SRC)fortran_kinds.o
Layer_pencil_Objects += $(MISC_SRC)read_command_line_args.o
Layer_pencil_Objects += $(FFTW_SRC)fftw3_wrap.o
Layer_pencil_Objects += $(TEXT_SRC)cwraps.o
Layer_pencil_Objects += $(TEXT_SRC)cfun_parse_text.o
Layer_pencil_Objects += $(PENCILS)LP_text_parsing.o
Layer_pencil_Objects += $(MPI_SRC)MPI_vars.o
Layer_pencil_Objects += $(PENCILS)LP_output.o          
Layer_pencil_Objects += $(PENCILS)LP_timings.o
Layer_pencil_Objects += $(MISC_SRC)time_mpi.o
Layer_pencil_Objects += $(TSTEP_SRC)IMEX_schemes.o
Layer_pencil_Objects += $(MISC_SRC)time_mpi.o
Layer_pencil_Objects += $(MISC_SRC)read_command_line_args.o
Layer_pencil_Objects += $(MISC_SRC)include_git_version.o
Layer_pencil_Objects += $(PENCILS)LP_wallclock.o
Layer_pencil_Objects += $(PENCILS)LP_domain_decomp.o
Layer_pencil_Objects += $(PENCILS)LP_geometry.o
Layer_pencil_Objects += $(PENCILS)LP_cheby_misc.o
Layer_pencil_Objects += $(MISC_SRC)include_git_version.o
Layer_pencil_Objects += $(CHEBY_SRC)lapack_module.o
Layer_pencil_Objects += $(PENCILS)LP_string_to_data.o
Layer_pencil_Objects += $(PENCILS)LP_equations.o             
Layer_pencil_Objects += $(SPARSE_SRC)sparse_formats_d.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_formats_z.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_formats.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_conversions.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_blas.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_manipulations.o
Layer_pencil_Objects += $(CHEBY_SRC)chebyshev_elementary.o
Layer_pencil_Objects += $(CHEBY_SRC)chebyshev_galerkin_2.o
Layer_pencil_Objects += $(PENCILS)LP_algebra_d_1D.o
Layer_pencil_Objects += $(PENCILS)LP_algebra_z_1D.o
Layer_pencil_Objects += $(PENCILS)LP_algebra_z_3D.o
Layer_pencil_Objects += $(PENCILS)LP_algebra.o             
Layer_pencil_Objects += $(PENCILS)LP_transforms.o
Layer_pencil_Objects += $(PENCILS)LP_IMEX_timestepping.o
Layer_pencil_Objects += $(PENCILS)layer_pencils_main.o

pencils: $(Layer_pencil_Objects)
	mkdir -p $(BUILD_DIR)
	$(MPIFC) $(MPIFLAGS) -o $(BUILD_DIR)coral_LP.exe $^ -lfftw3 -l2decomp_fft -L$(FFTW_LIB) -L$(DECOMP2D_ROOT)/lib $(MKL_LIB) -I$(MODDIR) -I$(DECOMP2D_ROOT)/include -I$(MKLROOT)/include  
	$(MAKE) clean


                        #####################
                        #####################
# ~~~~~~~~~~~~~~~~~~~   #  TRIPLY PERIODIC  #~~~~~~~~~~~~~~~~~~
                        #####################
                        #####################


Triply_periodic_Objects := $(MISC_SRC)fortran_kinds.o
Triply_periodic_Objects += $(CHEBY_SRC)lapack_module.o
Triply_periodic_Objects += $(TRIPLY)P3_lapack_wrappers.o
Triply_periodic_Objects += $(TRIPLY)P3_algebra.o
Triply_periodic_Objects += $(TRIPLY)algebra_driver.o

algebra: $(Triply_periodic_Objects)
	mkdir -p $(BUILD_DIR)
	$(MPIFC) $(MPIFLAGS) -o $(BUILD_DIR)algebra_driver.exe $^ $(MKL_LIB) -I$(MODDIR) -I$(MKLROOT)/include  
	$(MAKE) clean



                        ##################
                        ##################
# ~~~~~~~~~~~~~~~~~~~   #  SLABS  LAYER  #  ~~~~~~~~~~~~~~~~~~~
                        ##################
                        ##################


Slabs_layer_Objects := $(MISC_SRC)chdir_mod.o    
Slabs_layer_Objects += $(MISC_SRC)fortran_kinds.o
Slabs_layer_Objects += $(MISC_SRC)read_command_line_args.o
Slabs_layer_Objects += $(FFTW_SRC)fftw3_mpi.o
Slabs_layer_Objects += $(TEXT_SRC)cwraps.o
Slabs_layer_Objects += $(TEXT_SRC)cfun_parse_text.o
Slabs_layer_Objects += $(SLABS)LP_text_parsing.o
Slabs_layer_Objects += $(MPI_SRC)MPI_vars.o
Slabs_layer_Objects += $(SLABS)LP_output.o          
Slabs_layer_Objects += $(SLABS)LP_timings.o
Slabs_layer_Objects += $(MISC_SRC)time_mpi.o
Slabs_layer_Objects += $(TSTEP_SRC)IMEX_schemes.o
Slabs_layer_Objects += $(MISC_SRC)time_mpi.o
Slabs_layer_Objects += $(MISC_SRC)read_command_line_args.o
Slabs_layer_Objects += $(MISC_SRC)include_git_version.o
Slabs_layer_Objects += $(SLABS)LP_wallclock.o
Slabs_layer_Objects += $(SLABS)SL_domain_decomp.o
Slabs_layer_Objects += $(SLABS)SL_geometry.o
Slabs_layer_Objects += $(SLABS)SL_cheby_misc.o
Slabs_layer_Objects += $(MISC_SRC)include_git_version.o
Slabs_layer_Objects += $(CHEBY_SRC)lapack_module.o
Slabs_layer_Objects += $(SLABS)SL_string_to_data.o
Slabs_layer_Objects += $(SLABS)SL_equations.o             
Slabs_layer_Objects += $(SPARSE_SRC)sparse_formats_d.o
Slabs_layer_Objects += $(SPARSE_SRC)sparse_formats_z.o
Slabs_layer_Objects += $(SPARSE_SRC)sparse_formats.o
Slabs_layer_Objects += $(SPARSE_SRC)sparse_conversions.o
Slabs_layer_Objects += $(SPARSE_SRC)sparse_blas.o
Slabs_layer_Objects += $(SPARSE_SRC)sparse_manipulations.o
Slabs_layer_Objects += $(CHEBY_SRC)chebyshev_elementary.o
Slabs_layer_Objects += $(CHEBY_SRC)chebyshev_galerkin_2.o
Slabs_layer_Objects += $(SLABS)SL_algebra_d_1D.o
Slabs_layer_Objects += $(SLABS)SL_algebra_z_1D.o
Slabs_layer_Objects += $(SLABS)SL_algebra_z_3D.o
Slabs_layer_Objects += $(SLABS)SL_algebra.o             
Slabs_layer_Objects += $(SLABS)SL_transforms.o
#Slabs_layer_Objects += $(PENCILS)LP_IMEX_timestepping.o
Slabs_layer_Objects += $(SLABS)slabs_layer_main.o
#
slabs: $(Slabs_layer_Objects)
	mkdir -p $(BUILD_DIR)
	$(MPIFC) $(MPIFLAGS) -o $(BUILD_DIR)coral_SL.exe $^ $(MPI_FFTW_link) -L$(FFTW_LIB) -L$(DECOMP2D_ROOT)/lib $(MKL_LIB) -I$(MODDIR) -I$(MKLROOT)/include  
	$(MAKE) clean




