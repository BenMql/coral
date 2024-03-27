# Compiler Options: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MPIFC = mpif90 -cpp
MPICC = mpicc 


# Directories: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MODDIR     = $(CORAL_ROOT)/.mod/
BUILD_DIR  = $(CORAL_ROOT)/build/
SPARSE_SRC = $(CORAL_ROOT)/src/sparse_tools/
TEXT_SRC   = $(CORAL_ROOT)/src/text_parsing/
CHEBY_SRC  = $(CORAL_ROOT)/src/cheby_tools/
TSTEP_SRC  = $(CORAL_ROOT)/src/timesteppers/
FFTW_SRC   = $(CORAL_ROOT)/src/fftw_tools/
MISC_SRC   = $(CORAL_ROOT)/src/misc/
OUT_SRC    = $(CORAL_ROOT)/src/output_pack/
LAYER      = $(CORAL_ROOT)/src/plane_layer/
TRIPLY     = $(CORAL_ROOT)/src/triply_periodic/
SLABS      = $(CORAL_ROOT)/src/slab_layer/
MPI_SRC    = $(CORAL_ROOT)/src/MPI_tools/
MKL_LIB    = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64 -I${MKLROOT}/include
PENCILS_LAYER_2DECOMP  = $(CORAL_ROOT)/src/pencils.2decomp/
SLABS_LAYER_FFTW3MPI   = $(CORAL_ROOT)/src/paddedSlabs.fftw3mpi/
SLABS_49_FFTW3MPI   = $(CORAL_ROOT)/src/slabs.fftw3mpi/
FFTW_INC:=$(FFTW_ROOT)/include/
FFTW_LIB:=$(FFTW_ROOT)/lib/
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
	rm -f $(CORAL_ROOT)/src/*/*.o
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
Layer_pencil_Objects += $(TEXT_SRC)ftext_parsing.o
Layer_pencil_Objects += $(MPI_SRC)MPI_vars.o
Layer_pencil_Objects += $(OUT_SRC)output_misc.o          
Layer_pencil_Objects += $(MISC_SRC)timeKeeping.o
Layer_pencil_Objects += $(MISC_SRC)time_mpi.o
Layer_pencil_Objects += $(TSTEP_SRC)IMEX_schemes.o
Layer_pencil_Objects += $(MISC_SRC)read_command_line_args.o
Layer_pencil_Objects += $(MISC_SRC)include_git_version.o
Layer_pencil_Objects += $(MISC_SRC)wallclock.o
Layer_pencil_Objects += $(PENCILS_LAYER_2DECOMP)domain_decomposition.o
Layer_pencil_Objects += $(LAYER)PL_geometry.o
Layer_pencil_Objects += $(LAYER)PL_cheby_misc.o
Layer_pencil_Objects += $(MISC_SRC)include_git_version.o
Layer_pencil_Objects += $(CHEBY_SRC)lapack_module.o
Layer_pencil_Objects += $(LAYER)PL_string_to_data.o
Layer_pencil_Objects += $(LAYER)PL_equations.o             
Layer_pencil_Objects += $(SPARSE_SRC)sparse_formats_d.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_formats_z.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_formats.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_conversions.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_blas.o
Layer_pencil_Objects += $(SPARSE_SRC)sparse_manipulations.o
Layer_pencil_Objects += $(CHEBY_SRC)chebyshev_elementary.o
Layer_pencil_Objects += $(CHEBY_SRC)chebyshev_galerkin_2.o
Layer_pencil_Objects += $(LAYER)PL_algebra_d_1D.o
Layer_pencil_Objects += $(LAYER)PL_algebra_z_1D.o
Layer_pencil_Objects += $(LAYER)PL_algebra_z_3D.o
Layer_pencil_Objects += $(LAYER)PL_algebra.o             
Layer_pencil_Objects += $(PENCILS_LAYER_2DECOMP)transforms.o
Layer_pencil_Objects += $(LAYER)PL_IMEX_timestepping.o
Layer_pencil_Objects += $(LAYER)plane_layer_main.o

pencils: $(Layer_pencil_Objects)
	mkdir -p $(BUILD_DIR)
	$(MPIFC) $(MPIFLAGS) -o $(BUILD_DIR)coral.layer.pencils.2dcmp.x $^ -lfftw3 -l2decomp_fft -L$(FFTW_LIB) -L$(DECOMP2D_ROOT)/lib $(MKL_LIB) -I$(MODDIR) -I$(DECOMP2D_ROOT)/include -I$(MKLROOT)/include  
	cp $(BUILD_DIR)coral.layer.pencils.2dcmp.x $(BUILD_DIR)coral_LP.exe
	$(MAKE) clean


                        #############################
                        #############################
# ~~~~~~~~~~~~~~~~~~~   #  LAYER SLABS four ninths  #  ~~~~~~~~~~~~~~~~~~~
                        #############################
                        #############################



Slabs49_layer_Objects := $(MISC_SRC)chdir_mod.o    
Slabs49_layer_Objects += $(MISC_SRC)fortran_kinds.o
Slabs49_layer_Objects += $(MISC_SRC)read_command_line_args.o
Slabs49_layer_Objects += $(FFTW_SRC)fftw3_mpi.o
Slabs49_layer_Objects += $(TEXT_SRC)cwraps.o
Slabs49_layer_Objects += $(TEXT_SRC)cfun_parse_text.o
Slabs49_layer_Objects += $(TEXT_SRC)ftext_parsing.o
Slabs49_layer_Objects += $(MPI_SRC)MPI_vars.o
Slabs49_layer_Objects += $(OUT_SRC)output_misc.o          
Slabs49_layer_Objects += $(MISC_SRC)timeKeeping.o
Slabs49_layer_Objects += $(MISC_SRC)time_mpi.o
Slabs49_layer_Objects += $(TSTEP_SRC)IMEX_schemes.o
Slabs49_layer_Objects += $(MISC_SRC)read_command_line_args.o
Slabs49_layer_Objects += $(MISC_SRC)include_git_version.o
Slabs49_layer_Objects += $(MISC_SRC)wallclock.o
Slabs49_layer_Objects += $(SLABS_49_FFTW3MPI)domain_decomposition.o
Slabs49_layer_Objects += $(LAYER)PL_geometry.o
Slabs49_layer_Objects += $(LAYER)PL_cheby_misc.o
Slabs49_layer_Objects += $(MISC_SRC)include_git_version.o
Slabs49_layer_Objects += $(CHEBY_SRC)lapack_module.o
Slabs49_layer_Objects += $(LAYER)PL_string_to_data.o
Slabs49_layer_Objects += $(LAYER)PL_equations.o             
Slabs49_layer_Objects += $(SPARSE_SRC)sparse_formats_d.o
Slabs49_layer_Objects += $(SPARSE_SRC)sparse_formats_z.o
Slabs49_layer_Objects += $(SPARSE_SRC)sparse_formats.o
Slabs49_layer_Objects += $(SPARSE_SRC)sparse_conversions.o
Slabs49_layer_Objects += $(SPARSE_SRC)sparse_blas.o
Slabs49_layer_Objects += $(SPARSE_SRC)sparse_manipulations.o
Slabs49_layer_Objects += $(CHEBY_SRC)chebyshev_elementary.o
Slabs49_layer_Objects += $(CHEBY_SRC)chebyshev_galerkin_2.o
Slabs49_layer_Objects += $(LAYER)PL_algebra_d_1D.o
Slabs49_layer_Objects += $(LAYER)PL_algebra_z_1D.o
Slabs49_layer_Objects += $(LAYER)PL_algebra_z_3D.o
Slabs49_layer_Objects += $(LAYER)PL_algebra.o             
Slabs49_layer_Objects += $(SLABS_49_FFTW3MPI)transforms.o
Slabs49_layer_Objects += $(LAYER)PL_IMEX_timestepping.o
Slabs49_layer_Objects += $(LAYER)plane_layer_main.o


slabs: $(Slabs49_layer_Objects)
	mkdir -p $(BUILD_DIR)
	$(MPIFC) $(MPIFLAGS) -o $(BUILD_DIR)coral.layer.slabs.fftw3mpi.x $^ $(MPI_FFTW_link) -L$(FFTW_LIB) -I$(FFTW_INC) $(MKL_LIB) -I$(MODDIR) -I$(MKLROOT)/include  
	cp $(BUILD_DIR)coral.layer.slabs.fftw3mpi.x $(BUILD_DIR)coral_SL.exe
	$(MAKE) clean






                        #####################
                        #####################
# ~~~~~~~~~~~~~~~~~~~   #  TRIPLY PERIODIC  #~~~~~~~~~~~~~~~~~~
                        #####################
                        #####################


Triply_periodic_Objects := $(MISC_SRC)chdir_mod.o    
Triply_periodic_Objects += $(MISC_SRC)fortran_kinds.o
Triply_periodic_Objects += $(MISC_SRC)read_command_line_args.o
Triply_periodic_Objects += $(TEXT_SRC)cwraps.o
Triply_periodic_Objects += $(TEXT_SRC)cfun_parse_text.o
Triply_periodic_Objects += $(TEXT_SRC)ftext_parsing.o
Triply_periodic_Objects += $(MPI_SRC)MPI_vars.o
Triply_periodic_Objects += $(OUT_SRC)output_misc.o          
Triply_periodic_Objects += $(MISC_SRC)timeKeeping.o
Triply_periodic_Objects += $(TSTEP_SRC)IMEX_schemes.o
Triply_periodic_Objects += $(MISC_SRC)time_mpi.o
Triply_periodic_Objects += $(MISC_SRC)read_command_line_args.o
Triply_periodic_Objects += $(MISC_SRC)include_git_version.o
Triply_periodic_Objects += $(MISC_SRC)wallclock.o
Triply_periodic_Objects += $(TRIPLY)P3_domain_decomp.o
Triply_periodic_Objects += $(TRIPLY)P3_geometry.o
Triply_periodic_Objects += $(MISC_SRC)include_git_version.o
Triply_periodic_Objects += $(CHEBY_SRC)lapack_module.o
Triply_periodic_Objects += $(TRIPLY)P3_lapack_wrappers.o
Triply_periodic_Objects += $(TRIPLY)P3_algebra.o
Triply_periodic_Objects += $(TRIPLY)P3_string_to_data.o
Triply_periodic_Objects += $(TRIPLY)P3_equations.o             
Triply_periodic_Objects += $(TRIPLY)P3_transforms.o
Triply_periodic_Objects += $(TRIPLY)P3_IMEX_timestepping.o
Triply_periodic_Objects += $(TRIPLY)triply_periodic_main.o


triply: $(Triply_periodic_Objects)
	mkdir -p $(BUILD_DIR)
	$(MPIFC) $(MPIFLAGS) -o $(BUILD_DIR)coral_P3.exe $^ -l2decomp_fft -L$(DECOMP2D_ROOT)/lib $(MKL_LIB) -I$(MODDIR) -I$(DECOMP2D_ROOT)/include -I$(MKLROOT)/include  
	$(MAKE) clean





                        #############################
                        #############################
# ~~~~~~~~~~~~~~~~~~~   #      CHECK EQUATIONS      #  ~~~~~~~~~~~~~~~~~~~
                        #############################
                        #############################



Slabs49_layer_Objects := $(MISC_SRC)chdir_mod.o    
Slabs49_layer_Objects += $(MISC_SRC)fortran_kinds.o
Slabs49_layer_Objects += $(MISC_SRC)read_command_line_args.o
Slabs49_layer_Objects += $(TEXT_SRC)cwraps.o
Slabs49_layer_Objects += $(TEXT_SRC)cfun_parse_text.o
Slabs49_layer_Objects += $(TEXT_SRC)ftext_parsing.o
Slabs49_layer_Objects += $(MPI_SRC)MPI_vars.o
Slabs49_layer_Objects += $(MISC_SRC)read_command_line_args.o
Slabs49_layer_Objects += $(MISC_SRC)include_git_version.o
Slabs49_layer_Objects += $(LAYER)PL_string_to_data.o
Slabs49_layer_Objects += $(LAYER)PL_equations.o             
Slabs49_layer_Objects += $(LAYER)check_equations.o


checks: $(Slabs49_layer_Objects)
	mkdir -p $(BUILD_DIR)
	$(MPIFC) $(MPIFLAGS) -o $(BUILD_DIR)coral.checks.x $^ -I$(MODDIR)   
	$(MAKE) clean






