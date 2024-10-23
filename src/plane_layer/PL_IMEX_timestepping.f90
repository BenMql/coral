 !============================================================================= 
 !                            C O R A L 
 !============================================================================= 
 ! 
 ! MODULE: PL_IMEX_timestepping
 !  
 !> @author 
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com 
 ! 
 ! DESCRIPTION 
 !> core routines and data structure implementing the time-marching schemes
 !> of Ascher-Ruuth-Spiteri 1997 (Applied Numerical Mathematics 25, pp151-167)
 !> [denoted below as ARS97].
 ! 
 !============================================================================= 


module PL_IMEX_timestepping


 use PL_algebra
 use PL_cheby_misc
 use PL_equations
 use PL_geometry
 use output_misc
 use IMEX_schemes
 use sparse_manipulations
 use chebyshev_elementary
 use chebyshev_galerkin_2
 use transforms
 !use decomp_2d_io
 use timeKeeping, only: timings
 implicit none

 type :: intermediate_step_buffer_3d_T !< intermediate 3D fields (cf ARS97)
   complex(kind=dp), allocatable :: K_hat(:,:,:) !< eq.(2.1a)
   complex(kind=dp), allocatable :: K_std(:,:,:) !< eq.(2.1b)
 end type intermediate_step_buffer_3d_T


 type :: intermediate_step_buffer_1d_T !< intermediate 1D field (cf ARS97)
   real(kind=dp), allocatable :: K_hat(:) !< eq.(2.1a)
   real(kind=dp), allocatable :: K_std(:) !< eq.(2.1b)
 end type intermediate_step_buffer_1d_T

 type :: building_tools_T
    type(dcsr_matrix), allocatable :: cheb_IM(:,:) !< integration and mul. matrix 
 end type

 type :: system_shape_T !< some book-keeping
   integer, allocatable :: equation_firstIndex(:) !< starting row for each eqn
   integer, allocatable :: equation_lastIndex(:)  !< last row for each eqn
   integer, allocatable :: equation_size(:)       !< number of independant projections, for each eqn
   integer, allocatable :: variable_firstIndex(:) !< starting column for each var
   integer, allocatable :: variable_lastIndex(:)  !< last column for each eqn
   integer, allocatable :: variable_size(:)       !< number of independant coefs, for each var
   integer, allocatable :: spectral_local_NX      !< # of in-core kx modes
   integer, allocatable :: spectral_local_NY      !< # of in-core ky modes 
   integer, allocatable :: physical_local_NX      !< # of in-core kx grid points
   integer, allocatable :: physical_local_NY      !< # of in-core ky grid points
   integer :: total !< for each Fourier mode, # of Chebyshev pols. for *all* variables
 end type system_shape_T
 


 type :: coupled_set_1d_ops_vars_T
   real(kind=dp), allocatable :: field(:)
   real(kind=dp), allocatable :: rhs  (:)
   real(kind=dp), allocatable :: aux  (:)
   type(dOperator_1d_1coupled_T) :: mass
   type(dOperator_1d_1coupled_T) :: stif
   type(dOperator_1d_1coupled_T) :: evol
   ! type(sest)
   type(intermediate_step_buffer_1d_T), allocatable :: step(:)
   type(building_tools_T) :: building_tools
   type(dcsr_matrix), allocatable :: stencil(:)
   type(dOperator_1d_1coupled_T), allocatable :: square_stencil(:)
   type(system_shape_T) :: shape
   type(dcsr_matrix) :: shuffleP_polyDegree_to_var
   type(dcsr_matrix) :: shuffleQ_var_to_polyDegree
   type(dcsr_matrix) :: shuffleP_polyDegree_to_var_T
   type(dcsr_matrix) :: shuffleQ_var_to_polyDegree_T
   type(dcsr_matrix), allocatable :: padStencilExtractShuffle(:)
   type(dcsr_matrix), allocatable :: shuffleInsert(:)
   type(dcsr_matrix), allocatable :: shuffleTextractTtruncate(:)
 contains
   procedure :: add_block => add_dsca_times_MIblock_to_zero_operator
 end type coupled_set_1d_ops_vars_T

 type :: coupled_set_3d_ops_vars_T
   complex(kind=dp), allocatable :: field(:,:,:)
   complex(kind=dp), allocatable :: rhs  (:,:,:)
   complex(kind=dp), allocatable :: aux  (:,:,:)
   type(zOperator_3d_1coupled_T) :: mass
   type(zOperator_3d_1coupled_T) :: stif
   type(zOperator_3d_1coupled_T) :: evol
   ! type(sest)
   type(intermediate_step_buffer_3d_T), allocatable :: step(:)
   type(building_tools_T) :: building_tools
   type(dcsr_matrix), allocatable :: stencil(:)
   type(zOperator_1d_1coupled_T), allocatable :: square_stencil(:)
   type(dcsr_matrix) :: shuffleP_polyDegree_to_var
   type(dcsr_matrix) :: shuffleQ_var_to_polyDegree
   type(dcsr_matrix) :: shuffleP_polyDegree_to_var_T
   type(dcsr_matrix) :: shuffleQ_var_to_polyDegree_T
   type(system_shape_T) :: shape
   type(zcsr_matrix), allocatable :: padStencilExtractShuffle(:)
   type(zcsr_matrix), allocatable :: shuffleTextractTtruncate(:)
   type(zcsr_matrix), allocatable :: shuffleInsert(:)
 contains
   procedure :: add_block => add_zsca_times_MIblock_to_kxky_operators
 end type coupled_set_3d_ops_vars_T

 type :: output_bookkeeping_T
   integer :: dPosition = 1
   logical :: first_or_not = .true.
   character(len=:), allocatable :: output_directory
   integer :: output_dir_length
   real(dp) :: absolute_time
   integer :: rolling_integer = 1
   integer :: time_integer = 0
 end type
 
 type :: gauss_chebyshev_T
   real(kind=dp), allocatable :: weight1d(:)
   real(kind=dp), allocatable :: grid1d(:)
 end type

 type :: cargo_T
   real(kind=dp) :: KE_display
   real(kind=dp) :: cfl_based_DT
   real(kind=dp) :: cflFactor_along_z
   real(kind=dp) :: initialCondition_amp_kxky
   real(kind=dp) :: initialCondition_amp_zero
   real(kind=dp) :: smagorinsky_prefactor
 end type

 type :: sources_arrays_zero_T
   real(kind=dp), allocatable :: physical(:)
   real(kind=dp), allocatable :: spectral(:)
 end type sources_arrays_zero_T

 type :: full_problem_data_structure_T
   type(full_problem_recipe_T) :: recipe
   type(coupled_set_3d_ops_vars_T), allocatable :: coupled_kxky_set(:)
   type(coupled_set_1d_ops_vars_T), allocatable :: coupled_zero_set(:)
   type(scheme_handle) :: shandle
   type(geometry_vars_t) :: geometry
   type(PS_fields_T), allocatable :: quadratic_variables(:)
   type(PS_fields_T), allocatable :: linear_variables(:)
   type(output_bookkeeping_T) :: io_bookkeeping
   type(gauss_chebyshev_T) :: gauss_cheby
   type(cargo_T) :: cargo
   type(sources_arrays_zero_T), allocatable :: sources(:)
   integer :: numberOf_sources
   type(zOperator_1d_1coupled_T) :: Chebyshev_Integration_z
   type(dOperator_1d_1coupled_T) :: Chebyshev_Integration_d
   real(kind=dp), allocatable :: timeseries_buffer (:,:)
   integer :: buffer_length = 20 ! arbitrarily...
  contains 
   procedure :: add_K_std_to_rhs   => LPIMEX_add_K_std_to_rhs
   procedure :: add_K_hat_to_rhs   => LPIMEX_add_K_hat_to_rhs
   procedure :: allocate_all_buffers
   procedure :: backsolve_to_aux   => LPIMEX_backsolve_to_aux
   procedure :: build_operators    => build_all_matrices
   procedure :: build_kxky_mass_matrices_for_system     
   procedure :: build_kxky_stif_matrices_for_system     
   procedure :: build_kxky_shuffle_matrices_for_system     
   procedure :: build_zero_mass_matrices_for_system     
   procedure :: build_zero_stif_matrices_for_system     
   procedure :: build_zero_shuffle_matrices_for_system     
   procedure :: compute_full_variables_in_physical_space
   procedure :: compute_NL_terms_to_K_hat => LPIMEX_compute_NL_terms_to_K_hat
   procedure :: compute_systems_shapes
   procedure :: compute_cfl_based_timestep
   procedure :: copy_fields_to_aux => LPIMEX_copy_fields_to_aux
   procedure :: copy_aux_to_fields => LPIMEX_copy_aux_to_fields
   procedure :: differentiate_kxky, differentiate_zero
   generic   :: differentiate => differentiate_kxky, differentiate_zero
   procedure :: export_allPhys
   procedure :: export_allProfiles
   procedure :: export_slice
   procedure :: export_verticallyAvgedSlice => export_verticallyAvgedSlice_general
   procedure :: export_volume
   procedure :: export_Profile
   procedure :: export_Profile_fftw_mpi
   procedure :: export_Profile_decomp2d
   procedure :: export_CheckPoints
   procedure :: factorize_operators 
   procedure :: kxky_matrices_manyGalkn_to_uniqueCheby
   procedure :: import_quickSave
   procedure :: init => initialise_the_data_structure
   procedure :: initialise_bookkeeping_counters
   procedure :: march_forward => one_step_beyond
   procedure :: mass_field_to_rhs  => LPIMEX_mass_field_to_rhs
   procedure :: noisy_initial_conditions
   procedure :: output_global_quantities
   procedure :: output_slices_volumes_and_profiles
   procedure :: prepare_building_tools
   procedure :: prepare_stencils
   procedure :: prepare_chebyshev_integration
   procedure :: remove_vertical_mean_of_quadratic_var
   procedure :: remove_vertical_mean_of_linear_var
   procedure :: remove_vertical_fluctuations_of_linear_var
   procedure :: stif_aux_to_K_std  => LPIMEX_stif_aux_to_K_std
   procedure :: zero_matrices_manyGalkn_to_uniqueCheby
   procedure :: allocate_sources
   procedure :: add_source_array
   procedure :: add_source_dscalar
   generic   :: add_source => add_source_dscalar, &
                              add_source_array
   procedure :: transform_sources
   procedure :: sourceParams
   procedure :: init_all_sources
   procedure :: init_penalisation
   procedure :: horizontal_average_of_physical_quantity_inplace
 end type full_problem_data_structure_T

 contains

 subroutine remove_vertical_mean_of_quadratic_var( self, iQ)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iQ
   integer :: q
   self%quadratic_variables (iQ) %spec (1, :,:) = cmplx(0._dp, 0._dp, kind=dp)
   do q =1, self%geometry%NZ /2-1
      self%quadratic_variables(iQ) %spec (1,:,:)           &
       =   self%quadratic_variables(iQ) %spec (1,:,:)      &
       -   self%quadratic_variables(iQ) %spec (2*q+1,:,:)  &
           * (  0.5_dp  / (2._dp * q + 1._dp)              &
             -  0.5_dp  / (2._dp * q - 1._dp)  )
   end do
 end subroutine


 subroutine remove_vertical_mean_of_linear_var( self, iL)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iL
   integer :: q
   self%linear_variables (iL) %spec (1, :,:) = cmplx(0._dp, 0._dp, kind=dp)
   do q =1, self%geometry%NZ /2-1
      self%linear_variables (iL) %spec (1,:,:)           &
       =   self%linear_variables (iL) %spec (1,:,:)      &
       -   self%linear_variables (iL) %spec (2*q+1,:,:)  &
           * (  0.5_dp  / (2._dp * q + 1._dp)              &
             -  0.5_dp  / (2._dp * q - 1._dp)  )
   end do
 end subroutine

 subroutine remove_vertical_fluctuations_of_linear_var( self, iL)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iL
   integer :: q
   do q =1, self%geometry%NZ /2-1
      ! ...  
      ! instead of subtracting, we *add* the mean
      ! ...  
      self%linear_variables (iL) %spec (1,:,:)           &
       =   self%linear_variables (iL) %spec (1,:,:)      &
       +   self%linear_variables (iL) %spec (2*q+1,:,:)  &
           * (  0.5_dp  / (2._dp * q + 1._dp)              &
             -  0.5_dp  / (2._dp * q - 1._dp)  )
      ! ...  
      ! and these coefficients must then be zeroed
      ! warning: check that there is no ordering issue here 
      !          (the compiler must not re-order these instructions: 
      !           cancellation should occur *after* summation)
      ! ...  
      self%linear_variables (iL) %spec (2*q  ,:,:) = cmplx(0._dp, 0._dp, kind=dp)
      self%linear_variables (iL) %spec (2*q+1,:,:) = cmplx(0._dp, 0._dp, kind=dp)
   end do
 end subroutine

 subroutine differentiate_kxky ( self, field, dOrder)
   class(full_problem_data_structure_T), intent(inOut) :: self
   complex(kind=dp), allocatable, intent(inOut) :: field(:,:,:)
   complex(kind=dp), allocatable, target :: deriv(:,:,:)
   complex(kind=dp), pointer :: fDummyPtr(:)
   type(C_Ptr) :: dummyPtr
   integer :: ix,iy,iz
   integer, intent(in) :: dOrder
   integer :: iOrder
   if (dOrder.lt.1) then
     print *, 'subroutine differentiate exists only for dOrder > 0'
     error stop
   else 
   do iOrder = 1, dOrder
   allocate (deriv( domain_decomp% spec_iSize(1), &
                    domain_decomp% spec_iSize(2), &
                    domain_decomp% spec_iSize(3)))
      
   deriv = cmplx(0._dp, 0._dp, kind=dp)
   do ix = 1, domain_decomp% spec_iSize(3)
   do iy = 1, domain_decomp% spec_iSize(2)
   do iz = 1, self%geometry%NZ - 1
   deriv(iz,iy,ix) = field(iz+1, iy, ix)
   end do
   end do
   end do
   dummyPtr = C_loc(deriv(1,1,1))
   call C_F_pointer(dummyPtr, fDummyPtr, [domain_decomp% spec_iSize(1)*&
                                          domain_decomp% spec_iSize(2)*&
                                          domain_decomp% spec_iSize(3)])
   call self%Chebyshev_integration_z%backsolve( fDummyPtr,&
                                      domain_decomp% spec_iSize(3)*&
                                      domain_decomp% spec_iSize(2),&
                                      domain_decomp% spec_iSize(1))
              
   call move_alloc(from=deriv, to=field)
   end do
   do ix = 1, domain_decomp% spec_iSize(3)
   do iy = 1, domain_decomp% spec_iSize(2)
   field(self%geometry%NZ : domain_decomp% spec_iSize(1),iy,ix) = cmplx(0._dp, 0._dp, kind=dp)
   end do
   end do
   end if
 end subroutine
   
 subroutine differentiate_zero ( self, field, dOrder)
   class(full_problem_data_structure_T), intent(inOut) :: self
   real(kind=dp), allocatable, intent(inOut) :: field(:)
   real(kind=dp), allocatable, target :: deriv(:)
   real(kind=dp), pointer :: fDummyPtr(:)
   type(C_Ptr) :: dummyPtr
   integer :: iz
   integer, intent(in) :: dOrder
   integer :: iOrder
   if (dOrder.lt.1) then
     print *, 'subroutine differentiate exists only for dOrder > 0'
     error stop
   else 
   do iOrder = 1, dOrder
   allocate (deriv( domain_decomp% spec_iSize(1) ))
   deriv = 0._dp                             
   do iz = 1, self%geometry%NZ - 1
   deriv(iz) = field(iz+1)
   end do
   dummyPtr = C_loc(deriv(1))
   call C_F_pointer(dummyPtr, fDummyPtr, [domain_decomp% spec_iSize(1)])
   call self%Chebyshev_integration_d%backsolve( fDummyPtr,&
           1, domain_decomp% spec_iSize(1))
              
   call move_alloc(from=deriv, to=field)
   end do
   field(self%geometry%NZ : domain_decomp% spec_iSize(1)) = 0._dp
   end if
 end subroutine

 subroutine factorize_operators( self, dt_size, first_factor_bool )
   class(full_problem_data_structure_T), intent(inOut) :: self
   real(kind=dp), intent(in) :: dt_size
   logical, intent(in), optional :: first_factor_bool
   type(zcsr_matrix) :: aux_zCSR
   type(dcsr_matrix) :: aux_dCSR
   complex(kind=dp) :: zsca
   real(kind=dp) :: dsca
   integer :: i1, i2, iSys
   
   if (present(first_factor_bool)) then
   if (.not.first_factor_bool) then
      print *, 'By convention, first_factor_bool argument should be .True.'
      error stop
   else
      do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
         call self%coupled_kxky_set(iSys)%evol%set_size( &
                   self%coupled_kxky_set(iSys)%shape%spectral_local_NY, &
                   self%coupled_kxky_set(iSys)%shape%spectral_local_NX, &
                   self%coupled_kxky_set(iSys)%shape%total)
         call self%coupled_kxky_set(iSys)%evol%alloc_DIA_and_LU()
      end do
      if (self%geometry%this_core_has_zero_mode) then
      do iSys = 1, self%recipe%numberOf_coupled_zero_systems
         call self%coupled_zero_set(iSys)%evol%set_size( &
                   self%coupled_zero_set(iSys)%shape%total)
      end do
      end if
   end if
   end if
   ! in any case, we need to (re)factorize: 
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
      do i2 = 1, self%coupled_kxky_set(iSys)%evol%n2
      do i1 = 1, self%coupled_kxky_set(iSys)%evol%n1
      ! If the mode should be deAliased, the evolution operator is 
      ! set to the identity
      if ( self%geometry%deAliase_x (i2) .or. &
           self%geometry%deAliase_y (i1) ) then
           call aux_zCSR % identity( self%coupled_kxky_set(iSys)%shape%total )
      else 
         dsca = - dt_size * self%shandle%a_arr(1,1)    
         zsca = cmplx( dsca, 0._dp, kind=dp)
         call wrap_zcsradd(self%coupled_kxky_set(iSys)%mass%csr(i1,i2), &
                           self%coupled_kxky_set(iSys)%stif%csr(i1,i2), &
                           aux_zCSR, zsca)
      end if
      call aux_zCSR%to_dia( self%coupled_kxky_set(iSys)%evol%dia(i1,i2) ) 
      end do
      end do
      ! /////////////////////////////////////////////////////////
      ! Remember that we must set to zero the (kx,ky)=(0,0) modes
      ! ---------------------------------------------------------
      if (self%geometry%this_core_has_zero_mode) then
         call aux_zCSR % identity( self%coupled_kxky_set(iSys)%shape%total   )
         call aux_zCSR % to_dia  ( self%coupled_kxky_set(iSys)%evol%dia(1,1) )
      end if
      ! /////////////////////////////////////////////////////////

      call self%coupled_kxky_set(iSys)%evol%factorize()
   end do
 
   if (self%geometry%this_core_has_zero_mode) then
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
         dsca = - dt_size * self%shandle%a_arr(1,1)    
         call wrap_dcsradd(self%coupled_zero_set(iSys)%mass%csr, &
                           self%coupled_zero_set(iSys)%stif%csr, &
                           aux_dCSR, dsca)
         call aux_dCSR%to_dia( self%coupled_zero_set(iSys)%evol%dia ) 
         call self%coupled_zero_set(iSys)%evol%factorize()
   end do
   end if
 end subroutine 


 subroutine noisy_initial_conditions(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iSys, iVar
   integer :: ix, iy, iz
   integer :: nMode
   complex(kind=dp) :: zPhase
   real(kind=dp) :: random_dScalar
   real(kind=dp) :: amp
   real(kind=dp) :: amp0
   amp = self%cargo%initialCondition_amp_kxky
   amp0= self%cargo%initialCondition_amp_zero
   ! => linearly coupled systems, kxky modes 
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems 
     do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
     do ix = 1, self%coupled_kxky_set(isys)%shape%spectral_local_NX
     do iy = 1, self%coupled_kxky_set(isys)%shape%spectral_local_NY
     do iz = self%coupled_kxky_set(isys)%shape%variable_firstIndex(iVar), &
             self%coupled_kxky_set(isys)%shape%variable_lastIndex (iVar)
     nMode = iz - self%coupled_kxky_set(isys)%shape%variable_firstIndex(iVar)
     call random_number(random_dScalar)
     zPhase = cmplx(cos(random_dScalar*2.*3.14159),&
                    sin(random_dScalar*2.*3.14159),&
                    kind=dp)
     self%coupled_kxky_set(isys)%field(iz, iy, ix) = &
                   exp(-(nMode-3.5_dp)**2/2._dp/4._dp) *&
                   exp(-(sqrt(self%geometry%kx(ix)**2 + &
                              self%geometry%ky(iy)**2)- & 
                              self%geometry%kx_box*6 )**2/2._dp/4._dp ) &
                   * zPhase * amp
     call random_number(random_dScalar)
     zPhase = cmplx(cos(random_dScalar*2.*3.14159),&
                    sin(random_dScalar*2.*3.14159),&
                    kind=dp)
     self%coupled_kxky_set(isys)%field(iz, iy, ix) = &
                   self%coupled_kxky_set(isys)%field(iz, iy, ix) + &
                   exp(-(nMode-1.5_dp)**2/2._dp/4._dp) *&
                   exp(-(sqrt(self%geometry%kx(ix)**2 + &
                              self%geometry%ky(iy)**2)- & 
                              self%geometry%kx_box*3 )**2/2._dp/16._dp) &
                   * zPhase * 3._dp * amp
      end do
      end do
      end do
      end do
   end do
   ! // linearly coupled systems, kxky modes 
   ! => linearly coupled systems, zero modes 
   if (self%geometry%this_core_has_zero_mode) then
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems 
     do iVar = 1, self%recipe%zero_recipes(iSys)%n_coupled_vars
     do iz = self%coupled_zero_set(isys)%shape%variable_firstIndex(iVar), &
             self%coupled_zero_set(isys)%shape%variable_lastIndex (iVar)
     nMode = iz - self%coupled_zero_set(isys)%shape%variable_firstIndex(iVar)
     self%coupled_zero_set(isys)%field(iz) = &
                   exp(-(nMode-1.5_dp)**2/2._dp/4._dp)*cos(iz*1._dp)  * amp0
      end do
      end do
   end do
   end if
   ! // linearly coupled systems, zero modes 

 end subroutine

 subroutine prepare_chebyshev_integration(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, parameter :: padding = 10
   type(dcsr_matrix) :: CEM_I, CEM_M
   type(dcsr_matrix) :: auxCSR
   type(zcsr_matrix) :: auxZcsr
   call chebyshev_elementary_matrices(self%geometry%NZ, &
                                      self%geometry%center, &
                                      self%geometry%gap, &
                                      CEM_I, CEM_M, padding)
   call dCsr_truncate( CEM_I,  &
                       auxCSR, &
                       2, self%geometry%NZ, 1, self%geometry%NZ-1)
   call dcsr_convert_zcsr( auxCSR, auxZcsr)                        
   call self%chebyshev_integration_z%build_zOp_DIA_1d_1coupled_from_csr(auxZcsr)
   call self%chebyshev_integration_z%factorize()
   call self%chebyshev_integration_d%build_dOp_DIA_1d_1coupled_from_csr(auxCSR)
   call self%chebyshev_integration_d%factorize()
 end subroutine
   
   
   
 subroutine initialise_the_data_structure(self, scheme_id)
   class(full_problem_data_structure_T), intent(inOut) :: self
   character(len=7), intent(in) :: scheme_id
   integer :: ix, num_args, isys
   character(len=32), dimension(:), allocatable :: args
   logical :: output_matrices = .False.
   character(len=255) :: cwd
   character(len=:), allocatable :: matrices_path
   integer :: matrices_path_len
   integer :: istatus
   call self%shandle%init(scheme_id, (my_rank.eq.0))
   call self%allocate_all_buffers()
   call self%prepare_building_tools()
   call self%prepare_stencils()
   call self%build_operators()
   call self%prepare_chebyshev_integration()
   ! now handle the command line arguments
   self%cargo% initialCondition_amp_kxky = 1.d-16
   self%cargo% initialCondition_amp_zero = 1.d-16
   self%cargo% cflFactor_along_z = 1.d0             
   self%cargo% smagorinsky_prefactor = 1.d0
   num_args = command_argument_count()
   allocate(args(num_args))
   do ix = 1, num_args
      call get_command_argument(ix,args(ix))
   end do
   do ix = 1, num_args
      select case (args(ix)) 
             case ('--output-matrices')
                  output_matrices = .True.
             case ('--noises-amplitude')
                  read (args(ix+1),*) self%cargo%initialCondition_amp_kxky 
                  read (args(ix+1),*) self%cargo%initialCondition_amp_zero 
             case ('--noise-amplitude-k')
                  read (args(ix+1),*) self%cargo%initialCondition_amp_kxky 
             case ('--noise-amplitude-0')
                  read (args(ix+1),*) self%cargo%initialCondition_amp_zero 
             case ('--cflFactor-along-z')
                  read (args(ix+1),*) self%cargo%cflFactor_along_z
             case ('--smagorinsky-prefactor')
                  read (args(ix+1),*) self%cargo% smagorinsky_prefactor
      end select
   end do
   if (output_matrices) then
      if (my_rank.eq.0) then
      print *, '*** Per user request (flag ''--output-matrices'' is present), matrices are written to the disk'
      call getcwd(cwd)
      matrices_path_len = len(trim(cwd)) + 13
      allocate (character(len=matrices_path_len) :: matrices_path)
      matrices_path = trim(cwd) // '/../Matrices/'
      do isys = 1, self%recipe%numberOf_coupled_kxky_systems
         call self%coupled_kxky_set(isys)%mass%csr(2,2)%write2disk( &
                'mass', 4, my_rank, matrices_path, matrices_path_len)
         call self%coupled_kxky_set(isys)%stif%csr(2,2)%write2disk( &
                'stif', 4, my_rank, matrices_path, matrices_path_len)
      end do
      call chdir(trim(cwd), istatus)
      end if
   end if
 end subroutine


 subroutine allocate_all_buffers(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar, iStep, iSys, iObj
   integer :: N_timeseries_objects
   ! timeseries buffers:
   N_timeseries_objects = 0
   do iVar = 1, self%recipe%numberOf_linear_variables_full
      do iObj = 1, self%recipe%timeseries%numberOf_linearObjects( iVar )
         N_timeseries_objects = N_timeseries_objects + 1
      end do
   end do
   do iVar = 1, self%recipe%numberOf_quadratic_variables
      do iObj = 1, self%recipe%timeseries%numberOf_quadraObjects( iVar )
         N_timeseries_objects = N_timeseries_objects + 1
      end do
   end do
   allocate( self% timeseries_buffer(self%buffer_length, N_timeseries_objects))
   ! other buffers
   allocate( self%coupled_kxky_set( &
             self%recipe%numberOf_coupled_kxky_systems ))
   allocate( self%coupled_zero_set( &
               self%recipe%numberOf_coupled_zero_systems ))
   call self%compute_systems_shapes()
   ! => linearly coupled systems, kxky modes 
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems 
     allocate( self%coupled_kxky_set(iSys)%field(                   &
               self%coupled_kxky_set(iSys)%shape%total,             &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NY, &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NX) )
     allocate( self%coupled_kxky_set(iSys)%rhs  (                   &
               self%coupled_kxky_set(iSys)%shape%total,             &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NY, &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NX) )
     allocate( self%coupled_kxky_set(iSys)%aux  (                   &
               self%coupled_kxky_set(iSys)%shape%total,             &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NY, &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NX) )
     allocate( self%coupled_kxky_set(iSys)%step( self%shandle%imp_stages ))
     do iStep = 1, self%shandle%imp_stages
        allocate( self%coupled_kxky_set(iSys)%step(iStep)%K_hat(    &
               self%coupled_kxky_set(iSys)%shape%total,             &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NY, &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NX) )
        allocate( self%coupled_kxky_set(iSys)%step(iStep)%K_std(    &
               self%coupled_kxky_set(iSys)%shape%total,             &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NY, &
               self%coupled_kxky_set(iSys)%shape%spectral_local_NX) )
     end do
   end do
   ! // linearly coupled systems, kxky modes 
   ! => linearly coupled systems, zero modes 
   if (self%geometry%this_core_has_zero_mode) then
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
     allocate( self%coupled_zero_set(iSys)%field( &
               self%coupled_zero_set(iSys)%shape%total))
     allocate( self%coupled_zero_set(iSys)%rhs  ( &
               self%coupled_zero_set(iSys)%shape%total))
     allocate( self%coupled_zero_set(iSys)%aux  ( &
               self%coupled_zero_set(iSys)%shape%total))
     allocate( self%coupled_zero_set(iSys)%step( self%shandle%imp_stages ))
     do iStep = 1, self%shandle%imp_stages
        allocate( self%coupled_zero_set(iSys)%step(iStep)%K_hat(&
                  self%coupled_zero_set(iSys)%shape%total ))
        allocate( self%coupled_zero_set(iSys)%step(iStep)%K_std(&
                  self%coupled_zero_set(iSys)%shape%total )) 
     end do
   end do
   end if
   ! // linearly coupled systems, zero modes 
   ! => linear (full) variables buffer allocation
   allocate (self%linear_variables( self%recipe%numberOf_linear_variables_full ))
   do iVar = 1, self%recipe%numberOf_linear_variables_full
      call self%linear_variables(iVar)%alloc()
   end do
   ! // linear (full) variables buffer allocation
   ! => quadratic variables buffer allocation
   allocate (self%quadratic_variables( self%recipe%numberOf_quadratic_variables ))
   do iVar = 1, self%recipe%numberOf_quadratic_variables
      call self%quadratic_variables(iVar)%alloc()
   end do
   ! // quadratic variables buffer allocation
 end subroutine


 include "PL_IMEX_timestepping_marching.f90"
 include 'PL_IMEX_timestepping_output.f90'

 subroutine prepare_stencils(self)
   class(full_problem_data_structure_T) :: self
   integer :: iSys
   integer :: iVar
   type(zcsr_matrix) :: zcsrMat
   type(zcsr_matrix) :: zcsrMat1
   type(dcsr_matrix) :: dcsrMat

   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
      allocate( self%coupled_kxky_set(iSys)%stencil( &
                     self%recipe%kxky_recipes(iSys)%n_coupled_vars ) )
      allocate( self%coupled_kxky_set(iSys)%square_stencil( &
                     self%recipe%kxky_recipes(iSys)%n_coupled_vars ) )
      allocate( self%coupled_kxky_set(iSys)%padStencilExtractShuffle( &
                     self%recipe%kxky_recipes(iSys)%n_coupled_vars ) )
      allocate( self%coupled_kxky_set(iSys)%shuffleTextractTtruncate( &
                     self%recipe%kxky_recipes(iSys)%n_coupled_vars ) )
      allocate( self%coupled_kxky_set(iSys)%shuffleInsert( &
                     self%recipe%kxky_recipes(iSys)%n_coupled_vars ) )
      do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
         call chebyshev_galerkin_unique_stencil( self%geometry%NZ, &
                      self%recipe%kxky_recipes(iSys)%vars(iVar)%bc_code, &
                      self%coupled_kxky_set(iSys)%stencil(iVar) )
         call self%coupled_kxky_set(iSys)%stencil(iVar)%sort()
         call dcsr_convert_zcsr( self%coupled_kxky_set(iSys)%stencil(iVar), zcsrMat1)
         call zcsrMat1 % truncate_rows (&
                       zcsrMat, &
                       1, self%coupled_kxky_set(iSys)%stencil(iVar)%ncol)
         call self%coupled_kxky_set(iSys)%square_stencil(iVar)%constructor(zcsrMat)
         call self%coupled_kxky_set(iSys)%square_stencil(iVar)%factorize
      end do
   end do

   if (self%geometry%this_core_has_zero_mode) then
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
      allocate( self%coupled_zero_set(iSys)%stencil( &
                     self%recipe%zero_recipes(iSys)%n_coupled_vars ) )
      allocate( self%coupled_zero_set(iSys)%square_stencil( &
                     self%recipe%zero_recipes(iSys)%n_coupled_vars ) )
      allocate( self%coupled_zero_set(iSys)%shuffleInsert( &
                     self%recipe%zero_recipes(iSys)%n_coupled_vars ) )
      allocate( self%coupled_zero_set(iSys)%padStencilExtractShuffle( &
                     self%recipe%zero_recipes(iSys)%n_coupled_vars ) )
      allocate( self%coupled_zero_set(iSys)%shuffleTextractTTruncate( &
                     self%recipe%zero_recipes(iSys)%n_coupled_vars ) )
      do iVar = 1, self%recipe%zero_recipes(iSys)%n_coupled_vars
         call chebyshev_galerkin_unique_stencil( self%geometry%NZ, &
                      self%recipe%zero_recipes(iSys)%vars(iVar)%bc_code, &
                      self%coupled_zero_set(iSys)%stencil(iVar) )
         call self%coupled_zero_set(iSys)%stencil(iVar)%sort()
         call self%coupled_zero_set(iSys)%stencil(iVar) % truncate_rows(&
                       dcsrMat, &                                                
                       1, self%coupled_zero_set(iSys)%stencil(iVar)%ncol)
         call self%coupled_zero_set(iSys)%square_stencil(iVar)%constructor(dcsrMat)
         call self%coupled_zero_set(iSys)%square_stencil(iVar)%factorize()
 
      end do
   end do
   end if
  
 end subroutine 

 subroutine prepare_building_tools(self)
   class(full_problem_data_structure_T) :: self
   integer :: iSys
   integer :: iOrder, order, iPiece
   integer :: zPowerMax, iPower
   integer :: padding = 15
   type(dcsr_matrix) :: CEM_I, CEM_M, auxCSR
   !type(dcsr_matrix), allocatable :: untruncated_cheb_I(:)

   call chebyshev_elementary_matrices(self%geometry%NZ, &
                                      self%geometry%center, &
                                      self%geometry%gap, &
                                      CEM_I, CEM_M, padding)
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
   order = maxval(self%recipe%kxky_recipes(iSys)%eqn_order)
   !/////////////////////////////////////////
   ! what is the maximum z-power?
   zPowerMax = 0
   do iPiece = 1, self%recipe%kxky_recipes(iSys)%mass%n_pieces
      zPowerMax= max0 (zPowerMax,&
                 self%recipe%kxky_recipes(iSys)%mass%z_multiply(iPiece) )
   end do
   do iPiece = 1, self%recipe%kxky_recipes(iSys)%stif%n_pieces
      zPowerMax= max0 (zPowerMax,&
                 self%recipe%kxky_recipes(iSys)%stif%z_multiply(iPiece) )
   end do
   !/////////////////////////////////////////

   allocate( self%coupled_kxky_set(iSys)%building_tools%Cheb_IM(order+1, zPowerMax+1) )
   do iOrder = 1, Order+1
   do iPower = 1, zPowerMax + 1
   call dCsr_expExp(CEM_I, iOrder-1, CEM_M, iPower-1, auxCSR)
   call dCsr_truncate( auxCSR, &
                       self%coupled_kxky_set(iSys)%building_tools%Cheb_IM(iOrder, iPower), &
                       1, self%geometry%NZ, 1, self%geometry%NZ)
   end do
   end do
   end do
   if (self%geometry%this_core_has_zero_mode) then
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
   order = maxval(self%recipe%zero_recipes(iSys)%eqn_order)
   !/////////////////////////////////////////
   ! what is the maximum z-power?
   zPowerMax = 0
   do iPiece = 1, self%recipe%zero_recipes(iSys)%mass%n_pieces
      zPowerMax= max0 (zPowerMax,&
                 self%recipe%zero_recipes(iSys)%mass%z_multiply(iPiece) )
   end do
   do iPiece = 1, self%recipe%zero_recipes(iSys)%stif%n_pieces
      zPowerMax= max0 (zPowerMax,&
                 self%recipe%zero_recipes(iSys)%stif%z_multiply(iPiece) )
   end do
   !/////////////////////////////////////////
   allocate( self%coupled_zero_set(iSys)%building_tools%Cheb_IM(order+1, zPowerMax+1) )
   do iOrder = 1, Order+1
   do iPower = 1, zPowerMax + 1
   call dCsr_expExp(CEM_I, iOrder-1, CEM_M, iPower-1, auxCSR)
   call dCsr_truncate( auxCSR, &
                       self%coupled_zero_set(iSys)%building_tools%Cheb_IM(iOrder, iPower), &
                       1, self%geometry%NZ, 1, self%geometry%NZ)
   end do
   end do
   end do
   end if
   
 end subroutine prepare_building_tools

 subroutine add_dsca_times_MIblock_to_zero_operator (self, which, dsca, &
                                                         iz_expo, mz_power,&
                                                         ieqn, ivar, eq_order)
   class(coupled_set_1d_ops_vars_T), intent(inOut) :: self
   character(len=4), intent(in) :: which
   real(kind=dp) :: dsca
   integer, intent(in) :: iz_expo, ieqn, ivar, eq_order, mz_power
   type(dcsr_matrix) :: block, aux2
   call dcsr_truncate(self%building_tools%cheb_IM( iz_expo+1, mz_power+1 ), block, &
                                   eq_order+1, self%building_tools%cheb_IM(1,1)%nrow, &
                                   1,          self%building_tools%cheb_IM(1,1)%ncol)
   call wrap_dcsrmultcsr(block, self%stencil(ivar), aux2)
   !call dcsr_truncate(aux1, aux2, eq_order+1, aux1%nrow, 1, aux1%ncol)
   select case (which)
          case ('mass')
          call dcsr_add_block_dcsr_d(self%mass%csr, aux2, dsca, &
                                     self%shape%equation_firstIndex(iEqn), &
                                     self%shape%variable_firstIndex(iVar))
          case ('stif')
          call dcsr_add_block_dcsr_d(self%stif%csr, aux2, dsca, &
                                     self%shape%equation_firstIndex(iEqn), &
                                     self%shape%variable_firstIndex(iVar))
          case default
          stop 'bad character argument in add_zsca_times_MIblock...'
   end select
 end subroutine

 subroutine add_zsca_times_MIblock_to_kxky_operators (self, which, zsca, &
                                                         iz_expo, mz_power, &
                                                         ieqn, ivar, eq_order)
   class(coupled_set_3d_ops_vars_T), intent(inOut) :: self
   character(len=4), intent(in) :: which
   complex(kind=dp), allocatable :: zsca(:,:)
   integer, intent(in) :: iz_expo, ieqn, ivar, eq_order, mz_power
   type(dcsr_matrix) :: block, aux2
   integer :: i1, i2
   if (abs(iz_expo).gt.100) then
      print *, '/!\ /!\ routine add_zsca_times_MIblock_to_kxky_operators has been called with'
      print *, '/!\ /!\ unreallistic parameters. Problably an unintended bug.'
      print *, '/!\ /!\ other parameters [which, zsca(3,4), iz_expo, ieqn, ivar] are:'
      print *, which, zsca(3,4), iz_expo, ieqn, ivar
      stop
   end if
   call dcsr_truncate(self%building_tools%cheb_IM( iz_expo+1, mz_power+1 ), block, &
                                   eq_order+1, self%building_tools%cheb_IM(1,1)%nrow, &
                                   1,          self%building_tools%cheb_IM(1,1)%ncol)
   call wrap_dcsrmultcsr(block, self%stencil(ivar), aux2)
   do i2 = 1, self%mass%n2                                                         
   do i1 = 1, self%mass%n1                                                                
      select case (which)
             case ('mass')
             call zcsr_add_block_dcsr_z(self%mass%csr(i1,i2), aux2, zsca(i1,i2), &
                                        self%shape%equation_firstIndex(iEqn), &
                                        self%shape%variable_firstIndex(iVar))
             case ('stif')
             call zcsr_add_block_dcsr_z(self%stif%csr(i1,i2), aux2, zsca(i1,i2), &
                                        self%shape%equation_firstIndex(iEqn), &
                                        self%shape%variable_firstIndex(iVar))
             case default
             stop 'bad character argument in add_zsca_times_MIblock...'
      end select
   end do
   end do
 end subroutine

 subroutine compute_systems_shapes(self)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer :: iSys, iVar
   integer :: NZ
   integer :: total_dof1
   integer :: total_dof2
   integer :: integer_result
   integer :: nEqs, nVars
   NZ = self%geometry%NZ ! merely a shorthand notation
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
      Neqs = self%recipe%kxky_recipes(isys)%n_equations
      Nvars= self%recipe%kxky_recipes(isys)%n_coupled_vars
      if (Neqs.ne.NVars) Then              
         print *, '/!\ /!\ Incompatible number of coupled equations and variables'
         print *, '/!\ /!\ Violation for coupled kxky system #', iSys
         print *, Neqs, Nvars
         stop 
      else
         self%coupled_kxky_set(isys)%shape%spectral_local_NX = domain_decomp% spec_iSize(3)
         self%coupled_kxky_set(isys)%shape%spectral_local_NY = domain_decomp% spec_iSize(2)
         self%coupled_kxky_set(isys)%shape%physical_local_NX = domain_decomp% phys_iSize(3)
         self%coupled_kxky_set(isys)%shape%physical_local_NY = domain_decomp% phys_iSize(2)
         ! on one hand...
         total_dof1 = Neqs * NZ - &
                     sum(self%recipe%kxky_recipes(isys)%eqn_order)
         ! on the other hand
         total_dof2 = Nvars * NZ
         do iVar = 1, self%recipe%kxky_recipes(isys)%n_coupled_vars
           integer_result = self%recipe%kxky_recipes(isys)%vars(iVar)%bc_code / 10
           total_dof2 = total_dof2 - integer_result           
         end do
         if (total_dof1.ne.total_dof2) then
            print *, '/!\ /!\ Incompatible number of boundary conditions:', total_dof2
            print *, '/!\ /!\ Must match the order of the equations:', total_dof1       
            print *, '/!\ /!\ Violation for coupled kxky system #', iSys
            stop 
         else
            self%coupled_kxky_set(isys)%shape%total = total_dof1
            allocate ( self%coupled_kxky_set(isys)%shape%equation_firstIndex(nEqs) )
            allocate ( self%coupled_kxky_set(isys)%shape%equation_lastIndex (nEqs) )
            allocate ( self%coupled_kxky_set(isys)%shape%equation_size      (nEqs) )
            allocate ( self%coupled_kxky_set(isys)%shape%variable_firstIndex(nVars))
            allocate ( self%coupled_kxky_set(isys)%shape%variable_lastIndex (nVars))
            allocate ( self%coupled_kxky_set(isys)%shape%variable_size      (nVars))
            self%coupled_kxky_set(isys)%shape%equation_firstIndex (1) = 1
            self%coupled_kxky_set(isys)%shape%variable_firstIndex (1) = 1
            do iVar = 1, nVars-1
              self%coupled_kxky_set(isys)%shape%equation_size ( iVar ) = &
                   NZ - self%recipe%kxky_recipes(isys)%eqn_order ( iVar )
              self%coupled_kxky_set(isys)%shape%equation_lastIndex (iVar) = &
              self%coupled_kxky_set(isys)%shape%equation_firstIndex(iVar) + &
                   self%coupled_kxky_set(isys)%shape%equation_size( iVar) - 1 
              self%coupled_kxky_set(isys)%shape%equation_firstIndex(iVar+1) = &
              self%coupled_kxky_set(isys)%shape%equation_lastIndex (iVar) + 1 
              integer_result = self%recipe%kxky_recipes(isys)%vars (iVar)%bc_code / 10
              self%coupled_kxky_set(isys)%shape%variable_size (iVar) = &
                 NZ - integer_result
              self%coupled_kxky_set(isys)%shape%variable_lastIndex (iVar) = &
              self%coupled_kxky_set(isys)%shape%variable_firstIndex(iVar) - 1 &
                 + self%coupled_kxky_set(isys)%shape%variable_size (ivar) 
              self%coupled_kxky_set(isys)%shape%variable_firstIndex(iVar+1) = &
              self%coupled_kxky_set(isys)%shape%variable_lastIndex (iVar) + 1 
            end do
            self%coupled_kxky_set(isys)%shape%equation_size ( nVars) = &
                   NZ - self%recipe%kxky_recipes(isys)%eqn_order ( nVars)
            self%coupled_kxky_set(isys)%shape%equation_lastIndex(nVars)  = &
            self%coupled_kxky_set(isys)%shape%equation_firstIndex(nVars) + &
                   self%coupled_kxky_set(isys)%shape%equation_size(nVars)- 1 
            integer_result = self%recipe%kxky_recipes(isys)%vars (nVars)%bc_code / 10
            self%coupled_kxky_set(isys)%shape%variable_size (nvars) = &
               NZ - integer_result
            self%coupled_kxky_set(isys)%shape%variable_lastIndex (nvars) = &
            self%coupled_kxky_set(isys)%shape%variable_firstIndex(nvars) - 1 &
               + self%coupled_kxky_set(isys)%shape%variable_size (nvars) 
            ! a final check on the shape of the system:
            if (self%coupled_kxky_set(isys)%shape%total &
                             .ne.&
                self%coupled_kxky_set(isys)%shape%equation_lastIndex (nEqs) ) then
                print *, '/!\ /!\ The number of d.o.f. ', self%coupled_kxky_set(isys)%shape%total
                print *, '/!\ /!\ must match the computed last eqn index ', &
                        self%coupled_kxky_set(isys)%shape%equation_lastIndex (nEqs)
                print *, '/!\ /!\ Violation occured in kxky system #', iSys
                print *, '/!\ /!\ Computations details --- eqn order:',&
                        self%recipe%kxky_recipes(isys)%eqn_order
                print *, '/!\ /!\ Computations details --- eqn 1st index:',&
                        self%coupled_kxky_set(isys)%shape%equation_firstindex
                print *, '/!\ /!\ Computations details --- eqn last index:',&
                        self%coupled_kxky_set(isys)%shape%equation_lastindex
                stop 
            end if
            if (self%coupled_kxky_set(isys)%shape%total &
                             .ne.&
                self%coupled_kxky_set(isys)%shape%variable_lastIndex (nvars)) then
                print *, '/!\ /!\ The number of d.o.f. ', self%coupled_kxky_set(isys)%shape%total
                print *, '/!\ /!\ must match the computed last var index ', &
                        self%coupled_kxky_set(isys)%shape%variable_lastIndex (nvars)
                print *, '/!\ /!\ Violation occured in kxky system #', iSys
                !print *, '/!\ /!\ Computations details --- b.c. order:',&
                !         self%recipe%kxky_recipes(isys)%eqn_order
                print *, '/!\ /!\ Computations details --- var 1st index:',&
                        self%coupled_kxky_set(isys)%shape%variable_firstindex
                print *, '/!\ /!\ Computations details --- var last index:',&
                        self%coupled_kxky_set(isys)%shape%variable_lastindex
                stop 
            end if
         end if
      end if
   end do
      
   ! now we do the same for the zero-modes
   if (self%geometry%this_core_has_zero_mode) then
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
      Neqs = self%recipe%zero_recipes(isys)%n_equations
      Nvars= self%recipe%zero_recipes(isys)%n_coupled_vars
      if (Neqs.ne.NVars) Then              
         print *, '/!\ /!\ Incompatible number of coupled equations and variables'
         print *, '/!\ /!\ Violation for coupled zero system #', iSys
         stop 
      else
         ! on one hand...
         total_dof1 = Neqs * NZ - &
                     sum(self%recipe%zero_recipes(isys)%eqn_order)
         ! on the other hand
         total_dof2 = Nvars * NZ
         do iVar = 1, self%recipe%zero_recipes(isys)%n_coupled_vars
           integer_result = self%recipe%zero_recipes(isys)%vars(iVar)%bc_code / 10
           total_dof2 = total_dof2 - integer_result           
         end do
         if (total_dof1.ne.total_dof2) then
            print *, '/!\ /!\ Incompatible number of boundary conditions.'
            print *, '/!\ /!\ Must match the order of the equations.'       
            print *, '/!\ /!\ Violation for coupled zero system #', iSys
            stop 
         else
            self%coupled_zero_set(isys)%shape%total = total_dof1
            allocate ( self%coupled_zero_set(isys)%shape%equation_firstIndex(nEqs) )
            allocate ( self%coupled_zero_set(isys)%shape%equation_lastIndex (nEqs) )
            allocate ( self%coupled_zero_set(isys)%shape%equation_size      (nEqs) )
            allocate ( self%coupled_zero_set(isys)%shape%variable_firstIndex(nVars))
            allocate ( self%coupled_zero_set(isys)%shape%variable_lastIndex (nVars))
            allocate ( self%coupled_zero_set(isys)%shape%variable_size      (nVars))
            self%coupled_zero_set(isys)%shape%equation_firstIndex (1) = 1
            self%coupled_zero_set(isys)%shape%variable_firstIndex (1) = 1
            do iVar = 1, nVars-1
              self%coupled_zero_set(isys)%shape%equation_size ( iVar ) = &
                   NZ - self%recipe%zero_recipes(isys)%eqn_order ( iVar )
              self%coupled_zero_set(isys)%shape%equation_lastIndex (iVar) = &
              self%coupled_zero_set(isys)%shape%equation_firstIndex(iVar) + &
                   self%coupled_zero_set(isys)%shape%equation_size( iVar) - 1 
              self%coupled_zero_set(isys)%shape%equation_firstIndex(iVar+1) = &
              self%coupled_zero_set(isys)%shape%equation_lastIndex (iVar) + 1 
              integer_result = self%recipe%zero_recipes(isys)%vars (iVar)%bc_code / 10
              self%coupled_zero_set(isys)%shape%variable_size (iVar) = &
                 NZ - integer_result
              self%coupled_zero_set(isys)%shape%variable_lastIndex (iVar) = &
              self%coupled_zero_set(isys)%shape%variable_firstIndex(iVar) - 1 &
                 + self%coupled_zero_set(isys)%shape%variable_size (ivar) 
              self%coupled_zero_set(isys)%shape%variable_firstIndex(iVar+1) = &
              self%coupled_zero_set(isys)%shape%variable_lastIndex (iVar) + 1 
            end do
            self%coupled_zero_set(isys)%shape%equation_size ( nVars) = &
                   NZ - self%recipe%zero_recipes(isys)%eqn_order ( nVars)
            self%coupled_zero_set(isys)%shape%equation_lastIndex(nVars)  = &
            self%coupled_zero_set(isys)%shape%equation_firstIndex(nVars) + &
                   self%coupled_zero_set(isys)%shape%equation_size(nVars) - 1 
            integer_result = self%recipe%zero_recipes(isys)%vars (nVars)%bc_code / 10
            self%coupled_zero_set(isys)%shape%variable_size (nvars) = &
               NZ - integer_result
            self%coupled_zero_set(isys)%shape%variable_lastIndex (nvars) = &
            self%coupled_zero_set(isys)%shape%variable_firstIndex(nvars) - 1 &
               + self%coupled_zero_set(isys)%shape%variable_size (nvars) 
         end if
      end if
   end do
   end if
 end subroutine

 subroutine build_all_matrices(self)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer :: isys
   
   if (my_rank.eq.0) print *, '================================================================='
   if (my_rank.eq.0) print *, '>>> Assembling coupling matrices for (kx,ky)-modes systems.'          
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      call self%build_kxky_shuffle_matrices_for_system(isys)
      call self%kxky_matrices_manyGalkn_to_uniqueCheby(isys)
      call self%build_kxky_mass_matrices_for_system(isys)
      call self%build_kxky_stif_matrices_for_system(isys)
      call self%coupled_kxky_set(isys)%mass%change_of_basis_CSR(  &
           self%coupled_kxky_set(isys)%shuffleQ_var_to_polyDegree,&
           self%coupled_kxky_set(isys)%shuffleP_polyDegree_to_var) 
      call self%coupled_kxky_set(isys)%stif%change_of_basis_CSR(  &
           self%coupled_kxky_set(isys)%shuffleQ_var_to_polyDegree,&
           self%coupled_kxky_set(isys)%shuffleP_polyDegree_to_var) 
      ! the kxky=(0,0) is part of the operators and needs special care
      ! For this mode, we set M = 0 and L = id. 
      if (self%geometry%this_core_has_zero_mode) then
       call self%coupled_kxky_set(iSys)%mass%csr(1,1) % empty(&
                            self%coupled_kxky_set(iSys)%shape%total )
       call self%coupled_kxky_set(iSys)%stif%csr(1,1) % identity(&
                            self%coupled_kxky_set(iSys)%shape%total )
      end if
      ! at this point the change of basis is done. need to convert to DIA.
   end do
   if (my_rank.eq.0) print *, '/// Done with coupling matrices for (kx,ky)-modes systems.'
   if (my_rank.eq.0) print *, '================================================================='
   
   if (self%geometry%this_core_has_zero_mode) then
   print *, '================================================================='
   print *, '>>> Assembling coupling matrices for zero-mode systems.'          
   do isys = 1, self%recipe%numberOf_coupled_zero_systems
      call self%build_zero_shuffle_matrices_for_system(isys)
      call self%zero_matrices_manyGalkn_to_uniqueCheby(isys)
      call self%build_zero_mass_matrices_for_system(isys)
      call self%build_zero_stif_matrices_for_system(isys)
      call self%coupled_zero_set(isys)%mass%change_of_basis_CSR(  &
           self%coupled_zero_set(isys)%shuffleQ_var_to_polyDegree,&
           self%coupled_zero_set(isys)%shuffleP_polyDegree_to_var) 
      call self%coupled_zero_set(isys)%stif%change_of_basis_CSR(  &
           self%coupled_zero_set(isys)%shuffleQ_var_to_polyDegree,&
           self%coupled_zero_set(isys)%shuffleP_polyDegree_to_var) 
   end do
   print *, '/// Done with coupling matrices for zero-mode systems.'
   print *, '================================================================='
   end if

 end subroutine

 subroutine kxky_matrices_manyGalkn_to_uniqueCheby(self, isys)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iSys
   type(dcsr_matrix) ::    extractShuffle        !< some redundant auxilliary matrices
   type(dcsr_matrix) ::    shuffleInsert         !< some redundant auxilliary matrices
   type(dcsr_matrix) ::    StencilExtractShuffle !< that help a lot improve the legibility
   type(dcsr_matrix) :: padStencilExtractShuffle !< and understanding of the code.         
   type(dcsr_matrix) :: shuffleTextractTtruncate
   type(dcsr_matrix) :: shuffleTextractT_csr
   type(dcsr_matrix) :: extract_csr_T           
   type(dcsr_matrix) :: truncate_csr          
   type(dcoo_matrix) :: extract_coo_T           
   type(dcoo_matrix) :: truncate_coo          
   type(dcoo_matrix) :: extract_coo     
   type(dcsr_matrix) :: extract_csr
   type(dcoo_matrix) :: insert_coo     
   type(dcsr_matrix) :: insert_csr
   type(dcoo_matrix) :: padding_coo
   type(dcsr_matrix) :: padding_csr
   integer :: iVar, iElem
   ! => build the "padding matrix"   
   padding_coo%nelems = self%geometry%NZ                           
   padding_coo%nrow   = domain_decomp% spec_iSize(1)             
   padding_coo%ncol   = self%geometry%NZ
   allocate(padding_coo%row (padding_coo%nelems))
   allocate(padding_coo%col (padding_coo%nelems))
   allocate(padding_coo%dat (padding_coo%nelems))
   do iElem = 1, self%geometry%NZ       
      padding_coo%row(iElem) = iElem
      padding_coo%col(iElem) = iElem 
      padding_coo%dat(iElem) = 1._dp
   end do
   call d_coo2csr(padding_coo, padding_csr)
   ! // build the "padding matrix"
   do iVar = 1, self%recipe%kxky_recipes(isys)%n_coupled_vars
      ! => build the "truncation matrix"
      truncate_coo%nelems = self%coupled_kxky_set(isys)%shape%variable_size(iVar)
      truncate_coo%ncol   = domain_decomp% spec_iSize(1)             
      truncate_coo%nrow   = self%coupled_kxky_set(isys)%shape%variable_size(iVar)
      if (allocated( truncate_coo%row)) then
         deAllocate(truncate_coo%row)
         deAllocate(truncate_coo%col)
         deAllocate(truncate_coo%dat)
      end if
      allocate(truncate_coo%row (truncate_coo%nelems))
      allocate(truncate_coo%col (truncate_coo%nelems))
      allocate(truncate_coo%dat (truncate_coo%nelems))
      do iElem = 1, self%coupled_kxky_set(isys)%shape%variable_size(iVar)
         truncate_coo%row(iElem) = iElem
         truncate_coo%col(iElem) = iElem 
         truncate_coo%dat(iElem) = 1._dp
      end do
      call d_coo2csr(truncate_coo, truncate_csr)
   ! // build the "truncation matrix"
      !===========================================================
      !=> build the "extraction matrix" for variables
      extract_coo%nelems  = self%coupled_kxky_set(isys)%shape%variable_size(iVar)
      extract_coo%nrow    = self%coupled_kxky_set(isys)%shape%variable_size(iVar)
      extract_coo%ncol    = self%coupled_kxky_set(isys)%shape%total 
      if (allocated( extract_coo%row )) then
         deAllocate( extract_coo%col )
         deAllocate( extract_coo%dat )
         deAllocate( extract_coo%row )
      end if 
      allocate(extract_coo%row (extract_coo%nelems))
      allocate(extract_coo%col (extract_coo%nelems))
      allocate(extract_coo%dat (extract_coo%nelems))
      do iElem = 1, self%coupled_kxky_set(isys)%shape%variable_size(iVar)
         extract_coo%row(iElem) = iElem
         extract_coo%col(iElem) = iElem &
                 + self%coupled_kxky_set(isys)%shape%variable_firstIndex(iVar) -1
         extract_coo%dat(iElem) = 1._dp
      end do
      call d_coo2csr(extract_coo, extract_csr)
      call extract_coo%transpose( extract_coo_T)
      call d_coo2csr( extract_coo_T, extract_csr_T)
      !// build the "extraction matrix" for variables
      !===========================================================
      !=> build the "insertion matrix" for equations
      extract_coo%nelems  = self%coupled_kxky_set(isys)%shape%equation_size(iVar)
      extract_coo%nrow    = self%coupled_kxky_set(isys)%shape%equation_size(iVar)
      extract_coo%ncol    = self%coupled_kxky_set(isys)%shape%total 
      if (allocated( extract_coo%row )) then
         deAllocate( extract_coo%col )
         deAllocate( extract_coo%dat )
         deAllocate( extract_coo%row )
      end if 
      allocate(extract_coo%row (extract_coo%nelems))
      allocate(extract_coo%col (extract_coo%nelems))
      allocate(extract_coo%dat (extract_coo%nelems))
      do iElem = 1, self%coupled_kxky_set(isys)%shape%equation_size(iVar)
         extract_coo%row(iElem) = iElem
         extract_coo%col(iElem) = iElem &
                 + self%coupled_kxky_set(isys)%shape%equation_firstIndex(iVar) -1
         extract_coo%dat(iElem) = 1._dp
      end do
      call extract_coo%transpose(insert_coo)
      call d_coo2csr(insert_coo,  insert_csr)
      !// build the "insertion matrix" for equations
      !===========================================================

      call self%coupled_kxky_set(isys)%shuffleP_polyDegree_to_var_T%dot(&
           extract_csr_T,&
           shuffleTextractT_csr)
      call shuffleTextractT_csr%dot(&
           truncate_csr,&
           shuffleTextractTtruncate)
      call dcsr_convert_zcsr(           shuffleTextractTtruncate, &
            self%coupled_kxky_set(isys)%shuffleTextractTtruncate(iVar))
           
           


      call extract_csr%dot(&
           self%coupled_kxky_set(isys)%shuffleP_polyDegree_to_var, &
           ExtractShuffle)
      call self%coupled_kxky_set(isys)%stencil(iVar)%dot(&
                                       extractShuffle,   &
                                       stencilExtractShuffle)
      call padding_csr%dot(&
           stencilExtractShuffle,&
           padStencilExtractShuffle)                      
      call dcsr_convert_zcsr( padStencilExtractShuffle, &
            self%coupled_kxky_set(isys)%padStencilExtractShuffle(iVar))

      call self%coupled_kxky_set(isys)%shuffleQ_var_to_polyDegree%dot( &
                insert_csr, &
                shuffleInsert)
      call dcsr_convert_zcsr( shuffleInsert, &
            self%coupled_kxky_set(isys)%shuffleInsert(iVar))
      
   end do
 end subroutine

 subroutine zero_matrices_manyGalkn_to_uniqueCheby(self, isys)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iSys
   type(dcsr_matrix) ::    shuffleInsert         !< some redundant auxilliary matrices
   type(dcsr_matrix) ::    extractShuffle        !< some redundant auxilliary matrices
   type(dcsr_matrix) ::    stencilExtractShuffle !< some redundant auxilliary matrices
   type(dcsr_matrix) :: padStencilExtractShuffle !< that help a lot improve the legibility
   type(dcoo_matrix) :: extract_coo       !< and understanding of the code.           
   type(dcsr_matrix) :: extract_csr
   type(dcoo_matrix) :: insert_coo        !< and understanding of the code.           
   type(dcsr_matrix) :: insert_csr
   type(dcoo_matrix) :: padding_coo
   type(dcsr_matrix) :: padding_csr
   type(dcsr_matrix) :: shuffleTextractT_csr
   type(dcsr_matrix) :: extract_csr_T           
   type(dcsr_matrix) :: truncate_csr          
   type(dcoo_matrix) :: extract_coo_T           
   type(dcoo_matrix) :: truncate_coo          
   integer :: iVar, iElem
   ! => build the "padding matrix"   
   padding_coo%nelems = self%geometry%NZ                           
   padding_coo%nrow   = domain_decomp% spec_iSize(1)
   padding_coo%ncol   = self%geometry%NZ
   allocate(padding_coo%row (padding_coo%nelems))
   allocate(padding_coo%col (padding_coo%nelems))
   allocate(padding_coo%dat (padding_coo%nelems))
   do iElem = 1, self%geometry%NZ       
      padding_coo%row(iElem) = iElem
      padding_coo%col(iElem) = iElem 
      padding_coo%dat(iElem) = 1._dp
   end do
   call d_coo2csr(padding_coo, padding_csr)
   ! // build the "padding matrix"
   do iVar = 1, self%recipe%zero_recipes(isys)%n_coupled_vars
      ! => build the "truncation matrix"
      truncate_coo%nelems = self%coupled_zero_set(isys)%shape%variable_size(iVar)
      truncate_coo%ncol   = domain_decomp% spec_iSize(1)             
      truncate_coo%nrow   = self%coupled_zero_set(isys)%shape%variable_size(iVar)
      if (allocated( truncate_coo%row)) then
         deAllocate(truncate_coo%row)
         deAllocate(truncate_coo%col)
         deAllocate(truncate_coo%dat)
      end if
      allocate(truncate_coo%row (truncate_coo%nelems))
      allocate(truncate_coo%col (truncate_coo%nelems))
      allocate(truncate_coo%dat (truncate_coo%nelems))
      do iElem = 1, self%coupled_zero_set(isys)%shape%variable_size(iVar)
         truncate_coo%row(iElem) = iElem
         truncate_coo%col(iElem) = iElem 
         truncate_coo%dat(iElem) = 1._dp
      end do
      call d_coo2csr(truncate_coo, truncate_csr)
      !===========================================================
      !=> build the "extraction matrix" for variables
      extract_coo%nelems = self%coupled_zero_set(isys)%shape%variable_size(iVar)
      extract_coo%nrow   = self%coupled_zero_set(isys)%shape%variable_size(iVar)
      extract_coo%ncol   = self%coupled_zero_set(isys)%shape%total 
      if (allocated( extract_coo%row )) then
         deAllocate( extract_coo%col )
         deAllocate( extract_coo%dat )
         deAllocate( extract_coo%row )
      end if 
      allocate(extract_coo%row (extract_coo%nelems))
      allocate(extract_coo%col (extract_coo%nelems))
      allocate(extract_coo%dat (extract_coo%nelems))
      do iElem = 1, self%coupled_zero_set(isys)%shape%variable_size(iVar)
         extract_coo%row(iElem) = iElem
         extract_coo%col(iElem) = iElem + self%coupled_zero_set(isys)%shape%variable_firstIndex(iVar) -1
         extract_coo%dat(iElem) = 1._dp
      end do
      call d_coo2csr(extract_coo, extract_csr)
      call extract_coo%transpose( extract_coo_T)
      call d_coo2csr( extract_coo_T, extract_csr_T)
      ! // build the "extraction matrix"
      !===========================================================
      !=> build the "insertion matrix" for equations
      extract_coo%nelems = self%coupled_zero_set(isys)%shape%equation_size(iVar)
      extract_coo%nrow   = self%coupled_zero_set(isys)%shape%equation_size(iVar)
      extract_coo%ncol   = self%coupled_zero_set(isys)%shape%total 
      if (allocated( extract_coo%row )) then
         deAllocate( extract_coo%col )
         deAllocate( extract_coo%dat )
         deAllocate( extract_coo%row )
      end if 
      allocate(extract_coo%row (extract_coo%nelems))
      allocate(extract_coo%col (extract_coo%nelems))
      allocate(extract_coo%dat (extract_coo%nelems))
      do iElem = 1, self%coupled_zero_set(isys)%shape%equation_size(iVar)
         extract_coo%row(iElem) = iElem
         extract_coo%col(iElem) = iElem + self%coupled_zero_set(isys)%shape%equation_firstIndex(iVar) -1
         extract_coo%dat(iElem) = 1._dp
      end do
      call extract_coo%transpose(insert_coo)
      call d_coo2csr(insert_coo, insert_csr)
      !// build the "insertion matrix" for equations
      !===========================================================

      call self%coupled_zero_set(isys)%shuffleP_polyDegree_to_var_T%dot(&
           extract_csr_T,&
           shuffleTextractT_csr)
      call shuffleTextractT_csr%dot(&
           truncate_csr,&
            self%coupled_zero_set(isys)%shuffleTextractTtruncate(iVar))
           
      call extract_csr%dot(&
           self%coupled_zero_set(isys)%shuffleP_polyDegree_to_var, &
           ExtractShuffle)
      call self%coupled_zero_set(isys)%stencil(iVar)%dot(&
                                       extractShuffle,   &
                                       stencilExtractShuffle)
      call padding_csr%dot(&
           stencilExtractShuffle,&
           padStencilExtractShuffle)                      
      call dcsr_copy( padStencilExtractShuffle, &
            self%coupled_zero_set(isys)%padStencilExtractShuffle(iVar))
      ! for the nonlinearity
      call self%coupled_zero_set(isys)%shuffleQ_var_to_polyDegree%dot( &
                insert_csr, &
                shuffleInsert)
      call dcsr_copy ( shuffleInsert, &
            self%coupled_zero_set(isys)%shuffleInsert(iVar))
   end do
 end subroutine

 subroutine build_kxky_shuffle_matrices_for_system(self, isys)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: isys
   type(dcoo_matrix) :: ACoo
   integer :: iCounter, nCoefsMin
   integer :: iCoef, iRemainder, iVar
   integer :: offset
   aCoo%nElems = self%coupled_kxky_set(isys)%shape%total
   aCoo%nRow   = self%coupled_kxky_set(isys)%shape%total
   aCoo%nCol   = self%coupled_kxky_set(isys)%shape%total
   allocate( aCoo%dat (aCoo%nElems) )
   allocate( aCoo%row (aCoo%nElems) )
   allocate( aCoo%col (aCoo%nElems) )
   nCoefsMin = minval(self%coupled_kxky_set(isys)%shape%variable_size)
   iCounter = 0
   offset   = 0
   do iVar = 1, self%recipe%kxky_recipes(isys)%n_coupled_vars
   do iCoef = 1, nCoefsMin
      iCounter = iCounter + 1
      aCoo%row(iCounter) = iCounter
      aCoo%dat(iCounter) = 1._dp
      aCoo%col(iCounter) = iVar + (iCoef-1) * self%recipe%kxky_recipes(iSys)%n_coupled_vars
   end do
   do iRemainder = 1, self%coupled_kxky_set(isys)%shape%variable_size(iVar) - nCoefsMin
      iCounter = iCounter + 1
      offset = offset + 1
      aCoo%row(icounter) = iCounter
      aCoo%dat(icounter) = 1._dp
      aCoo%col(icounter) = nCoefsMin * self%recipe%kxky_recipes(isys)%n_coupled_vars + offset
   end do
   end do
   
   call d_coo2csr( aCoo, self%coupled_kxky_set(isys)%shuffleP_polyDegree_to_var )
   call aCoo%transpose()
   call d_coo2csr( aCoo, self%coupled_kxky_set(isys)%shuffleP_polyDegree_to_var_T) 

   !///////// code below need to be checked ///////////////////////
   ! do the same with equations
   aCoo%nElems = self%coupled_kxky_set(isys)%shape%total
   aCoo%nRow   = self%coupled_kxky_set(isys)%shape%total
   aCoo%nCol   = self%coupled_kxky_set(isys)%shape%total
   deAllocate( aCoo%dat )
   deAllocate( aCoo%row ) 
   deAllocate( aCoo%col )
   allocate( aCoo%dat (aCoo%nElems) )
   allocate( aCoo%row (aCoo%nElems) )
   allocate( aCoo%col (aCoo%nElems) )
   nCoefsMin = minval(self%coupled_kxky_set(isys)%shape%equation_size)
   icounter = 0
   offset   = 0
   do iVar = 1, self%recipe%kxky_recipes(isys)%n_coupled_vars
   do iCoef = 1, nCoefsMin
      icounter = icounter + 1
      aCoo%row(icounter) = iCounter
      aCoo%dat(icounter) = 1._dp
      aCoo%col(iCounter) = iVar + (iCoef-1) * self%recipe%kxky_recipes(iSys)%n_coupled_vars
   end do
   do iRemainder = 1, self%coupled_kxky_set(isys)%shape%equation_size(iVar) - nCoefsMin
      icounter = icounter + 1
      offset = offset + 1
      aCoo%row(icounter) = iCounter
      aCoo%dat(icounter) = 1._dp
      aCoo%col(icounter) = nCoefsMin * self%recipe%kxky_recipes(isys)%n_coupled_vars + offset
   end do
   end do
   
   call d_coo2csr( aCoo, self%coupled_kxky_set(isys)%shuffleQ_var_to_polyDegree_T)
   call aCoo%transpose()
   call d_coo2csr( aCoo, self%coupled_kxky_set(isys)%shuffleQ_var_to_polyDegree )
 
 end subroutine

 subroutine build_zero_shuffle_matrices_for_system(self, isys)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: isys
   type(dcoo_matrix) :: ACoo
   integer :: iCounter, nCoefsMin
   integer :: iCoef, iRemainder, iVar
   integer :: offset
   aCoo%nElems = self%coupled_zero_set(isys)%shape%total
   aCoo%nRow   = self%coupled_zero_set(isys)%shape%total
   aCoo%nCol   = self%coupled_zero_set(isys)%shape%total
   allocate( aCoo%dat (aCoo%nElems) )
   allocate( aCoo%row (aCoo%nElems) )
   allocate( aCoo%col (aCoo%nElems) )
   nCoefsMin = minval(self%coupled_zero_set(isys)%shape%variable_size)
   icounter = 0
   offset   = 0
   do iVar = 1, self%recipe%zero_recipes(isys)%n_coupled_vars
   do iCoef = 1, nCoefsMin
      icounter = icounter + 1
      aCoo%row(icounter) = iCounter
      aCoo%dat(icounter) = 1._dp
      aCoo%col(iCounter) = iVar + (iCoef-1) * self%recipe%zero_recipes(iSys)%n_coupled_vars
      !aCoo%col(icounter) = self%coupled_zero_set(isys)%shape%variable_firstIndex(iVar) - 1 + iCoef
   end do
   do iRemainder = 1, self%coupled_zero_set(isys)%shape%variable_size(iVar) - nCoefsMin
      icounter = icounter + 1
      offset = offset + 1
      aCoo%row(icounter) = iCounter
      aCoo%dat(icounter) = 1._dp
      aCoo%col(icounter) = nCoefsMin * self%recipe%zero_recipes(isys)%n_coupled_vars + offset
   end do
   end do
   
   call d_coo2csr( aCoo, self%coupled_zero_set(isys)%shuffleP_polyDegree_to_var   )
   call aCoo%transpose()
   call d_coo2csr( aCoo, self%coupled_zero_set(isys)%shuffleP_polyDegree_to_var_T )

   !///////// code below need to be checked ///////////////////////
   ! do the same with equations
   aCoo%nElems = self%coupled_zero_set(isys)%shape%total
   aCoo%nRow   = self%coupled_zero_set(isys)%shape%total
   aCoo%nCol   = self%coupled_zero_set(isys)%shape%total
   deAllocate( aCoo%dat )
   deAllocate( aCoo%row ) 
   deAllocate( aCoo%col )
   allocate( aCoo%dat (aCoo%nElems) )
   allocate( aCoo%row (aCoo%nElems) )
   allocate( aCoo%col (aCoo%nElems) )
   nCoefsMin = minval(self%coupled_zero_set(isys)%shape%equation_size)
   icounter = 0
   offset   = 0
   do iVar = 1, self%recipe%zero_recipes(isys)%n_coupled_vars
   do iCoef = 1, nCoefsMin
      icounter = icounter + 1
      aCoo%row(icounter) = iCounter
      aCoo%dat(icounter) = 1._dp
      !aCoo%col(icounter) = self%coupled_zero_set(isys)%shape%equation_firstIndex(iVar) - 1 + iCoef
      aCoo%col(iCounter) = iVar + (iCoef-1) * self%recipe%zero_recipes(iSys)%n_coupled_vars
   end do
   do iRemainder = 1, self%coupled_zero_set(isys)%shape%equation_size(iVar) - nCoefsMin
      icounter = icounter + 1
      offset = offset + 1
      aCoo%row(icounter) = iCounter
      aCoo%dat(icounter) = 1._dp
      aCoo%col(icounter) = nCoefsMin * self%recipe%zero_recipes(isys)%n_coupled_vars + offset
   end do
   end do
   
   call d_coo2csr( aCoo, self%coupled_zero_set(isys)%shuffleQ_var_to_polyDegree_T)
   call aCoo%transpose()
   call d_coo2csr( aCoo, self%coupled_zero_set(isys)%shuffleQ_var_to_polyDegree )
 
 end subroutine

 subroutine build_kxky_mass_matrices_for_system(self, isys)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer, intent(in) :: isys
   type(linear_operator_recipe_T), pointer :: my_recipe
   complex(kind=dp), allocatable :: numerical_prefactors_array(:,:)
   integer :: iPiece
   logical :: verbose = .False.
   if (my_rank.eq.0) verbose = .True.
   1987 Format (' >>> --- Mass matrix, system #',I2.2)
   if (verbose) write (*,1987) iSys
   1988 Format (' >>> --- ... coupling for var. #',I2.2,' in eqn. #',I2.2)
   my_recipe => self%recipe%kxky_recipes(isys)%mass 
   ! set the size of the zOperator_3d_1coupled_T object:
   call self%coupled_kxky_set(isys)%mass%set_size(  &
                     self%coupled_kxky_set(isys)%shape%spectral_local_NY, &
                     self%coupled_kxky_set(isys)%shape%spectral_local_NX, &
                     self%coupled_kxky_set(isys)%shape%total)
   call self%coupled_kxky_set(isys)%mass%init_CSR()
   do iPiece = 1, my_recipe%n_pieces
     if (verbose) write (*,1988) my_recipe%jVar(iPiece), my_recipe%iEqn(iPiece)
     ! check that the recipe is sensible before calling building function
     ! compute the numerical prefactor
     numerical_prefactors_array = recipe_to_prefactors_array(my_recipe, iPiece, self%geometry)
     call self%coupled_kxky_set(isys)%add_block(&
                                'mass', &
                                numerical_prefactors_array, &
                                my_recipe%Iz_exponent(iPiece), &
                                my_recipe%z_multiply (iPiece), &
                                my_recipe%iEqn(iPiece), &
                                my_recipe%jVar(iPiece), &
                                self%recipe%kxky_recipes(isys)%eqn_order( my_recipe%iEqn(iPiece) ))
   end do
 end subroutine

 subroutine build_kxky_stif_matrices_for_system(self, isys)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer, intent(in) :: isys
   type(linear_operator_recipe_T), pointer :: my_recipe
   complex(kind=dp), allocatable :: numerical_prefactors_array(:,:)
   integer :: iPiece
   logical :: verbose = .False.
   if (my_rank.eq.0) verbose = .True.
   1989 Format (' >>> --- Stiffness matrix, system #',I2.2)
   if (verbose) write (*,1989) iSys
   1990 Format (' >>> --- ... coupling for var. #',I2.2,' in eqn. #',I2.2)
   my_recipe => self%recipe%kxky_recipes(isys)%stif
   ! set the size of the zOperator_3d_1coupled_T object:
   call self%coupled_kxky_set(isys)%stif%set_size(  &
                     self%coupled_kxky_set(isys)%shape%spectral_local_NY, &
                     self%coupled_kxky_set(isys)%shape%spectral_local_NX, &
                     self%coupled_kxky_set(isys)%shape%total)
   call self%coupled_kxky_set(isys)%stif%init_CSR()
   do iPiece = 1, my_recipe%n_pieces
     if (verbose) write (*,1990) my_recipe%jVar(iPiece), my_recipe%iEqn(iPiece)
     ! compute the numerical prefactor
     numerical_prefactors_array = recipe_to_prefactors_array(my_recipe, iPiece, self%geometry)
     call self%coupled_kxky_set(isys)%add_block(&
                                'stif', &
                                numerical_prefactors_array, &
                                my_recipe%Iz_exponent(iPiece), &
                                my_recipe%z_multiply (iPiece), &
                                my_recipe%iEqn(iPiece), &
                                my_recipe%jVar(iPiece), &
                                self%recipe%kxky_recipes(isys)%eqn_order( my_recipe%iEqn(iPiece) ))
   end do
 end subroutine


pure function recipe_to_prefactors_array(a_recipe, iPiece, geo)
   integer, intent(in) :: iPiece
   Type(linear_operator_recipe_T), pointer, intent(in) :: a_recipe
   complex(kind=dp), allocatable :: recipe_to_prefactors_array (:,:)
   complex(kind=dp), allocatable :: sh(:,:)
   type(geometry_vars_T), intent(in) :: geo
   integer :: i1, i2
   allocate( sh (geo%spec%local_NY, geo%spec%local_NX) )
   do i2 = 1, geo%spec%local_NX ! X is the slow index in spectral!
   do i1 = 1, geo%spec%local_NY
      sh(i1,i2) = a_recipe%dsca(iPiece) &
                * ((geo%py(i1))**a_recipe%ky_exponent(iPiece))&
                * ((geo%px(i2))**a_recipe%kx_exponent(iPiece))
   end do
   end do
   call move_alloc(from=sh, to=recipe_to_prefactors_array)
 end function
   

 subroutine build_zero_stif_matrices_for_system(self, isys)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer, intent(in) :: isys
   type(linear_operator_recipe_T), pointer :: my_recipe
   real(kind=dp) :: numerical_prefactor
   integer :: iPiece
   logical :: verbose = .False.
   if (my_rank.eq.0) verbose = .True.
   1991 Format (' >>> --- Stiffness matrix, system #',I2.2)
   if (verbose) write (*,1991) iSys
   1992 Format (' >>> --- ... coupling for var. #',I2.2,' in eqn. #',I2.2)
   my_recipe => self%recipe%zero_recipes(isys)%stif
   call self%coupled_zero_set(isys)%stif%set_size(  &
                     self%coupled_zero_set(isys)%shape%total)
   call self%coupled_zero_set(isys)%stif%init_CSR()
   do iPiece = 1, my_recipe%n_pieces
     if (verbose) write (*,1992) my_recipe%jVar(iPiece), my_recipe%iEqn(iPiece)
     ! compute the numerical prefactor
     numerical_prefactor = my_recipe%dsca(iPiece) 
     call self%coupled_zero_set(isys)%add_block(&
                                'stif', &
                                numerical_prefactor, &
                                my_recipe%Iz_exponent(iPiece), &
                                my_recipe%z_multiply (iPiece), &
                                my_recipe%iEqn(iPiece), &
                                my_recipe%jVar(iPiece), &
                                self%recipe%zero_recipes(isys)%eqn_order( my_recipe%iEqn(iPiece) ))
   end do
 end subroutine
   
 subroutine build_zero_mass_matrices_for_system(self, isys)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer, intent(in) :: isys
   type(linear_operator_recipe_T), pointer :: my_recipe
   real(kind=dp) :: numerical_prefactor
   integer :: iPiece
   logical :: verbose = .False.
   if (my_rank.eq.0) verbose = .True.
   1993 Format (' >>> --- Mass matrix, system #',I2.2)
   if (verbose) write (*,1993) iSys
   1994 Format (' >>> --- ... coupling for var. #',I2.2,' in eqn. #',I2.2)
   my_recipe => self%recipe%zero_recipes(isys)%mass
   call self%coupled_zero_set(isys)%mass%set_size(  &
                     self%coupled_zero_set(isys)%shape%total)
   call self%coupled_zero_set(isys)%mass%init_CSR()
   do iPiece = 1, my_recipe%n_pieces
     if (verbose) write (*,1994) my_recipe%jVar(iPiece), my_recipe%iEqn(iPiece)
     ! compute the numerical prefactor
     numerical_prefactor = my_recipe%dsca(iPiece) 
     call self%coupled_zero_set(isys)%add_block(&
                                'mass', &
                                numerical_prefactor, &
                                my_recipe%Iz_exponent(iPiece), &
                                my_recipe%z_multiply (iPiece), &
                                my_recipe%iEqn(iPiece), &
                                my_recipe%jVar(iPiece), &
                                self%recipe%zero_recipes(isys)%eqn_order( my_recipe%iEqn(iPiece) ))
   end do
 end subroutine
 
end module PL_IMEX_timestepping
