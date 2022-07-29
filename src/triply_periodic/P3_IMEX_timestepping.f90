 !============================================================================= 
 !                            C O R A L 
 !============================================================================= 
 ! 
 ! MODULE: P3_IMEX_timestepping
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


module P3_IMEX_timestepping


 use P3_algebra
 use P3_equations
 use P3_geometry
 use output_misc
 use IMEX_schemes
 use P3_transforms
 use decomp_2d_io
 use timeKeeping, only: timings
 implicit none

 type :: intermediate_step_buffer_3d_T !< intermediate 3D fields (cf ARS97)
   complex(kind=dp), allocatable :: K_hat(:,:,:,:) !< eq.(2.1a)
   complex(kind=dp), allocatable :: K_std(:,:,:,:) !< eq.(2.1b)
 end type intermediate_step_buffer_3d_T

 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !type :: intermediate_step_buffer_1d_T !< intermediate 1D field (cf ARS97)
 !  complex(kind=dp), allocatable :: K_hat(:,:) !< eq.(2.1a)
 !  complex(kind=dp), allocatable :: K_std(:,:) !< eq.(2.1b)
 !end type intermediate_step_buffer_1d_T
 !///////////////////// ZERO MODE ////////////////////////////////////////////////


 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !type :: coupled_set_1d_ops_vars_T
 !  complex(kind=dp), allocatable :: field(:,:)
 !  complex(kind=dp), allocatable :: rhs  (:,:)
 !  complex(kind=dp), allocatable :: aux  (:,:)
 !  type(operator_1d_T) :: mass
 !  type(operator_1d_T) :: stif
 !  type(operator_1d_T) :: evol
 !  ! type(sest)
 !  type(intermediate_step_buffer_1d_T), allocatable :: step(:)
 !contains
 !  procedure :: add_block => add_zsca_times_MIblock_to_kxky_operators ! WRITEME
 !end type coupled_set_1d_ops_vars_T
 !///////////////////// ZERO MODE ////////////////////////////////////////////////

 type :: coupled_set_3d_ops_vars_T
   complex(kind=dp), allocatable :: field(:,:,:,:)
   complex(kind=dp), allocatable :: rhs  (:,:,:,:)
   complex(kind=dp), allocatable :: aux  (:,:,:,:)
   type(Operator_3d_T) :: mass
   type(Operator_3d_T) :: stif
   type(Operator_3d_T) :: evol
   ! type(sest)
   type(intermediate_step_buffer_3d_T), allocatable :: step(:)
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
 

 type :: cargo_T
   real(kind=dp) :: KE_display
   real(kind=dp) :: cfl_based_DT
   real(kind=dp) :: cflFactor_along_z
   real(kind=dp) :: initialCondition_amp_kxky
   real(kind=dp) :: initialCondition_amp_zero
 end type

 type :: full_problem_data_structure_T
   type(full_problem_recipe_T) :: recipe
   type(coupled_set_3d_ops_vars_T), allocatable :: coupled_kxky_set(:)
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !  type(coupled_set_1d_ops_vars_T), allocatable :: coupled_zero_set(:)
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
   type(scheme_handle) :: shandle
   type(geometry_vars_t) :: geometry
   type(PS_fields_T), allocatable :: quadratic_variables(:)
   type(PS_fields_T), allocatable :: linear_variables(:)
   type(output_bookkeeping_T) :: io_bookkeeping
   type(cargo_T) :: cargo
  contains 
   procedure :: add_K_std_to_rhs   => LPIMEX_add_K_std_to_rhs
   procedure :: add_K_hat_to_rhs   => LPIMEX_add_K_hat_to_rhs
   procedure :: allocate_all_buffers
   procedure :: backsolve_to_aux   => LPIMEX_backsolve_to_aux
   procedure :: build_operators    => build_all_matrices
   procedure :: build_kxky_mass_matrices_for_system     
   procedure :: build_kxky_stif_matrices_for_system     
   procedure :: compute_full_variables_in_physical_space
   procedure :: compute_NL_terms_to_K_hat => LPIMEX_compute_NL_terms_to_K_hat
   procedure :: compute_cfl_based_timestep
   procedure :: copy_fields_to_aux => LPIMEX_copy_fields_to_aux
   procedure :: copy_aux_to_fields => LPIMEX_copy_aux_to_fields
   procedure :: dealiase_the_rhs
   procedure :: export_allPhys
   procedure :: export_allProfiles
   procedure :: export_slice
   procedure :: export_verticallyAvgedSlice
   procedure :: export_volume
   procedure :: export_Profile
   procedure :: export_CheckPoints
   procedure :: invert_operators 
   procedure :: import_quickSave
   procedure :: init => initialise_the_data_structure
   procedure :: initialise_bookkeeping_counters
   procedure :: march_forward => one_step_beyond
   procedure :: mass_field_to_rhs  => LPIMEX_mass_field_to_rhs
   procedure :: noisy_initial_conditions
   procedure :: output_global_quantities
   procedure :: output_slices_volumes_and_profiles
   procedure :: stif_aux_to_K_std  => LPIMEX_stif_aux_to_K_std
   procedure :: horizontal_average_of_physical_quantity_inplace
 end type full_problem_data_structure_T

 contains

 subroutine invert_operators( self, dt_size, first_factor_bool )
   class(full_problem_data_structure_T), intent(inOut) :: self
   real(kind=dp), intent(in) :: dt_size
   logical, intent(in), optional :: first_factor_bool
   complex(kind=dp) :: zsca
   real   (kind=dp) :: dsca
   integer :: i1, i2, iSys
   
   if (present(first_factor_bool)) then
   if (.not.first_factor_bool) then
      print *, 'By convention, first_factor_bool argument should be .True.'
      error stop
   else
      do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
         call self%coupled_kxky_set(iSys)%evol%initialise( &
                   self% recipe% kxky_recipes(iSys)% n_coupled_vars, &
                   self% geometry% NZAA,                            &
                   self% geometry% spec% local_NY,                  &
                   self% geometry% spec% local_NX)
      end do
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !  some initialisation needed here
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
   end if
   end if
   ! in any case, we need to (re)build and invert:    
   dsca = - dt_size * self%shandle%a_arr(1,1)    
   zsca = cmplx( dsca, 0._dp, kind=dp)
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
      call self% coupled_kxky_set(iSys)% evol% & !<   evol ...
                     equals_sum_of(             & !<     =  ...
           self% coupled_kxky_set(iSys)% mass, & !<   mass ...
           self% coupled_kxky_set(iSys)% stif, & !< + stif ...
                     zsca)                       !< * dt*(scheme-dependent coef.)
      call self% coupled_kxky_set(iSys)% evol% invert() 
      !delme ! /pouet
      !delme call self% coupled_kxky_set(iSys)% evol% write2disk('evol')
      !delme ! pouet/
   end do
 
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !  zero-mode counterpart here
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
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
     do ix = 1, self%geometry% spec% local_NX
     do iy = 1, self%geometry% spec% local_NY
     do iz = 1, self%geometry% NZAA
     nMode = iz 
     call random_number(random_dScalar)
     zPhase = cmplx(cos(random_dScalar*2.*3.14159),&
                    sin(random_dScalar*2.*3.14159),&
                    kind=dp)
     self%coupled_kxky_set(isys)%field(iz, iy, ix, iVar) = &
                   exp(-(nMode-3.5_dp)**2/2._dp/4._dp) *&
                   exp(-(sqrt(self%geometry%kx(ix)**2 + &
                              self%geometry%ky(iy)**2)- & 
                              self%geometry%kx_box*6 )**2/2._dp/4._dp ) &
                   * zPhase * amp
     call random_number(random_dScalar)
     zPhase = cmplx(cos(random_dScalar*2.*3.14159),&
                    sin(random_dScalar*2.*3.14159),&
                    kind=dp)
     self%coupled_kxky_set(isys)%field(iz, iy, ix, iVar) = &
                   self%coupled_kxky_set(isys)%field(iz, iy, ix, iVar) + &
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
 
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !  zero-mode counterpart here
 !///////////////////// ZERO MODE ////////////////////////////////////////////////

 end subroutine

   
   
   
 subroutine initialise_the_data_structure(self, scheme_id)
   class(full_problem_data_structure_T), intent(inOut) :: self
   character(len=7), intent(in) :: scheme_id
   integer :: ix, num_args
   character(len=32), dimension(:), allocatable :: args
   call self%shandle%init(scheme_id, (my_rank.eq.0))
   call self%allocate_all_buffers()
   call self%build_operators()
   ! now handle the command line arguments
   self%cargo%initialCondition_amp_kxky = 1.d-16
   self%cargo%initialCondition_amp_zero = 1.d-16
   self%cargo%cflFactor_along_z = 1.d0             
   num_args = command_argument_count()
   allocate(args(num_args))
   do ix = 1, num_args
      call get_command_argument(ix,args(ix))
   end do
   do ix = 1, num_args
      select case (args(ix)) 
             case ('--noises-amplitude')
                  read (args(ix+1),*) self%cargo%initialCondition_amp_kxky 
                  read (args(ix+1),*) self%cargo%initialCondition_amp_zero 
             case ('--noise-amplitude-k')
                  read (args(ix+1),*) self%cargo%initialCondition_amp_kxky 
             case ('--noise-amplitude-0')
                  read (args(ix+1),*) self%cargo%initialCondition_amp_zero 
             case ('--cflFactor-along-z')
                  read (args(ix+1),*) self%cargo%cflFactor_along_z
      end select
   end do
 end subroutine

 subroutine allocate_all_buffers(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar, iStep, iSys
   allocate( self%coupled_kxky_set( &
             self%recipe%numberOf_coupled_kxky_systems ))
   ! => linearly coupled systems, kxky modes 
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems 
     allocate(                                          &
         self%coupled_kxky_set(iSys)%  field (          &
         self%geometry% NZAA,                           &
         self%geometry% spec% local_NY,                 &
         self%geometry% spec% local_NX,                 &
         self%recipe%kxky_recipes(iSys)%n_coupled_vars) &
            )
     allocate(                                          &
         self%coupled_kxky_set(iSys)%  rhs   (          &
         self%geometry% NZAA,                           &
         self%geometry% spec% local_NY,                 &
         self%geometry% spec% local_NX,                 &
         self%recipe%kxky_recipes(iSys)%n_coupled_vars) &
            )
     allocate(                                          &
         self%coupled_kxky_set(iSys)%  aux   (          &
         self%geometry% NZAA,                           &
         self%geometry% spec% local_NY,                 &
         self%geometry% spec% local_NX,                 &
         self%recipe%kxky_recipes(iSys)%n_coupled_vars) &
            )
     allocate( self%coupled_kxky_set(iSys)%step( self%shandle%imp_stages ))
     do iStep = 1, self%shandle%imp_stages
        allocate( self%coupled_kxky_set(iSys)%step(iStep)%K_hat( &
                  self%geometry% NZAA,                           &
                  self%geometry% spec% local_NY,                 &
                  self%geometry% spec% local_NX,                 &
                  self%recipe%kxky_recipes(iSys)%n_coupled_vars) &
                )
        allocate( self%coupled_kxky_set(iSys)%step(iStep)%K_std( &
                  self%geometry% NZAA,                           &
                  self%geometry% spec% local_NY,                 &
                  self%geometry% spec% local_NX,                 &
                  self%recipe%kxky_recipes(iSys)%n_coupled_vars) &
                )
     end do
   end do
   ! // linearly coupled systems, kxky modes 
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !  zero-mode counterpart here
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
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


 include "P3_IMEX_timestepping_marching.f90"
 include 'P3_IMEX_timestepping_output.f90'

 subroutine build_all_matrices(self)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer :: iSys, iVar
   
   if (my_rank.eq.0) print *, '================================================================='
   if (my_rank.eq.0) print *, '>>> Assembling coupling matrices for (kx,ky)-modes systems.'          
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
      call self%build_kxky_mass_matrices_for_system(iSys)
      call self%build_kxky_stif_matrices_for_system(iSys)
   end do
   if (my_rank.eq.0) print *, '/// Done with coupling matrices for (kx,ky)-modes systems.'
   if (my_rank.eq.0) print *, '================================================================='
   
 !delme call self%coupled_kxky_set(1)%mass% write2disk('mass')
 !delme call self%coupled_kxky_set(1)%stif% write2disk('stif')
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !  zero-mode counterpart here
 !///////////////////// ZERO MODE ////////////////////////////////////////////////

 end subroutine

 subroutine build_kxky_mass_matrices_for_system(self, isys)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer, intent(in) :: isys
   type(linear_operator_recipe_T), pointer :: my_recipe
   complex(kind=dp), allocatable :: numerical_prefactors_array(:,:,:)
   integer :: iPiece
   integer :: ix, iy, iz
   logical :: verbose = .False.
   if (my_rank.eq.0) verbose = .True.
   1987 Format (' >>> --- Mass matrix, system #',I2.2)
   if (verbose) write (*,1987) iSys
   1988 Format (' >>> --- ... coupling for var. #',I2.2,' in eqn. #',I2.2)
   my_recipe => self%recipe%kxky_recipes(isys)%mass 
   ! set the size of the zOperator_3d_1coupled_T object:
   call self%coupled_kxky_set(iSys)% mass% initialise( &
                   self% recipe% kxky_recipes(iSys)% n_coupled_vars, &
                   self% geometry% NZAA              , &
                   self% geometry% spec% local_NY    , &
                   self% geometry% spec% local_NX    )
   do iPiece = 1, my_recipe%n_pieces
     if (verbose) write (*,1988) my_recipe%jVar(iPiece), my_recipe%iEqn(iPiece)
     ! check that the recipe is sensible before calling building function
     ! compute the numerical prefactor
     numerical_prefactors_array = recipe_to_prefactors_array(my_recipe, iPiece, self%geometry)
     call self% coupled_kxky_set(iSys)% mass% fill_in( &
                  numerical_prefactors_array         , &
                  my_recipe% jVar (iPiece)           , &
                  my_recipe% iEqn (iPiece)           )
   end do

   !>>>>> 3/2 dealiasing
   !..... Dealiased wave numbers are set to zero
   !..... For these modes, we set M = 0 and L = id.
   do ix = 1, self%geometry% spec% local_NX
   do iy = 1, self%geometry% spec% local_NY
   do iz = 1, self%geometry% NZAA            
   if ( self%geometry% dealiase_x (ix)   .or. & 
        self%geometry% dealiase_y (iy)   .or. & 
        self%geometry% dealiase_z (iz) ) then
        self%coupled_kxky_set(iSys)% mass% mat (iz, iy, ix ,:,:) = cmplx(0._dp, 0._dp, kind=dp)
   end if
   end do
   end do
   end do
   !..... Treat the k=(0,0,0) mode similarly
   if (self%geometry%this_core_has_zero_mode) then
       self%coupled_kxky_set(iSys)% mass% mat (1,1,1,:,:) = cmplx(0._dp, 0._dp, kind=dp)
   end if
 end subroutine

 subroutine build_kxky_stif_matrices_for_system(self, isys)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer, intent(in) :: isys
   type(linear_operator_recipe_T), pointer :: my_recipe
   complex(kind=dp), allocatable :: numerical_prefactors_array(:,:,:)
   integer :: iPiece, iVar
   integer :: ix, iy, iz
   logical :: verbose = .False.
   if (my_rank.eq.0) verbose = .True.
   1987 Format (' >>> --- Mass matrix, system #',I2.2)
   if (verbose) write (*,1987) iSys
   1988 Format (' >>> --- ... coupling for var. #',I2.2,' in eqn. #',I2.2)
   my_recipe => self%recipe%kxky_recipes(isys)% stif
   ! set the size of the zOperator_3d_1coupled_T object:
   call self%coupled_kxky_set(iSys)% stif% initialise( &
                   self% recipe% kxky_recipes(iSys)% n_coupled_vars, &
                   self% geometry% NZAA              , &
                   self% geometry% spec% local_NY    , &
                   self% geometry% spec% local_NX    )
   do iPiece = 1, my_recipe%n_pieces
     if (verbose) write (*,1988) my_recipe%jVar(iPiece), my_recipe%iEqn(iPiece)
     ! check that the recipe is sensible before calling building function
     ! compute the numerical prefactor
     numerical_prefactors_array = recipe_to_prefactors_array(my_recipe, iPiece, self%geometry)
     call self% coupled_kxky_set(iSys)% stif% fill_in( &
                  numerical_prefactors_array         , &
                  my_recipe% jVar (iPiece)           , &
                  my_recipe% iEqn (iPiece)           )
   end do
   !>>>>> 3/2 dealiasing
   !..... Dealiased wave numbers are set to zero
   !..... For these modes, we set M = 0 and L = id.
   do ix = 1, self%geometry% spec% local_NX
   do iy = 1, self%geometry% spec% local_NY
   do iz = 1, self%geometry% NZAA            
   if ( self%geometry% dealiase_x (ix)   .or. & 
        self%geometry% dealiase_y (iy)   .or. & 
        self%geometry% dealiase_z (iz) ) then
        self%coupled_kxky_set(iSys)% stif% mat (iz, iy, ix, :,:) = cmplx(0._dp, 0._dp, kind=dp)
        do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
        self%coupled_kxky_set(iSys)% stif% mat (iz, iy, ix, iVar,iVar) = cmplx(1._dp, 0._dp, kind=dp)
        end do
   end if
   end do
   end do
   end do
   !..... Treat the k=(0,0,0) mode similarly
   if (self%geometry%this_core_has_zero_mode) then
       self%coupled_kxky_set(iSys)% stif% mat (1,1,1,:,:) = cmplx(0._dp, 0._dp, kind=dp)
       do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
       self%coupled_kxky_set(iSys)% stif% mat (1,1,1,iVar,iVar) = cmplx(1._dp, 0._dp, kind=dp)
       end do
   end if
 end subroutine


pure function recipe_to_prefactors_array(a_recipe, iPiece, geo)
   integer, intent(in) :: iPiece
   Type(linear_operator_recipe_T), pointer, intent(in) :: a_recipe
   complex(kind=dp), allocatable :: recipe_to_prefactors_array (:,:,:)
   complex(kind=dp), allocatable :: sh(:,:,:)
   type(geometry_vars_T), intent(in) :: geo
   integer :: ix, iy, iz
   allocate( sh (geo%NZAA, geo%spec%local_NY, geo%spec%local_NX) )
   do ix = 1, geo%spec%local_NX ! X is the slow index in spectral!
   do iy = 1, geo%spec%local_NY
   do iz = 1, geo%NZAA            
      sh(iz,iy,ix) = a_recipe%dsca(iPiece) &
                * ((geo%pz(iz))**a_recipe%kz_exponent(iPiece))&
                * ((geo%py(iy))**a_recipe%ky_exponent(iPiece))&
                * ((geo%px(ix))**a_recipe%kx_exponent(iPiece))
   end do
   end do
   end do
   call move_alloc(from=sh, to=recipe_to_prefactors_array)
 end function
 
end module P3_IMEX_timestepping
