 subroutine one_step_beyond(self, my_dt)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer :: i, j
   real(dp), intent(in) :: my_dt
   integer :: iSys, ix, iy

   call self%copy_fields_to_aux()

   do i = 1, self%shandle%imp_stages
      call self%compute_NL_terms_to_K_hat(i) !< stores directly the result in step(i)%K_hat
      ! ...
      ! building the right-hand side
      ! ...
      call self%mass_field_to_rhs()
      do j=1, i-1
         call self%add_K_std_to_rhs(my_dt*self%shandle%a_arr(i  ,j),j)
         call self%add_K_hat_to_rhs(my_dt*self%shandle%a_hat(i+1,j),j)
      end do
      call self%add_K_hat_to_rhs(my_dt*self%shandle%a_hat(i+1,i),i)
      call self%dealiase_the_rhs()
      ! ...
      ! implicit step (solve)
      ! ...
      call self%backsolve_to_aux()
      ! ...
      ! compute the next k_std
      ! ...
      call self%stif_aux_to_K_std(i)
   end do

   call self%copy_aux_to_fields()
   
 end subroutine one_step_beyond

 subroutine dealiase_the_rhs(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iSys, ix, iy, iz
   do ix = 1, self%geometry% spec% local_NX
   do iy = 1, self%geometry% spec% local_NY
   do iz = 1, self%geometry% NZAA            
   if ( self%geometry% dealiase_x (ix)   .or. & 
        self%geometry% dealiase_y (iy)   .or. & 
        self%geometry% dealiase_z (iz) ) then
        do isys = 1, self%recipe%numberOf_coupled_kxky_systems
        self%coupled_kxky_set(iSys)% rhs (iz, iy, ix, :) = cmplx(0._dp, 0._dp, kind=dp)
        end do
   end if
   end do
   end do
   end do
 end subroutine dealiase_the_rhs


 subroutine LPIMEX_backsolve_to_aux(self)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      call self%coupled_kxky_set(isys)%evol%backsolve( &
           self%coupled_kxky_set(isys)%rhs, &
           self%coupled_kxky_set(isys)%aux)
   end do
 end subroutine 


 subroutine LPIMEX_add_K_hat_to_rhs(self, dsca, i)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: i
   real(kind=dp), intent(in) :: dsca
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      self%coupled_kxky_set(isys)%rhs = self%coupled_kxky_set(isys)%rhs &
                                      + dsca*self%coupled_kxky_set(isys)%step(i)%K_hat
   end do
 end subroutine
   
 subroutine LPIMEX_add_K_std_to_rhs(self, dsca, i)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: i
   real(kind=dp), intent(in) :: dsca
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      self%coupled_kxky_set(isys)%rhs = self%coupled_kxky_set(isys)%rhs &
                                      + dsca*self%coupled_kxky_set(isys)%step(i)%K_std
   end do
 end subroutine


 subroutine LPIMEX_stif_aux_to_K_std(self, i)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: i
   integer :: iSys
   integer :: ix, iy
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
      call self%coupled_kxky_set(isys)%stif%dot(self%coupled_kxky_set(isys)%aux,         &
                                                self%coupled_kxky_set(isys)%step(i)%K_std)
   end do
 end subroutine

 subroutine LPIMEX_mass_field_to_rhs(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      call self%coupled_kxky_set(isys)%mass%dot(self%coupled_kxky_set(isys)%field, &
                                                self%coupled_kxky_set(isys)%rhs    )
   end do
 end subroutine

 subroutine LPIMEX_copy_aux_to_fields(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      self%coupled_kxky_set(isys)%field = self%coupled_kxky_set(isys)%aux
   end do
 end subroutine

 subroutine LPIMEX_copy_fields_to_aux(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      self%coupled_kxky_set(isys)%aux = self%coupled_kxky_set(isys)%field
   end do
 end subroutine


 subroutine compute_cfl_based_timestep(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   real(kind=dp) :: u_over_dx_max
   real(kind=dp) :: v_over_dy_max
   real(kind=dp) :: w_over_dz_max
   integer :: ix, iy, iz
   ! beware: potential stack problem by computing the max of the abs of an array?
   u_over_dx_max = maxVal(dAbs(self%linear_variables(1)%phys)) 
   u_over_dx_max = u_over_dx_max * self%geometry%NXAA / self%geometry%Lx
   v_over_dy_max = maxVal(dAbs(self%linear_variables(2)%phys))
   v_over_dy_max = v_over_dy_max * self%geometry%NYAA / self%geometry%Ly
   w_over_dz_max = maxVal(dAbs(self%linear_variables(3)%phys))
   w_over_dz_max = w_over_dz_max * self%geometry%NZAA / self%geometry%Lz

   self%cargo%cfl_based_DT = 1._dp/ dMax1 ( u_over_dx_max,&
                                            v_over_dy_max,&
                                            w_over_dz_max)

 end subroutine


 subroutine LPIMEX_compute_NL_terms_to_K_hat(self, i)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: i 
   integer :: isys
   integer :: iTerm
   integer :: iQvar, ix, iy, iz, iVar
   integer :: insert_in_equation
   integer :: insert_in_system
   complex(kind=dp), allocatable :: nl_buffer_deAliased(:,:,:)
   integer :: iSource

   call self%compute_full_variables_in_physical_space()
   
!
   call self%compute_cfl_based_timestep()

   do iQvar = 1, self%recipe%numberOf_quadratic_variables
      ! the recipe is self%recipe%nl_vars(iQvar
      ! largeArrayWarning 
      self%quadratic_variables(iQvar)%phys = &
              self%linear_variables( self%recipe%nl_vars(iQvar)%iVar1 )%phys *&
              self%linear_variables( self%recipe%nl_vars(iQvar)%iVar2 )%phys
      call self%quadratic_variables(iQvar)%phys_to_spec()
   end do

   if (i.eq.1) call self%output_global_quantities()
   if (i.eq.1) call self%output_slices_volumes_and_profiles()
  
   allocate( nl_buffer_deAliased( &
                        self%geometry%NZAA,&
                        self%geometry%spec%local_NY, &
                        self%geometry%spec%local_NX) )


   ! => for each set of coupled eqns: linear combinations of 
   ! self%quadratic_variables(iQvar)%spec are computed 
   ! and stored in the corresponding K_hat array
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
   self%coupled_kxky_set(isys)%step(i)%K_hat = cmplx(0._dp, 0._dp, kind=dp)
   do iTerm = 1, self%recipe%kxky_recipes(iSys)%NL%n_pieces
      nl_buffer_deAliased = cmplx(0._dp, 0._dp, kind=dp)
      do ix = 1, self%geometry%spec%local_NX
      do iy = 1, self%geometry%spec%local_NY
      do iz = 1, self%geometry%NZ
      nl_buffer_deAliased(iz,iy,ix) = self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dsca * &
             self%geometry%px(ix)**self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dx_exponent *&
             self%geometry%py(iy)**self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dy_exponent *&
             self%geometry%pz(iz)**self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dz_exponent *&
             self%quadratic_variables(&
                      self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%quad_var_index&
                                     )%spec(iz,iy,ix)
      end do
      end do
      end do

      insert_in_equation =  self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%eqn_index 

      self% coupled_kxky_set(iSys)% step(i)% K_hat(:,:,:,insert_in_equation) = &
      self% coupled_kxky_set(iSys)% step(i)% K_hat(:,:,:,insert_in_equation) + nl_buffer_deAliased

   end do
   end do
   
   !//////////////// DONE DIRECTLY FOR THE RHS ////////////////
   !deAliase K_hat
   !do isys = 1, self%recipe%numberOf_coupled_kxky_systems
   !do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
   !!do ix = 1, self%geometry%spec%local_NX
   !!do iy = 1, self%geometry%spec%local_NY
   !if (self%geometry%deAliase_x(ix) .or. &
       !self%geometry%deAliase_y(iy) ) then
       !self%coupled_kxky_set(isys)%step(i)%K_hat(:, iy, ix, iVar) = cmplx(0._dp, 0._dp, kind=dp)
   !end if
   !end do
   !end do
   !end do
   !end do
      
 end subroutine

 subroutine compute_full_variables_in_physical_space(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar, jVar
   integer :: iTerm
   integer :: ix, iy, iz
   integer :: system_of_interest, position_of_interest
   complex (kind=dp), allocatable :: datBuffer (:,:,:)

   allocate( datBuffer (self%geometry%NZAA, &
                        self%geometry%spec%local_NY, &
                        self%geometry%spec%local_NX) )
   ! ==================================================
   ! >>> spectral space computations
   do iVar = 1, self%recipe%numberOf_linear_variables_full
   self%linear_variables(iVar)%spec = cmplx(0._dp, 0._dp, kind=dp)
   do iTerm = 1, self%recipe%linear_vars_full( iVar )%n_terms
     select case (self%recipe%linear_vars_full( iVar)%Term(iTerm)%var_kind)
       case ('kxky')
          system_of_interest =               &
               self%recipe%linear_vars_kxky( &
               self%recipe%linear_vars_full( iVar )%term(iTerm)%var_index&
                                                  )%belongs_to_set
          position_of_interest =  &
               self%recipe%linear_vars_kxky( &
               self%recipe%linear_vars_full( iVar )%term(iTerm)%var_index&
                                           )%at_position
                   
          datBuffer = self% recipe% linear_vars_full(iVar)% term(iterm)% dsca * &
                      self% coupled_kxky_set(system_of_interest)% aux(:,:,:,position_of_interest)

          do ix = 1, self%geometry%spec%local_NX
          do iy = 1, self%geometry%spec%local_NY
          do iz = 1, self%geometry%Nz       
          self%linear_variables(iVar)%spec(iz,iy,ix) = &
                       self%linear_variables(iVar)%spec(iz,iy,ix) + &
                       datBuffer(iz,iy,ix)*&
                       (self%geometry%px(ix)**self%recipe%linear_vars_full(iVar)%Term(iTerm)%dx_exponent)*&
                       (self%geometry%pz(iz)**self%recipe%linear_vars_full(iVar)%Term(iTerm)%dz_exponent)*&
                       (self%geometry%py(iy)**self%recipe%linear_vars_full(iVar)%Term(iTerm)%dy_exponent)
          end do
          end do
          end do

       case ('zero')
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
 !  zero-mode counterpart here
 !///////////////////// ZERO MODE ////////////////////////////////////////////////
     end select
   end do
   ! /// done: spectral space computations
   ! ==================================================


   ! ==================================================
   ! >>> spectral to physical tranforms  
   call self%linear_variables(iVar)%spec_to_phys()
   ! /// done: transforms                     
   ! ==================================================
   end do
   deAllocate( datBuffer )
 end subroutine


 subroutine horizontal_average_of_physical_quantity_inPlace(self, iVar)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iVar
   real(dp), allocatable :: local_a1 (:)
   real(dp), allocatable :: global_a1(:)
   integer :: ix, iy
   ! stores the horizontal average of self%linear_variables(iVar)%phys in 
   ! the allocatable a1(:), and broadcast it back to the rank-3 array
   allocate (local_a1  (self%geometry%phys%local_NZ) )
   allocate (global_a1 (self%geometry%phys%local_NZ) )
   local_a1 = sum(sum(self%linear_variables(iVar)%phys, 3), 2) 
   local_a1 = local_a1 / self%geometry%NXAA 
   local_a1 = local_a1 / self%geometry%NYAA 
   ! averaged of distributed y
   call MPI_allReduce(local_a1,  & !send buffer
                      global_a1, & !recv buffer
                      self%geometry%phys%local_NZ, & !count
                      MPI_double, &     !dtype
                      MPI_sum, &
                      self%geometry%mpi_Zphys%comm, ierr)
   do ix = 1, self%geometry%phys%local_NX
   do iy = 1, self%geometry%phys%local_NY
      self%linear_variables(iVar)%phys(:, iy, ix) = global_a1
   end do
   end do
 end subroutine


 subroutine initialise_bookkeeping_counters (self, time_index)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar, iObj
   integer, intent(in) :: time_index
   self%io_bookkeeping%time_integer = time_index 
   do iVar = 1, self%recipe%numberOf_linear_variables_full
   do iObj = 1, self%recipe%output%numberOf_linearObjects( iVar )
   self%recipe%output%linear(iVar)%object(iObj)%counter = &
          modulo(time_index, self%recipe%output%linear(iVar)%object(iObj)%period)
   end do
   end do
   do iVar = 1, self%recipe%numberOf_quadratic_variables
   do iObj = 1, self%recipe%output%numberOf_quadraObjects( iVar )
   self%recipe%output%quadra(iVar)%object(iObj)%counter = &
          modulo(time_index, self%recipe%output%quadra(iVar)%object(iObj)%period)
   end do
   end do
 end subroutine

 subroutine output_slices_volumes_and_profiles(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar, iObj
   self%io_bookkeeping%time_integer = self%io_bookkeeping%time_integer + 1
   do iVar = 1, self%recipe%numberOf_linear_variables_full
   do iObj = 1, self%recipe%output%numberOf_linearObjects( iVar )
      self%recipe%output%linear(iVar)%object(iObj)%counter = self%recipe%output%linear(iVar)%object(iObj)%counter + 1
      if (self%recipe%output%linear(iVar)%object(iObj)%counter .eq. &
          self%recipe%output%linear(iVar)%object(iObj)%period) then
          self%recipe%output%linear(iVar)%object(iObj)%counter = 0 ! reset to 0, and output
          select case (self%recipe%output%linear(iVar)%object(iObj)%kind) 
                 case ('volume')
                 call self%export_volume(iVar, self%io_bookkeeping%time_integer, 'linear')
                 case ('zSlice')
                 call self%export_slice (iVar, self%io_bookkeeping%time_integer, &
                                        1, self%recipe%output%linear(iVar)%object(iObj)%slice_index, &
                                        'linear')
                 case ('ySlice')
                 call self%export_slice (iVar, self%io_bookkeeping%time_integer, &
                                        2, self%recipe%output%linear(iVar)%object(iObj)%slice_index, &
                                        'linear')
                 case ('xSlice')
                 call self%export_slice (iVar, self%io_bookkeeping%time_integer, &
                                        3, self%recipe%output%linear(iVar)%object(iObj)%slice_index, &
                                        'linear')
                 case ('zAvged')
                 call self%export_verticallyAvgedslice (iVar, self%io_bookkeeping%time_integer,&
                                        'linear')
                 case ('profil')
                 call self%export_profile (ivar, self%io_bookkeeping%time_integer, 'linear')
                 case default
                 print *, 'invalid self%recipe%output%linear(iVar)%object(iObj)%kind for iVar,iObj:', iVar, iObj
                 error stop
          end select
      end if
   end do
   end do
   do iVar = 1, self%recipe%numberOf_quadratic_Variables
   do iObj = 1, self%recipe%output%numberOf_quadraObjects( iVar )
      self%recipe%output%quadra(iVar)%object(iObj)%counter = self%recipe%output%quadra(iVar)%object(iObj)%counter + 1
      if (self%recipe%output%quadra(iVar)%object(iObj)%counter .eq. &
          self%recipe%output%quadra(iVar)%object(iObj)%period) then
          self%recipe%output%quadra(iVar)%object(iObj)%counter = 0 ! reset to 0, and output
          select case (self%recipe%output%quadra(iVar)%object(iObj)%kind) 
                 case ('volume')
                 call self%export_volume(iVar, self%io_bookkeeping%time_integer, 'quadra')
                 case ('zSlice')
                 call self%export_slice (iVar, self%io_bookkeeping%time_integer, &
                                        1, self%recipe%output%quadra(iVar)%object(iObj)%slice_index, &
                                        'quadra')
                 case ('ySlice')
                 call self%export_slice (iVar, self%io_bookkeeping%time_integer, &
                                        2, self%recipe%output%quadra(iVar)%object(iObj)%slice_index, &
                                        'quadra')
                 case ('xSlice')
                 call self%export_slice (iVar, self%io_bookkeeping%time_integer, &
                                        3, self%recipe%output%quadra(iVar)%object(iObj)%slice_index, &
                                        'quadra')
                 case ('profil')
                 call self%export_profile (ivar, self%io_bookkeeping%time_integer, 'quadra')
                 case default
                 print *, 'invalid self%recipe%output%quadra(iVar)%object(iObj)%kind for iVar,iObj:', iVar, iObj
                 error stop
          end select
      end if
   end do
   end do

 end subroutine


