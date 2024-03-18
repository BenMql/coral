

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
   !deAliase fields
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
   do ix = 1, self%geometry%spec%local_NX
   do iy = 1, self%geometry%spec%local_NY
   if (self%geometry%deAliase_x(ix) .or. &
       self%geometry%deAliase_y(iy) ) then
       self%coupled_kxky_set(isys)%field(:,iy,ix) = cmplx(0._dp, 0._dp, kind=dp)
   end if
   end do
   end do
   if (self%geometry%this_core_has_zero_mode) then
       self%coupled_kxky_set(isys)%field(:, 1, 1) = cmplx(0._dp, 0._dp, kind=dp)
   end if
   end do
   
 end subroutine one_step_beyond




 subroutine LPIMEX_backsolve_to_aux(self)
   class(full_problem_data_structure_T), intent(inOut), target :: self
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      call self%coupled_kxky_set(isys)%evol%backsolve( &
           self%coupled_kxky_set(isys)%rhs)
      !print*, my_rank , 'rhs/aux shape', &
                    !shape( self% coupled_kxky_set(isys)% rhs ), &
                    !shape( self% coupled_kxky_set(isys)% aux )
      self%coupled_kxky_set(isys)%aux(&
           1: self% coupled_kxky_set(iSys)% shape% total, &
           1: self% coupled_kxky_set(iSys)% shape% spectral_local_NY, &
           1: self% coupled_kxky_set(iSys)% shape% spectral_local_NX&
                                     )       = & 
      self%coupled_kxky_set(isys)%rhs(&
           1: self% coupled_kxky_set(iSys)% shape% total, &
           1: self% coupled_kxky_set(iSys)% shape% spectral_local_NY, &
           1: self% coupled_kxky_set(iSys)% shape% spectral_local_NX&
                                     )
   !do ix = 1, self% coupled_kxky_set(iSys)% shape% spectral_local_NX
   !do iy = 1, self% coupled_kxky_set(iSys)% shape% spectral_local_NY
      !self%coupled_kxky_set(isys)%aux(&
           !1: self% coupled_kxky_set(iSys)% shape% total, iy, ix &
                                     !) = &
      !self%coupled_kxky_set(isys)%rhs(&
           !1: self% coupled_kxky_set(iSys)% shape% total, iy, ix &
                                     !) 
   !end do 
   !end do
   end do
   if (self%geometry%this_core_has_zero_mode) then
   do isys = 1, self%recipe%numberOf_coupled_zero_systems
      call self%coupled_zero_set(isys)%evol%backsolve( &
           self%coupled_zero_set(isys)%rhs)
      self%coupled_zero_set(isys)%aux = self%coupled_zero_set(isys)%rhs
   end do
   end if
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
   if (self%geometry%this_core_has_zero_mode) then
   do isys = 1, self%recipe%numberOf_coupled_zero_systems
      self%coupled_zero_set(isys)%rhs = self%coupled_zero_set(isys)%rhs &
                                      + dsca*self%coupled_zero_set(isys)%step(i)%K_hat
   end do
   end if
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
   if (self%geometry%this_core_has_zero_mode) then
   do isys = 1, self%recipe%numberOf_coupled_zero_systems
      self%coupled_zero_set(isys)%rhs = self%coupled_zero_set(isys)%rhs &
                                      + dsca*self%coupled_zero_set(isys)%step(i)%K_std
   end do
   end if
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
   if (self%geometry%this_core_has_zero_mode) then
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
      call self%coupled_zero_set(isys)%stif%dot(self%coupled_zero_set(isys)%aux,         &
                                                self%coupled_zero_set(isys)%step(i)%K_std)
   end do
   end if
   !deAliase K_std
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
   do ix = 1, self%geometry%spec%local_NX
   do iy = 1, self%geometry%spec%local_NY
   if (self%geometry%deAliase_x(ix) .or. &
       self%geometry%deAliase_y(iy) ) then
       self%coupled_kxky_set(isys)%step(i)%K_std(:,iy,ix) = cmplx(0._dp, 0._dp, kind=dp)
   end if
   end do
   end do
   if (self%geometry%this_core_has_zero_mode) then
       self%coupled_kxky_set(isys)%step(i)%K_std(:, 1, 1) = cmplx(0._dp, 0._dp, kind=dp)
   end if
   end do
 end subroutine

 subroutine LPIMEX_mass_field_to_rhs(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      call self%coupled_kxky_set(isys)%mass%dot(self%coupled_kxky_set(isys)%field, &
                                                self%coupled_kxky_set(isys)%rhs    )
   end do
   if (self%geometry%this_core_has_zero_mode) then
   do isys = 1, self%recipe%numberOf_coupled_zero_systems
      call self%coupled_zero_set(isys)%mass%dot(self%coupled_zero_set(isys)%field, &
                                                self%coupled_zero_set(isys)%rhs    )
   end do
   end if
 end subroutine

 subroutine LPIMEX_copy_aux_to_fields(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      self%coupled_kxky_set(isys)%field = self%coupled_kxky_set(isys)%aux
   end do
   if (self%geometry%this_core_has_zero_mode) then
   do isys = 1, self%recipe%numberOf_coupled_zero_systems
      self%coupled_zero_set(isys)%field = self%coupled_zero_set(isys)%aux
   end do
   end if
 end subroutine

 subroutine LPIMEX_copy_fields_to_aux(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: isys
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
      self%coupled_kxky_set(isys)%aux = self%coupled_kxky_set(isys)%field
   end do
   if (self%geometry%this_core_has_zero_mode) then
   do isys = 1, self%recipe%numberOf_coupled_zero_systems
      self%coupled_zero_set(isys)%aux = self%coupled_zero_set(isys)%field
   end do
   end if
 end subroutine


 subroutine compute_cfl_based_timestep(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   real(kind=dp) :: u_over_dx_max
   real(kind=dp) :: v_over_dy_max
   real(kind=dp) :: w_over_dz_max
   real(kind=dp), allocatable :: local_cheby_weight(:)
   integer :: ix, iy 
   call gauss_chebyshev_weight_1d( self%gauss_cheby%weight1d, &
                                   self%geometry%NZAA,&
                                   self%geometry%gap)
                        
   allocate( local_cheby_weight ( self%geometry%phys%local_NZ ))
   local_cheby_weight = self%gauss_cheby%weight1d (domain_decomp%phys_iStart(1) : &
                                                   domain_decomp%phys_iStart(1) + &
                                                   self%geometry%phys%local_NZ  -1)
   if (domain_decomp% y_is_empty_in_phys) then
   self%cargo%cfl_based_DT = 1._dp
   else
   ! beware: potential stack problem by computing the max of the abs of an array?
   u_over_dx_max = maxVal(dAbs(self%linear_variables(1)%phys)) 
   u_over_dx_max = u_over_dx_max * self%geometry%NXAA / self%geometry%Lx
   v_over_dy_max = maxVal(dAbs(self%linear_variables(2)%phys))
   v_over_dy_max = v_over_dy_max * self%geometry%NYAA / self%geometry%Ly
   w_over_dz_max = v_over_dy_max
   do ix = 1, self%geometry%phys%local_NX
   do iy = 1, self%geometry%phys%local_NY
   w_over_dz_max = dMax1(w_over_dz_max,&
                         maxVal( dAbs( self%linear_variables(3)%phys(:,iy,ix) / &
                                       local_cheby_weight ) * &
                                       self%cargo%cflFactor_along_z   ) )
   end do
   end do
   

   self%cargo%cfl_based_DT = 1._dp/ dMax1 ( u_over_dx_max,&
                                            v_over_dy_max,&
                                            w_over_dz_max)
   end if 

   if (self%recipe%smagorinsky_flag) then
      self%cargo%cfl_based_DT = dMin1( self%cargo%cfl_based_DT, &
                                0.95/ ( (2*3.141592)**2*&
            maxVal(self%linear_variables(&
                     self%recipe%numberof_linear_variables_full&
                                          )%phys)  &
                                      *( &
                           (self%geometry%NXAA/self%geometry%Lx)**2&
                         + (self%geometry%NYAA/self%geometry%Ly)**2&
                                       ) &
                                      )  & 
                                     )  
   end if 
 end subroutine


 subroutine LPIMEX_compute_NL_terms_to_K_hat(self, i)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: i 
   integer :: isys
   integer :: iTerm
   integer :: iQvar, ix, iy
   integer :: insert_in_equation
   integer :: insert_in_system
   type(zcsr_matrix) :: zShuffleQ_insert_Iz
   type(dcsr_matrix) :: dShuffleQ_insert_Iz
   type(zcsr_matrix) :: zTruncI
   type(dcsr_matrix) :: TruncI
   complex(kind=dp), allocatable :: nl_buffer(:,:,:)
   complex(kind=dp), allocatable :: nl_buffer_deAliased(:,:,:)
   real   (kind=dp), allocatable :: nl_buffer_zero(:)
   integer :: iSource

   call self%compute_full_variables_in_physical_space()
   
!
   call self%compute_cfl_based_timestep()

   do iQvar = 1, self%recipe%numberOf_quadratic_variables
      ! largeArrayWarning 
      if (self%recipe%nl_vars(iQvar)%iVar1 == -1) then
        self%quadratic_variables(iQvar)%phys = &
              self%linear_variables( self%recipe%nl_vars(iQvar)%iVar2 )%phys
      else if (self%recipe%nl_vars(iQvar)%iVar2 == -1) then
        self%quadratic_variables(iQvar)%phys = &
              self%linear_variables( self%recipe%nl_vars(iQvar)%iVar1 )%phys
      else
        self%quadratic_variables(iQvar)%phys = &
              self%linear_variables( self%recipe%nl_vars(iQvar)%iVar1 )%phys *&
              self%linear_variables( self%recipe%nl_vars(iQvar)%iVar2 )%phys
      end if
     
      call self%quadratic_variables(iQvar)%phys_to_spec()
      if (self%recipe%nl_vars(iQvar)% remove_z_integral) then
      ! first we filter the vertically invariant component in spectral space
      call self%remove_vertical_mean_of_quadratic_var (iQvar)
      ! and we propagate that back to physical space, for output purposes
      call self%quadratic_variables(iQvar)%spec_to_phys()
      end if
   end do


   if (i.eq.1) call self%output_global_quantities()
   if (i.eq.1) call self%output_slices_volumes_and_profiles()
  
   allocate( nl_buffer( domain_decomp% spec_iSize(1), &
                        domain_decomp% spec_iSize(2), &
                        domain_decomp% spec_iSize(3)) )

   allocate( nl_buffer_deAliased( &
                        self%geometry%NZ,&
                        domain_decomp% spec_iSize(2), &
                        domain_decomp% spec_iSize(3)) )


   ! => for each set of coupled eqns: linear combinations of 
   ! self%quadratic_variables(iQvar)%spec are computed 
   ! and stored in the corresponding K_hat array
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
   self%coupled_kxky_set(isys)%step(i)%K_hat = cmplx(0._dp, 0._dp, kind=dp)
   do iTerm = 1, self%recipe%kxky_recipes(iSys)%NL%n_pieces
      do ix = 1, domain_decomp% spec_iSize(3)
      do iy = 1, domain_decomp% spec_iSize(2)
      nl_buffer(:,iy,ix) = self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dsca * &
             self%geometry%px(ix)**self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dx_exponent *&
             self%geometry%py(iy)**self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dy_exponent *&
             self%quadratic_variables(&
                      self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%quad_var_index&
                                     )%spec(:,iy,ix)
      end do
      end do
      nl_buffer_deAliased = nl_buffer(1:self%geometry%NZ, :, :)
      !print *,isys, iTerm,  sum(abs(nl_buffer_deAliased))
      !=> project the Iz integration over the higher Tchebychev polynomials
      !   (ditching the lowest projections, polluted by integration constants)
      !   kind of wasteful to do that here, but probably not a biggy
      insert_in_equation =  self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%eqn_index 

      call self%coupled_kxky_set(iSys)%building_tools%cheb_IM ( &
                          self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%Iz_exponent+1, &
                          self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%z_multiply +1  &
                                                              ) &
                        % truncate(  truncI, &
                         self%recipe%kxky_recipes(iSys)%eqn_order( insert_in_equation )+1, &
                         self%coupled_kxky_set(iSys)%building_tools%cheb_IM(1,1)%nrow, &
                    1,   self%coupled_kxky_set(iSys)%building_tools%cheb_IM(1,1)%ncol)
      call dcsr_convert_zcsr(truncI, zTruncI)
      call self%coupled_kxky_set(iSys)%shuffleInsert(insert_in_equation)%dot(ztruncI, zShuffleQ_insert_Iz)
      call zShuffleQ_insert_Iz%dot( NL_buffer_deAliased, self%coupled_kxky_set(isys)%step(i)%K_hat, 'cumul')
      !print *, isys, iTerm, sum(abs(self%coupled_kxky_set(isys)%step(i)%K_hat))
   end do
   end do

   !deAliase K_hat
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
   do ix = 1, domain_decomp% spec_iSize(3)
   do iy = 1, domain_decomp% spec_iSize(2)
   if (self%geometry%deAliase_x(ix) .or. &
       self%geometry%deAliase_y(iy) ) then
       self%coupled_kxky_set(isys)%step(i)%K_hat(:,iy,ix) = cmplx(0._dp, 0._dp, kind=dp)
   end if
   end do
   end do
   if (self%geometry%this_core_has_zero_mode) then
       self%coupled_kxky_set(isys)%step(i)%K_hat(:, 1, 1) = cmplx(0._dp, 0._dp, kind=dp)
   end if
   end do
      
   ! => for each set of coupled eqns: linear combinations of 
   ! self%quadratic_variables(iQvar)%spec are computed 
   ! and stored in the corresponding K_hat array
   if (self%geometry%this_core_has_zero_mode) then
           allocate( nl_buffer_zero( domain_decomp% spec_iSize(1) )) 
   do isys = 1, self%recipe%numberOf_coupled_zero_systems
   self%coupled_zero_set(isys)%step(i)%K_hat = 0._dp
   do iTerm = 1, self%recipe%zero_recipes(iSys)%NL%n_pieces
     nl_buffer_zero(:) = self%recipe%zero_recipes(iSys)%NL%term(iTerm)%dsca* &
       real(self%quadratic_variables(&
            self%recipe%zero_recipes(iSys)%NL%term(iTerm)%quad_var_index)%spec(:,1,1),&
            kind=dp)

     !=> project the Iz integration over the higher Tchebychev polynomials
     !   (ditching the lowest projections, polluted by integration constants)
     !   kind of wasteful to do that here, but probably not a biggy
     insert_in_equation =  self%recipe%zero_recipes(iSys)%NL%term(iTerm)%eqn_index 
     call self%coupled_zero_set(iSys)%building_tools%cheb_IM ( &
                       self%recipe%zero_recipes(iSys)%NL%term(iTerm)%Iz_exponent+1, &
                       self%recipe%zero_recipes(iSys)%NL%term(iTerm)%z_multiply +1  &
                                                             ) &
                       % truncate(  truncI, &
                        self%recipe%zero_recipes(iSys)%eqn_order( insert_in_equation )+1, &
                        self%coupled_zero_set(iSys)%building_tools%cheb_IM(1,1)%nrow, &
                   1,   self%coupled_zero_set(iSys)%building_tools%cheb_IM(1,1)%ncol)
     call self%coupled_zero_set(iSys)%shuffleInsert(insert_in_equation)%dot(truncI, dShuffleQ_insert_Iz)
     call dShuffleQ_insert_Iz%dot( NL_buffer_zero, self%coupled_zero_set(isys)%step(i)%K_hat, 'cumul')
   end do
   end do
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! add sources!                              
   !====================================
   do iSource = 1, self%recipe%sources%n_sources
     nl_buffer_zero = 0._dp
     nl_buffer_zero(:) = self%sources( self%recipe%sources%term(iSource)%source_index)%spectral (1: domain_decomp% spec_iSize(1)) &
                       * self%recipe%sources%term(iSource)%dsca
     insert_in_equation =  self%recipe%sources%term(iSource)%eqn_index
     insert_in_system   =  self%recipe%sources%term(iSource)%sys_index
     call self%coupled_zero_set(insert_in_system)%building_tools%cheb_IM&
                        (  self%recipe%sources%term(iSource)%Iz_exponent+1,1 )&
                       % truncate(  truncI, &
                        self%recipe%zero_recipes(insert_in_system)%eqn_order( insert_in_equation )+1, &
                        self%coupled_zero_set(insert_in_system)%building_tools%cheb_IM(1,1)%nrow, &
                   1,   self%coupled_zero_set(insert_in_system)%building_tools%cheb_IM(1,1)%ncol)
     call self%coupled_zero_set(insert_in_system)%shuffleInsert(insert_in_equation)%dot(truncI, dShuffleQ_insert_Iz)
     call dShuffleQ_insert_Iz%dot( NL_buffer_zero, self%coupled_zero_set(insert_in_system)%step(i)%K_hat, 'cumul')
   end do
   end if 

 end subroutine

 subroutine compute_full_variables_in_physical_space(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar, jVar
   integer :: iTerm
   integer :: ix, iy, iz
   integer :: system_of_interest, position_of_interest
   complex (kind=dp), allocatable :: datBuffer (:,:,:)
   real    (kind=dp), allocatable :: datBufferZero (:)
   integer :: parity_factor

   allocate( datBuffer (domain_decomp% spec_iSize(1), &
                        domain_decomp% spec_iSize(2), &
                        domain_decomp% spec_iSize(3)) )
   allocate( datBufferZero (domain_decomp% spec_iSize(1)) )
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
                   
          call self%coupled_kxky_set(system_of_interest &
                      )%padStencilExtractShuffle(position_of_interest&
                      )%dot(self%coupled_kxky_set(system_of_interest)%aux,&
                        datBuffer, 'overW')
          datBuffer = self%recipe%linear_vars_full(iVar)%Term(iTerm)%dsca* datBuffer
          if ((self%recipe%linear_vars_full(iVar)%Term(iTerm)%dx_exponent .eq. 0) .and.&
              (self%recipe%linear_vars_full(iVar)%Term(iTerm)%dy_exponent .eq. 0) .and.&
              (self%recipe%linear_vars_full(iVar)%Term(iTerm)%dz_exponent .eq. 0)) then
               self%linear_variables(iVar)%spec =  self%linear_variables(iVar)%spec &
                      + datBuffer
          else
          if (self%recipe%linear_vars_full(iVar)%Term(iTerm)%dz_exponent.ne.0) then
          call self%differentiate( datBuffer, &
                                   self%recipe%linear_vars_full(iVar)%Term(iTerm)%dz_exponent)
          end if
          do ix = 1, self%geometry%spec%local_NX
          do iy = 1, self%geometry%spec%local_NY
          self%linear_variables(iVar)%spec(:,iy,ix) = &
                       self%linear_variables(iVar)%spec(:,iy,ix) + &
                       datBuffer(:,iy,ix)*&
                       (self%geometry%px(ix)**self%recipe%linear_vars_full(iVar)%Term(iTerm)%dx_exponent)*&
                       (self%geometry%py(iy)**self%recipe%linear_vars_full(iVar)%Term(iTerm)%dy_exponent)
          end do
          end do
          end if
       case ('zero')
          if (self%geometry%this_core_has_zero_mode) then
          system_of_interest =               &
               self%recipe%linear_vars_zero( &
               self%recipe%linear_vars_full( iVar )%term(iTerm)%var_index&
                                                  )%belongs_to_set
          position_of_interest =             &
               self%recipe%linear_vars_zero( &
               self%recipe%linear_vars_full( iVar )%term(iTerm)%var_index&
                                                  )%at_position
          ! =======
          !call self%coupled_zero_set(system_of_interest &
                           !)%padStencilExtractShuffle(position_of_interest&
                           !)%dot(self%coupled_zero_set(system_of_interest)%aux,&
                                 !self%linear_variables(iVar)%spec, 'cumul')
          call self%coupled_zero_set(system_of_interest &
                           )%padStencilExtractShuffle(position_of_interest&
                           )%dot(self%coupled_zero_set(system_of_interest)%aux,&
                        datBufferZero, 'overW')
          if (self%recipe%linear_vars_full(iVar)%Term(iTerm)%dz_exponent.ne.0) then
          call self%differentiate( datBufferZero, &
                                   self%recipe%linear_vars_full(iVar)%Term(iTerm)%dz_exponent)
          end if
          self%linear_variables(iVar)%spec(:,1,1) =  self%linear_variables(iVar)%spec(:,1,1) &
                      + datBufferZero* self%recipe%linear_vars_full(iVar)%Term(iTerm)%dsca
          
          end if
     end select
   end do
   ! /// done: spectral space computations
   ! ==================================================

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
   ! todo => For the zero mode, we have assumed
   ! that no vertical derivatives are present.
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

   ! ==================================================
   ! >>> spectral to physical transforms  
   ! *except* for nuSmago
   if ( .not. self%recipe%linear_vars_full(iVar)% penalisation) then
   if ( self%recipe%linear_vars_full(iVar)%str .eq. 'nuSmago' ) then
       ! fixSmago fixSmago: below is not the exact expression!
       ! for the time being, we define nuSmago as the horizontal average of the 
       ! horizontal kinetic energy. We assume that uFull and vFull, the horizontal
       ! velocities are defined as the first and second 'full variable' in the 
       ! coral.equations file
       do jVar = 1,self%recipe%numberOf_linear_variables_full
        if (self%recipe%linear_vars_full(jVar)%str .eq. 'dxUfull') then
        self%linear_variables(iVar)%phys = self%linear_variables(iVar)%phys &
                                         + self%linear_variables(jVar)%phys**2
        else if (self%recipe%linear_vars_full(jVar)%str .eq. 'dyUfull') then
        self%linear_variables(iVar)%phys = self%linear_variables(iVar)%phys &
                                         + self%linear_variables(jVar)%phys**2
        else if (self%recipe%linear_vars_full(jVar)%str .eq. 'dxVfull') then
        self%linear_variables(iVar)%phys = self%linear_variables(iVar)%phys &
                                         + self%linear_variables(jVar)%phys**2
        else if (self%recipe%linear_vars_full(jVar)%str .eq. 'dyVfull') then
        self%linear_variables(iVar)%phys = self%linear_variables(iVar)%phys &
                                         + self%linear_variables(jVar)%phys**2
        end if
       end do
        
       ! self%linear_variables(iVar)%phys = 0._dp
       ! average it horizontally
       call self%horizontal_average_of_physical_quantity_inPlace( iVar )
       self%linear_variables(iVar)%phys = self%cargo% smagorinsky_prefactor * sqrt( self%linear_variables(iVar)%phys )
   else 
       ! print *, self%recipe%linear_vars_full(iVar)%extract_value_at 
       if (self%recipe%linear_vars_full(iVar)%extract_value_at == 'BotSurf') then
           do iz =2, self%geometry%NZ
           parity_factor = (-1)**iz
           self%linear_variables(iVar)%spec(1,:,:) = &
                   self%linear_variables(iVar)%spec (1,:,:) &
                 - self%linear_variables(iVar)%spec(iz,:,:) * parity_factor
           self%linear_variables(iVar)%spec(iz,:,:)  = cmplx(0._dp, 0._dp, kind=dp)
           end do
       else if (self%recipe%linear_vars_full(iVar)%extract_value_at == 'TopSurf') then
           do iz =2, self%geometry%NZ
           !call c_f_pointer(p1, self%spec,  [domain_decomp% NZ_long,&
                                             !domain_decomp% NYAA_long,&
                                             !domain_decomp% local_NX_spec]) 
           self%linear_variables(iVar)%spec(1,:,:) = &
                   self%linear_variables(iVar)%spec (1,:,:) &
                 + self%linear_variables(iVar)%spec(iz,:,:) 
           self%linear_variables(iVar)%spec(iz,:,:)  = cmplx(0._dp, 0._dp, kind=dp)
           end do
       end if
       call self%linear_variables(iVar)%spec_to_phys()
   end if
   end if
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
   !print *, global_a1
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



 real(kind=dp) function sourceParams(self, paramString)
 class(full_problem_data_structure_T), intent(in) :: self
 character(len=6), intent(in) :: paramString
 integer :: iParam
 sourceParams = 0._dp
 !print *, "DELME searching,", paramString,","
 do iParam = 1, self%recipe%sources% n_sourceParams
 !print *, "DELME reading  ,", self%recipe%sources%sourceParams(iParam)%str,","
    if (paramString .eq. self%recipe%sources%sourceParams(iParam)%str) then
       sourceParams= self%recipe%sources%sourceParams(iParam)%dsca
       exit
    end if
    if (iParam .eq. self%recipe%sources% n_sourceParams) then
        if (my_rank.eq.0) then
          print *, " ----------------- WRONG CORAL.EQUATION FILE AT RUNTIME ----------------"
        print *, 'In user-defined sources, LP_user_sources_definitions.f90, '
        print *, 'the parameter name did not match any definition from input file coral.equations'
          print *, " ----------------- WRONG CORAL.EQUATION FILE AT RUNTIME ----------------"
        end if
        call MPI_finalize(ierr)
        stop
    end if
 end do
 end function

 subroutine allocate_sources(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   allocate( self%sources (self%numberOf_sources) )
 end subroutine

 subroutine add_source_dscalar(self, sourceIndex, definition)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: sourceIndex
   real(dp), intent(in) :: definition
   if (sourceIndex .gt. self%numberOf_sources) then
      print *, 'In PL_user_sources_definitions.f90, the position of the'
      print *, 'source is greater than the total number of sources!'
      error stop
   end if
   allocate( self%sources(sourceIndex)%physical ( self%geometry%NZAA) )
   self%sources(sourceIndex)%physical = definition
 end subroutine

 subroutine add_source_array(self, sourceIndex, definition)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: sourceIndex
   real(dp), allocatable, intent(in) :: definition(:)
   if (sourceIndex .gt. self%numberOf_sources) then
      print *, 'In PL_user_sources_definitions.f90, the position of the'
      print *, 'source is greater than the total number of sources!'
      error stop
   end if
   allocate( self%sources(sourceIndex)%physical ( self%geometry%NZAA) )
   self%sources(sourceIndex)%physical = definition
 end subroutine


 subroutine transform_sources(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iSource
   do iSource = 1, self%numberOf_sources
   call Compute_1dCheby_transform(&
                  self%sources(iSource)%physical,&
                  self%sources(iSource)%spectral, self%geometry%NZAA) 
   end do
 end subroutine

 subroutine init_penalisation(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: ix, iy, iz
   real(dp), allocatable, dimension(:) :: penalisation_mask
   integer :: N_masks
   integer :: i_mask, iVar
   real(dp) :: pm, pi
   real(dp) :: x,y,z, lx, ly, zgap, zcenter

   call gauss_chebyshev_grid_1d( self%geometry%zGrid, &
                                 self%geometry%NZAA,&
                                 self%geometry%center,&
                                 self%geometry%gap)

   lx = self%geometry% Lx
   ly = self%geometry% Ly
   zgap = self%geometry% gap
   zcenter = self%geometry% center

   do iVar = 1, self%recipe% numberOf_linear_variables_full
      if (self%recipe% linear_vars_full (iVar) % penalisation) then
         do iy = 1, size(self%geometry% yGrid, 1) 
         do ix = 1, size(self%geometry% xGrid, 1) 
         do iz = 1, self%geometry% NZAA
          x = self%geometry% xGrid(ix)
          y = self%geometry% yGrid(iy)
          z = self%geometry% zGrid(iz)

          i_mask = 0 
          include "PL_user_penalisation_mask.f90"

          pm = maxval(penalisation_mask)
          pm = dtanh ( self%recipe% linear_vars_full(iVar)% penalisation_width* pm )
          pm = ( 1._dp + pm ) / 2._dp
          pm = self% recipe% linear_vars_full(iVar)% penalisation_strength * pm
        
          self%linear_variables(iVar)% phys(iz,ix,iy) = pm
         end do
         end do
         end do
         ! call self%linear_variables(iVar)% phys_to_spec()
      end if
   end do

 end subroutine init_penalisation

 !> @brief  
 !> Initialisation of the spectral coefficients array of volumic heat source.
 !> @details 
 !> If necessary: allocation, definition on the Gauss-Chebyshev grid, then
 !> 2/3 dealiasing of heat_source_cheb 
 Subroutine init_all_sources(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   real(kind=dp), allocatable :: z(:)
   real(kind=dp), allocatable :: sourceDef(:)

   call gauss_chebyshev_grid_1d( z, &
                                 self%geometry%NZAA,&
                                 self%geometry%center,&
                                 self%geometry%gap)
   allocate (sourceDef (self%geometry%NZAA) )

   include "PL_user_sources_definitions.f90"

   call self%transform_sources()
 End Subroutine init_all_sources

  Subroutine Compute_1dCheby_transform(arr_signal, arr_Coefs, N)
   Real(Kind=dp), Dimension(:), Allocatable, Intent(InOut):: arr_signal
   Real(Kind=dp), Dimension(:), Allocatable, Intent(InOut):: arr_coefs
   Integer, Intent(In) :: N
   Real(Kind=dp), pointer :: coefs(:)
   Real(kind=dp), Pointer :: signal(:) !we need a copy, aligned in memory
   !FFTW variables
   Type(C_ptr) :: plan_real2cos
   Type(C_ptr) :: pdum1, pdum2
   pdum1 = fftw_alloc_real(int(N, kind=8))
   Call C_F_Pointer(pdum1, coefs, [N])
   pdum2 = fftw_alloc_real(int(N, kind=8))
   Call C_F_Pointer(pdum2, signal, [N])
   plan_real2cos = fftw_plan_r2r_1d(Int(N, kind=4), signal, coefs,&
                                    FFTW_REDFT10, FFTW_ESTIMATE)
   signal = arr_signal
   Call fftw_execute_r2r(plan_real2cos, signal, coefs)
   arr_coefs = coefs / Real(N, kind=dp)
   arr_coefs(1) = arr_coefs(1)*0.5
   Call fftw_free(pdum1)
   Call fftw_free(pdum2)
 End Subroutine 

 
