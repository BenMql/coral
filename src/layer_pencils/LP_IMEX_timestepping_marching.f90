

 subroutine output_zero_fie_sample(self, time_integer)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: time_integer
   real(kind=dp), allocatable :: chebyCoefs_zeroMode(:)
   character(len=46) :: fileName
   integer :: iSys, iVar
   if (my_rank.eq.0) then
   1125 format ('chebyCoefs_fie0Smpl_sys',(i2.2),'_var',(i2.2),'_t',(i8.8),'.full')
   allocate ( chebyCoefs_zeroMode (self%geometry%NZAA) )
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
   do iVar = 1, self%recipe%zero_recipes(iSys)%n_coupled_vars
   chebyCoefs_zeroMode =  0._dp
   call self%coupled_zero_set(iSys)%padStencilExtractShuffle(iVar)%dot(&
        self%coupled_zero_set(iSys)%field, chebyCoefs_zeroMode, 'cumul')
   write (fileName, 1125) iSys, iVar, time_integer
   Open (Unit=9, File=fileName, Status='replace', Access='stream')
   Write(9) chebyCoefs_zeroMode(:)
   Close (Unit=9)
   end do
   end do
   end if
 end subroutine


 subroutine output_zero_aux_sample(self, time_integer)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: time_integer
   real(kind=dp), allocatable :: chebyCoefs_zeroMode(:)
   character(len=46) :: fileName
   integer :: iSys, iVar
   if (my_rank.eq.0) then
   1124 format ('chebyCoefs_aux0Smpl_sys',(i2.2),'_var',(i2.2),'_t',(i8.8),'.full')
   allocate ( chebyCoefs_zeroMode (self%geometry%NZAA) )
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
   do iVar = 1, self%recipe%zero_recipes(iSys)%n_coupled_vars
   chebyCoefs_zeroMode =  0._dp
   call self%coupled_zero_set(iSys)%padStencilExtractShuffle(iVar)%dot(&
        self%coupled_zero_set(iSys)%aux,   chebyCoefs_zeroMode, 'cumul')
   write (fileName, 1124) iSys, iVar, time_integer
   Open (Unit=9, File=fileName, Status='replace', Access='stream')
   Write(9) chebyCoefs_zeroMode(:)
   Close (Unit=9)
   end do
   end do
   end if
 end subroutine



 subroutine output_zero_rhs_sample(self, time_integer)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: time_integer
   character(len=46) :: fileName
   integer :: iSys, iVar
   if (my_rank.eq.0) then
   1123 format ('chebyCoefs_rhs0Smpl_sys',(i2.2),'_var',(i2.2),'_t',(i8.8),'.full')
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
   do iVar = 1, self%recipe%zero_recipes(iSys)%n_coupled_vars
   write (fileName, 1123) iSys, iVar, time_integer
   Open (Unit=9, File=fileName, Status='replace', Access='stream')
   Write(9) self%coupled_zero_set(iSys)%rhs
   Close (Unit=9)
   end do
   end do
   end if
 end subroutine
 subroutine output_kxky_fie_sample(self, time_integer)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: time_integer
   complex(kind=dp), allocatable :: chebyCoefs_kxkyMode(:,:,:)
   character(len=46) :: fileName
   integer :: iSys, iVar
   if (my_rank.eq.0) then
   1115 format ('chebyCoefs_fie_Smpl_sys',(i2.2),'_var',(i2.2),'_t',(i8.8),'.full')
   allocate ( chebyCoefs_kxkyMode (self%geometry%NZAA, self%geometry%spec%local_NY, self%geometry%spec%local_NX) )
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
   do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
   chebyCoefs_kxkyMode = cmplx( 0._dp, 0._dp, kind=dp)
   call self%coupled_kxky_set(iSys)%padStencilExtractShuffle(iVar)%dot(&
        self%coupled_kxky_set(iSys)%field, chebyCoefs_kxkyMode, 'cumul')
   write (fileName, 1115) iSys, iVar, time_integer
   Open (Unit=9, File=fileName, Status='replace', Access='stream')
   Write(9) chebyCoefs_kxkyMode(:,1,1)
   Close (Unit=9)
   end do
   end do
   end if
 end subroutine
 subroutine output_kxky_aux_sample(self, time_integer)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: time_integer
   complex(kind=dp), allocatable :: chebyCoefs_kxkyMode(:,:,:)
   character(len=46) :: fileName
   integer :: iSys, iVar
   if (my_rank.eq.0) then
   1113 format ('chebyCoefs_aux_Smpl_sys',(i2.2),'_var',(i2.2),'_t',(i8.8),'.full')
   allocate ( chebyCoefs_kxkyMode (self%geometry%NZAA, self%geometry%spec%local_NY, self%geometry%spec%local_NX) )
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
   do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
   chebyCoefs_kxkyMode = cmplx( 0._dp, 0._dp, kind=dp)
   call self%coupled_kxky_set(iSys)%padStencilExtractShuffle(iVar)%dot(&
        self%coupled_kxky_set(iSys)%aux, chebyCoefs_kxkyMode, 'cumul')
   write (fileName, 1113) iSys, iVar, time_integer
   open (Unit=9, File=fileName, Status='replace', Access='stream')
   write(9) chebyCoefs_kxkyMode(:,1,1)
   close (Unit=9)
   end do
   end do
   end if
 end subroutine

 subroutine output_kxky_rhs_sample(self, time_integer)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: time_integer
   complex(kind=dp), allocatable :: chebyCoefs_kxkyMode(:,:,:)
   character(len=46) :: fileName
   integer :: iSys, iVar
   if (my_rank.eq.0) then
   1112 format ('chebyCoefs_rhs_Smpl_sys',(i2.2),'_var',(i2.2),'_t',(i8.8),'.full')
   allocate ( chebyCoefs_kxkyMode (self%geometry%NZAA, self%geometry%spec%local_NY, self%geometry%spec%local_NX) )
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
   do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
   chebyCoefs_kxkyMode = cmplx( 0._dp, 0._dp, kind=dp)
   call self%coupled_kxky_set(iSys)%padStencilExtractShuffle(iVar)%dot(&
        self%coupled_kxky_set(iSys)%rhs, chebyCoefs_kxkyMode, 'cumul')
   write (fileName, 1112) iSys, iVar, time_integer
   Open (Unit=9, File=fileName, Status='replace', Access='stream')
   Write(9) chebyCoefs_kxkyMode(:,1,1)
   Close (Unit=9)
   end do
   end do
   end if
 end subroutine


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
      self%coupled_kxky_set(isys)%aux = self%coupled_kxky_set(isys)%rhs
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
   u_over_dx_max = maxVal(self%linear_variables(1)%phys)
   u_over_dx_max = u_over_dx_max * self%geometry%NXAA / self%geometry%Lx
   v_over_dy_max = maxVal(self%linear_variables(2)%phys)
   v_over_dy_max = v_over_dy_max * self%geometry%NYAA / self%geometry%Ly
   w_over_dz_max = v_over_dy_max
   do ix = 1, self%geometry%phys%local_NX
   do iy = 1, self%geometry%phys%local_NY
   w_over_dz_max = dMax1(w_over_dz_max,&
                         maxVal( self%linear_variables(3)%phys(:,iy,ix) / &
                                 local_cheby_weight) )
   end do
   end do
   
   !if (my_rank.eq.0) print *,  u_over_dx_max,& !delme
    !                                        v_over_dy_max,& !delme
     !                                       w_over_dz_max !delme

   self%cargo%cfl_based_DT = 1._dp/ dMax1 ( u_over_dx_max,&
                                            v_over_dy_max,&
                                            w_over_dz_max)

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
      ! the recipe is self%recipe%nl_vars(iQvar
      ! largeArrayWarning 
      self%quadratic_variables(iQvar)%phys = &
              self%linear_variables( self%recipe%nl_vars(iQvar)%iVar1 )%phys *&
              self%linear_variables( self%recipe%nl_vars(iQvar)%iVar2 )%phys
      call self%quadratic_variables(iQvar)%phys_to_spec()
   end do

   if (i.eq.1) call self%output_global_quantities()
   if (i.eq.1) call self%output_slices_volumes_and_profiles()
  
   allocate( nl_buffer( self%geometry%NZAA,&
                        self%geometry%spec%local_NY, &
                        self%geometry%spec%local_NX) )

   allocate( nl_buffer_deAliased( &
                        self%geometry%NZ,&
                        self%geometry%spec%local_NY, &
                        self%geometry%spec%local_NX) )


   ! => for each set of coupled eqns: linear combinations of 
   ! self%quadratic_variables(iQvar)%spec are computed 
   ! and stored in the corresponding K_hat array
   do isys = 1, self%recipe%numberOf_coupled_kxky_systems
   self%coupled_kxky_set(isys)%step(i)%K_hat = cmplx(0._dp, 0._dp, kind=dp)
   do iTerm = 1, self%recipe%kxky_recipes(iSys)%NL%n_pieces
      do ix = 1, self%geometry%spec%local_NX
      do iy = 1, self%geometry%spec%local_NY
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
      !if (my_rank.eq.0) print *, self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dsca, &
              !self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dx_exponent, &
              !self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%dy_exponent, &
              !self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%Iz_exponent, &
              !self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%quad_var_index, &
              !self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%eqn_index, &
              !self%recipe%kxky_recipes(iSys)%eqn_order( insert_in_equation )+1
      call self%coupled_kxky_set(iSys)%building_tools%cheb_IM ( &
                        (self%recipe%kxky_recipes(iSys)%NL%term(iTerm)%Iz_exponent+1), 1 ) &
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
   do ix = 1, self%geometry%spec%local_NX
   do iy = 1, self%geometry%spec%local_NY
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
   allocate( nl_buffer_zero( self%geometry%NZAA )) 
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
                       (self%recipe%zero_recipes(iSys)%NL%term(iTerm)%Iz_exponent+1), 1 ) &
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
     nl_buffer_zero(:) = self%sources( self%recipe%sources%term(iSource)%source_index) %spectral
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
   integer :: iVar
   integer :: iTerm
   integer :: ix, iy
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
          call self%coupled_zero_set(system_of_interest &
                           )%padStencilExtractShuffle(position_of_interest&
                           )%dot(self%coupled_zero_set(system_of_interest)%aux,&
                                 self%linear_variables(iVar)%spec, 'cumul')
          end if
     end select
   end do
   ! /// done: spectral space computations
   ! ==================================================

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
   ! todo => take into account the recipe if necessary
   ! In particular, we have assumed that the prefactor is 1., and that 
   ! no derivatives are present.
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

   ! ==================================================
   ! >>> spectral to physical tranforms      
   call self%linear_variables(iVar)%spec_to_phys()
   ! /// done: transforms                     
   ! ==================================================
   end do
   deAllocate( datBuffer )
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
 do iParam = 1, size(self%recipe%sources%sourceParams)
    if (paramString .eq. self%recipe%sources%sourceParams(iParam)%str) then
       sourceParams= self%recipe%sources%sourceParams(iParam)%dsca
       exit
    end if
    if (iParam .eq. size(self%recipe%sources%sourceParams)) then
        print *, 'In user-defined sources, LP_user_sources_definitions.f90, '
        print *, 'the parameter name did not match any definition from input file coral.equations'
        error stop
    end if
 end do
 end function

 subroutine allocate_sources(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   allocate( self%sources (self%numberOf_sources) )
 end subroutine

 subroutine add_source(self, sourceIndex, definition)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: sourceIndex
   real(dp), allocatable, intent(in) :: definition(:)
   if (sourceIndex .gt. self%numberOf_sources) then
      print *, 'In LP_user_sources_definitions.f90, the position of the'
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

   include "LP_user_sources_definitions.f90"

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

 
