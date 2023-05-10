



 subroutine output_global_quantities(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   real(kind=dp) :: varAvg, aux_dsca
   integer :: position_backup
   integer :: current_position
   logical :: first_or_not
   integer :: iVar, iObj
   integer :: ix, iy, iz
   character(len=17) :: fileExtension
   305 format ('_XYavg_z',(i5.5),'.dat')
   if (self%io_bookkeeping%first_or_not) then
   call gauss_chebyshev_weight_1d( self%gauss_cheby%weight1d, &
                                   self%geometry%NZAA,&
                                   self%geometry%gap)
   end if
   
   position_backup = self%io_bookkeeping%dPosition
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   self%io_bookkeeping%dPosition = position_backup
   current_position = (timings%i_timestep-1)*8+1
   first_or_not = .False.
   if (current_position.eq.1) first_or_not = .True.
   call output_dsca_in_unique_timeserie(&
                 self%io_bookkeeping%output_directory,&
                 self%io_bookkeeping%output_dir_length,&
                 'time.dat', 8,& 
                 self%io_bookkeeping%absolute_time,&
                 first_or_not, &                      
                 current_position)
      
   do iVar = 1, self%recipe%numberOf_linear_variables_full
   do iObj = 1, self%recipe%timeseries%numberOf_linearObjects( iVar )
      select case (self%recipe%timeseries%linear(iVar)%object(iObj)%kind)
          case ('volume')
                varAvg  = 0._dp                                                       
                do ix = 1, self%geometry%phys%local_NX
                do iy = 1, self%geometry%phys%local_NY
                do iz = 1, self%geometry%phys%local_NZ
                   varAvg = varAvg + self%linear_variables(iVar)%phys(iz,iy,ix) &
                         * self%gauss_cheby%weight1d( domain_decomp%phys_iStart(1) -1 + iz)
                end do
                end do
                end do
                varAvg = varAvg / real(self%geometry%NXAA, kind=dp)
                varAvg = varAvg / real(self%geometry%NYAA, kind=dp)
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                self%io_bookkeeping%dPosition = position_backup
                current_position = (timings%i_timestep-1)*8+1
   first_or_not = .False.
   if (current_position.eq.1) first_or_not = .True.
                call output_cumul_d_in_timeserie(&
                              self%io_bookkeeping%output_directory,&
                              self%io_bookkeeping%output_dir_length,&
                              self%recipe%linear_vars_full(iVar)%str,&
                              len(self%recipe%linear_vars_full(iVar)%str),&
                              varAvg,&
                              first_or_not, &                      
                              current_position)            
          case ('zSlice')
               if ((self%recipe%timeseries%linear(iVar)%object(iObj)%slice_index &
                               .ge. domain_decomp%phys_iStart(1)) &
                                       .and. &
                   (self%recipe%timeseries%linear(iVar)%object(iObj)%slice_index &
                               .le. domain_decomp%phys_iEnd  (1))) then
                varAvg  = 0._dp                                                       
                do ix = 1, self%geometry%phys%local_NX
                do iy = 1, self%geometry%phys%local_NY
                iz = self%recipe%timeseries%linear(iVar)%object(iObj)%slice_index - domain_decomp%phys_iStart(1) + 1
                varAvg = varAvg + self%linear_variables(iVar)%phys(iz,iy,ix) 
                end do
                end do
                varAvg = varAvg / real(self%geometry%NXAA, kind=dp)
                varAvg = varAvg / real(self%geometry%NYAA, kind=dp)
                call MPI_reduce(varAvg, &  !send buffer
                   aux_dsca,   & !recv buffer
                   1,          & !count
                   MPI_double, & !dtype
                   MPI_sum,    &
                   0, &
                   self%geometry%mpi_Zphys%comm, ierr)
                 if (self%geometry%mpi_Zphys%rank.eq.0) then
                 write (fileExtension, 305) self%recipe%timeseries%linear(iVar)%object(iObj)%slice_index 
                 current_position = (timings%i_timestep-1)*8+1
   first_or_not = .False.
   if (current_position.eq.1) first_or_not = .True.
                 call output_dsca_in_unique_timeserie(&
                    self%io_bookkeeping%output_directory,&
                    self%io_bookkeeping%output_dir_length,&
                    self%recipe%linear_vars_full(iVar)%str//fileExtension,&
                    len(self%recipe%linear_vars_full(iVar)%str)+17,&
                    aux_dsca, &                                 
                    first_or_not, &                      
                    current_position)            
                 end if 
               end if
          case default
               write(*,*) 'bad (L) object kind', self%recipe%timeseries%linear(iVar)%object(iObj)%kind, iVar, iObj
               error stop
          end select
   end do
   end do
   ! initialise the KE
   self%cargo%KE_display = 0._dp
   do iVar = 1, self%recipe%numberOf_quadratic_variables
   do iObj = 1, self%recipe%timeseries%numberOf_quadraObjects( iVar )
      select case (self%recipe%timeseries%quadra(iVar)%object(iObj)%kind)
          case ('volume')
                varAvg  = 0._dp                                                       
                do ix = 1, self%geometry%phys%local_NX
                do iy = 1, self%geometry%phys%local_NY
                do iz = 1, self%geometry%phys%local_NZ
                   varAvg = varAvg + self%quadratic_variables(iVar)%phys(iz,iy,ix) &
                         * self%gauss_cheby%weight1d( domain_decomp%phys_iStart(1) -1 + iz)
                end do
                end do
                end do
                varAvg = varAvg / real(self%geometry%NXAA, kind=dp)
                varAvg = varAvg / real(self%geometry%NYAA, kind=dp)
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                self%io_bookkeeping%dPosition = position_backup
                current_position = (timings%i_timestep-1)*8+1
   first_or_not = .False.
   if (current_position.eq.1) first_or_not = .True.
                call output_cumul_d_in_timeserie(&
                              self%io_bookkeeping%output_directory,&
                              self%io_bookkeeping%output_dir_length,&
                              self%recipe%nl_vars(iVar)%str,&
                              len(self%recipe%nl_vars(iVar)%str),&
                              varAvg,&
                              first_or_not, &                      
                              current_position)            
                if ((self%recipe%nl_vars(iVar)%str .eq. 'uu')   .or. &
                    (self%recipe%nl_vars(iVar)%str .eq. 'vv')   .or. &
                    (self%recipe%nl_vars(iVar)%str .eq. 'ww') ) then  
                        self%cargo%KE_display = self%cargo%KE_display + varAvg
                end if
          case ('zSlice')
               if ((self%recipe%timeseries%quadra(iVar)%object(iObj)%slice_index &
                               .ge. domain_decomp%phys_iStart(1)) &
                                       .and. &
                   (self%recipe%timeseries%quadra(iVar)%object(iObj)%slice_index &
                               .le. domain_decomp%phys_iEnd  (1))) then
                varAvg  = 0._dp                                                       
                do ix = 1, self%geometry%phys%local_NX
                do iy = 1, self%geometry%phys%local_NY
                iz = self%recipe%timeseries%quadra(iVar)%object(iObj)%slice_index - domain_decomp%phys_iStart(1) + 1
                varAvg = varAvg + self%quadratic_variables(iVar)%phys(iz,iy,ix) 
                end do
                end do
                varAvg = varAvg / real(self%geometry%NXAA, kind=dp)
                varAvg = varAvg / real(self%geometry%NYAA, kind=dp)
                call MPI_reduce(varAvg, &  !send buffer
                   aux_dsca,   & !recv buffer
                   1,          & !count
                   MPI_double, & !dtype
                   MPI_sum,    &
                   0, &
                   self%geometry%mpi_Zphys%comm, ierr)
                 if (self%geometry%mpi_Zphys%rank.eq.0) then
                 write (fileExtension, 305) self%recipe%timeseries%quadra(iVar)%object(iObj)%slice_index 
                 current_position = (timings%i_timestep-1)*8+1
   first_or_not = .False.
   if (current_position.eq.1) first_or_not = .True.
                 call output_dsca_in_unique_timeserie(&
                    self%io_bookkeeping%output_directory,&
                    self%io_bookkeeping%output_dir_length,&
                    self%recipe%nl_vars(iVar)%str//fileExtension,&
                    len(self%recipe%nl_vars(iVar)%str)+17,&
                    aux_dsca, &                                 
                    first_or_not, &                      
                    current_position)            
                 end if 
               end if
          case default
               print *,'bad (Q) object kind', self%recipe%timeseries%quadra(iVar)%object(iObj)%kind, iVar, iObj
               error stop
          end select
   end do
   end do
   self%io_bookkeeping%first_or_not=.False.

 end subroutine

   
 subroutine export_allPhys(self, iTime)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar
   integer, optional, intent(in) :: iTime
   integer :: myTime
   if (.not.present(iTime)) then
           myTime = 0
   else 
           myTime=iTime
   end if

   do iVar = 1, size( self%linear_variables)
      call self%export_volume(iVar, myTime, 'linear')
   end do
 end subroutine

 subroutine export_volume(self, ivar, iTime, linear_or_quadra)
   class(full_problem_data_structure_T), intent(inOut) :: self
   character(len=6), intent(in) :: linear_or_quadra
   integer, intent(in) :: iVar
   integer, intent(in) :: iTime
   character(len=44) :: fileName
   call chdir(self%io_bookkeeping%output_directory)
   405 format ('./Volumes/linear_var',(i2.2),'_time',(i8.8),'_full.dat')
   415 format ('./Volumes/quadra_var',(i2.2),'_time',(i8.8),'_full.dat')
   select case(linear_or_quadra)
          case ('linear')
          write (fileName,405) iVar, iTime
          call self%linear_variables(iVar)% write_phys_to_disk (fileName)
          case ('quadra')
          write (fileName,415) iVar, iTime
          call self%quadratic_variables(iVar)% write_phys_to_disk (fileName)
   end select
 end subroutine

 subroutine export_slice (self, ivar, iTime, slice_kind, slice_index, linear_or_quadra)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: slice_kind
   integer :: slice_index
   character(len=6), intent(in) :: linear_or_quadra
   integer, intent(in) :: iVar
   integer, intent(in) :: iTime
   character(len=50) :: fileName
   call chdir(self%io_bookkeeping%output_directory)
   445 format ('./Slices/linear_var',(i2.2),'_z',(i5.5),'_time',(i8.8),'_full.dat')
   446 format ('./Slices/linear_var',(i2.2),'_y',(i5.5),'_time',(i8.8),'_full.dat')
   447 format ('./Slices/linear_var',(i2.2),'_x',(i5.5),'_time',(i8.8),'_full.dat')
   455 format ('./Slices/quadra_var',(i2.2),'_z',(i5.5),'_time',(i8.8),'_full.dat')
   456 format ('./Slices/quadra_var',(i2.2),'_y',(i5.5),'_time',(i8.8),'_full.dat')
   457 format ('./Slices/quadra_var',(i2.2),'_x',(i5.5),'_time',(i8.8),'_full.dat')
   select case(linear_or_quadra)
          case ('linear')
          select case (slice_kind)
                 case (1)
                 write (fileName, 445) iVar, slice_index, iTime
                 case (2)
                 write (fileName, 446) iVar, slice_index, iTime
                 case (3)
                 write (fileName, 447) iVar, slice_index, iTime
          end select
          call self%linear_variables(iVar)% slice_phys_to_disk (slice_kind, slice_index, fileName)
          case ('quadra')
          select case (slice_kind)
                 case (1)
                 write (fileName, 455) iVar, slice_index, iTime
                 case (2)
                 write (fileName, 456) iVar, slice_index, iTime
                 case (3)
                 write (fileName, 457) iVar, slice_index, iTime
          end select
          call self%quadratic_variables(iVar)% slice_phys_to_disk (slice_kind, slice_index, fileName)
   end select
 end subroutine


 subroutine import_quicksave(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar, iSys
   character(len=58) :: fileName
   integer :: rolInt
   logical :: file_exists
   real(kind=dp), allocatable :: buff (:)
   if (my_rank.eq.0) print*, " ================================================================="
   if (my_rank.eq.0) print*, " >>> "
   if (my_rank.eq.0) print*, " >>> Reading a quicksave "
   if (my_rank.eq.0) print*, " >>> "
   if (my_rank.eq.0) print*, " ================================================================="
   call chdir(self%io_bookkeeping%output_directory)
   !check whether rolInt=1 or 2
   inquire (file="./Restart/QuickSave.kxky.sys001.var001.rolling1", exist=file_exists)
   if (file_exists) then
      rolInt=1
   else
      inquire (file="./Restart/QuickSave.kxky.sys001.var001.rolling2", exist=file_exists)
      if (file_exists) then
         rolInt=2
      else
         print *, "No QuickSave was found in the ./Restart directory"
         error stop
      end if
   end if
   506 format ('./Restart/QuickSave.kxky.sys',(i3.3),'.var',(i3.3),'.rolling',(i1.1))
   507 format ('./Restart/QuickSave.zero.sys',(i3.3),'.var',(i3.3),'.rolling',(i1.1))
   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
   self%coupled_kxky_set(isys)%field = cmplx(0._dp, 0._dp, kind=dp) ! it has been backed-up before
   do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
   ! copy spectral coefficients into the first array we can think of
   self%linear_variables(1)%spec = cmplx(0._dp, 0._dp, kind=dp)
   write (fileName,506) iSys, iVar, rolInt
   call self%linear_variables(1)% read_phys_from_disk (fileName)
   ! now tranform back to Galerkin coefficients, stored in 'field'
   call self%linear_variables(1)%phys_to_spec()
   ! backsolve for galerkin
   call self%coupled_kxky_set(iSys)%square_stencil(iVar)%backsolve(&
                  self%linear_variables(1)%spec, &
                  self%geometry%spec%local_NY*self%geometry%spec%local_NX, &
                  domain_decomp% spec_iSize(1))
   ! truncate, shuffle and insert at the right location
   call self%coupled_kxky_set(iSys)%shuffleTextractTtruncate(iVar)%dot(&
                  self%linear_variables(1)%spec, &
                  self%coupled_kxky_set(isys)%field,&
                  'cumul') 
   end do
   end do

   if (self%geometry%this_core_has_zero_mode) then
   allocate(buff( self%geometry%NZAA))
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
   self%coupled_zero_set(isys)%field = 0._dp
   do iVar = 1, self%recipe%zero_recipes(iSys)%n_Coupled_vars
      !call buff to phys()
      write (fileName,507) iSys, iVar, rolInt
      open (unit=9, file=fileName, status='old', access='stream')
      read(9) buff
      close(unit=9)
      call r2r_forward(buff)
      ! backsolve for galerkin
      call self%coupled_zero_set(iSys)%square_stencil(iVar)%backsolve(buff)
      ! truncate, shuffle and insert at the right location
      call self%coupled_zero_set(iSys)%shuffleTextractTtruncate(iVar)%dot(&
                  buff, &                                 
                  self%coupled_zero_set(isys)%field,&
                  'cumul')
      
   end do
   end do
   if (allocated(buff)) deAllocate(buff)
   end if

 end subroutine

 subroutine export_checkpoints(self)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar, iSys
   character(len=58) :: fileName
   real(kind=dp), allocatable :: buff (:)

   ! back-up fields
   call self%copy_fields_to_aux()

   call chdir(self%io_bookkeeping%output_directory)
   406 format ('./CheckPoints/QuickSave.kxky.sys',(i3.3),'.var',(i3.3),'.rolling',(i1.1))
   407 format ('./CheckPoints/QuickSave.zero.sys',(i3.3),'.var',(i3.3),'.rolling',(i1.1))

   do iSys = 1, self%recipe%numberOf_coupled_kxky_systems
   self%coupled_kxky_set(isys)%field = cmplx(0._dp, 0._dp, kind=dp) ! it has been backed-up before
   do iVar = 1, self%recipe%kxky_recipes(iSys)%n_coupled_vars
   ! compute Chebyshev coefficients into self%linear_variables(1)%spec
   self%linear_variables(1)%spec = cmplx(0._dp, 0._dp, kind=dp)
   call self%coupled_kxky_set(iSys) &
                      %padStencilExtractShuffle(iVar) &
                      %dot(self%coupled_kxky_set(iSys)%aux,&
                        self%linear_variables(1)%spec, 'overW')
   write (fileName,406) iSys, iVar, self%io_bookkeeping%rolling_integer
   call self%linear_variables(1)%spec_to_phys()
   call self%linear_variables(1)% write_phys_to_disk( filename )
   ! now tranform back to Galerkin coefficients, stored in 'field'
   call self%linear_variables(1)%phys_to_spec()
   ! backsolve for galerkin
   call self%coupled_kxky_set(iSys)%square_stencil(iVar)%backsolve(&
                  self%linear_variables(1)%spec, &
                  self%geometry%spec%local_NY*self%geometry%spec%local_NX, &
                  domain_decomp% spec_iSize(1))
   ! truncate, shuffle and insert at the right location
   call self%coupled_kxky_set(iSys)%shuffleTextractTtruncate(iVar)%dot(&
                  self%linear_variables(1)%spec, &
                  self%coupled_kxky_set(isys)%field,&
                  'cumul')
   end do
   end do

   if (self%geometry%this_core_has_zero_mode) then
           allocate(buff( domain_decomp% NZAA))
   !self%coupled_zero_set(isys)%field = 0._dp ! it has been backed-up before
   do iSys = 1, self%recipe%numberOf_coupled_zero_systems
   do iVar = 1, self%recipe%zero_recipes(iSys)%n_Coupled_vars
      ! compute Chebyshev coefficients into self%linear_variables(1)%spec
      buff = 0._dp
      call self%coupled_zero_set(iSys) &
                      %padStencilExtractShuffle(iVar) &
                      %dot(self%coupled_zero_set(iSys)%aux,&
                        buff, 'overW')
      call r2r_backward(buff)
      write (fileName,407) iSys, iVar, self%io_bookkeeping%rolling_integer
      open (unit=9, file=fileName, status='replace', access='stream')
      write(9) buff
      close(unit=9)

   end do
   end do
   deAllocate(buff)
   end if

   if (misc_cargo%rolling_qsaves) then
   if (self%io_bookkeeping%rolling_integer.eq.1) then
       self%io_bookkeeping%rolling_integer = 2
   else                                                
       self%io_bookkeeping%rolling_integer = 1
   end if
   end if
        
   
 end subroutine
   
 subroutine export_allProfiles(self, iTime)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer :: iVar
   integer, optional, intent(in) :: iTime
   integer :: myTime
   if (.not.present(iTime)) then
           myTime = 0
   else 
           myTime=iTime
   end if

   do iVar = 1, size( self%linear_variables)
      call self%export_profile(iVar, myTime, 'linear')
   end do
 end subroutine

 subroutine export_Profile(self, ivar, iTime, linear_or_quadra)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iVar
   integer, intent(in) :: iTime
   character(len=6), intent(in) :: linear_or_quadra
   select case (parallel_transforms_library)
          case ('fftw-mpi')
          call export_Profile_fftw_mpi( &
                               self, ivar, iTime, linear_or_quadra)
          case ('decomp2d')
          call export_Profile_decomp2d( &
                               self, ivar, iTime, linear_or_quadra)
   end select
 end subroutine

 subroutine export_Profile_fftw_mpi(self, ivar, iTime, linear_or_quadra)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iVar
   integer, intent(in) :: iTime
   character(len=6), intent(in) :: linear_or_quadra
   character(len=45) :: fileName
   real(kind=dp), allocatable :: local_profile(:)
   real(kind=dp), allocatable :: global_profile(:)
   !----------------------
   ! Internal variables
   !----------------------
   Integer(kind=MPI_Offset_Kind) :: displacement
                     !< Convert \p cumul_nElems to bytes
   Integer :: size_of_dble
                     !< size of a double precision number in bytes
   Integer :: file_id
                     !< mpiIO keeps tracks of what it's doing
   Integer, Dimension(:), Allocatable :: array_of_nelems
   Integer :: my_cumul
   Integer :: icore
   integer :: my_nElems

   allocate( local_profile  (self%geometry%phys%local_NZ) )
   allocate( global_profile (self%geometry%phys%local_NZ) )

   call chdir(self%io_bookkeeping%output_directory)
   409 format ('./Profiles/linear_var',(i2.2),'_time',(i8.8),'_full.dat')
   410 format ('./Profiles/quadra_var',(i2.2),'_time',(i8.8),'_full.dat')
   ! first, each process averages over x and y
   select case (linear_or_quadra)
          case ('linear')
          local_profile = sum(sum(self%linear_variables(iVar)%phys, 3),2)
          write (fileName,409) iVar, iTime
          case ('quadra')
          write (fileName,410) iVar, iTime
          local_profile = sum(sum(self%quadratic_variables(iVar)%phys, 3),2)
   end select
   local_profile = local_profile / self%geometry%NXAA                   
   local_profile = local_profile / self%geometry%NYAA                   
   ! averaged of distributed y
   call MPI_reduce(local_profile, &  !send buffer
                   global_profile, & !recv buffer
                   self%geometry%phys%local_NZ, & !count
                   MPI_double, &     !dtype
                   MPI_sum, &
                   0, &
                   self%geometry%mpi_Zphys%comm, ierr)
   ! now the z-distributed (but x,y-averaged) global_profiles are collectively written
   ! by the 0-rank processes of self%mpi_Yphys%comm
   ! // ugly, fixme one day
   if (my_rank .eq. 0 ) then 
      my_nElems    = self%geometry%phys%local_NZ
   else 
      my_nElems = 0
   end if
   !    ugly, fixme one day // 
   ! ~~~~~~~~~~~~~~~~~~~~
   ! compute cumul_nElems is missing, we compute it
   ! ~~~~~~~~~~~~~~~~~~~~
   Allocate(array_of_nElems(0:self%geometry%MPI_yphys%size-1))
   Call MPI_Gather(my_nElems, 1, MPI_Integer,&
          array_of_nElems, 1, MPI_integer, 0, self%geometry%MPI_yPhys%comm, ierr)
   Call MPI_Bcast(array_of_nElems, self%geometry%MPI_yphys%size, MPI_integer, 0, self%geometry%MPI_yPhys%comm, ierr)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! compute the number of elements owned by predecessors
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   my_cumul = 0
   Do icore = 0, self%geometry%MPI_yPhys%rank-1
     my_cumul = my_cumul + array_of_nelems(icore)
   End Do
   DeAllocate(array_of_nElems)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! compute the displacement
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Call MPI_Type_Size(MPI_Double, size_of_dble, ierr)
   displacement = my_cumul*size_of_dble

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
   !      Write the Profile           !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
   call MPI_File_Open(self%geometry%MPI_yPhys%comm, fileName,&
                      MPI_Mode_WROnly + MPI_Mode_Create,&
                      MPI_Info_Null, file_id, ierr)
   call MPI_File_set_view(file_id, displacement, MPI_Double, &
                          MPI_Double, 'native', &
                          MPI_Info_Null, ierr)
   call MPI_File_Write(file_id, global_profile, my_nElems, MPI_Double, &
                          MPI_Status_Ignore, Ierr)
   call MPI_File_Close(file_id, ierr)
 end subroutine




 subroutine export_Profile_decomp2d(self, iVar, iTime, linear_or_quadra)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iVar
   integer, intent(in) :: iTime
   character(len=6), intent(in) :: linear_or_quadra
   character(len=45) :: fileName
   real(kind=dp), allocatable :: local_profile(:)
   real(kind=dp), allocatable :: global_profile(:)
   !----------------------
   ! Internal variables
   !----------------------
   Integer(kind=MPI_Offset_Kind) :: displacement
                     !< Convert \p cumul_nElems to bytes
   Integer :: size_of_dble
                     !< size of a double precision number in bytes
   Integer :: file_id
                     !< mpiIO keeps tracks of what it's doing
   Integer, Dimension(:), Allocatable :: array_of_nelems
   Integer :: my_cumul
   Integer :: icore
   integer :: my_nElems

   allocate( local_profile  (self%geometry%phys%local_NZ) )
   allocate( global_profile (self%geometry%phys%local_NZ) )

   call chdir(self%io_bookkeeping%output_directory)
   409 format ('./Profiles/linear_var',(i2.2),'_time',(i8.8),'_full.dat')
   410 format ('./Profiles/quadra_var',(i2.2),'_time',(i8.8),'_full.dat')
   ! first, each process averages over x and y
   select case (linear_or_quadra)
          case ('linear')
          local_profile = sum(sum(self%linear_variables(iVar)%phys, 3),2)
          write (fileName,409) iVar, iTime
          case ('quadra')
          write (fileName,410) iVar, iTime
          local_profile = sum(sum(self%quadratic_variables(iVar)%phys, 3),2)
   end select
   local_profile = local_profile / self%geometry%NXAA                   
   local_profile = local_profile / self%geometry%NYAA                   
   ! averaged of distributed y
   call MPI_reduce(local_profile, &  !send buffer
                   global_profile, & !recv buffer
                   self%geometry%phys%local_NZ, & !count
                   MPI_double, &     !dtype
                   MPI_sum, &
                   0, &
                   self%geometry%mpi_Zphys%comm, ierr)
   ! now the z-distributed (but x,y-averaged) global_profiles are collectively written
   ! by the 0-rank processes of self%mpi_Yphys%comm
   my_nElems    = self%geometry%phys%local_NZ
   ! ~~~~~~~~~~~~~~~~~~~~
   ! compute cumul_nElems is missing, we compute it
   ! ~~~~~~~~~~~~~~~~~~~~
   Allocate(array_of_nElems(0:self%geometry%MPI_yphys%size-1))
   Call MPI_Gather(my_nElems, 1, MPI_Integer,&
          array_of_nElems, 1, MPI_integer, 0, self%geometry%MPI_yPhys%comm, ierr)
   Call MPI_Bcast(array_of_nElems, self%geometry%MPI_yphys%size, MPI_integer, 0, self%geometry%MPI_yPhys%comm, ierr)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! compute the number of elements owned by predecessors
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   my_cumul = 0
   Do icore = 0, self%geometry%MPI_yPhys%rank-1
     my_cumul = my_cumul + array_of_nelems(icore)
   End Do
   DeAllocate(array_of_nElems)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! compute the displacement
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Call MPI_Type_Size(MPI_Double, size_of_dble, ierr)
   displacement = my_cumul*size_of_dble

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
   !      Write the Profile           !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
   call MPI_File_Open(self%geometry%MPI_yPhys%comm, fileName,&
                      MPI_Mode_WROnly + MPI_Mode_Create,&
                      MPI_Info_Null, file_id, ierr)
   call MPI_File_set_view(file_id, displacement, MPI_Double, &
                          MPI_Double, 'native', &
                          MPI_Info_Null, ierr)
   call MPI_File_Write(file_id, global_profile, my_nElems, MPI_Double, &
                          MPI_Status_Ignore, Ierr)
   call MPI_File_Close(file_id, ierr)

 end subroutine
   
 subroutine export_verticallyAvgedSlice_general(self, ivar, iTime, linear_or_quadra)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iVar
   integer, intent(in) :: iTime
   character(len=6), intent(in) :: linear_or_quadra
   select case (parallel_transforms_library)
          case ('fftw-mpi')
          call export_verticallyAvgedSlice_fftw_mpi( &
                               self, ivar, iTime, linear_or_quadra)
          case ('decomp2d')
          call export_verticallyAvgedSlice_decomp2d( &
                               self, ivar, iTime, linear_or_quadra)
   end select
 end subroutine


 subroutine export_verticallyAvgedSlice_fftw_mpi(self, ivar, iTime, linear_or_quadra)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iVar
   integer, intent(in) :: iTime
   character(len=6), intent(in) :: linear_or_quadra
   character(len=50) :: fileName
   !----------------------
   !  Internal variables  
   !----------------------
   real(kind=dp), allocatable :: local_cheby_weight(:)
   ! ======================
   call gauss_chebyshev_weight_1d( self%gauss_cheby%weight1d, &
                                   self%geometry%NZAA,&
                                   self%geometry%gap)
   allocate( local_cheby_weight ( self%geometry%NZAA ))     ! superfluous? 
       !< could we use self%gauss_cheby% weight1d directly instead?        
   local_cheby_weight = self%gauss_cheby%weight1d (1 : self%geometry%NZAA) 
   call chdir(self%io_bookkeeping%output_directory)
   419 format ('./Slices/linear_var',(i2.2),'_zAvged_time',(i8.8),'_full.dat')
   418 format ('./Slices/quadra_var',(i2.2),'_zAvged_time',(i8.8),'_full.dat')
   select case(linear_or_quadra)
          case ('linear')
          write (fileName, 419) iVar, iTime
          call self%linear_variables(iVar)% zSummed_slice_phys_to_disk (local_cheby_weight, fileName)
          case ('quadra')
          write (fileName, 418) iVar, iTime
          call self%quadratic_variables(iVar)% zSummed_slice_phys_to_disk (local_cheby_weight, fileName)
   end select
 end subroutine
   


 subroutine export_verticallyAvgedSlice_decomp2d(self, ivar, iTime, linear_or_quadra)
   class(full_problem_data_structure_T), intent(inOut) :: self
   integer, intent(in) :: iVar
   integer, intent(in) :: iTime
   character(len=6), intent(in) :: linear_or_quadra
   character(len=41) :: fileName
   real(kind=dp), allocatable :: local_tile (:,:)
   real(kind=dp), allocatable :: global_tile(:,:)
   real(kind=dp), allocatable, target :: trposd_tile(:,:)
   !----------------------
   !  Internal variables  
   !----------------------
   Integer(kind=MPI_Offset_Kind) :: displacement
                     !< Convert \p cumul_nElems to bytes
   Integer :: size_of_dble
                     !< size of a double precision number in bytes
   Integer :: file_id
                     !< mpiIO keeps tracks of what it's doing
   Integer, Dimension(:), Allocatable :: array_of_nelems
   Integer :: my_cumul
   Integer :: icore
   integer :: my_nElems
   real(kind=dp), pointer :: fptr(:)        
   type(C_ptr) :: cptr
   real(kind=dp), allocatable :: local_cheby_weight(:)
   integer :: ix, iy, iz
   


   call gauss_chebyshev_weight_1d( self%gauss_cheby%weight1d, &
                                   self%geometry%NZAA,&
                                   self%geometry%gap)
   allocate( local_cheby_weight ( self%geometry%phys%local_NZ ))
   local_cheby_weight = self%gauss_cheby%weight1d (domain_decomp%phys_iStart(1) : &
                                                   domain_decomp%phys_iStart(1) + &
                                                   self%geometry%phys%local_NZ  -1)



   allocate( local_tile  (self%geometry%phys%local_NY, self%geometry%phys%local_NX) )
   allocate( global_tile (self%geometry%phys%local_NY, self%geometry%phys%local_NX) )
   allocate( trposd_tile (self%geometry%phys%local_NX, self%geometry%phys%local_NY) )

   241 format ('linear_var',(i2.2),'_zAvged_time',(i8.8),'_full.dat')
   242 format ('linear_var',(i2.2),'_zAvged_time',(i8.8),'_full.dat')
   select case (linear_or_quadra)
          case ('linear')
          write (fileName,241) iVar, iTime
          local_tile = 0._dp
          do ix = 1, self%geometry%phys%local_NX
          do iy = 1, self%geometry%phys%local_NY
          do iz = 1, self%geometry%phys%local_NZ
             local_tile (iy, ix) = local_tile(iy, ix) &
                                 + self%linear_variables(iVar)%phys(iz, iy, ix) &
                                 * local_cheby_weight(iz)
          end do 
          end do 
          end do 

          case ('quadra')
          write (fileName,242) iVar, iTime
          local_tile = 0._dp
          do ix = 1, self%geometry%phys%local_NX
          do iy = 1, self%geometry%phys%local_NY
          do iz = 1, self%geometry%phys%local_NZ
             local_tile (iy, ix) = local_tile(iy, ix) &
                                 + self%quadratic_variables(iVar)%phys(iz, iy, ix) &
                                 * local_cheby_weight(iz)
          end do 
          end do 
          end do 
   end select
   my_nElems    = self%geometry%phys%local_NY * self%geometry%phys%local_NX

   call mpi_reduce( local_tile, &  !send buffer
                    global_tile,&  !recv buffer
                    my_nElems,  &
                    mpi_double, &
                    mpi_sum,    &
                    0,          &
                    self%geometry%mpi_Yphys%comm, ierr )

   trposd_tile = transpose(global_tile)

   cptr = c_loc(trposd_tile)
   call C_F_pointer(cptr, fptr, [self%geometry%phys%local_NY * self%geometry%phys%local_NX])

   call chdir(self%io_bookkeeping%output_directory)
   call chdir('./Slices')
  
   ! ~~~~~~~~~~~~~~~~~~~~
   ! compute cumul_nElems is missing, we compute it
   ! ~~~~~~~~~~~~~~~~~~~~
   Allocate(array_of_nElems(0:self%geometry%mpi_Zphys%size-1))
   call MPI_allGather(self%geometry%phys%local_NY * self%geometry%phys%local_NX, &
                     1, MPI_Integer,&
          array_of_nElems, 1, MPI_integer,  self%geometry%mpi_Zphys%Comm, ierr)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! compute the number of elements owned by predecessors
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   my_cumul = 0
   Do icore = 0, self%geometry%mpi_ZPhys%rank-1
     my_cumul = my_cumul + array_of_nelems(icore)
   End Do
   DeAllocate(array_of_nElems)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! compute the displacement
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call MPI_Type_Size(MPI_Double, size_of_dble, ierr)
   displacement = my_cumul*size_of_dble


 
   if (self%geometry%mpi_Yphys%rank .eq. 0 ) then
      Call MPI_File_Open(self%geometry%mpi_zphys%comm, filename,&
                         MPI_Mode_WROnly + MPI_Mode_Create,&
                         MPI_Info_Null, file_id, ierr)
      Call MPI_File_set_view(file_id, displacement, MPI_Double, &
                             MPI_Double, 'native', &
                             MPI_Info_Null, ierr)
      Call MPI_File_Write(file_id, fptr, my_nElems, MPI_Double, &
                             MPI_Status_Ignore, Ierr)
      Call MPI_File_Close(file_id, ierr)
   end if

   

 end subroutine
