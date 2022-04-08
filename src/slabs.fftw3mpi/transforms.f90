module transforms

 use domain_decomposition
 use read_command_line_args
 use MPI_vars
 use fftw3mpi

 implicit none

 type :: PS_fields_T
   complex(C_double), pointer :: spec (:,:,:)
   real   (C_double), pointer :: phys (:,:,:)
   real   (C_double), pointer :: transposed_phys (:,:,:)
   type(c_ptr) :: trans_ptr
   logical :: has_been_transposed
 contains
   procedure :: transpose => transpose_zxy_to_zyx
   procedure :: alloc => allocate_fields
   procedure :: spec_to_phys => c2r_transform
   procedure :: phys_to_spec => r2c_transform
   procedure :: write_phys_to_disk
   procedure :: slice_phys_to_disk
   procedure :: zSummed_slice_phys_to_disk
   procedure :: read_phys_from_disk
 end type PS_Fields_T
 
 type(C_ptr), private, save :: p_c2c_backward_y
 type(C_ptr), private, save :: p_c2c_forward_y
 type(C_ptr), private, save :: p_c2r_x
 type(C_ptr), private, save :: p_r2c_x
 type(C_ptr), private, save :: p_r2r_backward_z
 type(C_ptr), private, save :: p_r2r_forward_z
 type(C_ptr), private, save :: dct_plan_backward_meanFields 
 type(C_ptr), private, save :: dct_plan_forward_meanFields 


 logical :: dct_includes_endpoints = .False.
 integer :: logical_NZ 

 complex(C_double), pointer, private :: buf_1c (:,:,:)
 complex(C_double), pointer, private :: buf_2c (:,:,:)
 complex(C_double), pointer, private :: buf_3c (:,:,:)
 real   (C_double), pointer, private :: buf_1r (:,:,:)
 real   (C_double), pointer, private :: buf_2r (:,:,:)
 real   (C_double), pointer, private :: buf_4r (:,:,:)
 real   (C_double), pointer, private :: buf_5r (:,:,:)
 !> buffers that are only used for plans creation
 complex(C_double), pointer, private :: buf_spec_cplx (:,:,:)
 real   (C_double), pointer, private :: buf_phys_real (:,:,:)

 contains

 subroutine write_phys_to_disk (self, filename)
   character(len=*), intent(in) :: fileName
   class(PS_fields_T), intent(inOut) :: self
   real(dp), pointer :: dum_ptr

   if (.not. self% has_been_transposed) call self% transpose()
   dum_ptr => self% transposed_phys(1,1,1)
   call DiskWrite_contiguous_memory_chunk(dum_ptr, &
                                          domain_decomp% zyx_iSize(1) * &
                                          domain_decomp% zyx_iSize(2) * &
                                          domain_decomp% zyx_iSize(3) , &
                                          filename,&
                                      int(domain_decomp% zyx_iSize(1) * &
                                          domain_decomp% zyx_iSize(2) * &
                                          domain_decomp% zyx_iSize(3) * &
                                          my_rank, kind= mpi_offset_kind))

   if (my_rank.eq.0) Print *, "... writing ", filename
 end subroutine

 subroutine read_phys_from_disk (self, filename)
   class(PS_fields_T), intent(inOut) :: self
   character(len=*), intent(in) :: fileName
   real(dp), pointer :: dum_ptr

   call fftw_free(self% trans_ptr)

   self% trans_ptr = &
        fftw_alloc_real   ( domain_decomp% NXAA_long * &
                            domain_decomp% NZAA_long * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(self% trans_ptr, self%transposed_phys, &
                                   [domain_decomp% NZAA_long,&
                                    domain_decomp% NYAA_long,&
                                    domain_decomp% NXAA_long / world_size])
   self% has_been_transposed = .True.
   ! read in a buffer
   dum_ptr => self % transposed_phys(1,1,1)
   call DiskRead_contiguous_memory_chunk(dum_ptr, &
                                          domain_decomp% zyx_iSize(1) * &
                                          domain_decomp% zyx_iSize(2) * &
                                          domain_decomp% zyx_iSize(3) , &
                                          filename,&
                                      int(domain_decomp% zyx_iSize(1) * &
                                          domain_decomp% zyx_iSize(2) * &
                                          domain_decomp% zyx_iSize(3) * &
                                          my_rank, kind= mpi_offset_kind))

   call fftw_execute_r2r(domain_decomp% plan_full_fields_transpose_zyx_to_zxy, &
                         self% transposed_phys, self% phys)
 end subroutine

 subroutine allocate_fields(self)
   class(PS_fields_T), intent(inOut) :: self
   type(c_ptr) :: p1, p7
   p1 = fftw_alloc_real   ((domain_decomp% NX_long)   * &
                            domain_decomp% NZ_long   * &
                            domain_decomp% NYAA_long / &
                            world_size )


   call c_f_pointer(p1, self%spec,  [domain_decomp% NZ_long,&
                                     domain_decomp% NYAA_long,&
                                     domain_decomp% local_NX_spec]) 

   p7 = fftw_alloc_real   ( domain_decomp% NXAA_long * &
                            domain_decomp% NZAA_long * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(p7, self%phys, [domain_decomp% NZAA_long,&
                                    domain_decomp% NXAA_long,&
                                    domain_decomp% local_NY_phys]) 

   self% trans_ptr = fftw_alloc_real( int(1, kind= C_intPtr_T) ) ! just because we need to free it later

 end subroutine allocate_fields

 subroutine transpose_zxy_to_zyx(self) 
   class(PS_fields_T), intent(inOut) :: self

   call fftw_free(self% trans_ptr)

   self% trans_ptr = &
        fftw_alloc_real   ( domain_decomp% NXAA_long * &
                            domain_decomp% NZAA_long * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(self% trans_ptr, self%transposed_phys, &
                                   [domain_decomp% NZAA_long,&
                                    domain_decomp% NYAA_long,&
                                    domain_decomp% NXAA_long / world_size])
   call fftw_execute_r2r(domain_decomp% plan_full_fields_transpose_zxy_to_zyx, &
                         self% phys, self% transposed_phys)

   self% has_been_transposed = .True.
 end subroutine transpose_zxy_to_zyx


 subroutine c2r_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   integer :: ix, iy
   type(c_ptr) :: dummy_p1, dummy_p2
   complex(c_double), pointer :: spec_ix(:)
   complex(c_double), pointer :: buf_1c_ix(:)
   complex(c_double), pointer :: buf_3c_iy(:)
   real   (c_double), pointer :: buf_4r_iy(:)
   

   self% has_been_transposed = .False.

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the y transform, from spec to buf_1
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do ix = 1, domain_decomp% local_NX_spec
   dummy_p1 = C_loc(self% spec(1,1,ix))
   dummy_p2 = C_loc(buf_1c    (1,1,ix))
   call C_F_pointer(dummy_p1, spec_ix,   [1])
   call C_F_pointer(dummy_p2, buf_1c_ix, [1])
   call fftw_execute_dft(p_c2c_backward_y, spec_ix, buf_1c_ix)                        
   end do
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> remember the cosine > Chebyshev normalisation
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   buf_1c(1,:,:) = buf_1c(1,:,:) * 2._dp
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the transpose from buf_1 to buf_2
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call fftw_execute_r2r(domain_decomp% plan_forward_transpose, &
                         buf_1r, buf_2r)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> padd buf_2 into buf_3 
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   buf_3c = cmplx(0._dp, 0._dp, kind=dp)
   buf_3c(:,1:domain_decomp% NX/2, :) = buf_2c*0.5_dp 
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the c2r-transform in x from buf_3 to buf_4        
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do iy = 1, domain_decomp% NYAA/world_size
   dummy_p1 = C_loc(buf_3c(1,1,iy))        
   dummy_p2 = C_loc(buf_4r(1,1,iy))        
   call C_F_pointer(dummy_p1, buf_3c_iy, [1])
   call C_F_pointer(dummy_p2, buf_4r_iy, [1])
   call fftw_execute_dft_c2r(p_c2r_x, buf_3c_iy, buf_4r_iy)                   
   end do
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> padd buf_4 into buf_5 
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   buf_5r = 0._dp
   buf_5r(1:domain_decomp% NZ, :,:) = buf_4r 
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the r2r-transform in z from buf_5 to phys_real    
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call fftw_execute_r2r (p_r2r_backward_z, buf_5r, self% phys)

 end subroutine c2r_transform

 subroutine r2c_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   real(dp) :: normalization_factor
   integer :: ix, iy
   type(c_ptr) :: dummy_p1, dummy_p2
   complex(c_double), pointer :: spec_ix(:)
   complex(c_double), pointer :: buf_1c_ix(:)
   complex(c_double), pointer :: buf_3c_iy(:)
   real   (c_double), pointer :: buf_4r_iy(:)

   self% has_been_transposed = .False.
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the r2r-transform in z from phys_real to buf_5    
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call fftw_execute_r2r (p_r2r_forward_z, self% phys, buf_5r)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> shrink buf_5 to buf_4
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   buf_4r = buf_5r(1:domain_decomp% NZ, :,:) 
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the r2c-transform in x from buf_4 to buf_3        
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do iy = 1, domain_decomp% NYAA/world_size
   dummy_p1 = C_loc(buf_3c(1,1,iy))        
   dummy_p2 = C_loc(buf_4r(1,1,iy))        
   call C_F_pointer(dummy_p1, buf_3c_iy, [1])
   call C_F_pointer(dummy_p2, buf_4r_iy, [1])
   call fftw_execute_dft_r2c(p_r2c_x, buf_4r_iy, buf_3c_iy)                   
   end do
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> shrink buf_3 into buf_2
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   buf_2c = buf_3c(:,1:domain_decomp% NX/2, :)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the transpose from buf_2 to buf_1
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call fftw_execute_r2r(domain_decomp% plan_inverse_transpose, &
                         buf_2r, buf_1r)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> remember the cosine > Chebyshev normalisation
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   buf_1c(1,:,:) = buf_1c(1,:,:) * 0.5_dp
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the y transform, from buf_1 to spec
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do ix = 1, domain_decomp% local_NX_spec
   dummy_p1 = C_loc(self% spec(1,1,ix))
   dummy_p2 = C_loc(buf_1c    (1,1,ix))
   call C_F_pointer(dummy_p1, spec_ix,   [1])
   call C_F_pointer(dummy_p2, buf_1c_ix, [1])
   call fftw_execute_dft(p_c2c_forward_y, buf_1c_ix, spec_ix)
   end do
   normalization_factor = 1._dp / (domain_decomp%NXAA * &
                                   domain_decomp%NYAA * &
                                   logical_NZ )
   self% spec = self% spec * normalization_factor




 end subroutine r2c_transform


 subroutine dct_planner()
   integer(kind=C_intPtr_T) :: n_fft(2)
   type(c_ptr) :: buf_A, buf_B
   real(dp), pointer :: x_real(:,:,:)
   complex(dp), pointer :: x_cplx(:,:,:)
   type(c_ptr) :: p1, p2
   real(kind=dp), pointer :: meanField_physical(:)
   real(kind=dp), pointer :: meanField_spectral(:)
   integer(C_fftw_r2r_kind) :: r2r_kind_forward (1), r2r_kind_backward(1)
   if (is_string_in_the_list('--grid-with-endpoints', 21)  &
       .or.                                                &
       is_string_in_the_list('--gauss-lobatto-grid',  20)) then
       
       dct_includes_endpoints = .True.
   end if

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! 
   ! first the 1d chebyshev transforms
   ! 
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   p1 = fftw_alloc_real(int( domain_decomp% NZAA, C_size_T))
   call c_f_pointer(p1, meanField_physical, [domain_decomp% NZAA])
   p2 = fftw_alloc_real(int( domain_decomp% NZAA, C_size_T))
   call c_f_pointer(p2, meanField_spectral, [domain_decomp% NZAA])
   if (dct_includes_endpoints) then
   ! adjust the logical size of the transforms
   logical_NZ = domain_decomp%NZAA - 1
   r2r_kind_backward = [FFTW_REDFT00]
   r2r_kind_forward  = [FFTW_REDFT00]
   else
   logical_NZ = domain_decomp% NZAA
   r2r_kind_backward = [FFTW_REDFT01]
   r2r_kind_forward  = [FFTW_REDFT10]
   end if

   dct_plan_backward_meanFields = fftw_plan_r2r_1d( domain_decomp% NZAA, &
                                                   meanField_spectral, &
                                                   meanField_physical, &
                                                   r2r_kind_backward(1), FFTW_PATIENT)
   dct_plan_forward_meanFields = fftw_plan_r2r_1d( domain_decomp% NZAA, &
                                                   meanField_physical, &
                                                   meanField_spectral, &
                                                   r2r_kind_forward(1), FFTW_PATIENT)
   call fftw_free(p1)
   call fftw_free(p2)

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   ! allocate all necessary buffers
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   buf_A = fftw_alloc_real   ((domain_decomp% NXAA_long+2)   * &
                               domain_decomp% NZAA_long   * &
                               domain_decomp% NYAA_long / &
                               world_size )
   buf_B = fftw_alloc_real   ((domain_decomp% NXAA_long+2)   * &
                               domain_decomp% NZAA_long   * &
                               domain_decomp% NYAA_long / &
                               world_size )

   call c_f_pointer(buf_B, buf_spec_cplx, [domain_decomp% NZ_long,&
                                           domain_decomp% NYAA_long,&
                                           domain_decomp% local_NX_spec]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> stores the result of the y transform
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call c_f_pointer(buf_A, buf_1c, &                         
                                  [domain_decomp% NZ_long,&
                                   domain_decomp% NYAA_long,&
                                   domain_decomp% local_NX_spec]) 
   call c_f_pointer(buf_A, buf_1r, &                         
                                  [domain_decomp% NZ_long,&
                                   domain_decomp% NYAA_long,&
                                   domain_decomp% local_NX_spec*2]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> stores the transposed x <=> y       
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call c_f_pointer(buf_B, buf_2c, &                         
                                  [domain_decomp% NZ_long,  &
                                   domain_decomp% NX_long/2,&
                                   domain_decomp% NYAA_long/&
                                   world_size]) 
   call c_f_pointer(buf_B, buf_2r, &                         
                                  [domain_decomp% NZ_long,&
                                   domain_decomp% NX_long,&
                                   domain_decomp% NYAA_long/&
                                   world_size]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> padds the x direction               
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call c_f_pointer(buf_A, buf_3c, &                         
                                  [domain_decomp% NZ_long,  &
                                   domain_decomp% NXAA_long/2 + 1,&
                                   domain_decomp% NYAA_long/&
                                   world_size]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> stores the result of the c2r-transform in x
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call c_f_pointer(buf_B, buf_4r, &                         
                                  [domain_decomp% NZ_long,  &
                                   domain_decomp% NXAA_long,&
                                   domain_decomp% NYAA_long/&
                                   world_size]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> padds the z direction                       
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call c_f_pointer(buf_A, buf_5r, &                         
                                  [domain_decomp% NZAA_long,  &
                                   domain_decomp% NXAA_long,&
                                   domain_decomp% NYAA_long/&
                                   world_size]) 


   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> stores the result of the r2r-transform in z 
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call c_f_pointer(buf_B, buf_phys_real, &                         
                                  [domain_decomp% NZAA_long,  &
                                   domain_decomp% NXAA_long,&
                                   domain_decomp% NYAA_long/&
                                   world_size]) 


   ! /////////////////////////////////
   ! /////////////////////////////////
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> plans the y transform, from spec to buf_1
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   p_c2c_backward_y = fftw_plan_many_dft( &
             1,&
            [domain_decomp% spec_iSize(2)], & !< size of the transform
             domain_decomp% spec_iSize(1)  & ! howmany
             ,&
             buf_spec_cplx, &
            [domain_decomp% spec_iSize(2)],& !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&                             !< idist                        
             buf_1c, &
            [domain_decomp% spec_iSize(2)],& !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&
             FFTW_BACKWARD, FFTW_PATIENT &
             )
   ! >>>>>   and its inverse
   p_c2c_forward_y = fftw_plan_many_dft( &
             1,&
            [domain_decomp% spec_iSize(2)], & !< size of the transform
             domain_decomp% spec_iSize(1)  & ! howmany
             ,&
             buf_1c, &
            [domain_decomp% spec_iSize(2)],& !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&                             !< idist                        
             buf_spec_cplx, &
            [domain_decomp% spec_iSize(2)],& !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&                             !< idist                        
             FFTW_FORWARD, FFTW_PATIENT &
             )
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> the transpose from buf_1 to buf_2 is already planned
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> padding from buf_2 to buf_3 is taken care of elsewhere
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> plan the c2r-transform in x from buf_3 to buf_4        
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   p_c2r_x = fftw_plan_many_dft_c2r( &
             1,&
            [domain_decomp% NXAA], & !< size of the transform
             domain_decomp% NZ & ! howmany
             ,&
             buf_3c, &
            [domain_decomp% NXAA/2 + 1],   & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&
             buf_4r, &
            [domain_decomp% NXAA],         & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&
             FFTW_PATIENT &
             )
   ! >>>>>   and its inverse
   p_r2c_x = fftw_plan_many_dft_r2c( &
             1,&
            [domain_decomp% NXAA], & !< size of the transform
             domain_decomp% NZ & ! howmany
             ,&
             buf_4r, &
            [domain_decomp% NXAA],         & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&
             buf_3c, &
            [domain_decomp% NXAA/2 + 1],   & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&
             FFTW_PATIENT &
             )
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> padding from buf_4 to buf_5 is taken care of elsewhere
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> plan the r2r-transform in z from buf_5 to phys_real    
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   p_r2r_backward_z = fftw_plan_many_r2r( &
             1,&
            [domain_decomp% NZAA], & !< size of the transform
             domain_decomp% NYAA * domain_decomp% NXAA/ world_size & ! howmany
             ,&
             buf_5r, &
            [domain_decomp% NZAA],   & !< inembed
             1, & !< istride
             domain_decomp% NZAA &
             ,&
             buf_phys_real, &
            [domain_decomp% NZAA],   & !< inembed
             1, & !< istride
             domain_decomp% NZAA &
             ,&
             r2r_kind_backward,FFTW_PATIENT &
             )
   ! >>>>>   and its inverse
   p_r2r_forward_z = fftw_plan_many_r2r( &
             1,&
            [domain_decomp% NZAA], & !< size of the transform
             domain_decomp% NYAA * domain_decomp% NXAA/ world_size & ! howmany
             ,&
             buf_phys_real, &
            [domain_decomp% NZAA],   & !< inembed
             1, & !< istride
             domain_decomp% NZAA &
             ,&
             buf_5r, &
            [domain_decomp% NZAA],   & !< inembed
             1, & !< istride
             domain_decomp% NZAA &
             ,&
             r2r_kind_forward,FFTW_PATIENT &
             )

 end subroutine dct_planner 
 
 subroutine r2r_backward( arr )
   real(kind=dp), allocatable, intent(inOut) :: arr(:)
   type(c_ptr) :: p1, p2
   real(kind=dp), pointer :: meanField_physical(:)
   real(kind=dp), pointer :: meanField_spectral(:)
   
   p1 = fftw_alloc_real(int( domain_decomp% NZAA, C_size_T))
   call c_f_pointer(p1, meanField_physical, [domain_decomp% NZAA])
   p2 = fftw_alloc_real(int( domain_decomp% NZAA, C_size_T))
   call c_f_pointer(p2, meanField_spectral, [domain_decomp% NZAA])
   
   meanField_spectral = arr
   
   meanfield_spectral(1) = meanField_spectral(1) *2._dp
   call fftw_execute_r2r( dct_plan_backward_meanFields, &
                          meanField_spectral, & 
                          meanField_physical)
   arr = meanField_physical/2._dp
   call fftw_free(p1)
   call fftw_free(p2)
 end subroutine

 subroutine r2r_forward( arr )
   real(kind=dp), allocatable, intent(inOut) :: arr(:)
   type(c_ptr) :: p1, p2
   real(kind=dp), pointer :: meanField_physical(:)
   real(kind=dp), pointer :: meanField_spectral(:)
   
   p1 = fftw_alloc_real(int( domain_decomp%NZAA, C_size_T))
   call c_f_pointer(p1, meanField_physical, [domain_decomp%NZAA])
   p2 = fftw_alloc_real(int( domain_decomp%NZAA, C_size_T))
   call c_f_pointer(p2, meanField_spectral, [domain_decomp%NZAA])

   meanField_physical = arr
   call fftw_execute_r2r( dct_plan_forward_meanFields, &
                          meanField_physical, & 
                          meanField_spectral)
   meanfield_spectral(1) = meanField_spectral(1) *0.5_dp
   arr                   = meanField_spectral / domain_decomp%NZAA 
   
   call fftw_free(p1)
   call fftw_free(p2)
 end subroutine

 !> @brief  
 !> Exports a field distributed over processes in contiguous chunks of
 !! memory. 
 !> @details 
 !> Detailled explanation (to be written).                                               
 Subroutine DiskWrite_contiguous_memory_chunk(vector_of_data,  my_nElems,&
                        file_name_str, cumul_nelems)
   !-------------
   ! Arguments
   !-------------
   Real(dp), Pointer :: vector_of_data
                     !< pointer to the data 
                     !! (N.B. this pointer may very well points to an
                     !! allocatable with an arbitrary rank, e.g. 2 or 3).
   Integer, Intent(In) :: my_nElems
                     !< Represents the size (in units of dble) of the slice of
                     !data
                     !! allocated to the current process 
   Integer(kind=mpi_offset_kind), Intent(In), Optional :: cumul_nElems
                     !< Represents the position of the slice of data 
                     !! treated by the current process, in relation to 
                     !! the full array.
                     !! In other terms, what is the displacement (i.e., how many
                     !! elements we need to jump in total) from the beginning of 
                     !! the file.
   Character (len=*), Intent(In) :: file_name_str
                     !< File name 
   !----------------------
   ! Internal variables
   !----------------------
   Integer(kind=MPI_Offset_Kind) :: displacement
                     !< Convert \p cumul_nElems to bytes
   Integer :: size_of_dble
                     !< size of a double precision number in bytes
   Integer :: file_id
                     !< mpiIO keeps tracks of what it's doing
   Integer(kind=MPI_offset_kind) :: my_cumul


      ! ~~~~~~~~~~~~~~~~~~~~
      ! if cumul_nElems is missing, we compute it
      ! ~~~~~~~~~~~~~~~~~~~~
      if (present(cumul_nElems)) Then
         my_cumul = cumul_nElems
      else
         my_cumul = compute_cumul_nElems(my_nElems)
      End If
      Call MPI_Type_Size(MPI_Double, size_of_dble, ierr)
      displacement = my_cumul*size_of_dble

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
      !      Write the Quicksave         !                         
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
      Call MPI_File_Open(MPI_Comm_world, file_name_str,&
                         MPI_Mode_WROnly + MPI_Mode_Create,&
                         MPI_Info_Null, file_id, ierr)
      Call MPI_File_set_view(file_id, displacement, MPI_Double, &
                             MPI_Double, 'native', &
                             MPI_Info_Null, ierr)
      Call MPI_File_Write(file_id, vector_of_data, my_nElems, MPI_Double, &
                             MPI_Status_Ignore, Ierr)
      Call MPI_File_Close(file_id, ierr)
 End Subroutine DiskWrite_contiguous_memory_chunk

 integer(kind=mpi_offset_kind) function compute_cumul_nElems(my_nElems)
 integer, intent(in) :: my_nElems
 integer, dimension(:), allocatable :: array_of_nelems
 integer :: icore
 !integer :: ierr

         Allocate(array_of_nElems(0:world_size-1))
         Call MPI_allGather(my_nElems, 1, MPI_Integer,&
                array_of_nElems, 1, MPI_integer,  MPI_Comm_World, ierr)
         !Call MPI_Gather(my_nElems, 1, MPI_Integer,&
                !array_of_nElems, 1, MPI_integer, 0, MPI_Comm_World, ierr)
         !Call MPI_Bcast(array_of_nElems, world_size, MPI_integer, 0, MPI_Comm_world, ierr)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! compute the number of elements owned by predecessors
         compute_cumul_nElems = 0
         Do icore = 0, my_rank-1
           compute_cumul_nElems = compute_cumul_nElems + array_of_nelems(icore)
         End Do
         deAllocate(array_of_nElems)
 end function

 !> @brief  
 !> Reads a files and distributes the data to processes
 !> @details 
 !> Detailled explanation (to be written).                                               
 Subroutine DiskRead_contiguous_memory_chunk(vector_of_data,  my_nElems,&
                        file_name_str,  cumul_nelems)
   !-------------
   ! Arguments
   !-------------
   Real(dp), Pointer :: vector_of_data
                     !< pointer to the data 
                     !! (N.B. this pointer may very well points to an
                     !! allocatable with an arbitrary rank, e.g. 2 or 3).
   Integer, Intent(In) :: my_nElems
                     !< Represents the size (in units of dble) of the slice of data
                     !! allocated to the current process 
   Integer(kind=mpi_offset_kind), Intent(In), Optional :: cumul_nElems
                     !< Represents the position of the slice of data 
                     !! treated by the current process, in relation to 
                     !! the full array.
                     !! In other terms, what is the displacement (i.e., how many
                     !! elements we need to jump in total) from the beginning of 
                     !! the file.
   Character (len=*), Intent(In) :: file_name_str
                     !< File name 
   !----------------------
   ! Internal variables
   !----------------------
   Integer(kind=MPI_Offset_Kind) :: displacement
                     !< Convert \p cumul_nElems to bytes
   Integer :: size_of_dble
                     !< size of a double precision number in bytes
   Integer :: file_id
                     !< mpiIO keeps tracks of what it's doing
   Integer(kind=mpi_offset_kind) :: my_cumul


      ! ~~~~~~~~~~~~~~~~~~~~
      ! if cumul_nElems is missing, we compute it
      ! ~~~~~~~~~~~~~~~~~~~~
      if (present(cumul_nElems)) Then 
         my_cumul = cumul_nElems
      else 
         my_cumul = compute_cumul_nElems(my_nElems)
      End If
      
      Call MPI_Type_Size(MPI_Double, size_of_dble, ierr)
      displacement = my_cumul*size_of_dble
      

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
      !      Write the Quicksave         !                         
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
      Call MPI_File_Open(MPI_Comm_world, file_name_str,&
                         MPI_Mode_RdOnly,&
                         MPI_Info_Null, file_id, ierr)
      Call MPI_File_set_view(file_id, displacement, MPI_Double, &
                             MPI_Double, 'native', &
                             MPI_Info_Null, ierr)
      Call MPI_File_Read(file_id, vector_of_data, my_nElems, MPI_Double, &
                             MPI_Status_Ignore, Ierr)
      Call MPI_File_Close(file_id, ierr)
 End Subroutine DiskRead_contiguous_memory_chunk

 subroutine zSummed_slice_phys_to_disk (self, cheby_weight, fileName)
   class(PS_fields_T), intent(inOut) :: self
   character(len=*), intent(in) :: fileName
   integer :: n1, n2
   real(dp), allocatable, target :: aux2d(:,:)
   real(dp), pointer :: dum_ptr
   real(dp), allocatable, dimension(:), intent(in):: cheby_weight
   integer :: ix, iy, iz

   if (.not. self% has_been_transposed) call self% transpose()

   n1 = domain_decomp% NYAA
   n2 = domain_decomp% NXAA / world_size
   allocate (aux2d( n1, n2 ))                                            
   aux2d = 0._dp
   do ix = 1, domain_decomp% NXAA / world_size
   do iy = 1, domain_decomp% NYAA 
   do iz = 1, domain_decomp% NZAA
   aux2d(iy,ix) = aux2d(iy,ix) + &
                     self%transposed_phys(iz,iy,ix) * cheby_weight(iz)
   end do
   end do
   end do
   dum_ptr =>  aux2d(1,1)
   call diskWrite_contiguous_memory_chunk(dum_ptr, &
                                            n1*n2,   &
                                         fileName,&
                                int(n1*n2*my_rank, &
                             kind=mpi_offset_kind))
   deallocate (aux2d)
   if (my_rank.eq.0) Print *, "... writing ", filename


 end subroutine

 subroutine slice_phys_to_disk (self, slice_kind, slice_index, fileName)
   class(PS_fields_T), intent(inOut) :: self
   character(len=*), intent(in) :: fileName
   integer, intent(in) :: slice_kind
   integer, intent(in) :: slice_index
   integer :: n1, n2
   real(dp), allocatable, target :: aux2d(:,:)
   real(dp), pointer :: dum_ptr


   if (.not. self% has_been_transposed) call self% transpose()

   select case (slice_kind)
          case (1)
             n1 = domain_decomp% NYAA
             n2 = domain_decomp% NXAA / world_size
             allocate (aux2d( n1, n2 ))                                            
             aux2d = self%transposed_phys( slice_index, :, :)
             dum_ptr =>  aux2d(1,1)
             call diskWrite_contiguous_memory_chunk(dum_ptr, &
                                                    n1*n2,   &
                                                    fileName,&
                                                int(n1*n2*my_rank, &
                                                kind=mpi_offset_kind))
             deallocate (aux2d)
          case (2)
             n1 = domain_decomp% NZAA
             n2 = domain_decomp% NXAA / world_size
             allocate (aux2d( n1, n2 ))                                            
             aux2d = self%transposed_phys(:, slice_index,  :)
             dum_ptr =>  aux2d(1,1)
             call diskWrite_contiguous_memory_chunk(dum_ptr, &
                                                    n1*n2,   &
                                                    fileName,&
                                                int(n1*n2*my_rank, &
                                                kind=mpi_offset_kind))
             deallocate (aux2d)
          case (3)
            if     ( (slice_index .gt. domain_decomp%NXAA / world_size *  my_rank) &
               .and. (slice_index .le. domain_decomp%NXAA / world_size * (my_rank+1)) ) then
               n1 = domain_decomp% NZAA
               n2 = domain_decomp% NYAA 
               allocate (aux2d( n1, n2 ))                                            
               aux2d = self%transposed_phys(:, :, slice_index - &
                                   domain_decomp%NXAA / world_size )
               open (unit=9, file=fileName, status='replace', access='stream')
               write(9) aux2d
               close(unit=9)
               deallocate (aux2d)
            end if
    end select

   if (my_rank.eq.0) Print *, "... writing ", filename

 end subroutine


end module transforms

