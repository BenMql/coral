module transforms

 use domain_decomposition
 use read_command_line_args
 use MPI_vars
 use fftw3mpi

 implicit none

 type :: PS_fields_T
   complex(C_double), pointer :: spec (:,:,:)
   real   (C_double), pointer :: phys (:,:,:)
 contains
   procedure :: alloc => allocate_fields
   procedure :: spec_to_phys => c2r_transform
   procedure :: phys_to_spec => r2c_transform
   procedure :: write_phys_to_disk
   procedure :: slice_phys_to_disk
   procedure :: read_phys_from_disk
 end type PS_Fields_T
 
 type(C_ptr), private, save :: p_c2c_backward_y
 type(C_ptr), private, save :: p_c2c_forward_y
 type(C_ptr), private, save :: p_c2r_x
 type(C_ptr), private, save :: p_r2c_x
 type(C_ptr), private, save :: p_r2r_backward_z
 type(C_ptr), private, save :: p_r2r_forward_z


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

 subroutine slice_phys_to_disk (self, slice_kind, slice_index, fileName)
   class(PS_fields_T), intent(inOut) :: self
   character(len=*), intent(in) :: fileName
   integer, intent(in) :: slice_kind
   integer, intent(in) :: slice_index
   print *, 'todo'
 end subroutine

 subroutine write_phys_to_disk (self, filename)
   character(len=*), intent(in) :: fileName
   class(PS_fields_T), intent(inOut) :: self
   print *, 'todo'
 end subroutine

 subroutine read_phys_from_disk (self, filename)
   class(PS_fields_T), intent(inOut) :: self
   character(len=*), intent(in) :: fileName
   print *, 'todo'
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

 end subroutine allocate_fields


 subroutine c2r_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   integer :: u 
   
   self%spec(1,:,:) = self%spec(1,:,:) * 2._dp
   !call fftw_mpi_execute_dft_c2r( p_c2r, self%spec, phys_padded)          
   !phys_buffer(:,:,:) = phys_padded(:, 1:domain_decomp% NXAA, :)  * 0.5_dp 
   !call fftw_execute_r2r        ( p_r2r_i, phys_buffer, self%phys)

 end subroutine c2r_transform

 subroutine r2c_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   real(dp) :: normalization_factor
   real(kind=dp), pointer :: recast_spec_real(:,:,:)

   !integer :: u 
!
   !call fftw_execute_r2r        ( p_r2r_f, self%phys, phys_buffer)
   !phys_padded = 0._dp
   !phys_padded(:, 1:domain_decomp% NXAA, :) = phys_buffer(:,:,:)

   !call fftw_mpi_execute_dft_r2c( p_r2c, phys_padded, self%spec )

   self% spec(1,:,:) = self% spec(1,:,:) * 0.5_dp
   !self% spec = self% spec  /  Real(domain_decomp% NXAA,  Kind = dp)
   !self% spec = self% spec  /  Real(domain_decomp% NYAA,  Kind = dp)
   !self% spec = self% spec  /  Real(domain_decomp% NZAA,  Kind = dp)
 end subroutine r2c_transform


 subroutine dct_planner()
   integer(kind=C_intPtr_T) :: n_fft(2)
   type(c_ptr) :: buf_A, buf_B
   real(dp), pointer :: x_real(:,:,:)
   complex(dp), pointer :: x_cplx(:,:,:)
   real(kind=dp), pointer :: a1(:,:,:)
   if (is_string_in_the_list('--grid-with-endpoints', 21)  &
       .or.                                                &
       is_string_in_the_list('--gauss-lobatto-grid',  20)) then
       
       dct_includes_endpoints = .True.
   end if

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
             domain_decomp% spec_iSize(1) * domain_decomp% spec_iSize(3) & ! howmany
             ,&
             buf_spec_cplx, &
             domain_decomp% spec_iSize(2), & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             domain_decomp% spec_iSize(1)* domain_decomp% spec_iSize(2) &
             ,&
             buf_1c, &
             domain_decomp% spec_iSize(2), & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             domain_decomp% spec_iSize(1)* domain_decomp% spec_iSize(2) &
             ,&
             FFTW_BACKWARD, FFTW_PATIENT &
             )
   ! >>>>>   and its inverse
   p_c2c_forward_y = fftw_plan_many_dft( &
             1,&
            [domain_decomp% spec_iSize(2)], & !< size of the transform
             domain_decomp% spec_iSize(1) * domain_decomp% spec_iSize(3) & ! howmany
             ,&
             buf_1c, &
             domain_decomp% spec_iSize(2), & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             domain_decomp% spec_iSize(1)* domain_decomp% spec_iSize(2) &
             ,&
             buf_spec_cplx, &
             domain_decomp% spec_iSize(2), & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             domain_decomp% spec_iSize(1)* domain_decomp% spec_iSize(2) &
             ,&
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
             domain_decomp% NYAA * domain_decomp% NZ/ world_size & ! howmany
             ,&
             buf_3c, &
            [domain_decomp% NXAA/2 + 1],   & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             domain_decomp% spec_iSize(1)*(domain_decomp% NXAA/2 + 1) &
             ,&
             buf_4r, &
            [domain_decomp% NXAA],         & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             domain_decomp% spec_iSize(1)* domain_decomp% NXAA & 
             ,&
             FFTW_PATIENT &
             )
   ! >>>>>   and its inverse
   p_r2c_x = fftw_plan_many_dft_r2c( &
             1,&
            [domain_decomp% NXAA], & !< size of the transform
             domain_decomp% NYAA * domain_decomp% NZ / world_size& ! howmany
             ,&
             buf_4r, &
            [domain_decomp% NXAA],         & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             domain_decomp% spec_iSize(1)* domain_decomp% NXAA & 
             ,&
             buf_3c, &
            [domain_decomp% NXAA/2 + 1],   & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             domain_decomp% spec_iSize(1)*(domain_decomp% NXAA/2 + 1) &
             ,&
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
   if (dct_includes_endpoints) then
     ! adjust the logical size of the transforms
     logical_NZ = domain_decomp%NZAA - 1
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
            [FFTW_REDFT00],FFTW_PATIENT &
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
            [FFTW_REDFT00],FFTW_PATIENT &
             )
   else
     ! adjust the logical size of the transforms
     logical_NZ = domain_decomp%NZAA
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
            [FFTW_REDFT01],FFTW_PATIENT &
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
            [FFTW_REDFT10],FFTW_PATIENT &
             )
   end if

 end subroutine dct_planner 
 
 subroutine r2r_backward( arr )
   real(kind=dp), allocatable, intent(inOut) :: arr(:)
   type(c_ptr) :: p1, p2
   real(kind=dp), pointer :: meanField_physical(:)
   real(kind=dp), pointer :: meanField_spectral(:)
   type(C_ptr) :: dct_plan_backward_meanFields
   
   p1 = fftw_alloc_real(int( domain_decomp%spec_iSize(1), C_size_T))
   call c_f_pointer(p1, meanField_physical, [domain_decomp%spec_iSize(1)])
   p2 = fftw_alloc_real(int( domain_decomp%spec_iSize(1), C_size_T))
   call c_f_pointer(p2, meanField_spectral, [domain_decomp%spec_iSize(1)])
   
   dct_plan_backward_meanFields = fftw_plan_r2r_1d( domain_decomp%spec_iSize(1), &
                                                   meanField_spectral, &
                                                   meanField_physical, &
                                                   FFTW_REDFT01, FFTW_MEASURE)
   meanField_spectral = arr
   
   meanfield_spectral(1) = meanField_spectral(1) *2._dp
   call fftw_execute_r2r( dct_plan_backward_meanFields, &
                          meanField_spectral, & 
                          meanField_physical)
   arr = meanField_physical/2._dp
   call fftw_destroy_plan( dct_plan_backward_meanFields)
   call fftw_free(p1)
   call fftw_free(p2)
 end subroutine

 subroutine r2r_forward( arr )
   real(kind=dp), allocatable, intent(inOut) :: arr(:)
   type(c_ptr) :: p1, p2
   real(kind=dp), pointer :: meanField_physical(:)
   real(kind=dp), pointer :: meanField_spectral(:)
   type(C_ptr) :: dct_plan_forward_meanFields
   
   p1 = fftw_alloc_real(int( domain_decomp%spec_iSize(1), C_size_T))
   call c_f_pointer(p1, meanField_physical, [domain_decomp%spec_iSize(1)])
   p2 = fftw_alloc_real(int( domain_decomp%spec_iSize(1), C_size_T))
   call c_f_pointer(p2, meanField_spectral, [domain_decomp%spec_iSize(1)])

   dct_plan_forward_meanFields = fftw_plan_r2r_1d( domain_decomp%spec_iSize(1), &
                                                   meanField_physical, &
                                                   meanField_spectral, &
                                                   FFTW_REDFT10, FFTW_MEASURE)
   meanField_physical = arr
   call fftw_execute_r2r( dct_plan_forward_meanFields, &
                          meanField_physical, & 
                          meanField_spectral)
   meanfield_spectral(1) = meanField_spectral(1) *0.5_dp
   arr                   = meanField_spectral / domain_decomp%NZAA 
   
   call fftw_destroy_plan( dct_plan_forward_meanFields)
   call fftw_free(p1)
   call fftw_free(p2)
 end subroutine


end module transforms

