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
   integer :: ix, iy
   type(c_ptr) :: dummy_p1, dummy_p2
   complex(c_double), pointer :: spec_ix(:)
   complex(c_double), pointer :: buf_1c_ix(:)
   complex(c_double), pointer :: buf_3c_iy(:)
   real   (c_double), pointer :: buf_4r_iy(:)
   
   !!pouet
   !!do ix = 1, domain_decomp% local_NX_spec
      !print *, '////////////////////////////////////////////'
      !print *, self%spec (:,:,1)
      !print *, '////////////////////////////////////////////'
   !!end do
   !!pouet


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
   call fftw_execute_r2r (p_r2r_forward_z, buf_5r, self% phys)

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

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !
   !> execute the r2r-transform in z from phys_real to buf_5    
   !
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call fftw_execute_r2r (p_r2r_backward_z, self% phys, buf_5r)
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
   !call fftw_execute_dft(p_c2c_backward_y, self% spec(1,1,ix), buf_1c(1,1,ix))
   call fftw_execute_dft(p_c2c_forward_y, buf_1c_ix, spec_ix)
   end do
   normalization_factor = 1._dp / (domain_decomp%NXAA * &
                                   domain_decomp%NYAA * &
                                   logical_NZ )
   self% spec = self% spec * normalization_factor


   !!pouet
   !!do ix = 1, domain_decomp% local_NX_spec
      !print *, '////////////////////////////////////////////'
      !print *, self%spec (:,:,1)
      !print *, '////////////////////////////////////////////'
   !!end do
   !!pouet


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

   !print *, "pouet BUF_A,B have size:", (domain_decomp% NXAA_long+2)   * &
                               !domain_decomp% NZAA_long   * &
                               !domain_decomp% NYAA_long / &
                               !world_size
   call c_f_pointer(buf_B, buf_spec_cplx, [domain_decomp% NZ_long,&
                                           domain_decomp% NYAA_long,&
                                           domain_decomp% local_NX_spec]) 
   !print *, "pouet spec has size:", shape(buf_spec_cplx)

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
   !print *, "pouet buf_1c has size:", shape(buf_1c)
   !print *, "pouet buf_1r has size:", shape(buf_1r)             

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
   !print *, "pouet p_c2c has size            :", domain_decomp% spec_iSize(2)
   !print *, "pouet p_c2c has howmany         :",  domain_decomp% spec_iSize(1) * domain_decomp% spec_iSize(3)
   !print *, "pouet p_c2c has inembed         :", domain_decomp% spec_iSize(2)
   !print *, "pouet p_c2c has istride         :", domain_decomp% spec_iSize(1)
   !print *, "pouet p_c2c has iDist           :", domain_decomp% spec_iSize(1)* domain_decomp% spec_iSize(2) 

   p_c2c_backward_y = fftw_plan_many_dft( &
             1,&
            [domain_decomp% spec_iSize(2)], & !< size of the transform
             domain_decomp% spec_iSize(1)  & ! howmany
             ,&
             buf_spec_cplx, &
             domain_decomp% spec_iSize(2), & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&                             !< idist                        
             buf_1c, &
             domain_decomp% spec_iSize(2), & !< inembed
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
             domain_decomp% spec_iSize(2), & !< inembed
             domain_decomp% spec_iSize(1), & !< istride
             1,&                             !< idist                        
             buf_spec_cplx, &
             domain_decomp% spec_iSize(2), & !< inembed
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
   
   p1 = fftw_alloc_real(int( domain_decomp% NZAA, C_size_T))
   call c_f_pointer(p1, meanField_physical, [domain_decomp% NZAA])
   p2 = fftw_alloc_real(int( domain_decomp%spec_iSize(1), C_size_T))
   call c_f_pointer(p2, meanField_spectral, [domain_decomp% NZAA])
   
   dct_plan_backward_meanFields = fftw_plan_r2r_1d( domain_decomp% NZAA, &
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
   
   p1 = fftw_alloc_real(int( domain_decomp%NZAA, C_size_T))
   call c_f_pointer(p1, meanField_physical, [domain_decomp%NZAA])
   p2 = fftw_alloc_real(int( domain_decomp%spec_iSize(1), C_size_T))
   call c_f_pointer(p2, meanField_spectral, [domain_decomp%NZAA])

   dct_plan_forward_meanFields = fftw_plan_r2r_1d( domain_decomp%NZAA, &
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

