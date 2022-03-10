module transforms

 use decomp_2d
 use decomp_2d_fft
 use decomp_2d_io
 use MPI_vars
 use domain_decomposition
 use fftw3_wrap
 use read_command_line_args

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
 
 type :: specField_T
   complex(C_double), pointer :: spec (:,:,:)
 contains
   procedure :: alloc => allocate_specField
 end type specField_T

 type(C_ptr), private, save :: dct_plan_forward
 type(C_ptr), private, save :: dct_plan_backward

 logical :: dct_includes_endpoints = .False.
 integer :: logical_NZ 

 contains

 subroutine slice_phys_to_disk (self, slice_kind, slice_index, fileName)
   class(PS_fields_T), intent(inOut) :: self
   character(len=*), intent(in) :: fileName
   integer, intent(in) :: slice_kind
   integer, intent(in) :: slice_index
   call decomp_2d_write_plane(3, self%phys, slice_kind, &
                                     slice_index, fileName)              
 end subroutine

 subroutine write_phys_to_disk (self, filename)
   character(len=*), intent(in) :: fileName
   class(PS_fields_T), intent(inOut) :: self
   call decomp_2d_write_one(3, self%phys, fileName)              
 end subroutine

 subroutine read_phys_from_disk (self, filename)
   class(PS_fields_T), intent(inOut) :: self
   character(len=*), intent(in) :: fileName
   call decomp_2d_read_one(3, self%phys, fileName)              
 end subroutine

 subroutine allocate_fields(self)
   class(PS_fields_T), intent(inOut) :: self
   type(c_ptr) :: p1, p2
   202 format ('allocating spectral fields of shape ', (i4.4), &
                                                 ', ', (i4.4), &
                                                 ', ', (i4.4))

   201 format ('allocating physical fields of shape ', (i4.4), &
                                                 ', ', (i4.4), &
                                                 ', ', (i4.4))

   if (my_rank.eq.0) write (*,202) domain_decomp%spec_iSize(1),&
                                   domain_decomp%spec_iSize(2),&
                                   domain_decomp%spec_iSize(3)

   if (my_rank.eq.0) write (*,201) domain_decomp%phys_iSize(1),&
                                   domain_decomp%phys_iSize(2),&
                                   domain_decomp%phys_iSize(3)

   p1 = fftw_alloc_complex(int( domain_decomp%spec_iSize(1) * &
                                domain_decomp%spec_iSize(2) * &
                                domain_decomp%spec_iSize(3) , &
                                C_size_T))

   call c_f_pointer(p1, self%spec, [domain_decomp%spec_iSize(1),&
                                    domain_decomp%spec_iSize(2),&
                                    domain_decomp%spec_iSize(3)]) 

   p2 = fftw_alloc_real   (int( domain_decomp%phys_iSize(1) * &
                                domain_decomp%phys_iSize(2) * &
                                domain_decomp%phys_iSize(3) , &
                                C_size_T))

   call c_f_pointer(p2, self%phys, [domain_decomp%phys_iSize(1),&
                                    domain_decomp%phys_iSize(2),&
                                    domain_decomp%phys_iSize(3)]) 

 end subroutine allocate_fields

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


 subroutine c2r_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   type(c_ptr) :: dummy_ptr
   real(kind=dp), pointer :: recast_spec_real(:,:,:)
   
   self%spec(1,:,:) = self%spec(1,:,:) * 2._dp
   dummy_ptr = c_loc(self%spec)
   call C_F_pointer (dummy_ptr, recast_spec_real, [2*domain_decomp%spec_iSize(1),&
                                                     domain_decomp%spec_iSize(2),&
                                                     domain_decomp%spec_iSize(3)])
   ! Inverse DCT, real part
   call fftw_execute_r2r ( dct_plan_backward, recast_spec_real, recast_spec_real )
   dummy_ptr = transfer(transfer( dummy_ptr, 1_c_intptr_T) &
                                + c_sizeOf(1._c_double), dummy_ptr) 
   call C_F_pointer (dummy_ptr, recast_spec_real, [2*domain_decomp%spec_iSize(1),&
                                                     domain_decomp%spec_iSize(2),&
                                                     domain_decomp%spec_iSize(3)])
   ! Inverse DCT, imaginary part
   call fftw_execute_r2r ( dct_plan_backward, recast_spec_real, recast_spec_real )
   call decomp_2d_fft_3d(self%spec, self%phys)
   self%phys = self%phys*0.5_dp
 end subroutine c2r_transform

 subroutine r2c_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   real(dp) :: normalization_factor
   type(c_ptr) :: dummy_ptr
   real(kind=dp), pointer :: recast_spec_real(:,:,:)
   Call decomp_2d_fft_3d(self%phys, self%spec)
   dummy_ptr = c_loc(self%spec)
   call C_F_pointer (dummy_ptr, recast_spec_real, [2*domain_decomp%spec_iSize(1),&
                                                     domain_decomp%spec_iSize(2),&
                                                     domain_decomp%spec_iSize(3)])
   ! Forward DCT, real part
   call fftw_execute_r2r ( dct_plan_forward,  recast_spec_real, recast_spec_real )
   dummy_ptr = transfer(transfer( dummy_ptr, 1_c_intptr_T) &
                                + c_sizeOf(1._C_double), dummy_ptr) 
   call C_F_pointer (dummy_ptr, recast_spec_real, [2*domain_decomp%spec_iSize(1),&
                                                     domain_decomp%spec_iSize(2),&
                                                     domain_decomp%spec_iSize(3)])
   ! Forward DCT, imaginary part
   call fftw_execute_r2r ( dct_plan_forward,  recast_spec_real, recast_spec_real )
   normalization_factor = 1._dp / (domain_decomp%NXAA * &
                                   domain_decomp%NYAA * &
                                   logical_NZ )

   self%spec(1,:,:) = self%spec(1,:,:) * 0.5_dp
   self%spec = self%spec *normalization_factor
 end subroutine r2c_transform

 subroutine allocate_specfield(self)
   class(specfield_T), intent(inOut) :: self
   allocate( self%spec (domain_decomp%spec_iSize(1),&
                        domain_decomp%spec_iSize(2),&
                        domain_decomp%spec_iSize(3)))
 end subroutine allocate_specfield


 subroutine dct_planner()
   integer(kind=C_intPtr_T) :: sz
   type(c_ptr) :: a1_p
   real(kind=dp), pointer :: a1(:,:,:)
   if (is_string_in_the_list('--grid-with-endpoints', 21)  &
       .or.                                                &
       is_string_in_the_list('--gauss-lobatto-grid',  20)) then
       
       dct_includes_endpoints = .True.
         
   end if
   sz = domain_decomp%spec_iSize(1) *&
        domain_decomp%spec_iSize(2) *&
        domain_decomp%spec_iSize(3)
   a1_p = fftw_alloc_complex(int(sz, C_intptr_T))
   call c_f_pointer(a1_p,a1,[2*domain_decomp%spec_iSize(1), &
                               domain_decomp%spec_iSize(2), &
                               domain_decomp%spec_iSize(3)])
   if (dct_includes_endpoints) then
     ! adjust the logical size of the transforms
     logical_NZ = domain_decomp%NZAA - 1
     ! plan the transforms
     dct_plan_forward =  fftw_plan_many_r2r(1, [domain_decomp%spec_iSize(1)], &
         domain_decomp%spec_iSize(2)*domain_decomp%spec_iSize(3), &
         a1, [2*domain_decomp%spec_iSize(1)], 2, 2*domain_decomp%spec_iSize(1), &
         a1, [2*domain_decomp%spec_iSize(1)], 2, 2*domain_decomp%spec_iSize(1), &
         [FFTW_REDFT00], FFTW_MEASURE)
     dct_plan_backward=  fftw_plan_many_r2r(1, [domain_decomp%spec_iSize(1)], &
         domain_decomp%spec_iSize(2)*domain_decomp%spec_iSize(3), &
         a1, [2*domain_decomp%spec_iSize(1)], 2, 2*domain_decomp%spec_iSize(1), &
         a1, [2*domain_decomp%spec_iSize(1)], 2, 2*domain_decomp%spec_iSize(1), &
         [FFTW_REDFT00], FFTW_MEASURE)
   else
     ! adjust the logical size of the transforms
     logical_NZ = domain_decomp%NZAA
     ! plan the transforms
     dct_plan_forward =  fftw_plan_many_r2r(1, [domain_decomp%spec_iSize(1)], &
         domain_decomp%spec_iSize(2)*domain_decomp%spec_iSize(3), &
         a1, [2*domain_decomp%spec_iSize(1)], 2, 2*domain_decomp%spec_iSize(1), &
         a1, [2*domain_decomp%spec_iSize(1)], 2, 2*domain_decomp%spec_iSize(1), &
         [FFTW_REDFT10], FFTW_MEASURE)
     dct_plan_backward=  fftw_plan_many_r2r(1, [domain_decomp%spec_iSize(1)], &
         domain_decomp%spec_iSize(2)*domain_decomp%spec_iSize(3), &
         a1, [2*domain_decomp%spec_iSize(1)], 2, 2*domain_decomp%spec_iSize(1), &
         a1, [2*domain_decomp%spec_iSize(1)], 2, 2*domain_decomp%spec_iSize(1), &
         [FFTW_REDFT01], FFTW_MEASURE)
   end if

   
   call fftw_free(a1_p)   

 end subroutine dct_planner 
 

end module transforms

