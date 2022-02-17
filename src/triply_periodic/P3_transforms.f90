module P3_transforms

 use decomp_2d
 use decomp_2d_fft
 use MPI_vars
 use P3_domain_decomp
 use read_command_line_args

 implicit none

 type :: PS_fields_T
   complex(kind=dp), allocatable, dimension(:,:,:) :: spec
   real   (kind=dp), allocatable, dimension(:,:,:) :: phys
 contains
   procedure :: alloc => allocate_fields
   procedure :: spec_to_phys => c2r_transform
   procedure :: phys_to_spec => r2c_transform
 end type PS_Fields_T
 
 type :: specField_T
   complex(kind=dp), allocatable, dimension(:,:,:) :: spec
 contains
   procedure :: alloc => allocate_specField
 end type specField_T

 contains



 subroutine allocate_fields(self)
   class(PS_fields_T), intent(inOut) :: self
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
   allocate( self%spec (domain_decomp%spec_iSize(1),&
                        domain_decomp%spec_iSize(2),&
                        domain_decomp%spec_iSize(3)))
   allocate( self%phys (domain_decomp%phys_iSize(1),&
                        domain_decomp%phys_iSize(2),&
                        domain_decomp%phys_iSize(3)))
 end subroutine allocate_fields


 subroutine c2r_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   
   ! Inverse DCT, imaginary part
   call decomp_2d_fft_3d(self%spec, self%phys)
 end subroutine c2r_transform

 subroutine r2c_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   real(dp) :: normalization_factor
   Call decomp_2d_fft_3d(self%phys, self%spec)
   normalization_factor = 1._dp / (domain_decomp%NXAA * &
                                   domain_decomp%NYAA * &
                                   domain_decomp%NZAA )
   self%spec = self%spec *normalization_factor
 end subroutine r2c_transform

 subroutine allocate_specfield(self)
   class(specfield_T), intent(inOut) :: self
   allocate( self%spec (domain_decomp%spec_iSize(1),&
                        domain_decomp%spec_iSize(2),&
                        domain_decomp%spec_iSize(3)))
 end subroutine allocate_specfield



end module P3_transforms

