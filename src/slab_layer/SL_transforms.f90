module SL_transforms

 use SL_geometry

 implicit none

 type :: PS_fields_T
   complex(kind=dp), pointer, dimension(:,:,:) :: spec
   real   (kind=dp), pointer, dimension(:,:,:) :: phys
 contains
   procedure :: alloc => allocate_fields
   !procedure :: spec_to_phys => c2r_transform
   !procedure :: phys_to_spec => r2c_transform
 end type PS_Fields_T
 
 type :: specField_T
   complex(kind=dp), pointer, dimension(:,:,:) :: spec
 contains
   procedure :: alloc => allocate_specField
 end type specField_T


 contains

subroutine DFTs_plans_creation(p_r2c, p_c2r, p_r2r_f, p_r2r_i, X_padd, X_mixt, X_phys)
   Complex(dp), Pointer :: x_padd(:,:,:)
   Real(dp), Pointer :: x_mixt (:,:,:)
   Real(dp), Pointer :: x_phys (:,:,:)
   Type(C_Ptr) :: p_r2c
   Type(C_Ptr) :: p_c2r
   Type(C_Ptr) :: p_r2r_f
   Type(C_Ptr) :: p_r2r_i
   Integer (C_IntPtr_t) :: N_fft(2)

   N_fft = [geom%NXAA, geom%NYAA]

   p_c2r = fftw_mpi_plan_many_dft_c2r( 2,& ! rank
                                   N_fft,& ! logical size
                                   int(geom%NZAA, kind=C_intPtr_T),& ! how many transforms
                  FFTW_MPI_DEFAULT_BLOCK,& ! let fftw decide
                  FFTW_MPI_DEFAULT_BLOCK,& ! let fftw decide
                                  X_padd,& ! complex input (TRANSPOSED)
                                  X_mixt,& ! real output
                          MPI_COMM_WORLD,&
                          ior(FFTW_estimate, FFTW_MPI_TRANSPOSED_In))

   p_r2c = fftw_mpi_plan_many_dft_r2c( 2,& ! rank
                                   N_fft,& ! logical size
                                   int(geom%NZAA, kind=C_intPtr_T),& ! how many transforms
                  FFTW_MPI_DEFAULT_BLOCK,& ! let fftw decide
                  FFTW_MPI_DEFAULT_BLOCK,& ! let fftw decide
                                  X_mixt,& ! Real input
                                  X_padd,& ! Complex output (TRANSPOSED)
                          MPI_COMM_WORLD,&
                          ior(FFTW_estimate, FFTW_MPI_TRANSPOSED_Out))

   p_r2r_i = fftw_plan_many_r2r(1, int([geom%NZAA],kind=4), int(2*(geom%NYAA/2+1)*geom%localSize_phys, kind=4), &
                                X_mixt, [int(geom%NZAA, kind=4)], 1, int(geom%NZAA, kind=4), &
                                X_phys, [int(geom%NZAA, kind=4)], 1, int(geom%NZAA, kind=4), &
                                [FFTW_REDFT01], FFTW_estimate)

   p_r2r_f = fftw_plan_many_r2r(1, int([geom%NZAA],kind=4), int(2*(geom%NYAA/2+1)*geom%localSize_phys, kind=4), &
                                X_mixt, [int(geom%NZAA, kind=4)], 1, int(geom%NZAA, kind=4), &
                                X_phys, [int(geom%NZAA, kind=4)], 1, int(geom%NZAA, kind=4), &
                                [FFTW_REDFT10], FFTW_estimate)

 end subroutine DFTs_plans_creation




 subroutine allocate_fields(self)
   class(PS_fields_T), intent(inOut) :: self
   type(c_ptr) :: pdum
   202 format ('allocating spectral fields of shape ', (i4.4), &
                                                 ', ', (i4.4), &
                                                 ', ', (i4.4))
   201 format ('allocating physical fields of shape ', (i4.4), &
                                                 ', ', (i4.4), &
                                                 ', ', (i4.4))
   if (my_rank.eq.0) write (*,202) domain_decomp%NZAA,&
                        domain_decomp%NXAA,&
                        domain_decomp%localSize_spec
   if (my_rank.eq.0) write (*,201) domain_decomp%NZAA,&
                        domain_decomp%NYAA + 2,&
                        domain_decomp%localSize_phys
   pdum = fftw_alloc_complex(geom%alloc_local)
   call C_F_pointer( pdum, self%spec, &
            [geom%NZAA, geom%NXAA, &
             int(geom%localSize_spec, kind=C_int)])
   pdum = fftw_alloc_real( 2*geom%alloc_local)
   call C_F_pointer( pdum, self%phys, &
            [geom%NZAA, geom%NXAA+2, &
             int(geom%localSize_phys, kind=C_int)])
 end subroutine allocate_fields

 subroutine allocate_specfield(self)
   class(specfield_T), intent(inOut) :: self
   type(c_ptr) :: pdum
   pdum = fftw_alloc_complex(geom%alloc_local)
   call C_F_pointer( pdum, self%spec, &
            [geom%NZAA, geom%NXAA, &
             int(geom%localSize_spec, kind=C_int)])
 end subroutine allocate_specfield



end module SL_transforms

