module transforms

 use domain_decomposition
 use read_command_line_args
 use MPI_vars
 use fftw3mpi

 implicit none

 type :: PS_fields_T
   real   (C_double), pointer :: spec_real (:,:,:)
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
 
 type(C_ptr), private, save :: p_c2c_y

 type(C_ptr), private, save :: p_r2c
 type(C_ptr), private, save :: p_c2r
 type(C_ptr), private, save :: p_r2r_f           
 type(C_ptr), private, save :: p_r2r_i           

 logical :: dct_includes_endpoints = .False.
 integer :: logical_NZ 

 complex(C_double), pointer, private :: buf_spec_cplx (:,:,:)
 complex(C_double), pointer, private :: buf_YphysPadded_xDist_cplx (:,:,:)
 real   (C_double), pointer, private :: buf_YphysPadded_xDist_real (:,:,:)
 complex(C_double), pointer, private :: buf_YphysPadded_yDist_cplx (:,:,:)
 real   (C_double), pointer, private :: buf_YphysPadded_yDist_real (:,:,:)
 complex(C_double), pointer, private :: buf_XspecPadded_cplx (:,:,:)
 real   (C_double), pointer, private :: buf_XphysPadded_real (:,:,:)
 real   (C_double), pointer, private :: buf_ZspecPadded_real (:,:,:)
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

   p1 = fftw_alloc_real   ( domain_decomp% NX_long   * &
                            domain_decomp% NZ_long   * &
                            domain_decomp% NYAA_long / &
                            world_size )


   call c_f_pointer(p1, self%spec,      [domain_decomp% NZ_long,&
                                         domain_decomp% NYAA_long,&
                                         domain_decomp% local_NX_spec]) 

   call c_f_pointer(p1, self%spec_real, [domain_decomp% NZ_long,&
                                         domain_decomp% NYAA_long,&
                                         domain_decomp% local_NX_spec*2]) 

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
   call fftw_mpi_execute_dft_c2r( p_c2r, self%spec, phys_padded)          
   phys_buffer(:,:,:) = phys_padded(:, 1:domain_decomp% NXAA, :)  * 0.5_dp 
   call fftw_execute_r2r        ( p_r2r_i, phys_buffer, self%phys)

 end subroutine c2r_transform

 subroutine r2c_transform(self)
   class(PS_fields_T), intent(inOut), target :: self
   real(dp) :: normalization_factor
   real(kind=dp), pointer :: recast_spec_real(:,:,:)

   integer :: u 

   call fftw_execute_r2r        ( p_r2r_f, self%phys, phys_buffer)
   phys_padded = 0._dp
   phys_padded(:, 1:domain_decomp% NXAA, :) = phys_buffer(:,:,:)

   call fftw_mpi_execute_dft_r2c( p_r2c, phys_padded, self%spec )

   self% spec(1,:,:) = self% spec(1,:,:) * 0.5_dp
   self% spec = self% spec  /  Real(domain_decomp% NXAA,  Kind = dp)
   self% spec = self% spec  /  Real(domain_decomp% NYAA,  Kind = dp)
   self% spec = self% spec  /  Real(domain_decomp% NZAA,  Kind = dp)
 end subroutine r2c_transform


 subroutine dct_planner()
   integer(kind=C_intPtr_T) :: n_fft(2)
   type(c_ptr) :: x_real_ptr, x_cplx_ptr, p2, p3
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

   p1 = fftw_alloc_real   ( domain_decomp% NX_long   * &
                            domain_decomp% NZ_long   * &
                            domain_decomp% NYAA_long / &
                            world_size )


   call c_f_pointer(p1, buf_spec_cplx, [domain_decomp% NZ_long,&
                                        domain_decomp% NYAA_long,&
                                        domain_decomp% local_NX_spec]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   p2 = fftw_alloc_real   ( domain_decomp% NX_long   * &
                            domain_decomp% NZ_long   * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(p2, buf_YphysPadded_xDist_cplx, [domain_decomp% NZ_long,&
                                   domain_decomp% NYAA_long,&
                                   domain_decomp% local_NX_spec]) 
   call c_f_pointer(p2, buf_YphysPadded_xDist_real, [domain_decomp% NZ_long*2,&
                                   domain_decomp% NYAA_long,&
                                   domain_decomp% local_NX_spec]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   p3 = fftw_alloc_real   ( domain_decomp% NX_long   * &
                            domain_decomp% NZ_long   * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(p3, buf_YphysPadded_yDist_cplx, [domain_decomp% NZ_long,&
                                   domain_decomp% NX_long/2,&
                                   domain_decomp% local_NY_phys]) 
   call c_f_pointer(p3, buf_YphysPadded_yDist_real, [domain_decomp% NZ_long,&
                                   domain_decomp% NX_long,&
                                   domain_decomp% local_NY_phys]) 


   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   p4 = fftw_alloc_real   ((domain_decomp% NXAA_long +2) * &
                            domain_decomp% NZ_long   * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(p4, buf_XspecPadded_cplx, [domain_decomp% NZ_long,&
                                   domain_decomp% NXAA_long/2+1,&
                                   domain_decomp% local_NY_phys]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   p5 = fftw_alloc_real   ((domain_decomp% NXAA_long+2) * &
                            domain_decomp% NZ_long   * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(p5, buf_XphysPadded_real, [domain_decomp% NZ_long,&
                                   domain_decomp% NXAA_long + 2,&
                                   domain_decomp% local_NY_phys]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   p6 = fftw_alloc_real   ( domain_decomp% NXAA_long * &
                            domain_decomp% NZ_long   * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(p6, buf_ZspecPadded_real, [domain_decomp% NZ_long,&
                                   domain_decomp% NXAA_long,&
                                   domain_decomp% local_NY_phys]) 

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   p7 = fftw_alloc_real   ( domain_decomp% NXAA_long * &
                            domain_decomp% NZAA_long * &
                            domain_decomp% NYAA_long / &
                            world_size )

   call c_f_pointer(p7, buf_phys_real, [domain_decomp% NZAA_long,&
                                   domain_decomp% NXAA_long,&
                                   domain_decomp% local_NY_phys]) 

   ! /////////////////////////////////
   ! /////////////////////////////////
   p_c2c_y = fftw_plan_many_dft(1, [domain_decomp% spec_iSize(1)], &
                             domain_decomp% spec_iSize(1) * domain_decomp% local_NX, &
                             buf_spec_cplx, domain_decomp% NYAA, domain_decomp%NZ, &
                             domain_decomp% NZ * domain_decomp%NYAA, 
                             buf_yphys_padded_xDist, domain_decomp% NYAA, domain_decomp%NZ, &
                             domain_decomp% NZ * domain_decomp%NYAA, 
                             FFTW_BACKWARD, FFTW_PATIENT)
   !!! STOPPED HERE
   !!! STOPPED HERE
   !!! STOPPED HERE

   n_fft = [domain_decomp% NYAA_long, domain_decomp% NXAA_long]

   p_c2r = fftw_mpi_plan_many_dft_c2r( 2,& ! rank
                                   n_fft,& ! logical size
                domain_decomp% NZAA_long,& ! how many transforms
                  FFTW_MPI_DEFAULT_BLOCK,& ! let fftw decide
                  FFTW_MPI_DEFAULT_BLOCK,& ! let fftw decide
                                  X_cplx,& ! complex input (TRANSPOSED)
                                  X_real,& ! real output
                          MPI_COMM_WORLD,&
                          ior(FFTW_patient, FFTW_MPI_TRANSPOSED_In))

   p_r2c = fftw_mpi_plan_many_dft_r2c( 2,& ! rank
                                   N_fft,& ! logical size
                domain_decomp% NZAA_long,& ! how many transforms
                  FFTW_MPI_DEFAULT_BLOCK,& ! let fftw decide
                  FFTW_MPI_DEFAULT_BLOCK,& ! let fftw decide
                                  X_real,& ! Real input
                                  X_cplx,& ! Complex output (TRANSPOSED)
                          MPI_COMM_WORLD,&
                          ior(FFTW_patient, FFTW_MPI_TRANSPOSED_Out))



   if (dct_includes_endpoints) then
     ! adjust the logical size of the transforms
     logical_NZ = domain_decomp%NZAA - 1
     ! plan the transforms
     p_r2r_f =  fftw_plan_many_r2r(1, [domain_decomp%phys_iSize(1)], &
         domain_decomp%phys_iSize(2)*domain_decomp%phys_iSize(3), &
         x_real, [domain_decomp%phys_iSize(1)], 1, domain_decomp%phys_iSize(1), &
         x_real, [domain_decomp%phys_iSize(1)], 1, domain_decomp%phys_iSize(1), &
         [FFTW_REDFT00], FFTW_patient)
     p_r2r_i =  fftw_plan_many_r2r(1, [domain_decomp%phys_iSize(1)], &
         domain_decomp%phys_iSize(2)*domain_decomp%phys_iSize(3), &
         x_real, [domain_decomp%phys_iSize(1)], 1, domain_decomp%phys_iSize(1), &
         x_real, [domain_decomp%phys_iSize(1)], 1, domain_decomp%phys_iSize(1), &
         [FFTW_REDFT00], FFTW_patient)
   else
     ! adjust the logical size of the transforms
     logical_NZ = domain_decomp%NZAA
     ! plan the transforms
     p_r2r_f =  fftw_plan_many_r2r(1, [domain_decomp%phys_iSize(1)], &
         domain_decomp%phys_iSize(2)*domain_decomp%phys_iSize(3), &
         x_real, [domain_decomp%phys_iSize(1)], 1, domain_decomp%phys_iSize(1), &
         x_real, [domain_decomp%phys_iSize(1)], 1, domain_decomp%phys_iSize(1), &
         [FFTW_REDFT10], FFTW_patient)
     p_r2r_i =  fftw_plan_many_r2r(1, [domain_decomp%phys_iSize(1)], &
         domain_decomp%phys_iSize(2)*domain_decomp%phys_iSize(3), &
         x_real, [domain_decomp%phys_iSize(1)], 1, domain_decomp%phys_iSize(1), &
         x_real, [domain_decomp%phys_iSize(1)], 1, domain_decomp%phys_iSize(1), &
         [FFTW_REDFT01], FFTW_patient)
   end if

   
   call fftw_free(x_real_ptr)
   call fftw_free(x_cplx_ptr)

   p2 = fftw_alloc_real   ( 2*domain_decomp% alloc_local )

   call c_f_pointer(p2, phys_padded,[domain_decomp% NZAA_long,&
                                     domain_decomp% NXAA_long+2,&
                                     domain_decomp% local_NY_phys]) 

   p3 = fftw_alloc_real   ( 2*domain_decomp% alloc_local )

   call c_f_pointer(p3, phys_buffer,[domain_decomp% NZAA_long,&
                                     domain_decomp% NXAA_long,&
                                     domain_decomp% local_NY_phys]) 

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
