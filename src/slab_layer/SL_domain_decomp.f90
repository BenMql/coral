 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: SL_domain_decomp                               
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> TODO
 !>  
 !
 !=============================================================================
 Module SL_domain_decomp

   use fortran_kinds
   use fftw3mpi 
   use MPI_vars


   Implicit None

   type :: domain_decomposition_T
      integer(C_intPtr_T) :: alloc_local
      integer(C_intPtr_T) :: localSize_spec !formerly local_NX
      integer(C_intPtr_T) :: localSize_phys !formerly local_NY
      integer(C_intPtr_T) :: spec_offset    !formerly local_mx_offset
      integer(C_intPtr_T) :: phys_offset    !formerly local_ky_offset
      integer :: NXAA, NYAA, NZAA
    Contains
      procedure :: decomp_init => domain_decomp_init
   end type domain_decomposition_T            

   type(domain_decomposition_T), save :: domain_decomp

   

   Contains


   subroutine domain_decomp_init(self, n1, n2, n3)

      class(domain_decomposition_T) :: self
      integer, intent(in) :: n1, n2, n3
      integer(C_intPtr_T) :: n_fft(2)
      integer(C_intPtr_T) :: howmany

      self%NXAA = n3
      self%NYAA = n2
      self%NZAA = n1

      howmany = self%NZAA
      n_fft = [self%NXAA, self%NYAA/2+1]

      self%alloc_local = fftw_mpi_local_size_many_transposed(&
                       2, n_fft, howmany,& 
                       fftw_mpi_default_block,&
                       fftw_mpi_default_block,&
                       mpi_comm_world,&
                       self%localSize_phys, self%phys_offset,&
                       self%localSize_spec, self%spec_offset)

   end subroutine domain_decomp_init
      

 End Module SL_domain_decomp
