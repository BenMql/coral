 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: domain_decomposition
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> Domain decomposition for linking with 2decomp library
 !>  
 !
 !=============================================================================
 Module domain_decomposition

   use fortran_kinds
   use decomp_2d
   use decomp_2d_fft
   use decomp_2d_mpi
   use decomp_2d_constants
   use MPI_vars


   Implicit None

   Type :: domain_decomposition_T
      Integer :: p_row = 0
      Integer :: p_col = 0
      Integer, dimension(3) :: phys_iStart, phys_iEnd, phys_iSize
      Integer, dimension(3) :: spec_iStart, spec_iEnd, spec_iSize
      integer :: NXAA, NYAA, NZAA
      logical :: x_is_empty_in_spec = .false.
      logical :: y_is_empty_in_phys = .false.
    Contains
      Procedure :: init => domain_decomp_init
   End Type domain_decomposition_T            

   Type(domain_decomposition_T), save :: domain_decomp

   !type(decomp_info) :: decomplib_setup

   !> This module insulates the rest of the program from the
   !! variables defined in the decomp_2d library.
   !! Internal variables are accessed only through the
   !! variable domain_decomp

   Private

   Public :: domain_decomp
   character(len=8), Public :: parallel_transforms_library = 'decomp2d'

   Contains


   Subroutine domain_decomp_init(self, n1, n2, n3)

      Class(domain_decomposition_T) :: self
      Integer, intent(in) :: n1, n2, n3
      self%NXAA = n3
      self%NYAA = n2
      self%NZAA = n1
      self%p_row = 0
      self%p_col = 0
      call decomp_2d_init(n1, n2, n3, self%p_row, self%p_col)
      self%phys_istart = zstart !< copy 2decomp's internal var.
      self%phys_iEnd   = zend   !< copy 2decomp's internal var.
      self%phys_iSize  = zsize  !< copy 2decomp's internal var. 
      call decomp_2d_fft_init(PHYSICAL_IN_Z, &
                              opt_skip_XYZ_c2c=[.True., .False., .False.])
      call decomp_2d_fft_get_size(self%spec_iStart, &
                                  self%spec_iEnd,   &
                                  self%spec_iSize   )

   End Subroutine domain_decomp_init


 End Module domain_decomposition
