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
 !> Domain decomposition for linking with fftw3-mpi library.
 !>  
 !
 !=============================================================================
 Module domain_decomposition

   use fortran_kinds
   use MPI_vars
   use fftw3mpi


   Implicit None

   Type :: domain_decomposition_T
      Integer :: p_row = 0
      Integer :: p_col = 0
      Integer, dimension(3) :: phys_iStart, phys_iEnd, phys_iSize
      Integer, dimension(3) :: spec_iStart, spec_iEnd, spec_iSize
      integer :: NXAA, NYAA, NZAA
      integer(kind=C_intPtr_T) :: NXAA_long, NYAA_long, NZAA_long
      integer(kind=C_intPtr_T) :: alloc_local
      integer(kind=C_intPtr_T) :: local_mx_offset
      integer(kind=C_intPtr_T) :: local_ky_offset
      integer(kind=C_intPtr_T) :: local_NX_spec
      integer(kind=C_intPtr_T) :: local_NY_phys  
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

   Contains


   Subroutine domain_decomp_init(self, n1, n2, n3)

      class(domain_decomposition_T) :: self
      integer, intent(in) :: n1, n2, n3
      !-  internal vars
      integer (C_intPtr_T) :: n_fft(2)
      integer (C_intPtr_T) :: howmany

      self% NXAA = n3
      self% NYAA = n2
      self% NZAA = n1
      self% NXAA_long = int( self% NXAA, kind= C_intPtr_T)
      self% NYAA_long = int( self% NYAA, kind= C_intPtr_T)
      self% NZAA_long = int( self% NZAA, kind= C_intPtr_T)

      howmany =   self% NZAA_long
      n_fft   = [ self% NYAA_long, self% NXAA_long/2+1 ] 

      call fftw_mpi_init()

      self% alloc_local = fftw_mpi_local_size_many_transposed (2, n_fft, howmany, &
                                                  FFTW_MPI_DEFAULT_BLOCK, & 
                                                  FFTW_MPI_DEFAULT_BLOCK, & 
                                                  MPI_COMM_WORLD, &             
                                                  self% local_NY_phys,&
                                                  self% local_ky_offset,&
                                                  self% local_NX_spec,&
                                                  self% local_mx_offset )
      self% p_row = 1
      self% p_col = world_size
      !/////////////////////////////////////////////
      !
      !      Spectral space distributions 
      !
      !/////////////////////////////////////////////
      !-- In spectral space, the fast index z (1) is in-core
      self% spec_iStart(1) = 1 
      self% spec_iSize (1) = self% NZAA
      self% spec_iEnd  (1) = self% NZAA
      !-- In spectral space, the intermediate index y (2) is in-core
      self% spec_iStart(2) = 1 
      self% spec_iSize (2) = self% NYAA
      self% spec_iEnd  (2) = self% NYAA
      !-- In spectral space, the slow index x (3) is distributed    
      self% spec_iStart(3) = self% local_mx_offset + 1
      self% spec_iSize (3) = self% local_NX_spec
      self% spec_iEnd  (3) = self% local_NX_spec

      !/////////////////////////////////////////////
      !
      !      Physical space distributions 
      !
      !/////////////////////////////////////////////
      !-- In physical space, the fast index z (1) is in-core
      self% phys_iStart(1) = 1 
      self% phys_iSize (1) = self% NZAA
      self% phys_iEnd  (1) = self% NZAA
      !-- In physical space, the intermediate index x (2) is in-core
      self% phys_iStart(2) = 1 
      self% phys_iSize (2) = self% NXAA
      self% phys_iEnd  (2) = self% NXAA
      !-- In physical space, the slow index y (3) is distributed    
      self% phys_iStart(3) = self% local_ky_offset + 1
      self% phys_iSize (3) = self% local_NY_phys
      self% phys_iEnd  (3) = self% local_NY_phys

   202 format ('Core ', (i4.4), ' has ', (i4.4), ' x modes and ',(i4.4),' y gridpoints (offsets: ', (i4.4),',',(i4.4))

   write (*,202) my_rank, self% local_NX_spec, self% local_NY_phys, self% local_mx_offset, self% local_ky_offset

   End Subroutine domain_decomp_init


 End Module domain_decomposition
