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
 !>        
 !
 !=============================================================================
 Module domain_decomposition

   use fortran_kinds
   use MPI_vars
   use fftw3mpi


   Implicit None

   Type :: domain_decomposition_T
      Integer, dimension(2) :: phys_iStart, phys_iEnd, phys_iSize
      Integer, dimension(2) :: spec_iStart, spec_iEnd, spec_iSize
      integer :: NXAA, NZAA
      integer :: NX, NZ
      integer(kind=C_intPtr_T) :: NXAA_long, NZAA_long
      integer(kind=C_intPtr_T) :: NX_long, NZ_long
      integer(kind=C_intPtr_T) :: alloc_local
      integer(kind=C_intPtr_T) :: local_mx_offset
      integer(kind=C_intPtr_T) :: local_Tn_offset
      integer(kind=C_intPtr_T) :: local_NX_spec
      integer(kind=C_intPtr_T) :: local_NZ_phys  
      type(C_ptr) :: plan_forward_transpose 
      type(C_ptr) :: plan_inverse_transpose 
    Contains
      Procedure :: init => domain_decomp_init
   End Type domain_decomposition_T            

   Type(domain_decomposition_T), save :: domain_decomp


   Private

   Public :: domain_decomp
   character(len=8), Public :: parallel_transforms_library = 'fftw-mpi'

   Contains


   Subroutine domain_decomp_init(self, n1, n3)

      class(domain_decomposition_T) :: self
      integer, intent(in) :: n1, n3
      !-  internal vars
      type(c_ptr) :: cplx_ptr_1, cplx_ptr_2
      Real(C_double), pointer :: X_in_core(:,:), z_in_core(:,:)
      


      ! start by checking that the number of processes divides n3/3 and n2
      if ((mod(n3/3, world_size) .ne. 0 ) &
                                 .or.     & 
          (mod(  n1, world_size) .ne. 0 ) ) then
          call MPI_finalize(ierr)
          print *, " ----------------- WRONG PARAMETERS PROVIDED AT RUNTIME ----------------"
          print *, "Coral ribbon version: number of processes incompatible with the resolution"
          print *, "Please run with a number of cores that divides NX/2 *and* 3.NZ/2, ",&
                   "where NX and NZ are the values entered in coral.parameters.in"
          print *, "Terminating now. " 
          print *, " ----------------- WRONG PARAMETERS PROVIDED AT RUNTIME ----------------"
          stop
      end if
      self% NXAA = n3
      self% NZAA = n1
      self% NX = n3*2/3
      self% NZ = n1*2/3
      self% NXAA_long = int( self% NXAA, kind= C_intPtr_T)
      self% NZAA_long = int( self% NZAA, kind= C_intPtr_T)
      self% NX_long = int( self% NX, kind= C_intPtr_T)
      self% NZ_long = int( self% NZ, kind= C_intPtr_T)

      self% alloc_local = (self%NXAA_long + 2) * self%NZAA_long / world_size


      call fftw_mpi_init()

           
      cplx_ptr_1 = fftw_alloc_real( domain_decomp% alloc_local )
      cplx_ptr_2 = fftw_alloc_real( domain_decomp% alloc_local )

      call c_f_pointer(cplx_ptr_1, x_in_core, [domain_decomp% NZAA_long / world_size,&
                                    domain_decomp% NX_long]) 

      call c_f_pointer(cplx_ptr_2, z_in_core, &
                                   [domain_decomp% NX_long / world_size, &
                                    domain_decomp% NZAA_long ]) 
      
      
      self%plan_forward_transpose = fftw_mpi_plan_many_transpose( &
                                         self% NX_long/2, & !n0
                                         self% NZAA_long, & !n1
                                         int(2, kind=c_intPtr_T),  & !howmany
                                         self% NX_long/2/ world_size, &
                                         self% NZAA_long/ world_size, &
                                         z_in_core, &
                                         X_in_core, &
                                         MPI_COMM_WORLD, &             
                                         FFTW_PATIENT)
             
      self% plan_inverse_transpose = fftw_mpi_plan_many_transpose( &
                                         self% NzAA_long, & !n1
                                         self% NX_long/2, & !n0
                                         int(2, kind=c_intPtr_T),  & !howmany
                                         self% NZAA_long/ world_size, &
                                         self% NX_long/2/ world_size, &
                                         X_in_core, &
                                         Z_in_core, &
                                         MPI_COMM_WORLD, &             
                                         FFTW_PATIENT)


      call fftw_free(cplx_ptr_1)
      call fftw_free(cplx_ptr_2)


      !/////////////////////////////////////////////
      !
      !      Spectral space distributions 
      !
      !/////////////////////////////////////////////
      !-- In spectral space, the fast index z (1) is in-core
      self% spec_iStart(1) = 1 
      self% spec_iSize (1) = self% NZ
      self% spec_iEnd  (1) = self% NZ
      !-- In spectral space, the slow index x (2) is distributed    
      self% spec_iStart(2) = self% NX/2/world_size * my_rank + 1
      self% spec_iSize (2) = self% NX/2/world_size
      self% spec_iEnd  (2) = self% NX/2/world_size

      self% local_NX_spec   = self% NX/2/world_size
      self% local_mx_offset = self% NX/2/world_size * my_rank

      !/////////////////////////////////////////////
      !
      !      Physical space distributions 
      !
      !/////////////////////////////////////////////
      !-- In physical space, the fast index x (1) is in-core
      self% phys_iStart(1) = 1 
      self% phys_iSize (1) = self% NXAA
      self% phys_iEnd  (1) = self% NXAA
      !-- In physical space, the slow index z (2) is distributed    
      self% phys_iStart(2) =  self% NZAA / world_size * my_rank + 1
      self% phys_iSize (2) =  self% NZAA / world_size       
      self% phys_iEnd  (2) =  self% NZAA / world_size       

      self% local_NZ_phys   = self% NZAA / world_size
      self% local_TN_offset = self% NZAA / world_size * my_rank


   202 format ('Core ', (i4.4), ' has ', (i4.4), ' x modes and ',(i4.4),' y gridpoints (offsets: ', (i4.4),',',(i4.4))

   write (*,202) my_rank, self% local_NX_spec, self% local_NZ_phys, self% local_mx_offset, self% local_tn_offset


   End Subroutine domain_decomp_init


 End Module domain_decomposition
