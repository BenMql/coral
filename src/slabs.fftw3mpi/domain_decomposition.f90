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
 !> We hope to decrease communication costs by padding the directions just before
 !> the transform is computed. 
 !> 
 !> (NZ, 3/2 NY, NX/2 + 1), complex,   x-distributed
 !>        ||
 !>        ||   DFT along Y
 !>        ||    transpose                
 !>        \/
 !> (NZ, 3/2 NY, NX/2 + 1), complex,   y-distributed
 !>        ||
 !>        ||      pad X                    
 !>        \/
 !> (NZ, 3/2 NY, 3/4 NX + 1), complex, y-distributed
 !>        ||
 !>        ||  c2r DFT along X          
 !>        \/
 !> (NZ, 3/2 NY, 3/2 NX), real,        y-distributed
 !>        ||
 !>        ||      pad Z                    
 !>        \/
 !> (3/2 NZ, 3/2 NY, 3/2 NX), real,    y-distributed
 !>        ||
 !>        ||   DCT along Z              
 !>        \/
 !> (3/2 NZ, 3/2 NY, 3/2 NX), real,    y-distributed
 !>        
 !> for NX = NY = NZ, complexity is (3/2 + 9/4 + 27/8) N**3 log (N) 
 !>        
 !> Thus complexity is down by a factor 57/81 ~ 70%
 !>        
 !> In addition (and perhaps more importantly), a factor 9/4 is saved for
 !> communications!
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
      Integer, dimension(3) :: zyx_iStart, zyx_iEnd, zyx_iSize !< the transposed shapes
      Integer, dimension(3) :: phys_iStart, phys_iEnd, phys_iSize
      Integer, dimension(3) :: spec_iStart, spec_iEnd, spec_iSize
      integer :: NXAA, NYAA, NZAA
      integer :: NX, NY, NZ
      integer(kind=C_intPtr_T) :: NYAA_long, NXAA_long, NZAA_long
      integer(kind=C_intPtr_T) :: NX_long, NZ_long
      integer(kind=C_intPtr_T) :: alloc_local
      integer(kind=C_intPtr_T) :: local_mx_offset
      integer(kind=C_intPtr_T) :: local_ky_offset
      integer(kind=C_intPtr_T) :: local_NX_spec
      integer(kind=C_intPtr_T) :: local_NY_phys  
      logical :: x_is_empty_in_spec
      logical :: y_is_empty_in_phys
      type(C_ptr) :: plan_forward_transpose !< in:  third  dimension distributed 
                                            !< out: second dimension distributed
      type(C_ptr) :: plan_inverse_transpose !< in:  second dimension distributed 
                                            !< out: third  dimension distributed
      type(C_ptr) :: plan_full_fields_transpose_zxy_to_zyx
      type(C_ptr) :: plan_full_fields_transpose_zyx_to_zxy
                                                
    Contains
      Procedure :: init => domain_decomp_init
   End Type domain_decomposition_T            

   Type(domain_decomposition_T), save :: domain_decomp


   Private

   Public :: domain_decomp
   character(len=8), Public :: parallel_transforms_library = 'fftw-mpi'

   Contains


   Subroutine domain_decomp_init(self, n1, n2, n3)

      class(domain_decomposition_T) :: self
      integer, intent(in) :: n1, n2, n3
      !-  internal vars
      type(c_ptr) :: cplx_ptr_1, cplx_ptr_2
      Real(C_double), pointer :: X_in_core(:,:,:), Y_in_core(:,:,:)
      


      ! start by checking that the number of processes divides n3/3 and n2
      if ((mod(n3/3, world_size) .ne. 0 ) &
                                 .or.     & 
          (mod(  n2, world_size) .ne. 0 ) ) then
          call MPI_finalize(ierr)
          print *, " ----------------- WRONG PARAMETERS PROVIDED AT RUNTIME ----------------"
          print *, "Coral slab version: number of processes incompatible with the resolution"
          print *, "Please run with a number of cores that divides NX/2 *and* 3.NY/2, ",&
                   "where NX and NY are the values entered in coral.parameters.in"
          print *, "Terminating now. " 
          print *, " ----------------- WRONG PARAMETERS PROVIDED AT RUNTIME ----------------"
          stop
      end if
      self% NXAA = n3
      self% NYAA = n2
      self% NZAA = n1
      self% NX = n3*2/3
      self% NY = n2*2/3
      self% NZ = n1*2/3
      self% NXAA_long = int( self% NXAA, kind= C_intPtr_T)
      self% NYAA_long = int( self% NYAA, kind= C_intPtr_T)
      self% NZAA_long = int( self% NZAA, kind= C_intPtr_T)
      self% NX_long = int( self% NX, kind= C_intPtr_T)
      self% NZ_long = int( self% NZ, kind= C_intPtr_T)

      self% p_row = 1
      self% p_col = world_size

      self% alloc_local = self% NYAA * self%NXAA_long * self%NZAA_long / world_size


      call fftw_mpi_init()

           
      cplx_ptr_1 = fftw_alloc_real( domain_decomp% alloc_local )
      cplx_ptr_2 = fftw_alloc_real( domain_decomp% alloc_local )

      call c_f_pointer(cplx_ptr_1, x_in_core, [domain_decomp% NZAA_long,&
                                    domain_decomp% NYAA_long/ world_size,&
                                    domain_decomp% NXAA_long]) 

      call c_f_pointer(cplx_ptr_2, y_in_core, [domain_decomp% NZAA_long,&
                                    domain_decomp% NYAA_long,&
                                    domain_decomp% NXAA_long / world_size]) 
      
      
      self%plan_forward_transpose = fftw_mpi_plan_many_transpose( &
                                         self% NX_long/2, & !n0
                                         self% NYAA_long, & !n1
                                         2*self%NZ_long,  & !howmany
                                         self% NX_long/2/ world_size, &
                                         self% NYAA_long/ world_size, &
                                         Y_in_core, &
                                         X_in_core, &
                                         MPI_COMM_WORLD, &             
                                         FFTW_PATIENT)
             
      !transpose_plan = fftw_mpi_plan_many_transpose( &
                                         !self% NX_long/2, & !n0
                                         !self% NYAA_long, & !n1
                                         !2*self%NZ_long,  & !howmany
                                         !self% NX_long/2/ world_size, &
                                         !self% NYAA_long/ world_size, &
                                         !Y_in_core, &
                                         !X_in_core, &
                                         !MPI_COMM_WORLD, &             
                                        !ior( FFTW_PATIENT, FFTW_MPI_Transposed_out))
             
      self% plan_inverse_transpose = fftw_mpi_plan_many_transpose( &
                                         self% NYAA_long, & !n1
                                         self% NX_long/2, & !n0
                                         2*self%NZ_long,  & !howmany
                                         self% NYAA_long/ world_size, &
                                         self% NX_long/2/ world_size, &
                                         X_in_core, &
                                         Y_in_core, &
                                         MPI_COMM_WORLD, &             
                                         FFTW_PATIENT)

      self% plan_full_fields_transpose_zxy_to_zyx = fftw_mpi_plan_many_transpose( &
                                         self% NYAA_long, & !n1
                                         self% NXAA_long, & !n0
                                         self% NZAA_long,  & !howmany
                                         self% NYAA_long/ world_size, &
                                         self% NXAA_long/ world_size, &
                                         X_in_core, &
                                         Y_in_core, &
                                         MPI_COMM_WORLD, &             
                                         FFTW_PATIENT)

      self% plan_full_fields_transpose_zyx_to_zxy = fftw_mpi_plan_many_transpose( &
                                         self% NXAA_long, & !n1
                                         self% NYAA_long, & !n0
                                         self% NZAA_long,  & !howmany
                                         self% NXAA_long/ world_size, &
                                         self% NYAA_long/ world_size, &
                                         X_in_core, &
                                         Y_in_core, &
                                         MPI_COMM_WORLD, &             
                                         FFTW_PATIENT)

      call fftw_free(cplx_ptr_1)
      call fftw_free(cplx_ptr_2)

      !inverse_transpose_plan = fftw_mpi_plan_many_transpose( &
                                         !self% NYAA_long, & !n1
                                         !self% NX_long/2, & !n0
                                         !2*self%NZ_long,  & !howmany
                                         !self% NYAA_long/ world_size, &
                                         !self% NX_long/2/ world_size, &
                                         !X_in_core, &
                                         !Y_in_core, &
                                         !MPI_COMM_WORLD, &             
                                         !ior(FFTW_PATIENT, FFTW_MPI_transposed_in))
      !do ix = 1, self% NX_long/2/world_size
      !do iy = 1, self% NYAA
      !do iz = 1, self% NZ_long
      !    Y_in_core_complex(iz,iy,ix) = cmplx(iz + iy*100 + (ix+my_rank*self% NX_long/2/world_size)*10000,&
      !                                        iz + iy*100 + (ix+my_rank*self% NX_long/2/world_size)*10000, kind=C_double)
      !end do
      !end do
      !end do
!
      !print *, my_rank, Y_in_core_complex
!
!      call fftw_execute_r2r(self% plan_forward_transpose, Y_in_core, X_in_core!)

     ! print *, my_rank, X_in_core_complex

      !call fftw_execute_r2r(self% plan_inverse_transpose, X_in_core, Y_in_core)
      !print *, my_rank, Y_in_core_complex
!

      !/////////////////////////////////////////////
      !
      !      Spectral space distributions 
      !
      !/////////////////////////////////////////////
      !-- In spectral space, the fast index z (1) is in-core
      self% spec_iStart(1) = 1 
      self% spec_iSize (1) = self% NZ
      self% spec_iEnd  (1) = self% NZ
      !-- In spectral space, the intermediate index y (2) is in-core
      self% spec_iStart(2) = 1 
      self% spec_iSize (2) = self% NYAA
      self% spec_iEnd  (2) = self% NYAA
      !-- In spectral space, the slow index x (3) is distributed    
      self% spec_iStart(3) = self% NX/2/world_size * my_rank + 1
      self% spec_iSize (3) = self% NX/2/world_size
      self% spec_iEnd  (3) = self% NX/2/world_size

      self% local_NX_spec   = self% NX/2/world_size
      self% local_mx_offset = self% NX/2/world_size * my_rank

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
      self% phys_iStart(3) =  self% NYAA / world_size * my_rank + 1
      self% phys_iSize (3) =  self% NYAA / world_size       
      self% phys_iEnd  (3) =  self% NYAA / world_size       

      self% local_NY_phys   = self% NYAA / world_size
      self% local_ky_offset = self% NYAA / world_size * my_rank

      !/////////////////////////////////////////////
      !
      !      Transposed distributions for output 
      !
      !/////////////////////////////////////////////
      !-- For transposed fields, the fast index z (1) is in-core
      self% zyx_iStart(1) = 1 
      self% zyx_iSize (1) = self% NZAA
      self% zyx_iEnd  (1) = self% NZAA
      !-- For transposed fields, the intermediate index y (2) is in-core
      self% zyx_iStart(2) = 1        
      self% zyx_iSize (2) = self% NYAA 
      self% zyx_iEnd  (2) = self% NYAA 
      !-- For transposed fields, the slow index x (3) is distributed    
      self% zyx_iStart(3) =  self% NXAA / world_size * my_rank + 1
      self% zyx_iSize (3) =  self% NXAA / world_size       
      self% zyx_iEnd  (3) =  self% NXAA / world_size       

   202 format ('Core ', (i4.4), ' has ', (i4.4), ' x modes and ',(i4.4),' y gridpoints (offsets: ', (i4.4),',',(i4.4))

   write (*,202) my_rank, self% local_NX_spec, self% local_NY_phys, self% local_mx_offset, self% local_ky_offset

   if (self% local_NX_spec .eq. 0) self% x_is_empty_in_spec = .True.
   if (self% local_NY_phys .eq. 0) self% y_is_empty_in_phys = .True.

   End Subroutine domain_decomp_init


 End Module domain_decomposition
