 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: PL_geometry
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> TODO
 !>  
 !
 !=============================================================================
module PL_geometry

 use chdir_mod
 use fortran_kinds
 use domain_decomposition
 use MPI_vars
 use read_command_line_args
        
 implicit None
 
 type :: local_shapes_T
   integer :: local_NX, local_NY, local_NZ
 end type local_shapes_T

 type :: subset_communicator_T
    integer :: rank
    integer :: size
    integer :: comm
 end type

 type :: geometry_vars_T
   integer :: NX,  NXAA
   integer :: NY,  NYAA
   integer :: NZ,  NZAA
   type(local_shapes_T) :: phys, spec
   integer :: Nvar
   real(dp) :: Lx, Ly, Lz
   real(dp) :: kx_box
   real(dp) :: ky_box
   real(dp) :: gap, center
   real(dp), allocatable, dimension(:) :: kx
   real(dp), allocatable, dimension(:) :: ky
   complex(dp), allocatable, dimension(:) :: px
   complex(dp), allocatable, dimension(:) :: py
   logical, allocatable, dimension(:) :: deAliase_x
   logical, allocatable, dimension(:) :: deAliase_y
   logical :: this_core_has_zero_mode
   type(subset_communicator_T) :: mpi_zPhys
   type(subset_communicator_T) :: mpi_yPhys
  contains
   procedure :: init => initialize_geometry
   procedure :: read_params => read_input_file_parameters
   procedure :: export_to_disk => export_geometry_to_disk
 end type geometry_vars_T
 
 type :: miscellaneous_cargo_T
   character(len=7) :: scheme_id
   logical :: qsave_restart
   logical :: rolling_qsaves=.False.
   real(kind=dp) :: time_limit
   real(kind=dp) :: dt_max
   real(kind=dp) :: tolerance_factor
   real(kind=dp) :: CFL_C
   character(len=:), allocatable :: output_directory
   integer :: output_dir_length
 end type 

 type(miscellaneous_cargo_T) :: misc_cargo


 contains

 subroutine Initialize_geometry(self)
   class(geometry_vars_T) :: self
   integer :: k_index 
   real(dp) :: kcut
   real(dp), parameter :: pi = 4._dp*atan(1.0_dp)
   integer :: ix, iy

   !> Variables below have already been initialised 
   !! when calling read_input_file_parameters from 
   !! P3_read_parameters module:
   !!   - self%NX
   !!   - self%NY
   !!   - self%NZ
   !!   - self%Lx
   !!   - self%Ly
   !!   - self%Lz
   !!
   call self%read_params()

   self%NXAA = (3 * self%NX) / 2
   self%NYAA = (3 * self%NY) / 2
   self%NZAA = (3 * self%NZ) / 2

   ! if we include endpoints, we need an additional point to ensure that
   ! the logical size of the transforms has a simple prime factors decomposition
   if (is_string_in_the_list('--grid-with-endpoints', 21)  &
       .or.                                                &
       is_string_in_the_list('--gauss-lobatto-grid',  20)) then
       
      self%NZAA = (3*self%NZ) / 2 + 1 
   end if

   if (my_rank.eq.0) print*, '================================================================='
   if (my_rank.eq.0) print*, ' ...'
   if (my_rank.eq.0) print*, ' ... Now initialising fftw3-mpi.'
   if (my_rank.eq.0) print*, ' ... The optimal data distribution is being found'
   if (my_rank.eq.0) print*, ' ...'
   call domain_decomp%init( self%NZAA, &
                            self%NYAA, &
                            self%NXAA  )

   self%spec%local_NX = domain_decomp%spec_iSize(3)                 
   self%spec%local_NY = domain_decomp%spec_iSize(2)                 
   self%spec%local_NZ = domain_decomp%spec_iSize(1)                 

   self%phys%local_NX = domain_decomp%phys_iSize(3)
   self%phys%local_NY = domain_decomp%phys_iSize(2)                 
   self%phys%local_NZ = domain_decomp%phys_iSize(1)                 

   allocate (self%kx (self%spec%local_NX ))
   allocate (self%px (self%spec%local_NX ))
   allocate (self%ky (self%spec%local_NY ))
   allocate (self%py (self%spec%local_NY ))
   allocate (self%deAliase_x (self%spec%local_NX ))
   allocate (self%deAliase_y (self%spec%local_NY ))
  
   ! we determine which processes share a common z in physical space
   call MPI_comm_split(MPI_comm_world,&
                       domain_decomp%phys_iStart(1),& !color: mpi group gathers common z
                       domain_decomp%phys_iStart(2),& !key  : mpi group ranked by y
                       self%MPI_zPhys%comm,&
                       ierr)
   call MPI_Comm_Size(self%MPI_zPhys%comm, self%MPI_zPhys%size, ierr)
   call MPI_Comm_Rank(self%MPI_zPhys%comm, self%MPI_zPhys%rank, ierr) 
   ! we determine which processes share a common y in physical space
   call MPI_comm_split(MPI_comm_world,&
                       domain_decomp%phys_iStart(2),& !color: mpi group gathers common y
                       domain_decomp%phys_iStart(1),& !key  : mpi group ranked by z
                       self%MPI_yPhys%comm,&
                       ierr)
   call MPI_Comm_Size(self%MPI_yPhys%comm, self%MPI_yPhys%size, ierr)
   call MPI_Comm_Rank(self%MPI_yPhys%comm, self%MPI_yPhys%rank, ierr) 

   if ((domain_decomp%spec_iStart(1)+ &
        domain_decomp%spec_iStart(2)+ &
        domain_decomp%spec_iStart(3)) .eq. 3) then
        self%this_core_has_zero_mode = .true.
   else
        self%this_core_has_zero_mode = .false.
   end if

   
   !!
   !> next, initialize the y-wavenumbers
   !!
   self%ky_box = 2._dp * pi / self%Ly
   kcut = self%ky_box * self%NY / 2._dp
   do iy = 1, self%spec%local_NY
      k_index = domain_decomp%spec_iStart(2) - 1 + iy
      !> array of wave numbers
      if (k_index.ge.(self%NYAA/2+1)) then
         self%ky (iy) =-real( &
                        (self%NYAA + 1 - k_index) * self%ky_box ,& 
                        kind=dp )
      else
         self%ky (iy) = real( (k_index-1) * self%ky_box ,& 
                        kind=dp )
      end if
      !> dealiasing (set the mode to zero if True)
      if ( abs(self%ky (iy)) .ge. kcut) then
         self%deAliase_y (iy) = .True.
      else 
         self%deAliase_y (iy) = .False.
      end if
   end do

   !! 
   !> finally, initialize the x-wavenumbers
   !!
   self%kx_box = 2._dp * pi / self%Lx
   kcut = self%kx_box * self%NX / 2._dp
   do ix = 1, self%spec%local_NX
      k_index = domain_decomp%spec_iStart(3) - 1 + ix
      !> array of wave numbers
      if (k_index.ge.(self%NXAA/2+1)) then
         self%kx (ix) =-real( &
                        (self%NXAA + 1 - k_index) * self%kx_box ,& 
                        kind=dp )
      else
         self%kx (ix) = real( (k_index-1) * self%kx_box ,& 
                        kind=dp )
      end if
      !> dealiasing (set the mode to zero if True)
      if ( abs(self%kx (ix)) .ge. kcut) then
         self%deAliase_x (ix) = .True.
      else 
         self%deAliase_x (ix) = .False.
      end if
   end do
   
   self%px = Cmplx(0._dp, self%kx, kind=dp)
   self%py = Cmplx(0._dp, self%ky, kind=dp)
   if (my_rank.eq.0) then
   print*, '================================================================='
   print *, " Initializing geometry ... (PL_geometry.f90)"
   print *, 'gridSize before dealiasing :', self%NXAA, self%NYAA, self%NZAA
   print *, 'gridSize after  dealiasing :', self%NX, self%NY, self%NZ
   print *, "                       ... done. "
   print*, '================================================================='
   end if

   call self%export_to_disk()

 end subroutine Initialize_geometry





 subroutine export_geometry_to_disk(self)
   class(geometry_vars_T) :: self
   character (len = 4) :: corestr               
   character (len =17) :: file_str             
   
   call chdir('./Geometry/')
   write(corestr,"(i4.4)") my_rank
   file_str = "domDecmp.core"//corestr
   Open (Unit=9, File=file_str, Status='replace', Access='stream')
   Write(9) self%NX, self%NY, self%NZ, self%NXAA, self%NYAA, self%NZAA, &
            self%spec%local_NX, self%spec%local_NY, &!self%spec%local_NZ, 
            self%Nvar, &
            domain_decomp%phys_iStart, &
            domain_decomp%phys_iSize, &
            domain_decomp%phys_iEnd, &
            domain_decomp%spec_iStart, &
            domain_decomp%spec_iSize, &
            domain_decomp%spec_iEnd
   Close (Unit=9)
   file_str = "geometry.core"//corestr
   Open (Unit=9, File=file_str, Status='replace', Access='stream')
   Write(9) self%Lx,     self%Ly,     self%gap, self%center, &
            self%kx_box, self%ky_box, &
            self%kx,     self%ky
   Close (Unit=9)
 end subroutine export_geometry_to_disk

 subroutine read_input_file_parameters(geom)
   class(geometry_vars_T), intent(inOut) :: geom
   integer :: empty_entry
   character(len=:), Allocatable :: temporary_string
   
   misc_cargo%output_dir_length=200
   allocate (Character(misc_cargo%output_dir_length) :: temporary_string)
   open(3,File='coral.parameters.in',Status='unknown')
   read(3,*) empty_entry        ! ........... RESOLUTION ............ !
   read(3,*) geom%NX            ! 01 ! (INT) ... along X         ! 01 !
   read(3,*) geom%NY            ! 02 ! (INT) ... along Y         ! 02 !
   read(3,*) geom%NZ            ! 03 ! (INT) ... along Z         ! 03 !
   read(3,*) empty_entry        ! ............ GEOMETRY ............. !
   read(3,*) geom%Lx            ! 05 ! (DBL) x/z aspect ratio    ! 05 !
   read(3,*) geom%Ly            ! 06 ! (DBL) y/z aspect ratio    ! 06 !
   read(3,*) geom%gap           ! 07 ! (DBL) gap (should=1.0_dp) ! 07 !
   read(3,*) geom%center        ! 08 ! (DBL) center              ! 08 !
   read(3,*) empty_entry        ! .......... TIMESTEPPING ........... !
   read(3,*) misc_cargo%dt_max           ! 10 ! (DBL) maximum timestep    ! 10 !
   read(3,*) misc_cargo%tolerance_factor ! 11 ! (DBL) timestep tolerance  ! 11 !
   read(3,*) misc_cargo%CFL_C            ! 12 ! (DBL) CFL Factor          ! 12 !
   read(3,*) misc_cargo%scheme_id        ! 13 ! scheme                    ! 13 !
   read(3,*) empty_entry               ! ......... QUICKSAVE RESTART ....... !
   read(3,*) misc_cargo%qsave_restart  ! 15 ! (BOO) quicksave restart   ! 15 !
   read(3,*) empty_entry               ! ............ TIMER ................ !
   read(3,*) misc_cargo%time_limit     ! 17 ! (DBL) Time limit (in sec.)! 17 !
   close(3)      

   call get_environment_variable('PWD',temporary_string)

   misc_cargo%output_dir_length = len(trim(temporary_string))

   allocate (Character(misc_cargo%output_dir_length) :: misc_cargo%output_directory)
   misc_cargo%output_directory = trim(temporary_string)
   deAllocate(temporary_string)
 
   End Subroutine read_input_file_parameters

end module PL_geometry
