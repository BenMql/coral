 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: SL_geometry
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> TODO
 !>  
 !
 !=============================================================================
module SL_geometry

 use chdir_mod
 use SL_domain_decomp

 implicit None
 
 type :: local_shapes_T
   integer :: local_N3, local_N2, local_N1
 end type local_shapes_T

 type :: subset_communicator_T
    integer :: rank
    integer :: size
    integer :: comm
 end type

 type, extends(domain_decomposition_T) :: geometry_vars_T
   integer :: NX
   integer :: NY
   integer :: NZ
   type(local_shapes_T) :: phys, spec, linAlg
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

 type(geometry_vars_T) :: geom

 contains

 subroutine Initialize_geometry(self)
   class(geometry_vars_T) :: self
   real(dp) :: kcut
   real(dp), parameter :: pi = 4._dp*atan(1.0_dp)
   integer :: ix, iy, global_iy

   !> Variables below are initialised 
   !! when reading the parameters file:
   !!   - self%NX
   !!   - self%NY
   !!   - self%NZ
   !!   - self%Lx
   !!   - self%Ly
   !!   - self%Lz
   !!
   call self%read_params()

   if (my_rank.eq.0) print*, '================================================================='
   if (my_rank.eq.0) print*, ' ...'
   if (my_rank.eq.0) print*, ' ... Now calling fftw3_mpi.'   
   if (my_rank.eq.0) print*, ' ... Initializing domain decomposition.'
   if (my_rank.eq.0) print*, ' ...'
   call self%decomp_init( (3*self%NZ)/2, &
                          (3*self%NY)/2, &
                          (3*self%NX)/2  )

   self%linAlg%local_N1 = self%NZ
   self%linAlg%local_N2 = self%NX
   self%linAlg%local_N3 = self%localSize_spec
   self% spec %local_N1 = self%NZAA
   self% spec %local_N2 = self%NXAA
   self% spec %local_N3 = self%localSize_spec
   self% phys %local_N1 = self%NZAA
   self% phys %local_N2 = self%NYAA + 2
   self% phys %local_N3 = self%localSize_phys

   allocate (self%kx (self%NX ))                
   allocate (self%ky (self%localSize_spec ))                
   allocate (self%px (self%NX ))                
   allocate (self%py (self%localSize_spec ))                
   allocate (self%deAliase_x (self%NXAA ))          
   allocate (self%deAliase_y (self%localSize_spec ))
  

   if (self%spec_offset .eq. 0 ) then
        self%this_core_has_zero_mode = .true.
   else
        self%this_core_has_zero_mode = .false.
   end if

   self%deAliase_y = .False.
   self%deAliase_x = .true. 
   
                        !
                        !
   !##########################################!
   ! Wave Numbers Along y:                    !
   !..........................................!
   self%ky_box = 2._dp * pi / self%Ly
   kcut = self%ky_box * self%NY / 2._dp
   do iy = 1, self%localSize_spec
      global_iy = self%spec_offset + iy
      self%ky (iy) = self%ky_box*(global_iy-1)
      if (global_iy.gt.(self%NY/2+1)) Then
         self%deAliase_y(iy) = .True.
      end if
   end do
   !..........................................!
   !##########################################!
                        !
                        !
   !##########################################!
   !..........................................!
   ! Wave Numbers Along x:                    ! 
   !..........................................! 
   ! Constant component k=0:                  ! 
   !..........................................! 
   self%kx(1)  = 0._dp 
   !..........................................! 
   ! Nyquist along y:                         ! 
   !..........................................! 
   ix = self%NX/2+1 
   self%kx (ix) = Real(ix-1, Kind=dp)*self%kx_box 
   !..........................................! 
   ! Other components:                        ! 
   !..........................................! 
   do ix = 2,(self%Nx/2) 
      self%kx (      ix) = Real(ix-1, Kind=dp)*self%kx_box 
      self%kx (self%Nx+2-ix) =-Real(ix-1, Kind=dp)*self%kx_box 
   End do
   !..........................................!
   !##########################################!

   self%px = cmplx( 0._dp, self%kx, kind=dp)
   self%py = cmplx( 0._dp, self%ky, kind=dp)

   if (my_rank.eq.0) then
   print*, '================================================================='
   print *, " Initializing geometry ... (SL_geometry.f90)"
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
            self%localSize_spec,self%localSize_phys,&!self%spec%local_NZ, 
            self%Nvar
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

end module SL_geometry
