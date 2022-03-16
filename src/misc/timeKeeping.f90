 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: timeKeeping
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> Contains variables related to keeping track of DNS time (as opposed to
 !! wall-clock time. 
 !
 !=============================================================================
module timeKeeping
 use fortran_kinds
 use MPI_vars

implicit none

 type :: timings_T
 !> Contains variables related to keeping track of DNS time (as opposed to
 !! wall-clock time. 
   Integer :: i_timestep !< time-step index
   real(dp) :: dt !< target timestep
   real(dp) :: dt_based_on_CFL !< target timestep (difference with above?)
   real(dp) :: dt_max !< bounds the timestep in use
   real(dp) :: dt_current !< timestep actually used.
   real(dp) :: absolute_time !< time in whatever units have been used
                             !! for non-dimensionalization.
   real(dp) :: tolerance_factor !< do not refactorize matrices if the relative
                                !! timestep variation is less than that
   real(dp) :: cfl_C !< dt = CFL_C * min (dx/u, dy/v, dz/w)
   logical :: CFL_flag !< apply CFL criterion if .True.
   logical :: refactorization_is_due !< .True. : the timestep has been actualized.
                                     !! Matrices refactorizations are in order.
   Real(dp), Dimension(:), Allocatable :: array_of_dts !< MPI communications
   real(dp):: dt_min !< tolerance. Below that value, the code is very likely to
                     !! be blowing-up. Triggers STOP.
   integer :: checkpoint_period=200 !< how often we output checkpoints.
   integer :: checkpoint_counter   !< a dummy variable for checkpoint outputs.
 contains
   procedure :: init      => initialize_timings
   procedure :: increment => increment_timings 
   procedure :: apply_cfl => compare_dt_against_CFL
   procedure :: check_for_update => update_dt_if_needed   
   procedure :: fromfile => read_timings_from_file            
   procedure :: tofile   => write_timings_to_file            
 end type timings_T

 type(timings_T), save :: timings

contains

!> @brief
!> Initialises timing counters and time-step sizes.
subroutine initialize_timings(self, dt_max, tolerance_fac, CFL_C)
  class(timings_T) :: self
  real(kind=dp), intent(in) :: dt_max, tolerance_fac, CFL_C
  self%dt_max = dt_max
  self%tolerance_factor = tolerance_fac
  self%CFL_C = CFL_C
  self%i_timestep = 1
  self%checkpoint_counter = 1
  Allocate(self%array_of_dts(world_size))
  self%dt_current = self%dt_max
  self%dt         = self%dt_max
  self%dt_min = 0.000000001_dp
  self%CFL_flag = .True. 
  self%absolute_time = 0._dp
end subroutine initialize_timings
 
!> @brief
!> Computes CFL-based dt. Returned value is bound by dt_max.
!> The program's execution is halted whenever dt_min is reached.
!> @details
!> In practice, reaching dt_min happens during blows-up anyway...
subroutine compare_dt_against_CFL(self, dx_over_velocity_min)
  class(timings_T) :: self
  real(dp), intent(In) :: dx_over_velocity_min 
                              !< minimum of (dz/w, dy/v, dx/u)
  if (self%CFL_flag) self%dt = dx_over_velocity_min * self%CFL_C
  if (self%dt.gt.self%dt_max) self%dt = self%dt_max
  !if (my_rank.eq.0) print *, dx_over_velocity_min,  self%CFL_C delme
  !if (my_rank.eq.0) print *, self%dt, self%dt_min  delme
  if (self%dt.lt.self%dt_min) then
   stop 'minimal timestep reached'
  end if
end subroutine compare_dt_against_CFL


subroutine read_timings_from_file(self,fileStr, fileLen)
   class(timings_T) :: self
   integer, intent(In) :: fileLen
   character(len=fileLen), intent(In) :: fileStr
   character(len=fileLen+16) :: time_file_str
   logical :: time_file_exists
   integer :: file_size

   ! first assess what is in the file:
   inquire (file=fileStr, size=file_size)
   select case (file_size)
          case (16)
              if (my_rank.eq.0) then
                 open(unit=9, file=fileStr, status='old', access='stream')
                 read(9) self%dt, self%absolute_time
                 close(9)
                 self%dt = self%dt / 50.
              end if
              call MPI_Bcast(self%dt,            1, MPI_DOUBLE, 0, MPI_Comm_world, ierr)
              call MPI_Bcast(self%absolute_time, 1, MPI_DOUBLE, 0, MPI_Comm_world, ierr)
          case (20)
              if (my_rank.eq.0) then
                 open(unit=9, file=fileStr, status='old', access='stream')
                 read(9) self%dt, self%absolute_time, self%i_timestep
                 close(9)
                 self%dt = self%dt / 50.
              end if
              time_file_str = fileStr(1:fileLen-6)//'../Timeseries/time.dat'
              inquire( file=time_file_str, exist=time_file_exists)
              if (.not.(time_file_exists)) self%i_timestep = 1
              call MPI_Bcast(self%dt,            1, MPI_DOUBLE, 0, MPI_Comm_world, ierr)
              call MPI_Bcast(self%absolute_time, 1, MPI_DOUBLE, 0, MPI_Comm_world, ierr)
              call MPI_Bcast(self%i_timestep,    1, MPI_INT,    0, MPI_Comm_world, ierr)
   end select
              

   if (self%dt.le.self%dt_max) then
      self%dt_current = self%dt
   else                                
      self%dt_current = self%dt_max
      self%dt         = self%dt_max
   end if
end subroutine
 

subroutine write_timings_to_file(self,fileStr, fileLen)
   class(timings_T) :: self
   integer, intent(In) :: fileLen
   character(len=fileLen), intent(In) :: fileStr
   ! process 0 takes care of the writing for everybody.
   if (my_rank.eq.0) then
   open(unit=9, file=fileStr, status='replace', access='stream')
   write(9,POS=1) self%dt_current, self%absolute_time, self%i_timestep
   close(9)
   end if
end subroutine

!> @brief
!> Computes CFL-based dt. Returned value is bound by dt_max.
!> The program's execution is halted whenever dt_min is reached.
!> @details
!> In practice, reaching dt_min happens during blows-up anyway...
subroutine increment_timings(self)
  class(timings_T) :: self
  self%i_timestep = self%i_timestep + 1
  self%absolute_time = self%absolute_time + self%dt_current
  self%checkpoint_counter = self%checkpoint_counter+1
end subroutine increment_timings

!> @brief
!> All processes agree on the necessity to actualize the time-step, 
!> based on the CFL-based dt.
!> @details
!> Locally computed CFL-based dt are gathered by process 0, and the 
!> minimum is then broadcast back to all processes. If the relative
!> difference between the target timestep and the timestep currently 
!> used exceeds a user-defined tolerance, then the target timestep is
!> adopted and a flag is switched to indicate that refactorization
!> of the matrices is in order.
subroutine update_dt_if_needed(self)
  class(timings_T) :: self
  Call MPI_Gather(self%dt, 1, MPI_DOUBLE, self%array_of_dts, 1, &
                                 MPI_DOUBLE, 0, MPI_Comm_World, ierr)
  if (my_rank.eq.0) self%dt = minval(self%array_of_dts)
  Call MPI_Bcast(self%dt, 1, MPI_DOUBLE, 0, MPI_Comm_world, ierr)
  If (abs(self%dt-self%dt_current).gt.(self%tolerance_factor*self%dt_current)) Then
      1521 Format ('>>> REFACTORIZATION IS DUE ====================> dt =',ES10.3)
      if (my_rank.eq.0) write(*,1521) self%dt
      self%refactorization_is_due = .True.
      self%dt_current = self%dt
  End If
end subroutine update_dt_if_needed

end module timeKeeping
