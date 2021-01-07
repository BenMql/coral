 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: LP_wallclock  
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> Contains wallclock-time related variables, and manipulations thereof.    
 !
 !=============================================================================
 module LP_wallclock
  use Fortran_kinds
  use MPI_vars
 Implicit None

 !> @brief
 !> type for keeping track of wall-clock time
  type :: wallclock_T
    Logical :: time_is_elapsed !< timer. DNS stops if .True.
    real(dp) :: limit          !< timer limit
    Real(dp) :: start_time     !< when did start?
    Real(dp) :: time_now       !< what time is it?
    Integer :: days, Hours, minutes, seconds, msecs !< integers for human-readable display
   contains
    procedure :: init  => get_the_starting_time
    procedure :: check => check_the_time
    procedure :: convert => convert_to_human_format
  end type wallclock_T

  type(wallclock_T) :: Wclock !< type for keeping track of wall-clock time

 Private
 public :: wclock

 Contains

 Subroutine get_the_starting_time(self)
   class(wallclock_T), intent(inOut) :: self
   if (my_rank.eq.0) self%start_time = MPI_Wtime()
 End Subroutine get_the_starting_time
   


 Subroutine Check_the_time(self, verb)
   class(wallclock_T), intent(inOut) :: self
   Logical, Intent(In)  :: verb
   if (my_rank.eq.0) then
      self%Time_now = MPI_Wtime()
      if (verb) print *, "Wall clock, time elapsed:", self%time_now-self%start_time
      if ((self%Time_now-self%start_time).lt.self%limit) Then
         self%time_is_elapsed = .False.
      Else
         self%time_is_elapsed = .True.
      End If
   End If
   Call MPI_Bcast(self%time_is_elapsed, 1, MPI_Logical, 0, MPI_Comm_world, ierr)
   call self%convert()
 End Subroutine Check_the_time



 subroutine convert_to_human_format(self)
   class(wallclock_T), intent(inOut) :: self
   Real(dp) :: tsec
   tsec = self%time_now - self%start_time
   self%seconds = Int(tsec)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   self%msecs = Int(1000._dp*(tsec-self%seconds))
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   self%days = self%seconds / (3600*24)
   self%seconds = self%seconds - 3600*24 * self%days
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   self%hours = self%seconds / (3600)
   self%seconds = self%seconds - 3600* self%hours
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   self%minutes = self%seconds / (60)
   self%seconds = self%seconds - 60* self%minutes
 end subroutine convert_to_human_format


 end module LP_wallclock
