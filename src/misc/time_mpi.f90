 Module Time_MPI

 Use fortran_kinds
 Use MPI_vars

 Implicit None

 Contains

 Subroutine Check_the_time(T_start, T_limit, T_now, elapsed, verb)
   Real(dp), Intent(In) :: T_start
   Real(dp), Intent(In) :: T_Limit
   Logical, Intent(Out) :: elapsed
   Logical, Intent(In)  :: verb
   Real(dp), intent(out) :: T_now
  
   if (my_rank.eq.0) then 
      T_now = MPI_Wtime()
      if (verb) print *, "Wall clock, time elapsed:", T_now-T_start
      if ((T_now-T_start).lt.T_limit) Then
         elapsed = .False.
      Else 
         elapsed = .True.
      End If
   End If
   Call MPI_Bcast(elapsed, 1, MPI_Logical, 0, MPI_Comm_world, ierr)
 End Subroutine Check_the_time

 End Module Time_MPI
