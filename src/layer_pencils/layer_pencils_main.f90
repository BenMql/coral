Program layer_pencils_main

 use include_git_version 
 use MPI_vars
 use LP_wallclock
 use LP_text_parsing
 use LP_equations
 use LP_algebra
 use LP_IMEX_timestepping
 use LP_timings
 use LP_transforms
 use decomp_2d_io
 use fftw3_wrap
 use read_command_line_args
 use, intrinsic :: iso_c_binding
 Implicit None

 type(full_problem_data_structure_T) :: main

 Logical :: verbose

 integer :: i1, i2, i3
 real(kind=dp) :: first_moment
 real(kind=dp) :: second_moment
 real(kind=dp) :: third_moment


 call MPI_init(ierr)
 call MPI_Comm_Size(MPI_Comm_world, world_size, ierr)
 call MPI_Comm_Rank(MPI_Comm_world, My_rank,    ierr) 

 call wClock%init()

 !! only the first process displays stuff (usually...)
 if (my_rank.eq.0) then
    verbose = .True.
 else 
    verbose = .False.
 end if

    

 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) print *, '================================================================='
 if (my_rank.eq.0) print *, '================================================================='
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) print *, '             ,o888888o.        ,o888888o.    '
 if (my_rank.eq.0) print *, '            8888     `88.   . 8888     `88.  '
 if (my_rank.eq.0) print *, '         ,8 8888       `8. ,8 8888       `8b '
 if (my_rank.eq.0) print *, '         88 8888           88 8888        `8b '
 if (my_rank.eq.0) print *, '         88 8888           88 8888         88 '
 if (my_rank.eq.0) print *, '         88 8888           88 8888         88 '
 if (my_rank.eq.0) print *, '         88 8888           88 8888        ,8P '
 if (my_rank.eq.0) print *, '         `8 8888       .8'' `8 8888       ,8P '
 if (my_rank.eq.0) print *, '            8888     ,88''   ` 8888     ,88'' '
 if (my_rank.eq.0) print *, '            `8888888P''         `8888888P''   '
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) print *, '    8 888888888o.            .8.          8 8888         '
 if (my_rank.eq.0) print *, '    8 8888    `88.          .888.         8 8888         '
 if (my_rank.eq.0) print *, '    8 8888     `88         :88888.        8 8888         '
 if (my_rank.eq.0) print *, '    8 8888     ,88        . `88888.       8 8888         '
 if (my_rank.eq.0) print *, '    8 8888.   ,88''       .8. `88888.      8 8888         '
 if (my_rank.eq.0) print *, '    8 888888888P''       .8`8. `88888.     8 8888         '
 if (my_rank.eq.0) print *, '    8 8888`8b          .8'' `8. `88888.    8 8888         '
 if (my_rank.eq.0) print *, '    8 8888 `8b.       .8''   `8. `88888.   8 8888         '
 if (my_rank.eq.0) print *, '    8 8888   `8b.    .888888888. `88888.  8 8888         '
 if (my_rank.eq.0) print *, '    8 8888     `88. .8''       `8. `88888. 8 888888888888 '
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) print *, ' Module...... Pencils'
 if (my_rank.eq.0) print *, ' Version..... : '
 if (my_rank.eq.0) call display_git_version()
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) print *, '================================================================='
 if (my_rank.eq.0) print *, '================================================================='
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) write(*,*) New_line('a')

!
   call read_equations_text_file()
   call main%recipe%build (arr_of_lines_eqn)
   call read_usrOutput_text_file()
   call main%recipe%add_outputs (arr_of_lines_out)
   call read_timeseries_text_file()
   call main%recipe%add_timeseries (arr_of_lines_tms)
   deAllocate( arr_of_lines_out, arr_of_lines_eqn, arr_of_lines_tms)

   if (verbose) call main%recipe%summarize()

   call main%geometry%init()
   call timings%init(misc_cargo%dt_max,&
                     misc_cargo%tolerance_factor,&
                     misc_cargo%CFL_C) 

   main%io_bookkeeping%output_directory = misc_cargo%output_directory
   main%io_bookkeeping%output_dir_length= misc_cargo%output_dir_length

   call dct_planner()             

   call main%init(misc_cargo%scheme_id)

   if (misc_cargo%qsave_restart) then
      if (verbose) print *, ">> loading a quicksave"
      call main%import_quickSave()
      call timings%fromFile(main%io_bookkeeping%output_directory//'/Restart/dt.sav',&
                            main%io_bookkeeping%output_dir_length+15)
   else
      call main%noisy_initial_conditions()
   end if

   call main%Factorize_operators (timings%dt, .True.)

 if (is_string_in_the_list('--checks',8)) then

   include "LP_checks.f90"

 else
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !.                                                                .
   !.            ~~~~  TIME-STEPPING STARTS HERE  ~~~~               .
   !.            ~~~~  TIME-STEPPING STARTS HERE  ~~~~               .
   !.            ~~~~  TIME-STEPPING STARTS HERE  ~~~~               .
   !.            ~~~~  TIME-STEPPING STARTS HERE  ~~~~               .
   !.                                                                .
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 wClock%limit = misc_cargo%time_limit
 call main%init_all_sources()

 do while (.not.wClock%Time_is_elapsed)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !.                                                                .
   !.  IF NEEDED, CHANGE dt AND RECOMPUTE THE LU FACTORIZATIONS.     .
   !.                                                                .
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call timings%check_for_update()
   main%io_bookkeeping%absolute_time = timings%absolute_time
   if (timings%refactorization_is_due) then 
      call main%factorize_operators(timings%dt)
      timings%refactorization_is_due = .False.
   End If

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !.                                                                .
   !.               PROPAGATE OVER ONE TIME-STEP                     .
   !.                                                                .
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call main%march_forward ( timings%dt_current )


   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!.                                                                .
   !!.                  DISPLAY SOME INFORMATIONS                     .
   !!.                                                                .
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
   if (verbose) Then
   4242 Format ('>>> Timestep:',I7,'; Wallclock:', I2,'T',I2,':',I2.2,':',I2.2,'.',I3.3,'; dt =',ES10.3,'; U(rms)=',ES10.3)
   write (*,4242) timings%i_timestep, Wclock%days, Wclock%hours, Wclock%minutes, Wclock%seconds, &
                               Wclock%msecs, timings%dt_current, &
                               sqrt(main%cargo%KE_display)
   End if

   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!.                                                                .
   !!.                      CFL CONSIDERATIONS                        .
   !!.                                                                .
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call timings%apply_cfl( main%cargo%cfl_based_DT ) 
   call timings%increment()

   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!.                                                                .
   !!.                       SOME OUTPUTS                             .
   !!.                                                                .
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (timings%checkpoint_counter .eq. timings%checkpoint_period) then
      call main%export_checkpoints()
      call timings%tofile(main%io_bookkeeping%output_directory//'/CheckPoints/dt.sav',&
                          main%io_bookkeeping%output_dir_length+19)
      timings%checkpoint_counter = 0
   end if
!
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!.                                                                .
   !!.                    CHECK THE WALL-CLOCK                        .
   !!.                                                                .
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Call wclock%check(.False.) !< boolean argument: verbosity 
                              !!< (.False. is quiet)               



 end do
 
 end if

   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!.                                                                .
   !!.        WAIT FOR ALL THE PROCESSES TO COMPLETE THEIR TASK       .
   !!.                                                                .
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Call MPI_Barrier(MPI_Comm_world, ierr)
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 If (verbose) Print*, "The time limit is reached."
 If (verbose) Print*, "Going to sleep now."
 If (verbose) Print*, "Good night all the ships at sea ..."

 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !                                                                             .
 !                  ~~~~  TERMINATE. END OF MAIN PROGRAM  ~~~~                 .
 !                  ~~~~  TERMINATE. END OF MAIN PROGRAM  ~~~~                 .
 !                  ~~~~  TERMINATE. END OF MAIN PROGRAM  ~~~~                 .
 !                  ~~~~  TERMINATE. END OF MAIN PROGRAM  ~~~~                 .
 !                                                                             .
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  call MPI_finalize(ierr)


End Program layer_pencils_main




