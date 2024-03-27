Program check_equations

 use MPI_vars
 use ftext_parsing
 use PL_equations
 use, intrinsic :: iso_c_binding
 Implicit None

 type(full_problem_recipe_T) :: recipe


 call MPI_init(ierr)
 call MPI_Comm_Size(MPI_Comm_world, world_size, ierr)
 call MPI_Comm_Rank(MPI_Comm_world, My_rank,    ierr) 



    

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
 if (my_rank.eq.0) print *, '================================================================='
 if (my_rank.eq.0) print *, '================================================================='
 if (my_rank.eq.0) write(*,*) New_line('a')
 if (my_rank.eq.0) write(*,*) New_line('a')

!
   call read_equations_text_file()
   call recipe%build (arr_of_lines_eqn)
   !deAllocate( arr_of_lines_out, arr_of_lines_eqn, arr_of_lines_tms)

   call recipe%summarize()
      print *, ' * * * THE EQUATION FILE WAS SUCCESSFULLY READ AND INTERPRETED !! * * *'



  call MPI_finalize(ierr)


End Program check_equations




