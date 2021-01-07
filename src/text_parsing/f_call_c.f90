Program f_call_c
   use cwraps
   Implicit None
   Integer :: num
   Integer :: ierr
   integer, allocatable, dimension(:), target :: length
   character(len=:), allocatable , target :: arr_of_lines(:)
   integer :: iline, icar


   ierr = parse_text(num, 1, C_loc(length), C_loc(arr_of_lines))   
   print *, num   
   allocate(length(1:num))

   ierr = parse_text(num, 2, C_Loc(Length), C_loc(arr_of_lines))   
   Allocate(character(1024) :: arr_of_lines(num))
   
   Do iline = 1,num
   Do icar = 1,1024
        arr_of_lines(iline)(icar:icar) = '$'
   end Do
   end Do
   ierr = parse_text(num, 3, C_Loc(Length), C_loc(arr_of_lines))   
   print *, length(1:num)
   do iline=1,num
      write(*,*) arr_of_lines(iline)
      write(*,*) '\n'
   end do
End Program f_call_c

