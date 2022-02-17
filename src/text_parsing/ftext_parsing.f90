Module ftext_parsing
 
 use cwraps
 implicit none
 
 
 character(len=:), allocatable, target :: arr_of_lines_eqn(:)
 character(len=:), allocatable, target :: arr_of_lines_out(:)
 character(len=:), allocatable, target :: arr_of_lines_tms(:)
 

 contains
 

 subroutine read_equations_text_file()
   integer :: num
   integer :: ierr
   integer :: iline, icar
   integer, allocatable, dimension(:), target :: length

   ! call the c routines on temporary variables we can point to
   ! (Target attribute forbidden in a type definition) 
   ! temporary alloc to please the compiler
   allocate(character(1) :: arr_of_lines_eqn(1))
   allocate( length  (1) )
   ierr = parse_text(num, 1, C_loc(length), C_loc(arr_of_lines_eqn(1) ))
   deAllocate ( length )
   deAllocate(arr_of_lines_eqn)
   Allocate(length(1:num))
   Allocate(character(1024) :: arr_of_lines_eqn(num))
   Do iline = 1,num
   Do icar = 1,1024
       arr_of_lines_eqn(iline)(icar:icar) = '$'
   end Do
   end Do
   ierr = parse_text(num, 2, C_Loc(Length), C_loc(arr_of_lines_eqn(1)))   
 end subroutine 

 subroutine read_usrOutput_text_file()
   integer :: num
   integer :: ierr
   integer :: iline, icar
   integer, allocatable, dimension(:), target :: length

   ! call the c routines on temporary variables we can point to
   ! (Target attribute forbidden in a type definition) 
   ! temporary alloc to please the compiler
   allocate(character(1) :: arr_of_lines_out(1))
   allocate( length  (1) )
   ierr = parse_output(num, 1, C_loc(length), C_loc(arr_of_lines_out(1) ))
   deAllocate ( length )
   deAllocate(arr_of_lines_out)
   Allocate(length(1:num))
   Allocate(character(1024) :: arr_of_lines_out(num))
   Do iline = 1,num
   Do icar = 1,1024
       arr_of_lines_out(iline)(icar:icar) = '$'
   end Do
   end Do
   ierr = parse_output(num, 2, C_Loc(Length), C_loc(arr_of_lines_out(1)))   
 end subroutine 

 subroutine read_timeseries_text_file()
   integer :: num
   integer :: ierr
   integer :: iline, icar
   integer, allocatable, dimension(:), target :: length

   ! call the c routines on temporary variables we can point to
   ! (Target attribute forbidden in a type definition) 
   ! temporary alloc to please the compiler
   allocate(character(1) :: arr_of_lines_tms(1))
   allocate( length  (1) )
   ierr = parse_timeseries(num, 1, C_loc(length), C_loc(arr_of_lines_tms(1) ))
   deAllocate ( length )
   deAllocate(arr_of_lines_tms)
   Allocate(length(1:num))
   Allocate(character(1024) :: arr_of_lines_tms(num))
   Do iline = 1,num
   Do icar = 1,1024
       arr_of_lines_tms(iline)(icar:icar) = '$'
   end Do
   end Do
   ierr = parse_timeseries(num, 2, C_Loc(Length), C_loc(arr_of_lines_tms(1)))   
 end subroutine 

End module ftext_parsing


