Module read_command_line_args
Implicit none

 Contains


 Logical Function is_string_in_the_list(the_string, LenOfStr)
 Integer, Intent(In) :: LenOfStr
 Character(len=LenOfStr) :: the_string
 Integer :: iarg, my_length, num_of_commandLine_Args
 Character(len=:), Allocatable :: What_I_read

 Num_of_commandLine_args = Command_Argument_Count()

 is_string_in_the_list = .False.
 Do iarg =1 , Num_of_commandline_args
   Call Get_Command_Argument(iarg, Length=my_length)
   Allocate(Character(len=my_length) :: what_I_read)
   Call Get_Command_Argument(iarg, value = what_I_read)
   if (what_I_read.eq.the_string) is_string_in_the_list= .True.
   DeAllocate(What_I_read)
 End Do
  
 End Function is_String_in_the_list

End Module read_command_line_args
