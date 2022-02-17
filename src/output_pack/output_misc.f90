 Module output_misc

 use MPI_vars
 use chdir_mod
 use fortran_kinds, Only: dp
 Implicit None
 Private :: dp

 Contains

 Subroutine Output_cumul_d_in_timeserie(path2file, path2file_len, file_str, file_len, dsca, first_or_not, my_position)
   Integer, Intent(in) :: path2file_len
   Character(len=path2file_len), intent(in) :: path2file
   Integer, Intent(In) :: file_len
   Character(len=file_len), Intent(In) :: file_str
   Real(Kind=dp), intent(InOut) :: dsca
   Logical, Intent(In) :: first_or_not
   Integer, Intent(InOut) :: my_position
   Character(len=100) :: my_file
   Real(kind=dp), Allocatable, Dimension(:) :: array_of_dsca
   Real(Kind=dp) :: sum_of_dsca

   Allocate(array_of_dsca(world_size))

   my_file = file_str//"_volAvg.dat"
   Call MPI_Gather(dsca, 1, MPI_DOUBLE, array_of_dsca, 1, MPI_DOUBLE, 0, MPI_Comm_World, ierr)

   sum_of_dsca = sum(array_of_dsca)
   dsca = sum_of_dsca
   
   If (My_Rank.eq.0) Then
     call chdir(trim(path2file))
     call chdir('./Timeseries')
     If (first_or_not) Then
     Open (Unit=9, File=trim(my_file), Status='replace', access='stream')
     Else 
     Open (Unit=9, File=trim(my_file), Status='old',     access='stream')
     End If
     Write(9, POS = my_position) sum_of_dsca
     Inquire(Unit=9, POS = my_position)
     Close(Unit=9)
   End If
   Call MPI_Bcast(my_position, 1, MPI_INTEGER, 0, MPI_Comm_world, ierr)
 End Subroutine Output_cumul_d_in_timeserie

 Subroutine Output_dsca_in_timeserie(path2file, path2file_len, file_str, file_len, core_num, dsca, first_or_not, my_position)
   Integer, Intent(in) :: path2file_len
   Character(len=path2file_len), intent(in) :: path2file
   Integer, Intent(In) :: file_len
   Character(len=file_len), Intent(In) :: file_str
   Integer, Intent(In) :: core_num
   Real(Kind=dp), intent(In) :: dsca
   Logical, Intent(In) :: first_or_not
   Integer, Intent(InOut) :: my_position
   Character(len=100) :: my_file
   Character(len=4) :: core_str

   write(core_str, "(i4.4)") core_num
   my_file = file_str//"_core"//core_str//".dat"
   
   call chdir(trim(path2file))
   call chdir('./Timeseries')
   If (first_or_not) Then
   Open (Unit=9, File=trim(my_file), Status='replace', access='stream')
   Else 
   Open (Unit=9, File=trim(my_file), Status='old',     access='stream')
   End If
   Write(9, POS = my_position) dsca
   Inquire(Unit=9, POS = my_position)
   Close(Unit=9)
 End Subroutine Output_dsca_in_timeserie

 Subroutine Output_dsca_in_unique_timeserie(path2file, path2file_len, file_str, file_len, dsca, first_or_not, my_position)
   Integer, Intent(in) :: path2file_len
   Character(len=path2file_len), intent(in) :: path2file
   Integer, Intent(In) :: file_len
   Character(len=file_len), Intent(In) :: file_str
   Real(Kind=dp), intent(In) :: dsca
   Logical, Intent(In) :: first_or_not
   Integer, Intent(InOut) :: my_position
   
   call chdir(trim(path2file))
   call chdir('./Timeseries')
   If (first_or_not) Then
   Open (Unit=9, File=trim(file_str), Status='replace', access='stream')
   Else 
   Open (Unit=9, File=trim(file_str), Status='old',     access='stream')
   End If
   Write(9, POS = my_position) dsca
   Inquire(Unit=9, POS = my_position)
   Close(Unit=9)
 End Subroutine Output_dsca_in_unique_timeserie

 End Module output_misc
