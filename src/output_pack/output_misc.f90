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
   Real(Kind=dp) :: sum_of_dsca

   my_file = file_str//"_volAvg.dat"
   call MPI_reduce(&
       dsca, &
       sum_of_dsca, &
       1, MPI_DOUBLE, MPI_sum, 0, MPI_Comm_World, ierr)
   dsca = sum_of_dsca
   
   rank0: if (My_Rank.eq.0) then
     call chdir(trim(path2file))
     call chdir('./Timeseries')
     if (first_or_not) then
     Open (Unit=9, File=trim(my_file), Status='replace', access='stream')
     Else 
     Open (Unit=9, File=trim(my_file), Status='old',     access='stream')
     End If
     Write(9, POS = my_position) sum_of_dsca
     Inquire(Unit=9, POS = my_position)
     Close(Unit=9)
   end if rank0
   call MPI_Bcast(my_position, 1, MPI_INTEGER, 0, MPI_Comm_world, ierr)
 end subroutine Output_cumul_d_in_timeserie

 Subroutine Output_bundled_cumul_d_in_timeserie(path2file, path2file_len, file_str, file_len, &
                   arr_dsca, count, first_or_not, my_position)
   integer, intent(in) :: path2file_len
   integer, intent(in) :: count
   character(len=path2file_len), intent(in) :: path2file
   integer, intent(In) :: file_len
   character(len=file_len), intent(In) :: file_str
   real(Kind=dp), intent(inOut), allocatable :: arr_dsca(:)
   logical, intent(in) :: first_or_not
   integer, intent(inOut) :: my_position
   character(len=100) :: my_file
   real(kind=dp), allocatable :: sum_of_dsca(:)

   allocate(sum_of_dsca(count))
   my_file = file_str//"_volAvg.dat"
   call MPI_reduce(&
       arr_dsca, &
       sum_of_dsca, &
       count, MPI_DOUBLE, MPI_sum, 0, MPI_Comm_World, ierr)

   arr_dsca = sum_of_dsca
   
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
 End Subroutine Output_bundled_cumul_d_in_timeserie
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

 subroutine Output_dsca_in_unique_timeserie(path2file, path2file_len, file_str, file_len, dsca, first_or_not, my_position)
   integer, intent(in) :: path2file_len
   character(len=path2file_len), intent(in) :: path2file
   integer, intent(in) :: file_len
   character(len=file_len), intent(in) :: file_str
   real(Kind=dp), intent(in) :: dsca
   logical, intent(in) :: first_or_not
   integer, intent(inOut) :: my_position
   
   call chdir(trim(path2file))
   call chdir('./Timeseries')
   if (first_or_not) then
   open (unit=9, file=trim(file_str), status='replace', access='stream')
   else 
   open (unit=9, file=trim(file_str), status='old',     access='stream')
   end if
   write(9, pos = my_position) dsca
   inquire(unit=9, pos = my_position)
   close(unit=9)
 end subroutine Output_dsca_in_unique_timeserie

 subroutine output_bundled_dsca_in_unique_timeserie(path2file, path2file_len, &
                    file_str, file_len, arr_dsca, first_or_not, my_position)
   integer, intent(in) :: path2file_len
   character(len=path2file_len), intent(in) :: path2file
   integer, intent(in) :: file_len
   character(len=file_len), intent(in) :: file_str
   real(Kind=dp), intent(in), allocatable :: arr_dsca(:)
   logical, intent(in) :: first_or_not
   integer, intent(inOut) :: my_position
   
   call chdir(trim(path2file))
   call chdir('./Timeseries')
   if (first_or_not) Then
   open (unit=9, file=trim(file_str), status='replace', access='stream')
   else 
   open (unit=9, file=trim(file_str), status='old',     access='stream')
   end if
   write(9, pos = my_position) arr_dsca
   inquire(unit=9, pos = my_position)
   close(unit=9)
 end subroutine output_bundled_dsca_in_unique_timeserie

 !!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!
 !!!!                                                          !!!!
 !!!!   ..... THE TWO SUBROUTINES BELOW ARE NEVER CALLED ..... !!!!
 !!!!                                                          !!!!
 !!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!

 subroutine Output_2Dbundled_dsca_in_unique_timeserie(path2file, path2file_len, file_str, file_len, &
         arr_dsca, index2, first_or_not, my_position)
   integer, intent(in) :: path2file_len
   integer, intent(in) :: index2
   character(len=path2file_len), intent(in) :: path2file
   integer, intent(in) :: file_len
   character(len=file_len), intent(in) :: file_str
   real(Kind=dp), intent(in), allocatable :: arr_dsca(:,:)
   logical, intent(in) :: first_or_not
   integer, intent(inOut) :: my_position
   
   call chdir(trim(path2file))
   call chdir('./Timeseries')
   if (first_or_not) Then
   open (unit=9, file=trim(file_str), status='replace', access='stream')
   else 
   open (unit=9, file=trim(file_str), status='old',     access='stream')
   end if
   write(9, pos = my_position) arr_dsca(:,index2)
   inquire(unit=9, pos = my_position)
   close(unit=9)
 end subroutine Output_2Dbundled_dsca_in_unique_timeserie

 Subroutine Output_2Dbundled_cumul_d_in_timeserie(path2file, path2file_len, file_str, file_len, &
                 arr_dsca, index2, count, first_or_not, my_position)
   integer, intent(in) :: path2file_len
   integer, intent(in) :: index2
   integer, intent(in) :: count
   character(len=path2file_len), intent(in) :: path2file
   integer, intent(In) :: file_len
   character(len=file_len), intent(In) :: file_str
   real(Kind=dp), intent(inOut), allocatable :: arr_dsca(:,:)
   logical, intent(in) :: first_or_not
   integer, intent(inOut) :: my_position
   character(len=100) :: my_file
   real(kind=dp), allocatable :: sum_of_dsca(:)

   allocate(sum_of_dsca(count))
   my_file = file_str//"_volAvg.dat"
   call MPI_reduce(&
       arr_dsca(1,index2), &
       sum_of_dsca, &
       count, MPI_DOUBLE, MPI_sum, 0, MPI_Comm_World, ierr)

   arr_dsca(:,index2) = sum_of_dsca(:)
   
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
 End Subroutine Output_2Dbundled_cumul_d_in_timeserie


 end module output_misc
