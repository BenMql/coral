module PL_string_to_data
 use fortran_kinds
 implicit none

 type :: unsorted_T
    character(len=:), allocatable :: str
 end type unsorted_T

 type :: atom_linearOp_T
    character(len=:), allocatable :: str
    type(unsorted_T), dimension(5) :: unsorted
    integer :: nOf_unsortedStrings
    character(len=:), allocatable :: param_str
    character(len=:), allocatable :: var_str
    real(dp):: sgn           
    integer :: dx_exponent = 0 
    integer :: dy_exponent = 0 
    integer :: Iz_exponent = 0 
    integer :: z_multiply  = 0 
    real(dp) :: dsca
 end type atom_linearOp_T

 contains

 subroutine interpret_quadratic_variable_definition( myLine, qv, lv1, lv2)
   character(len=1024), intent(in) :: myLine
   character(len=:), allocatable, intent(out) :: qv, lv1, lv2
   character(len=:), allocatable :: val_str
   integer :: double_colon_index
   integer :: equal_sign_index
   integer :: dollar_sign_index
   integer :: dot_index
   double_colon_index = index(myLine,'::')+2
   equal_sign_index   = index(myLine,':=') -1
   dollar_sign_index  = index(myLine,'<<') -2
   if (Allocated(val_str)) deAllocate(val_str)
   Allocate(character(equal_sign_index-double_colon_index+1) :: val_str)
   val_str = myLine(double_colon_index:equal_sign_index)
   Call clean_str(val_str, qv)                   
   if (allocated(val_str)) DeAllocate(val_str)
   val_str = trim(adjustl(myLine(equal_sign_index+4:dollar_sign_index)))
   if (val_str(len(val_str):len(val_str))==char(0)) then
     val_str(len(val_str):len(val_str))=char(32)
   end if
   ! 
   dot_index = index(val_str,'.')  
   call clean_str_(val_str(1:dot_index-1), dot_index-1, lv1)                   
   call clean_str_(val_str(dot_index+1:len(val_str)), len(val_str) - dot_index, lv2)
    
 end subroutine interpret_quadratic_variable_definition


 subroutine interpret_full_variable_definition(myLine, atom_linear_ops, dzStr)
   character(len=1024), intent(in) :: myLine
   character(len=2), intent(In) :: dzStr
   character(len=:), allocatable :: my_str
   character(len=:), allocatable :: aux
   character(len=:), allocatable :: operator_chunk
   type(atom_linearOp_T), dimension(:), allocatable, intent(out) :: atom_linear_ops
   type(atom_linearOp_T), dimension(:), allocatable :: signs_and_strings
   integer :: next_dot
   integer :: unsorted_index
   integer :: iAtom
   integer :: n_atoms
  
   call break_linear_operators_in_pieces(myLine, signs_and_strings)
   n_atoms = size(signs_and_strings)
   allocate(atom_linear_ops(n_atoms))
   do iAtom = 1, n_atoms
      atom_linear_ops(iAtom)%sgn = signs_and_strings(iAtom)%sgn
      atom_linear_ops(iAtom)%dsca = signs_and_strings(iAtom)%sgn
      if (allocated(my_str)) deAllocate(my_str)
      allocate(character(len=len(signs_and_strings(iAtom)%str)) :: my_str )
      my_str = signs_and_strings(iAtom)%str
      unsorted_index = 1
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> Analyze each chunk of the linear operators
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 
        next_dot = index(my_str, '.')
        if (next_dot.gt.0) then
          if (allocated(aux)) DeAllocate(aux)
          allocate(character(next_dot-1) :: aux)
          aux = my_str(1:next_dot-1)
          Call clean_str(aux, operator_chunk)
          aux = my_str(next_dot+1:)
          deallocate(my_str)
          my_str=aux
          deallocate(aux)
        Else
          operator_chunk = adjustl(my_str)
        End if
        !! check whether this is a derivative
        if (operator_chunk=='dx') then
            atom_linear_ops(iAtom)%dx_exponent = atom_linear_ops(iAtom)%dx_exponent + 1
        else if (operator_chunk=='dy') then
            atom_linear_ops(iAtom)%dy_exponent = atom_linear_ops(iAtom)%dy_exponent + 1
        else if (operator_chunk==dzStr) then 
            atom_linear_ops(iAtom)%Iz_exponent = atom_linear_ops(iAtom)%Iz_exponent + 1
        else if (operator_chunk=='mz') then 
            atom_linear_ops(iAtom)%z_multiply  = atom_linear_ops(iAtom)%z_multiply  + 1
        else
            atom_linear_ops(iAtom)%unsorted(unsorted_index)%str = operator_chunk
            unsorted_index = unsorted_index + 1
        end if
      if (.not.(next_dot.gt.0)) exit
      if (unsorted_index.gt.5) then
         STOP 'too many unrecognized characters in full variable definition'                
      end if
      end do
      atom_linear_ops(iAtom)%nOf_unsortedStrings = unsorted_index-1
  

   end do

 end subroutine

 subroutine extract_one_float(myLine, myFloat)
   real(dp), intent(out) :: myFloat  
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=:), allocatable :: val_str
                                  !< raw value, read as a string               
   Integer :: double_colon_index, dollar_sign_index
                                  !< position of some delimiters

      double_colon_index = index(myLine,'::')+2
      dollar_sign_index  = index(myLine,'<<') -2
      if (Allocated(val_str)) deAllocate(val_str)
      Allocate(character(dollar_sign_index-double_colon_index+1) :: val_str)
      val_str = myLine(double_colon_index:dollar_sign_index)
      if (val_str(len(val_str):len(val_str))==char(0)) then
        val_str(len(val_str):len(val_str))=char(32)
      end if
      read(val_str,*) myFloat  
 end subroutine

 subroutine extract_one_integer(myLine, myInteger)
   integer, intent(out) :: myInteger
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=:), allocatable :: val_str
                                  !< raw value, read as a string               
   Integer :: double_colon_index, dollar_sign_index
                                  !< position of some delimiters

      double_colon_index = index(myLine,'::')+2
      dollar_sign_index  = index(myLine,'<<') -2
      if (Allocated(val_str)) deAllocate(val_str)
      Allocate(character(dollar_sign_index-double_colon_index+1) :: val_str)
      val_str = myLine(double_colon_index:dollar_sign_index)
      if (val_str(len(val_str):len(val_str))==char(0)) then
        val_str(len(val_str):len(val_str))=char(32)
      end if
      read(val_str,*) myInteger
 end subroutine
   
 subroutine get_BC_code(myLine, bc_code)
   integer, intent(out) :: bc_code
   character(len=1024), intent(In) :: myLine
   call extract_one_integer(myLine, bc_Code)
 end subroutine

 subroutine get_eqn_order(myLine, eq_order)
   integer, intent(out) :: eq_order
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   call extract_one_integer(myLine, eq_order)
 end subroutine

 
 subroutine get_parameter_name_and_value(myLine, param_values_str, param_values_val) 
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=:), allocatable, intent(out) :: param_values_str
                                  !< list of parameter names and values (double)
   real(dp), intent(out) :: param_values_val
   character(len=:), allocatable :: val_str
                                  !< raw value, read as a string               
   Integer :: double_colon_index, equal_sign_index, dollar_sign_index
                                  !< position of some delimiters

      double_colon_index = index(myLine,'::')+2
      equal_sign_index   = index(myLine,':=') -1
      dollar_sign_index  = index(myLine,'<<') -2
      if (Allocated(val_str)) deAllocate(val_str)
      Allocate(character(equal_sign_index-double_colon_index+1) :: val_str)
      val_str = myLine(double_colon_index:equal_sign_index)
      Call clean_str(val_str, param_values_str)
      if (allocated(val_str)) DeAllocate(val_str)
      val_str = trim(adjustl(myLine(equal_sign_index+4:dollar_sign_index)))
      if (val_str(len(val_str):len(val_str))==char(0)) then
        val_str(len(val_str):len(val_str))=char(32)
      end if
      read(val_str,*) param_values_val  
 end subroutine

 subroutine get_timeseriesVarName_period(myLine, variableName_str)
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=:), allocatable, intent(out) :: variableName_str
                                  !< list of parameter names and values (double)
   character(len=:), allocatable :: val_str
                                  !< raw value, read as a string               
   Integer :: double_colon_index, equal_sign_index
                                  !< position of some delimiters

      double_colon_index = index(myLine,'::')+2
      equal_sign_index   = index(myLine,'<<') -1
      if (Allocated(val_str)) deAllocate(val_str)
      Allocate(character(equal_sign_index-double_colon_index+1) :: val_str)
      val_str = myLine(double_colon_index:equal_sign_index)
      Call clean_str(val_str, variableName_str)
 end subroutine

 subroutine get_timeseriesVarName_period_position(myLine, variableName_str, positionInt)
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=:), allocatable, intent(out) :: variableName_str
                                  !< list of parameter names and values (double)
   integer, intent(out) :: positionInt
   character(len=:), allocatable :: val_str
                                  !< raw value, read as a string               
   Integer :: double_colon_index, dollar_sign_index, position_index
                                  !< position of some delimiters

      double_colon_index = index(myLine,'::')+2
      position_index   = index(myLine,';; Position=') -1
      dollar_sign_index  = index(myLine,'<<') -2
      if (Allocated(val_str)) deAllocate(val_str)
      Allocate(character(position_index-double_colon_index+1) :: val_str)
      val_str = myLine(double_colon_index:position_index)
      Call clean_str(val_str, variableName_str)
      if (allocated(val_str)) DeAllocate(val_str)
      val_str = trim(adjustl(myLine(position_index+14:dollar_sign_index)))
      if (val_str(len(val_str):len(val_str))==char(0)) then
        val_str(len(val_str):len(val_str))=char(32)
      end if
      read(val_str,*) positionInt  
 end subroutine

 subroutine get_outputVarName_period(myLine, variableName_str, periodInt) 
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=:), allocatable, intent(out) :: variableName_str
                                  !< list of parameter names and values (double)
   integer, intent(out) :: periodInt
   character(len=:), allocatable :: val_str
                                  !< raw value, read as a string               
   Integer :: double_colon_index, equal_sign_index, dollar_sign_index
                                  !< position of some delimiters

      double_colon_index = index(myLine,'::')+2
      equal_sign_index   = index(myLine,';; Period=') -1
      dollar_sign_index  = index(myLine,'<<') -2
      if (Allocated(val_str)) deAllocate(val_str)
      Allocate(character(equal_sign_index-double_colon_index+1) :: val_str)
      val_str = myLine(double_colon_index:equal_sign_index)
      Call clean_str(val_str, variableName_str)
      if (allocated(val_str)) DeAllocate(val_str)
      val_str = trim(adjustl(myLine(equal_sign_index+12:dollar_sign_index)))
      if (val_str(len(val_str):len(val_str))==char(0)) then
        val_str(len(val_str):len(val_str))=char(32)
      end if
      read(val_str,*) periodInt  
 end subroutine

 subroutine get_outputVarName_period_position(myLine, variableName_str, periodInt, positionInt)
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=:), allocatable, intent(out) :: variableName_str
                                  !< list of parameter names and values (double)
   integer, intent(out) :: periodInt
   integer, intent(out) :: positionInt
   character(len=:), allocatable :: val_str
                                  !< raw value, read as a string               
   Integer :: double_colon_index, equal_sign_index, dollar_sign_index, position_index
                                  !< position of some delimiters

      double_colon_index = index(myLine,'::')+2
      equal_sign_index   = index(myLine,';; Period=') -1
      position_index   = index(myLine,';; Position=') -1
      dollar_sign_index  = index(myLine,'<<') -2
      if (Allocated(val_str)) deAllocate(val_str)
      Allocate(character(equal_sign_index-double_colon_index+1) :: val_str)
      val_str = myLine(double_colon_index:equal_sign_index)
      Call clean_str(val_str, variableName_str)
      if (allocated(val_str)) DeAllocate(val_str)
      val_str = trim(adjustl(myLine(equal_sign_index+12:position_index)))
      if (val_str(len(val_str):len(val_str))==char(0)) then
        val_str(len(val_str):len(val_str))=char(32)
      end if
      read(val_str,*) periodInt  
      if (allocated(val_str)) DeAllocate(val_str)
      val_str = trim(adjustl(myLine(position_index+14:dollar_sign_index)))
      if (val_str(len(val_str):len(val_str))==char(0)) then
        val_str(len(val_str):len(val_str))=char(32)
      end if
      read(val_str,*) positionInt  
 end subroutine

 subroutine get_linear_variable_name(myLine, param_values_str)
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=:), allocatable, intent(out) :: param_values_str
                                  !< list of parameter names and values (double)
   character(len=:), allocatable :: val_str
                                  !< raw value, read as a string               
   Integer :: double_colon_index, eol_index
                                  !< position of some delimiters

      double_colon_index = index(myLine,'::')+2
      eol_index   = index(myLine,'<<') -1
      if (Allocated(val_str)) deAllocate(val_str)
      Allocate(character(eol_index-double_colon_index+1) :: val_str)
      val_str = myLine(double_colon_index:eol_index)
      Call clean_str(val_str, param_values_str)
 end subroutine
   
 subroutine get_linear_variable_zDecomposition(myLine, baroclinic_or_barotropic)
   character(len=1024), intent(In) :: myLine
                                  !< input line, extracted from the equations file
   character(len=10), intent(out) :: baroclinic_or_barotropic
                                  !< list of parameter names and values (double)
   Integer :: double_colon_index, eol_index
                                  !< position of some delimiters

   double_colon_index = index(myLine,'::')+2
   eol_index   = index(myLine,'<<') -2
   ! this routine is supposed to read either 'baroclinic' or 'barotropic',
   ! both of which have length 10. This is being checked now. 
   if (.not.(eol_index-double_colon_index == 10)) then
      ! return an error message:
      print *, '.',myLine(double_colon_index+1:eol_index),'.'
      print *, eol_index-double_colon_index
      error stop "get_linear_variable_zDecomposition read a string with a length different from 10"
   else 
      baroclinic_or_barotropic = myLine(double_colon_index+1:eol_index)
   end if
 end subroutine
   
 Subroutine clean_str(polluted_string, clean_string)
  Character(len=:), allocatable, intent(InOut) :: polluted_string
  Character(len=:), allocatable, intent(Out) :: clean_string   
  Character(len=:), allocatable :: aux
  integer :: int_char
  integer :: keep_int
  Integer :: ascii_code
  
  keep_int=0
  allocate(character(len(polluted_string)) :: aux)
  Do int_char = 1,len(polluted_string)
   ascii_code= iachar(polluted_string(int_char:int_char))
   if ((ascii_code.ne.32).and.(ascii_code.ne.0)) then
   keep_int=keep_int+1
   aux(keep_int:keep_int) = polluted_string(int_char:int_char)
   end if
  end do
  Allocate(character(keep_int) :: clean_string)
  clean_string(1:keep_int) = aux(1:keep_int)
 End Subroutine clean_str

 Subroutine clean_str_(polluted_string, ps_len, clean_string)
  integer, intent(in) :: ps_len
  Character(len=ps_len), intent(In) :: polluted_string
  Character(len=:), allocatable, intent(Out) :: clean_string   
  Character(len=:), allocatable :: aux
  integer :: int_char
  integer :: keep_int
  Integer :: ascii_code
  
  keep_int=0
  allocate(character(len(polluted_string)) :: aux)
  Do int_char = 1,len(polluted_string)
   ascii_code= iachar(polluted_string(int_char:int_char))
   if ((ascii_code.ne.32).and.(ascii_code.ne.0)) then
   keep_int=keep_int+1
   aux(keep_int:keep_int) = polluted_string(int_char:int_char)
   end if
  end do
  Allocate(character(keep_int) :: clean_string)
  clean_string(1:keep_int) = aux(1:keep_int)
 End Subroutine clean_str_

 Subroutine break_linear_operators_in_pieces(myline, pieces_of_linearOps)
   Type(atom_linearOp_T), allocatable, dimension(:), intent(Out) :: pieces_of_linearOps
   character(len=1024), intent(in) :: myLine                                           
   character(len=:), allocatable :: aux
                                  !< line being scanned, content
   Integer :: counter_pieces = 0  !< stores the number of parameters as we read them
   Integer :: new_start           !< search the +/- character starting from this position
   Integer :: end_expression

   counter_pieces = 0
   aux = myLine
   !> We count the number of + signs ...
   do while (index(aux   , '+').gt.0) 
      new_start = index(aux   , '+')+1
      counter_pieces = counter_pieces + 1
      aux    = aux   (new_start:)
   end do
   aux = myLine
   !> ... and the number of - signs.
   do while (index(aux   , '-').gt.0) 
      new_start = index(aux   , '-')+1
      counter_pieces = counter_pieces + 1
      aux    = aux   (new_start:)
   end do

   !============================
   !
   ! Allocate, and do a second scan to extract       
   ! the strings.                                      
   !
   !============================
   Allocate(pieces_of_linearOps(counter_pieces))
   aux = myLine
   counter_pieces=0
   do while (index(aux   , '+').gt.0) 
      new_start = index(aux   , '+')+1
      counter_pieces = counter_pieces + 1
      aux    = aux   (new_start:)
      end_expression = find_first_occurence_of(aux   ,'+-<',3)
      pieces_of_linearOps(counter_pieces)%str = trim(aux   (:end_expression-1))
      pieces_of_linearOps(counter_pieces)%sgn = +1.d0                           
   end do
   aux = myLine
   do while (index(aux   , '-').gt.0) 
      new_start = index(aux, '-')+1
      counter_pieces = counter_pieces + 1
      aux    = aux   (new_start:)
      end_expression = find_first_occurence_of(aux   ,'+-<',3)
      pieces_of_linearOps(counter_pieces)%str = trim(aux   (:end_expression-1))
      pieces_of_linearOps(counter_pieces)%sgn = -1.d0                           
   end do

 End Subroutine break_linear_operators_in_pieces

 Integer Pure Function find_first_occurence_of(a_string, some_substrings, len_substr)
   Character(len=:), allocatable, intent(In) :: a_string
   Integer, intent(In) :: len_substr
   Character(len=len_substr), intent(In) :: some_substrings
   integer :: ind, this_ind
   integer :: ic
   logical :: Present_bool
  
   ind = 120000!< a ridiculously large value

   Present_bool = .False.
   Do ic = 1,len(some_substrings)
      this_ind = index(a_string, some_substrings(ic:ic))
      if (this_ind.ne.0) Present_bool=.True.
      if (this_ind.ne.0) ind = min(ind,this_ind)
   End Do

   if (.not.Present_bool) Then
     ind = 0       
   end if
 
   find_first_occurence_of =  ind


 End Function
end module PL_string_to_data
