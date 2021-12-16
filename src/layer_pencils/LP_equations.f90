!=============================================================================
!                            C O R A L
!=============================================================================
!
! MODULE: LP_equations
! 
!> @author
!> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
!
! DESCRIPTION
!> Gathers computer-usable information for building linear operators and 
!! non-linear terms from user-input instructions.
!
!=============================================================================

Module LP_equations
 use Fortran_kinds
 use MPI_vars, only: my_rank
 use LP_string_to_data
 implicit None

 type :: atom_NL_term_T   
    integer :: dx_exponent
    integer :: dy_exponent
    integer :: Iz_exponent
    integer :: z_multiply 
    integer :: quad_var_index
    integer :: eqn_index
    real(dp) :: dsca
 end type atom_NL_term_T   

 type :: nonlinear_rhs_recipe_T
    integer :: n_pieces
    type(atom_NL_term_T), dimension(:), allocatable :: term
 end type nonlinear_rhs_recipe_T

 type :: source_term_T   
    integer :: dx_exponent
    integer :: dy_exponent
    integer :: Iz_exponent
    integer :: source_index
    integer :: eqn_index
    integer :: sys_index
    real(dp) :: dsca
 end type source_term_T   

 type :: sourceParams_T
    character(len=6) :: str !< user-given name 
    real(dp) :: dsca        !< value
 end type sourceParams_T

 type :: source_rhs_recipe_T
    integer :: n_sources
    type(sourceParams_T), allocatable :: sourceParams(:)
    type(source_term_T), dimension(:), allocatable :: term
 end type source_rhs_recipe_T

 type :: linear_operator_recipe_T   
    real(dp), allocatable :: dsca(:)
    integer, allocatable :: kx_exponent(:)
    integer, allocatable :: ky_exponent (:)
    integer, allocatable :: Iz_exponent (:)
    integer, allocatable :: z_multiply  (:)
    integer, allocatable :: jVar (:)
    integer, allocatable :: iEqn (:)
    integer :: n_pieces
 end type linear_operator_recipe_T   

 
 type :: a_linear_contribution_T
   integer :: dx_exponent
   integer :: dy_exponent
   integer :: dz_exponent
   integer :: z_multiply 
   real(dp) :: dsca
   character(len=4) :: var_kind ! 'kxky' or 'zero'
   integer :: var_index
 end type  a_linear_contribution_T

 type :: list_of_linVars_T
    character(len=:), allocatable :: str
    integer :: belongs_to_set
    integer :: at_position    
 end type list_of_linVars_T

 type :: list_of_parameters_T
    character(len=:), allocatable :: str
    real(C_Double), allocatable :: val
 end type list_of_parameters_T

 type :: list_of_nlVars_T
    character(len=:), allocatable :: str
    integer :: ivar1
    integer :: ivar2
    character(len=:), allocatable :: svar1
    character(len=:), allocatable :: svar2
 end type list_of_nlVars_T


 type :: fullVars_recipe_T
   !< contains the recipe for building fields that will be computed
   !! in physical space. 
   !! For instance, to compute advection in flows, we need the horizontal velocity
   !! that can be expressed in terms of toroidal/poloidal potentials:
   !!    U = [d/dy] psi + [d/dx][d/dz] phi + u0
   !! In the example above, self%N_terms = 3 and each of these 3 terms is stored in
   !! self % term (1:3)
   character(len=:), allocatable :: str
   integer :: N_terms
   type(a_linear_contribution_T), dimension(:), allocatable :: term
 end type fullVars_recipe_T

 type :: coupled_system_vars_info_T
    character(len=:), allocatable :: name
    integer :: bc_code
 end type


 type :: coupled_system_recipe_T
    type(linear_operator_recipe_T) :: mass
    type(linear_operator_recipe_T) :: stif
    type(nonlinear_rhs_recipe_T) :: NL
    integer :: n_coupled_vars
    type(coupled_system_vars_info_T), dimension(:), allocatable :: vars
    integer :: n_equations
    integer, dimension(:), allocatable :: eqn_order
 end type coupled_system_recipe_T

 type :: one_output_object_T
    integer :: counter
    integer :: period
    character(len=6) :: kind
    integer :: slice_index
 end type one_output_object_T

 type :: list_of_output_objects_T
    type(one_output_object_T), allocatable :: object(:)
 end type list_of_output_objects_T

 type :: output_lists_T
    integer, allocatable :: numberOf_linearObjects(:)
    integer, allocatable :: numberOf_quadraObjects(:)
    type(list_of_output_objects_T), allocatable :: linear(:)
    type(list_of_output_objects_T), allocatable :: quadra(:)
 end type output_lists_T

 type :: full_problem_recipe_T 
   integer :: numberOf_coupled_kxky_systems
   integer :: numberOf_coupled_zero_systems
   integer :: numberOf_parameters                  
   integer :: numberOf_linear_variables_k
   integer :: numberOf_linear_variables_0
   integer :: numberOf_quadratic_variables         
   integer :: numberOf_linear_variables_full       
   type(coupled_system_recipe_T), dimension(:), allocatable :: kxky_recipes
   type(coupled_system_recipe_T), dimension(:), allocatable :: zero_recipes
   type(list_of_parameters_T), dimension(:), allocatable :: list_parameters
   type(list_of_linVars_T ), dimension(:), allocatable :: linear_vars_kxky
   type(list_of_linVars_T ), dimension(:), allocatable :: linear_vars_zero
   type(fullVars_recipe_T ), dimension(:), allocatable :: linear_vars_full
   type(list_of_NLVars_T  ), dimension(:), allocatable :: nl_vars
   type(output_lists_T) :: output
   type(output_lists_T) :: timeseries
   type(source_rhs_recipe_T) :: sources
   logical :: smagorinsky_flag
  contains
   procedure :: add_outputs
   procedure :: add_timeseries
   procedure :: copy_params => copy_list_of_parameters_from_recipe
   !procedure :: copy_linear_vars_k => copy_linear_vars_k_from_recipe
   !procedure :: copy_linear_vars_0 => copy_linear_vars_0_from_recipe
   procedure :: build => build_full_recipe_from_text_file
   procedure :: add_parameter => add_parameter_to_full_recipe
   procedure :: add_source_parameter => add_source_parameter_to_full_recipe
   procedure :: add_linear_var => add_linear_variable_to_full_recipe
   procedure :: add_quadratic_var => add_quadratic_variable_to_full_recipe
   procedure :: add_full_var => add_full_linear_variable_to_full_recipe
   procedure :: add_coupled_kxky_eqns => add_coupled_kxky_equations_to_full_recipe
   procedure :: add_coupled_zero_eqns => add_coupled_zero_equations_to_full_recipe
   procedure :: build_current_full_var => build_current_full_var_in_full_recipe
   procedure :: init => initialize_full_recipe
   procedure :: add_linear_var_to_kxky_system => add_lin_var_to_kxky_eqns_in_full_recipe
   procedure :: add_linear_var_to_zero_system => add_lin_var_to_zero_eqns_in_full_recipe
   procedure :: specify_BC_k => specify_boundaryConditions_in_kxky_eqns
   procedure :: specify_BC_0 => specify_boundaryConditions_in_zero_eqns
   procedure :: add_kxky_eqn => add_equation_to_kxky_coupled_set
   procedure :: add_zero_eqn => add_equation_to_zero_coupled_set
   procedure :: build_kxky_mass_matrix => build_kxky_mass_matrix_in_full_recipe
   procedure :: build_kxky_stif_matrix => build_kxky_stif_matrix_in_full_recipe
   procedure :: build_zero_mass_matrix => build_zero_mass_matrix_in_full_recipe
   procedure :: build_zero_stif_matrix => build_zero_stif_matrix_in_full_recipe
   procedure :: build_zero_NL_term     => build_zero_NL_term_in_full_recipe
   procedure :: build_zero_sources     => build_zero_sources_in_full_recipe
   procedure :: build_kxky_NL_term     => build_kxky_NL_term_in_full_recipe
   procedure :: summarize              => summarize_full_recipe
 end type full_problem_recipe_T

 type(full_problem_recipe_T) :: master_recipe

 contains 

 Subroutine add_timeseries( self, text_list)
   class(full_problem_recipe_T) :: self
   character(len=:), allocatable, intent(in) :: text_list(:)
   character(len=1024) :: restOfmyLine
   integer :: iLine, iVar
   integer :: numOfLines
   character(len=:), allocatable :: varNameStr
   integer :: positionInt
   character(len=6) :: kindStr
   type(one_output_object_T), allocatable :: larger_list_of_objects(:)
   iLine = 0
   do 
     iLine = iLine + 1
     if (text_list(iLine)(3:5).eq.'EOF') exit
   end do
   numOfLines = iLine -1
   allocate (self%timeseries%linear( self%numberOf_linear_variables_full ))
   allocate (self%timeseries%quadra( self%numberOf_quadratic_variables ))
   allocate (self%timeseries%numberOf_linearObjects( self%numberOf_linear_variables_full ))
   allocate (self%timeseries%numberOf_quadraObjects( self%numberOf_quadratic_variables ))
   self%timeseries%numberOf_linearObjects = 0
   self%timeseries%numberOf_quadraObjects = 0
   do iVar = 1,self%numberOf_linear_variables_full
      allocate (self%timeseries%linear( iVar )% object (0) )                           
   end do
   do iVar = 1,self%numberOf_quadratic_variables
      allocate (self%timeseries%quadra( iVar )% object (0) )                           
   end do

   do iLine = 1, numOfLines
      restOfmyLine = text_list(iLine)(1:1024)
      select case (text_list(iLine)(1:19))
             case ('>>Volume_average ::')
                  call get_timeseriesVarName_period( restOfMyLine, varNameStr)
                  kindStr = 'volume'
             case ('>>Horizontal_avg ::')
                  call get_timeseriesVarName_period_position( restOfMyLine, varNameStr, positionInt)
                  kindStr = 'zSlice'
             case default
                  print *, 'Error in coral.timeseries'
                  print *, 'First nineteen characters of each line should be:'
                  print *, '>>Volume_average ::, >>Horizontal_avg ::'
                  error stop
      end select
      !scan all the variables
      do iVar = 1,self%numberOf_linear_variables_full
      if (varNameStr.eq.self%linear_vars_full(iVar)%str) then
         self%timeseries%numberOf_linearObjects(iVar) = self%timeseries%numberOf_linearObjects(iVar) + 1
         allocate (larger_list_of_objects ( self%timeseries%numberOf_linearObjects(iVar) ))
         larger_list_of_objects(1:self%timeseries%numberOf_linearObjects(iVar)-1) = &
                 self%timeseries%linear(iVar)%object(1:self%timeseries%numberOf_linearObjects(iVar)-1)
         larger_list_of_objects(  self%timeseries%numberOf_linearObjects(iVar)) % counter = 0
         larger_list_of_objects(  self%timeseries%numberOf_linearObjects(iVar)) % period  = 0 
         larger_list_of_objects(  self%timeseries%numberOf_linearObjects(iVar)) % slice_index = positionInt
         larger_list_of_objects(  self%timeseries%numberOf_linearObjects(iVar)) % kind = kindStr
         call move_alloc( from = larger_list_of_objects, to = self%timeseries%linear(iVar)%object)
      end if
      end do
      do iVar = 1,self%numberOf_quadratic_variables  
      if (varNameStr.eq.self%nl_vars(iVar)%str) then
         self%timeseries%numberOf_quadraObjects(iVar) = self%timeseries%numberOf_quadraObjects(iVar) + 1
         allocate (larger_list_of_objects ( self%timeseries%numberOf_quadraObjects(iVar) ))
         larger_list_of_objects(1:self%timeseries%numberOf_quadraObjects(iVar)-1) = &
                 self%timeseries%quadra(iVar)%object(1:self%timeseries%numberOf_quadraObjects(iVar)-1)
         larger_list_of_objects(  self%timeseries%numberOf_quadraObjects(iVar)) % counter = 0
         larger_list_of_objects(  self%timeseries%numberOf_quadraObjects(iVar)) % period  = 0
         larger_list_of_objects(  self%timeseries%numberOf_quadraObjects(iVar)) % slice_index = positionInt
         larger_list_of_objects(  self%timeseries%numberOf_quadraObjects(iVar)) % kind = kindStr
         call move_alloc( from = larger_list_of_objects, to = self%timeseries%quadra(iVar)%object)
      end if
      end do
   end do
 end subroutine
 

 Subroutine add_outputs( self, text_list)
   class(full_problem_recipe_T) :: self
   character(len=:), allocatable, intent(in) :: text_list(:)
   character(len=1024) :: restOfmyLine
   integer :: iLine
   integer :: numOfLines
   character(len=:), allocatable :: varNameStr
   integer :: periodInt, positionInt, iVar
   character(len=6) :: kindStr
   type(one_output_object_T), allocatable :: larger_list_of_objects(:)
   iLine = 0
   do 
     iLine = iLine + 1
     if (text_list(iLine)(3:5).eq.'EOF') exit
   end do
   numOfLines = iLine -1
   
   allocate (self%output%linear( self%numberOf_linear_variables_full ))
   allocate (self%output%quadra( self%numberOf_quadratic_variables ))
   allocate (self%output%numberOf_linearObjects( self%numberOf_linear_variables_full ))
   allocate (self%output%numberOf_quadraObjects( self%numberOf_quadratic_variables ))
   self%output%numberOf_linearObjects = 0
   self%output%numberOf_quadraObjects = 0
   do iVar = 1,self%numberOf_linear_variables_full
      allocate (self%output%linear( iVar )% object (0) )                           
   end do
   do iVar = 1,self%numberOf_quadratic_variables
      allocate (self%output%quadra( iVar )% object (0) )                           
   end do

   do iLine = 1, numOfLines
      restOfmyLine = text_list(iLine)(1:1024)
      select case (text_list(iLine)(1:11))
             case ('>>Volume ::')
                  call get_outputVarName_period( restOfMyLine, varNameStr, periodint)
                  kindStr = 'volume'
             case ('>>XSlice ::')
                  call get_outputVarName_period_position( restOfMyLine, varNameStr, periodint, positionInt)
                  kindStr = 'xSlice'
             case ('>>YSlice ::')
                  call get_outputVarName_period_position( restOfMyLine, varNameStr, periodint, positionInt)
                  kindStr = 'ySlice'
             case ('>>ZSlice ::')
                  call get_outputVarName_period_position( restOfMyLine, varNameStr, periodint, positionInt)
                  kindStr = 'zSlice'
             case ('>>Zsummed::')
                  call get_outputVarName_period( restOfMyLine, varNameStr, periodint)
                  kindStr = 'zAvged'
             case ('>>Profile::')
                  call get_outputVarName_period( restOfMyLine, varNameStr, periodInt)
                  kindStr = 'profil'
             case default
                  print *, 'Error in coral.usrOuput'
                  print *, 'First eleven characters of each line should be:'
                  print *, '>>Volume ::, >>XSlice ::, >>YSlice ::, ZSlice ::,'
                  print *, '>>Profile::, >>Zsummed::'
                  print *, '(Common mistake: no space between Profile and ::'               
                  error stop
      end select
      !scan all the variables
      do iVar = 1,self%numberOf_linear_variables_full
      if (varNameStr.eq.self%linear_vars_full(iVar)%str) then
         self%output%numberOf_linearObjects(iVar) = self%output%numberOf_linearObjects(iVar) + 1
         allocate (larger_list_of_objects ( self%output%numberOf_linearObjects(iVar) ))
         larger_list_of_objects(1:self%output%numberOf_linearObjects(iVar)-1) = &
                 self%output%linear(iVar)%object(1:self%output%numberOf_linearObjects(iVar)-1)
         larger_list_of_objects(  self%output%numberOf_linearObjects(iVar)) % counter = 0
         larger_list_of_objects(  self%output%numberOf_linearObjects(iVar)) % period = periodInt
         larger_list_of_objects(  self%output%numberOf_linearObjects(iVar)) % slice_index = positionInt
         larger_list_of_objects(  self%output%numberOf_linearObjects(iVar)) % kind = kindStr
         call move_alloc( from = larger_list_of_objects, to = self%output%linear(iVar)%object)
      end if
      end do
      do iVar = 1,self%numberOf_quadratic_variables  
      if (varNameStr.eq.self%nl_vars(iVar)%str) then
         self%output%numberOf_quadraObjects(iVar) = self%output%numberOf_quadraObjects(iVar) + 1
         allocate (larger_list_of_objects ( self%output%numberOf_quadraObjects(iVar) ))
         larger_list_of_objects(1:self%output%numberOf_quadraObjects(iVar)-1) = &
                 self%output%quadra(iVar)%object(1:self%output%numberOf_quadraObjects(iVar)-1)
         larger_list_of_objects(  self%output%numberOf_quadraObjects(iVar)) % counter = 0
         larger_list_of_objects(  self%output%numberOf_quadraObjects(iVar)) % period = periodInt
         larger_list_of_objects(  self%output%numberOf_quadraObjects(iVar)) % slice_index = positionInt
         larger_list_of_objects(  self%output%numberOf_quadraObjects(iVar)) % kind = kindStr
         call move_alloc( from = larger_list_of_objects, to = self%output%quadra(iVar)%object)
      end if
      end do
   end do
 end subroutine

 Subroutine build_full_recipe_from_text_file( self, text_list)
   class(full_problem_recipe_T) :: self
   character(len=:), allocatable, intent(in) :: text_list(:)
   character(len=1024) :: myLine
   integer :: iLine !< scanning the lines of text_list, filling-up self.
   integer :: numOfLines                                                  
   integer :: building_step
   logical :: eof_signal, stop_signal
   character(len=:), allocatable :: param_name
   character(len=:), allocatable :: linear_var_name
   real(dp) :: dsca
   type(atom_linearOp_T), dimension(:), allocatable :: pieces_of_definition
   character(len=:), allocatable :: QVar_name, LVar1_name, LVar2_name, nuSmagoStr
   integer :: bc_code
   integer :: eqn_order
   logical :: have_not_considered_smagorinsky_yet = .True.
   building_step = 1
   call self%init()
   iLine = 0
   do 
     iLine = iLine + 1
     if (text_list(iLine)(1:3).eq.'>>=') then
        cycle
     end if
     myLine = text_list(iLine)(1:1024)
     call determine_kind_of_information(myLine, building_step, eof_signal, stop_signal)
     if (stop_signal) write (*,*) 'Error in coral.equations at line', iLine
     if (stop_signal) write (*,*) 'A different prefix was expected. Stopping now.'
     if (stop_signal) STOP
     if (eof_signal) exit
   end do
   numOfLines = iLine-1

   building_step = 1
   do iLine = 1, numOfLines
     if (my_rank.eq.0) print *, text_list(iLine)(1:78)
     ! if the line starts with '>>=', ignore and proceed to next line
     if (text_list(iLine)(1:3).eq.'>>=') then
        cycle
     end if
     myLine = text_list(iLine)(1:1024)
     call determine_kind_of_information(myLine, building_step, eof_signal, stop_signal)
     ! What kind of information is provided by the current line?
     Select case (building_step)
       case (1)
           call get_parameter_name_and_value(text_list(iLine)(1:1024), param_name, dsca)
           call self%add_parameter( param_name, dsca)
       case (11)
           call get_parameter_name_and_value(text_list(iLine)(1:1024), param_name, dsca)
           call self%add_source_parameter( param_name, dsca)
       case (2)
           call get_linear_variable_name(text_list(iLine)(1:1024), linear_var_name)
           call self%add_linear_var( linear_var_name, 'kxky')
       case (3)
           call get_linear_variable_name(text_list(iLine)(1:1024), linear_var_name)
           call self%add_linear_var( linear_var_name, 'zero')
       case (4)
           call get_linear_variable_name(text_list(iLine)(1:1024), linear_var_name)
           call self%add_full_var( linear_var_name)
       case (5) 
           call interpret_full_variable_definition(text_list(iLine)(1:1024), &
                                                   pieces_of_definition, 'dz')
           call self%build_current_full_var( pieces_of_definition)
       case (6)
           if (have_not_considered_smagorinsky_yet) then
              if (self%smagorinsky_flag) then
                  !allocate( character(len=7) :: nuSmagoStr, source='nuSmago')
                  allocate( nuSmagoStr, source='nuSmago')
                  call self%add_full_var( nuSmagoStr )     
              end if
              have_not_considered_smagorinsky_yet = .False.
           end if
           call interpret_quadratic_variable_definition(text_list(iLine)(1:1024), Qvar_name, &
                                                                        Lvar1_name, Lvar2_name)
           call self%add_quadratic_var(Qvar_name, Lvar1_name, Lvar2_name)
       case(7)
           call self%add_coupled_kxky_eqns()
       case(71)
           call get_linear_variable_name(text_list(iLine)(1:1024), linear_var_name)
           call self%add_linear_var_to_kxky_system (linear_var_name)
       case(72)
           call get_BC_code(text_list(iLine)(1:1024), bc_code)
           call self%specify_BC_k(bc_code)
       case(73)
           call get_eqn_order(text_list(iLine)(1:1024), eqn_order)
           call self%add_kxky_eqn(eqn_order)
       case(731)
           call interpret_full_variable_definition(text_list(iLine)(1:1024),&
                                                   pieces_of_definition, 'Iz')
           call self%build_kxky_mass_matrix( pieces_of_definition)
       case(732)
           call interpret_full_variable_definition(text_list(iLine)(1:1024),&
                                                   pieces_of_definition, 'Iz')
           call self%build_kxky_stif_matrix( pieces_of_definition)
       case(733)
           call interpret_full_variable_definition(text_list(iLine)(1:1024),&
                                                   pieces_of_definition, 'Iz')
           call self%build_kxky_NL_term( pieces_of_definition)
       case(8)
           call self%add_coupled_zero_eqns()
       case(81)
           call get_linear_variable_name(text_list(iLine)(1:1024), linear_var_name)
           call self%add_linear_var_to_zero_system (linear_var_name)
       case(82)
           call get_BC_code(text_list(iLine)(1:1024), bc_code)
           call self%specify_BC_0(bc_code)
       case(83)
           call get_eqn_order(text_list(iLine)(1:1024), eqn_order)
           call self%add_zero_eqn(eqn_order)
       case(831)
           call interpret_full_variable_definition(text_list(iLine)(1:1024),&
                                                   pieces_of_definition, 'Iz')
           call self%build_zero_mass_matrix( pieces_of_definition)
       case(832)
           call interpret_full_variable_definition(text_list(iLine)(1:1024),&
                                                   pieces_of_definition, 'Iz')
           call self%build_zero_stif_matrix( pieces_of_definition)
       case(833)
           call interpret_full_variable_definition(text_list(iLine)(1:1024),&
                                                   pieces_of_definition, 'Iz')
           call self%build_zero_NL_term( pieces_of_definition)
       case(834)
           call interpret_full_variable_definition(text_list(iLine)(1:1024),&
                                                   pieces_of_definition, 'Iz')
           call self%build_zero_sources( pieces_of_definition)
       case default
           if (my_rank.eq.0) print *, text_list(iLine)(1:78)
           stop 'the line above could not be read'
     end select

   end do 


   
 end subroutine build_full_recipe_from_text_file


 subroutine summarize_full_recipe(self)
   class(full_problem_recipe_T) :: self
   integer :: iline, iSys, iPiece, iEq, iPieceStop
   print *, '================================================================='
   print *, ' summary of the parameters and equations, as read from'
   print *, ' the text input:'
   print *, '.................................................................'
   do iline =1, self%numberOf_parameters
     print *, '::', self%list_parameters(iline)%str,':=',self%list_parameters(iline)%val
   end do
   print *, ':: linear variables -- fluctuations'
   do iline = 1, self%numberOf_linear_variables_k
     print *, self%linear_vars_kxky(iline)%str,','
   end do
   print *, ':: linear variables -- mean fields'  
   do iline = 1, self%numberOf_linear_variables_0
     print *, self%linear_vars_zero(iline)%str,','
   end do
   print *, ':: linear variables -- full fields for NL computations'
   do iline = 1, self%numberOf_linear_variables_full
     print *, self%linear_vars_full(iline)%str,','
   end do
   print *, ':: Finite (kx, ky) modes:'
   print *, '~~~~~~~~~~~~~~~~~~~~~~~~~'
   746 format ('System #', I2,' couples ',I2,' n_coupled_vars')
   do isys = 1, self%numberOf_coupled_kxky_systems
     write (*,746) iSys, self%kxky_recipes(iSys)%n_coupled_vars
     do iEq = 1, self%kxky_recipes(iSys)%n_equations
        748 format ('Equation #', I2,' is of order ',I2)
        750 format ('In equation #', I2,' there is a [mass] contribution of order',I2)
        752 format ('In equation #', I2,' there is a [stif] contribution of order',I2)
        751 format ('mass operator #',I2,' equation ',I2,', variable ',I2,', dx:',I2,', dy:',I2,', Iz:',I2)
        write (*,748) iEq, self%kxky_recipes(iSys)%eqn_order(iEq)
        !check that all operators contains no primitive of higher order
        !start with mass operator
        do iPiece = 1,self%kxky_recipes(iSys)%mass%n_pieces
         if (self%kxky_recipes(iSys)%mass%iEqn(iPiece).eq.iEq) then
           if (self%kxky_recipes(iSys)%mass%Iz_exponent(iPiece).gt. self%kxky_recipes(iSys)%eqn_order(iEq)) then
           write (*,750) iEq, self%kxky_recipes(iSys)%mass%Iz_exponent(iPiece)
           do iPieceStop = 1,self%kxky_recipes(iSys)%mass%n_pieces
           write (*,751) iPieceStop, self%kxky_recipes(iSys)%mass%iEqn(iPieceStop),&
                                 self%kxky_recipes(iSys)%mass%jVar(iPieceStop),&
                                 self%kxky_recipes(iSys)%mass%kx_exponent(iPieceStop),&
                                 self%kxky_recipes(iSys)%mass%ky_exponent(iPieceStop),&
                                 self%kxky_recipes(iSys)%mass%Iz_exponent(iPieceStop)
           end do
           error stop
           end if
         end if
        end do
        !now the stiffness operator
        do iPiece = 1,self%kxky_recipes(iSys)%stif%n_pieces
         if (self%kxky_recipes(iSys)%stif%iEqn(iPiece).eq.iEq) then
           if (self%kxky_recipes(iSys)%stif%Iz_exponent(iPiece).gt. self%kxky_recipes(iSys)%eqn_order(iEq)) then
           write (*,752) iEq, self%kxky_recipes(iSys)%stif%Iz_exponent(iPiece)
           error stop
           end if
         end if
        end do
     end do
   end do
   print *, ':: zero modes:'
   print *, '~~~~~~~~~~~~~~~~~~~~~~~~~'
   do isys = 1, self%numberOf_coupled_zero_systems
     write (*,746) iSys, self%zero_recipes(iSys)%n_coupled_vars
   end do
   print *, '================================================================='
   print *, '================================================================='
 end subroutine

 subroutine initialize_full_recipe(self)
   class(full_problem_recipe_T) :: self
   self%numberOf_coupled_kxky_systems = 0 
   self%numberOf_coupled_zero_systems = 0 
   self%numberOf_parameters           = 0 
   self%numberOf_linear_variables_k   = 0 
   self%numberOf_linear_variables_0   = 0 
   self%numberOf_linear_variables_full= 0 
   self%numberOf_quadratic_variables  = 0 
   allocate (self%kxky_recipes( 0 ) )
   allocate (self%zero_recipes( 0 ) )
   allocate (self%list_parameters  ( 0 ) )
   allocate (self%linear_vars_kxky ( 0 ) )
   allocate (self%linear_vars_zero ( 0 ) )
   allocate (self%linear_vars_full ( 0 ) )
   allocate (self%nl_vars     ( 0 ) )
   self%sources%n_sources = 0
   allocate(self%sources%term(0))
   allocate(self%sources%sourceParams(0))                           
 End Subroutine initialize_full_recipe

 subroutine specify_boundaryConditions_in_kxky_eqns(self, bc_code)
    class(full_problem_recipe_T) :: self
    integer :: this_recipe, this_var
    integer, intent(in) :: bc_code
    this_recipe = self%numberOf_coupled_kxky_systems 
    this_var = self%kxky_recipes(this_recipe)%n_coupled_vars
    self%kxky_recipes(this_recipe)%vars(this_var)%bc_code = bc_code
 end subroutine

 subroutine specify_boundaryConditions_in_zero_eqns(self, bc_code)
    class(full_problem_recipe_T) :: self
    integer :: this_recipe, this_var
    integer, intent(in) :: bc_code
    this_recipe = self%numberOf_coupled_zero_systems 
    this_var = self%zero_recipes(this_recipe)%n_coupled_vars
    self%zero_recipes(this_recipe)%vars(this_var)%bc_code = bc_code
 end subroutine

 subroutine add_lin_var_to_kxky_eqns_in_full_recipe(self, linvarname)
    class(full_problem_recipe_T) :: self
    type(coupled_system_vars_info_T), dimension(:), allocatable :: aux
    character(len=:), allocatable, intent(in) :: linvarName
    integer :: this_recipe, this_var
    integer :: ivar
   
    this_recipe = self%numberOf_coupled_kxky_systems 
    self%kxky_recipes(this_recipe)%n_coupled_vars = self%kxky_recipes(this_recipe)%n_coupled_vars + 1
    allocate (aux(self%kxky_recipes(this_recipe)%n_coupled_vars))
    aux(1:self%kxky_recipes(this_recipe)%n_coupled_vars -1) = self%kxky_recipes(this_recipe)%vars
    deAllocate(self%kxky_recipes(this_recipe)%vars)
    call move_alloc(from=aux, to=self%kxky_recipes(this_recipe)%vars)
    this_var = self%kxky_recipes(this_recipe)%n_coupled_vars
    self%kxky_recipes(this_recipe)%vars(this_var)%name = linvarname
    self%kxky_recipes(this_recipe)%vars(this_var)%bc_code = 0
    !we need to keep track of which set of equation this variable belongs to
    do ivar = 1, self%numberOf_linear_variables_k
        if (linvarname==self%linear_vars_kxky(ivar)%str) then
           if (self%linear_vars_kxky(ivar)%belongs_to_set.ne.0) then
              print *, self%linear_vars_kxky(ivar)%str,' belongs to set of equations', &
                       self%linear_vars_kxky(ivar)%belongs_to_set
              stop 'a linear variable already belongs to a set'
           end if
           if (self%linear_vars_kxky(ivar)%at_position.ne.0) then
              stop 'a linear variable already has a position'  
           end if
           self%linear_vars_kxky(ivar)%belongs_to_set = this_recipe
           self%linear_vars_kxky(ivar)%at_position    = this_var
           exit
        end if
     end do
 end subroutine


 subroutine add_lin_var_to_zero_eqns_in_full_recipe(self, linvarname)
    class(full_problem_recipe_T) :: self
    type(coupled_system_vars_info_T), dimension(:), allocatable :: aux
    character(len=:), allocatable, intent(in) :: linvarName
    integer :: this_recipe, this_var
    integer :: ivar
    logical :: foundIt
   
    this_recipe = self%numberOf_coupled_zero_systems 
    self%zero_recipes(this_recipe)%n_coupled_vars = self%zero_recipes(this_recipe)%n_coupled_vars + 1
    allocate (aux(self%zero_recipes(this_recipe)%n_coupled_vars))
    aux(1:self%zero_recipes(this_recipe)%n_coupled_vars -1) = self%zero_recipes(this_recipe)%vars
    deAllocate(self%zero_recipes(this_recipe)%vars)
    call move_alloc(from=aux, to=self%zero_recipes(this_recipe)%vars)
    this_var = self%zero_recipes(this_recipe)%n_coupled_vars
    self%zero_recipes(this_recipe)%vars(this_var)%name = linvarname
    self%zero_recipes(this_recipe)%vars(this_var)%bc_code = 0
    !we need to keep track of which set of equation this variable belongs to
    foundIt = .False.
    do ivar = 1, self%numberOf_linear_variables_0
        if (linvarname==self%linear_vars_zero(ivar)%str) then
           if (self%linear_vars_zero(ivar)%belongs_to_set.ne.0) then
              print *, self%linear_vars_zero(ivar)%str,' belongs to set of equations', &
                       self%linear_vars_zero(ivar)%belongs_to_set
              stop 'a linear variable already belongs to a set'
           end if
           if (self%linear_vars_zero(ivar)%at_position.ne.0) then
              stop 'a linear variable already has a position'  
           end if
           self%linear_vars_zero(ivar)%belongs_to_set = this_recipe
           self%linear_vars_zero(ivar)%at_position    = this_var
           foundIt = .True.
           exit
        end if
     end do
     if (.not.foundIt) then
        print *, 'Below the >>add_set_of_coupled_zero_equations << tag,'
        print *, 'the following variable was entered:'
        print *, linvarname
        print *, 'This name does not correspond to any variable declared' 
        print *, 'in the >>linear_variable mean section, as it should.'
        stop     'Please fix coral.equations'
     end if
 end subroutine

 subroutine add_equation_to_zero_coupled_set(self, eqn_order)
   class(full_problem_recipe_T) :: self
   integer, intent(in) :: eqn_order
   integer, allocatable :: aux(:)
   self%zero_recipes(self%numberOf_coupled_zero_systems)%n_equations = &
   self%zero_recipes(self%numberOf_coupled_zero_systems)%n_equations + 1
   allocate(aux( self%zero_recipes(self%numberOf_coupled_zero_systems)%n_equations))
   aux(1:size(aux,1)-1) = self%zero_recipes(self%numberOf_coupled_zero_systems)%eqn_order
   aux(size(aux,1)) = eqn_order
   deAllocate( self%zero_recipes(self%numberOf_coupled_zero_systems)%eqn_order)
   call move_alloc(from=aux, to =self%zero_recipes(self%numberOf_coupled_zero_systems)%eqn_order)  
 end subroutine
    

 subroutine add_equation_to_kxky_coupled_set(self, eqn_order)
   class(full_problem_recipe_T) :: self
   integer, intent(in) :: eqn_order
   integer, allocatable :: aux(:)
   self%kxky_recipes(self%numberOf_coupled_kxky_systems)%n_equations = &
   self%kxky_recipes(self%numberOf_coupled_kxky_systems)%n_equations + 1
   allocate(aux( self%kxky_recipes(self%numberOf_coupled_kxky_systems)%n_equations))
   aux(1:size(aux,1)-1) = self%kxky_recipes(self%numberOf_coupled_kxky_systems)%eqn_order
   aux(size(aux,1)) = eqn_order
   deAllocate( self%kxky_recipes(self%numberOf_coupled_kxky_systems)%eqn_order)
   call move_alloc(from=aux, to =self%kxky_recipes(self%numberOf_coupled_kxky_systems)%eqn_order)  
 end subroutine
    
 subroutine add_coupled_zero_equations_to_full_recipe(self)
    class(full_problem_recipe_T) :: self
    type(coupled_system_recipe_T), dimension(:), allocatable :: aux
    integer :: this_recipe
    self%numberOf_coupled_zero_systems = self%numberOf_coupled_zero_systems + 1 
    allocate(aux(self%numberOf_coupled_zero_systems))
    aux(1:self%numberOf_coupled_zero_systems-1) = self%zero_recipes
    deAllocate(self%zero_recipes)
    call move_alloc(from=aux, to=self%zero_recipes)
   
    this_recipe = self%numberOf_coupled_zero_systems 
    self%zero_recipes(this_recipe)%n_coupled_vars = 0 
    allocate(self%zero_recipes(this_recipe)%vars (0))
    allocate(self%zero_recipes(this_recipe)%eqn_order(0))
    self%zero_recipes(this_recipe)%n_equations = 0
    self%zero_recipes(this_recipe)%stif%n_pieces = 0
    self%zero_recipes(this_recipe)%mass%n_pieces = 0
    self%zero_recipes(this_recipe)%NL%n_pieces = 0
    allocate(self%zero_recipes(this_recipe)%mass%dsca (0))
    allocate(self%zero_recipes(this_recipe)%mass%jVar (0))
    allocate(self%zero_recipes(this_recipe)%mass%iEqn (0))
    allocate(self%zero_recipes(this_recipe)%mass%kx_exponent (0))
    allocate(self%zero_recipes(this_recipe)%mass%ky_exponent (0))
    allocate(self%zero_recipes(this_recipe)%mass%Iz_exponent (0))
    allocate(self%zero_recipes(this_recipe)%mass%z_multiply  (0))
    allocate(self%zero_recipes(this_recipe)%stif%dsca (0))
    allocate(self%zero_recipes(this_recipe)%stif%jVar (0))
    allocate(self%zero_recipes(this_recipe)%stif%iEqn (0))
    allocate(self%zero_recipes(this_recipe)%stif%kx_exponent (0))
    allocate(self%zero_recipes(this_recipe)%stif%ky_exponent (0))
    allocate(self%zero_recipes(this_recipe)%stif%Iz_exponent (0))
    allocate(self%zero_recipes(this_recipe)%stif%z_multiply  (0))
    allocate(self%zero_recipes(this_recipe)%NL%term(0))
   
 end subroutine

 subroutine add_coupled_kxky_equations_to_full_recipe(self)
    class(full_problem_recipe_T) :: self
    type(coupled_system_recipe_T), dimension(:), allocatable :: aux
    integer :: this_recipe
    self%numberOf_coupled_kxky_systems = self%numberOf_coupled_kxky_systems + 1 
    allocate(aux(self%numberOf_coupled_kxky_systems))
    aux(1:self%numberOf_coupled_kxky_systems-1) = self%kxky_recipes
    deAllocate(self%kxky_recipes)
    call move_alloc(from=aux, to=self%kxky_recipes)
   
    this_recipe = self%numberOf_coupled_kxky_systems 
    self%kxky_recipes(this_recipe)%n_coupled_vars = 0 
    allocate(self%kxky_recipes(this_recipe)%vars (0))
    allocate(self%kxky_recipes(this_recipe)%eqn_order(0))
    self%kxky_recipes(this_recipe)%n_equations = 0
    self%kxky_recipes(this_recipe)%mass%n_pieces = 0
    self%kxky_recipes(this_recipe)%stif%n_pieces = 0
    self%kxky_recipes(this_recipe)%NL%n_pieces = 0
    allocate(self%kxky_recipes(this_recipe)%mass%dsca (0))
    allocate(self%kxky_recipes(this_recipe)%mass%jVar (0))
    allocate(self%kxky_recipes(this_recipe)%mass%iEqn (0))
    allocate(self%kxky_recipes(this_recipe)%mass%kx_exponent (0))
    allocate(self%kxky_recipes(this_recipe)%mass%ky_exponent (0))
    allocate(self%kxky_recipes(this_recipe)%mass%Iz_exponent (0))
    allocate(self%kxky_recipes(this_recipe)%mass%z_multiply  (0))
    allocate(self%kxky_recipes(this_recipe)%stif%dsca (0))
    allocate(self%kxky_recipes(this_recipe)%stif%jVar (0))
    allocate(self%kxky_recipes(this_recipe)%stif%iEqn (0))
    allocate(self%kxky_recipes(this_recipe)%stif%kx_exponent (0))
    allocate(self%kxky_recipes(this_recipe)%stif%ky_exponent (0))
    allocate(self%kxky_recipes(this_recipe)%stif%Iz_exponent (0))
    allocate(self%kxky_recipes(this_recipe)%stif%z_multiply  (0))
    allocate(self%kxky_recipes(this_recipe)%NL%term(0))
   
 end subroutine

 Subroutine add_quadratic_variable_to_full_recipe(self, qvname, lv1name, lv2name)
   class(full_problem_recipe_T) :: self
   type(list_of_NLVars_T), dimension(:), allocatable :: temporary_list
   character(len=:), intent(in), allocatable :: qvName 
   character(len=:), intent(in), allocatable :: lv1Name, lv2Name
   integer :: ivar, entry_index
   ! first initialize the list if it was empty
   if (self%numberOf_quadratic_variables.eq.0) then
       self%numberOf_quadratic_variables = 1
       deAllocate(self%nl_vars)
       allocate  (self%nl_vars(1))
   ! or extend the existing allocatable                           
   else
      self%numberOf_quadratic_variables = self%numberOf_quadratic_variables + 1 
      allocate(temporary_list(self%numberOf_quadratic_variables))
      temporary_list(1:self%numberOf_quadratic_variables-1) = self%nl_vars
      deAllocate(self%nl_vars)
      Allocate(self%nl_vars ( self%numberOf_quadratic_variables))
      call move_alloc(from=temporary_list, to=self%nl_vars)
   end if 
   ! now fill-in the last entry
   entry_index = self%numberOf_quadratic_variables 
   self%nl_vars(entry_index)%str = qvname
   self%nl_vars(entry_index)%svar1 = lv1name
   self%nl_vars(entry_index)%svar2 = lv2name
   self%nl_vars(entry_index)%ivar1 = 0
   self%nl_vars(entry_index)%ivar2 = 0
   do ivar = 1, self%numberOf_Linear_variables_full
     if (lv1name==self%linear_vars_full(ivar)%str) then 
        self%nl_vars(entry_index)%ivar1 = ivar
     end if
     if (lv2name==self%linear_vars_full(ivar)%str) then 
        self%nl_vars(entry_index)%ivar2 = ivar
     end if
   end do
   if ((self%NL_vars(entry_index)%ivar1==0) .or. &
       (self%NL_vars(entry_index)%ivar2==0) ) then
       STOP 'bad NL var definition (wrong reference to linear var.)'
   end if 
 End Subroutine add_quadratic_variable_to_full_recipe

 Subroutine add_linear_variable_to_full_recipe(self, varName, varKind)
   class(full_problem_recipe_T) :: self
   type(list_of_linVars_T), dimension(:), allocatable :: temporary_list
   character(len=:), intent(in), allocatable :: varName
   character(len=4), intent(in) :: varKind
   if (varKind.eq.'kxky') then
      self%numberOf_linear_variables_k = self%numberOf_linear_variables_k + 1 
      allocate(temporary_list(self%numberOf_linear_variables_k))
      temporary_list(1:self%numberOf_linear_variables_k-1) = self%linear_vars_kxky
      call move_alloc(from=temporary_list, to=self%linear_vars_kxky)
      allocate(character(len(varName)) :: self%linear_vars_kxky(self%numberOf_linear_variables_k)%str)
      self%linear_vars_kxky(self%numberOf_linear_variables_k)%str = varName
      self%linear_vars_kxky(self%numberOf_linear_variables_k)%belongs_to_set = 0
      self%linear_vars_kxky(self%numberOf_linear_variables_k)%at_position   = 0 
   else if (varKind.eq.'zero') then
      self%numberOf_linear_variables_0 = self%numberOf_linear_variables_0 + 1 
      allocate(temporary_list(self%numberOf_linear_variables_0))
      temporary_list(1:self%numberOf_linear_variables_0-1) = self%linear_vars_zero
      call move_alloc(from=temporary_list, to=self%linear_vars_zero)
      allocate(character(len(varName)) :: self%linear_vars_zero(self%numberOf_linear_variables_0)%str)
      self%linear_vars_zero(self%numberOf_linear_variables_0)%str = varName
      self%linear_vars_zero(self%numberOf_linear_variables_0)%belongs_to_set = 0
      self%linear_vars_zero(self%numberOf_linear_variables_0)%at_position   = 0
   else
      Stop 'bad varKind value in add_linear_variable_to_full_recipe'
   end if
 End Subroutine add_linear_variable_to_full_recipe


 Subroutine add_source_parameter_to_full_recipe(self, string, dsca)
   class(full_problem_recipe_T) :: self
   type(sourceParams_T), allocatable :: temporary_list(:)
   character(len=:), intent(in), allocatable :: string
   real(dp), intent(in) :: dsca
   integer :: numberOf_source_params
   numberOf_source_params = size(self%sources%sourceParams) + 1 
   allocate (temporary_list( numberOf_source_params ) )
   temporary_list (1: numberOf_source_params-1) = self%sources%sourceParams(1:numberOf_source_params-1)
   temporary_list (   numberOf_source_params  ) % str = string
   temporary_list (   numberOf_source_params  ) % dsca= dsca
   call move_alloc (from= temporary_list, to= self%sources%sourceParams )
 End Subroutine add_source_parameter_to_full_recipe

 Subroutine add_parameter_to_full_recipe(self, string, dsca)
   class(full_problem_recipe_T) :: self
   type(list_of_parameters_T), dimension(:), allocatable :: temporary_list
   character(len=:), intent(in), allocatable :: string
   real(dp), intent(in) :: dsca
   self%numberOf_parameters = self%numberOf_parameters + 1 
   call self%copy_params(temporary_list)
   deAllocate(self%list_parameters)
   Allocate(self%list_parameters( self%numberOf_parameters))
   call copy_list_of_parameters( temporary_list, self%list_parameters)
   allocate(character(len(string)) :: self%list_parameters(self%numberOf_parameters)%str)
   self%list_parameters(self%numberOf_parameters)%str = string
   self%list_parameters(self%numberOf_parameters)%val = dsca   
 End Subroutine add_parameter_to_full_recipe


 Subroutine copy_list_of_parameters_from_recipe(self, list_of_p)
   class(full_problem_recipe_T) :: self
   type(list_of_parameters_T), dimension(:), allocatable, intent(inOut) :: list_of_p
   call copy_list_of_parameters(self%list_parameters, list_of_p)
 end subroutine copy_list_of_parameters_from_recipe

 Subroutine copy_list_of_parameters(self, other)
   type(list_of_parameters_T), dimension(:), allocatable :: self
   type(list_of_parameters_T), dimension(:), allocatable, intent(inOut) :: other
   integer :: i 
   if (.not.allocated(other)) then
     allocate(other(size(self)))
   end if
   do i = 1,size(self)
      other(i)%val = self(i)%val
      allocate(character(len(self(i)%str)) :: other(i)%str)
      other(i)%str = self(i)%str
   end do
 End Subroutine copy_list_of_parameters

   

 elemental subroutine move_full_recipes(mySource, myTarget)
   type(fullVars_recipe_T), intent(inOut) :: mySource, myTarget
   call move_alloc( from=mySource%str, to=myTarget%str)
   myTarget%N_terms = mySource%N_terms
   call move_alloc( from=mySource%term, to=myTarget%term)
 end subroutine
 
 subroutine add_full_linear_variable_to_full_recipe(self, linear_var_name)
   class(full_problem_recipe_T) :: self
   character(len=:), allocatable, intent(in) :: linear_var_name
   type(fullVars_recipe_T), dimension(:), allocatable :: tmp
   ! increment the total number of "full variables"
   self%numberOf_linear_variables_full = self%numberOf_linear_variables_full + 1
   ! now copy the existing recipes in a bigger array
   allocate(tmp(self%numberOf_linear_variables_full))
   call move_full_recipes(self%linear_vars_full, tmp(1:self%numberOf_linear_variables_full -1))
   tmp(self%numberof_linear_variables_full)%N_terms = 0
   allocate(character(len=len(linear_var_name)) :: tmp(self%numberOf_linear_variables_full)%str)
   tmp(self%numberOf_linear_variables_full)%str = linear_var_name
   allocate(tmp(self%numberOf_linear_variables_full)%term(0))
   call move_alloc(tmp, self%linear_vars_full)
 end subroutine add_full_linear_variable_to_full_recipe

 subroutine build_kxky_NL_term_in_full_recipe( self, atom_linear_ops)
  class(full_problem_recipe_T) :: self
  type(atom_linearOp_T), dimension(:), allocatable, intent(in) :: atom_linear_ops
  type(atom_NL_term_T), allocatable :: buffer_NL(:)
  integer :: n_old 
  integer :: n_new_pieces, counter_eqn
  integer :: ivar, i_str
  integer :: iAtom, this_recipe
  logical :: this_string_is_identified
  ! first expand the list of terms, and copy the old terms
  this_recipe = self%numberOf_coupled_kxky_systems
  counter_eqn = self%kxky_recipes(this_recipe)%n_equations
  n_new_pieces = size(atom_linear_ops,1)
  
  n_old = self%kxky_recipes(this_recipe)%NL%n_pieces
  self%kxky_recipes(this_recipe)%NL%n_pieces = &
  self%kxky_recipes(this_recipe)%NL%n_pieces + n_new_pieces

  allocate( buffer_NL( n_old + n_new_pieces))
  buffer_NL(1:n_old) = self%kxky_recipes(this_recipe)%NL%term(:)
  call move_alloc(from=buffer_NL, to=self%kxky_recipes(this_recipe)%NL%term)  

  ! now add the contributions encoded in atom_linear_ops
  do iAtom = 1, size(atom_linear_ops)
   self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%dx_exponent = atom_linear_ops(iAtom)%dx_exponent
   self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%dy_exponent = atom_linear_ops(iAtom)%dy_exponent
   self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%Iz_exponent = atom_linear_ops(iAtom)%Iz_exponent
   self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%z_multiply  = atom_linear_ops(iAtom)%z_multiply 
   self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%dsca = atom_linear_ops(iAtom)%dsca         
   self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%quad_var_index = 0     
   self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%eqn_index = counter_eqn
   do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
       this_string_is_identified = .False.
       ! perhaps the unsorted string is a field?
       do ivar = 1, self%numberOf_quadratic_variables               
         if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%nl_vars(ivar)%str) then 
          self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%quad_var_index = ivar
          this_string_is_identified = .True.
          exit
         end if
       end do
       ! or perhaps this is a parameter ?   
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
          do ivar = 1, self%numberOf_parameters                  
            if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%list_parameters(ivar)%str) then 
             self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%dsca = self%list_parameters(ivar)%val&
            *self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%dsca
            this_string_is_identified = .True.
            exit ! we identified the string. move on to the next
            end if
          end do
          if (.not.this_string_is_identified) then
          print *, 'unknown string in equations:', atom_linear_ops(iAtom)%unsorted(i_str)%str
          error stop                                                                          
          end if
       end if

   end do
   if (self%kxky_recipes(this_recipe)%NL%term(n_old+iAtom)%quad_var_index==0) then
      stop 'invalid quad_var_index'
   end if

  end do
 end subroutine

 subroutine build_zero_sources_in_full_recipe( self, atom_linear_ops)
  class(full_problem_recipe_T) :: self
  type(atom_linearOp_T), dimension(:), allocatable, intent(in) :: atom_linear_ops
  type(source_term_T), allocatable :: buffer_source(:)
  integer :: n_old 
  integer :: n_new_pieces, counter_eqn
  integer :: ivar, i_str
  integer :: iAtom, this_recipe  
  integer :: stat
   logical :: this_string_is_identified
  ! first expand the list of terms, and copy the old terms
  this_recipe = self%numberOf_coupled_zero_systems
  counter_eqn = self%zero_recipes(this_recipe)%n_equations
  n_new_pieces = size(atom_linear_ops,1)
  
  n_old = self%sources%n_sources
  self%sources%n_sources = &
  self%sources%n_sources + n_new_pieces

  allocate( buffer_source( n_old + n_new_pieces))
  buffer_source(1:n_old) = self%sources%term(:)
  call move_alloc(from=buffer_source, to=self%sources%term)  

  ! now add the contributions encoded in atom_linear_ops
  do iAtom = 1, size(atom_linear_ops)
   if (atom_linear_ops(iAtom)%dy_exponent.ne.0) print *, 'error, y-derivative found in 0-recipe!!'
   if (atom_linear_ops(iAtom)%dy_exponent.ne.0) error stop
   if (atom_linear_ops(iAtom)%dx_exponent.ne.0) print *, 'error, x-derivative found in 0-recipe!!'
   if (atom_linear_ops(iAtom)%dx_exponent.ne.0) error stop
   self%sources%term(n_old+iAtom)%dx_exponent = 0
   self%sources%term(n_old+iAtom)%dy_exponent = 0
   self%sources%term(n_old+iAtom)%Iz_exponent = atom_linear_ops(iAtom)%Iz_exponent
   self%sources%term(n_old+iAtom)%dsca = atom_linear_ops(iAtom)%dsca
   self%sources%term(n_old+iAtom)%source_index = 0     
   self%sources%term(n_old+iAtom)%eqn_index = counter_eqn
   self%sources%term(n_old+iAtom)%sys_index = this_recipe
   do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
       this_string_is_identified = .False.
       ! perhaps the unsorted string is a source?
       if (len(atom_linear_ops(iAtom)%unsorted(i_str)%str).eq.8) then
       if     (atom_linear_ops(iAtom)%unsorted(i_str)%str(1:6) .eq.'source' ) then
          read(atom_linear_ops(iAtom)%unsorted(i_str)%str(7:8),*,iostat=stat)&
               self%sources%term(n_old+iAtom)%source_index
          if (stat.ne.0) print *, 'Source terms must be of the form sourceNN, where NN are 2 integers'
          if (stat.ne.0) error stop
          this_string_is_identified = .True.
          exit
       end if
       end if
       ! or perhaps this is a parameter ?   
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
       do ivar = 1, self%numberOf_parameters                  
         if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%list_parameters(ivar)%str) then 
           self%sources%term(n_old+iAtom)%dsca = self%list_parameters(ivar)%val&
          *self%sources%term(n_old+iAtom)%dsca
            this_string_is_identified = .True.
            exit ! we identified the string. move on to the next
         end if
       end do
          if (.not.this_string_is_identified) then
          print *, 'unknown string in equations:', atom_linear_ops(iAtom)%unsorted(i_str)%str
          error stop                                                                          
          end if
       end if

   end do
   if (self%sources%term(n_old+iAtom)%source_index==0) then
      stop 'invalid source_index'
   end if

  end do
 end subroutine


 subroutine build_zero_NL_term_in_full_recipe( self, atom_linear_ops)
  class(full_problem_recipe_T) :: self
  type(atom_linearOp_T), dimension(:), allocatable, intent(in) :: atom_linear_ops
  type(atom_NL_term_T), allocatable :: buffer_NL(:)
  integer :: n_old 
  integer :: n_new_pieces, counter_eqn
  integer :: ivar, i_str
  integer :: iAtom, this_recipe
   logical :: this_string_is_identified
  ! first expand the list of terms, and copy the old terms
  this_recipe = self%numberOf_coupled_zero_systems
  counter_eqn = self%zero_recipes(this_recipe)%n_equations
  n_new_pieces = size(atom_linear_ops,1)
  
  n_old = self%zero_recipes(this_recipe)%NL%n_pieces
  self%zero_recipes(this_recipe)%NL%n_pieces = &
  self%zero_recipes(this_recipe)%NL%n_pieces + n_new_pieces

  allocate( buffer_NL( n_old + n_new_pieces))
  buffer_NL(1:n_old) = self%zero_recipes(this_recipe)%NL%term(:)
  call move_alloc(from=buffer_NL, to=self%zero_recipes(this_recipe)%NL%term)  

  ! now add the contributions encoded in atom_linear_ops
  do iAtom = 1, size(atom_linear_ops)
   if (atom_linear_ops(iAtom)%dy_exponent.ne.0) print *, 'error, y-derivative found in 0-recipe!!'
   if (atom_linear_ops(iAtom)%dy_exponent.ne.0) error stop
   if (atom_linear_ops(iAtom)%dx_exponent.ne.0) print *, 'error, x-derivative found in 0-recipe!!'
   if (atom_linear_ops(iAtom)%dx_exponent.ne.0) error stop
   self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%dx_exponent = 0
   self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%dy_exponent = 0
   self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%Iz_exponent = atom_linear_ops(iAtom)%Iz_exponent
   self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%z_multiply  = atom_linear_ops(iAtom)%z_multiply 
   self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%dsca = atom_linear_ops(iAtom)%dsca         
   self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%quad_var_index = 0     
   self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%eqn_index = counter_eqn
   do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
       this_string_is_identified = .False.
       ! perhaps the unsorted string is a field?
       do ivar = 1, self%numberOf_quadratic_variables               
         if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%nl_vars(ivar)%str) then 
          self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%quad_var_index = ivar
          this_string_is_identified = .True.
          exit
         end if
       end do
       ! or perhaps this is a parameter ?   
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
       do ivar = 1, self%numberOf_parameters                  
         if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%list_parameters(ivar)%str) then 
           self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%dsca = self%list_parameters(ivar)%val&
          *self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%dsca
          this_string_is_identified = .True.
         end if
       end do
          if (.not.this_string_is_identified) then
          print *, 'unknown string in equations:', atom_linear_ops(iAtom)%unsorted(i_str)%str
          error stop                                                                          
          end if
       end if

   end do
   if (self%zero_recipes(this_recipe)%NL%term(n_old+iAtom)%quad_var_index==0) then
      stop 'invalid quad_var_index'
   end if

  end do
 end subroutine

 subroutine build_current_full_var_in_full_recipe( self, atom_linear_ops)
   class(full_problem_recipe_T) :: self
   type(atom_linearOp_T), dimension(:), allocatable, intent(in) :: atom_linear_ops
   type(a_linear_contribution_T), dimension(:), allocatable :: buf_term
   integer :: n_old 
   integer :: i_str
   integer :: ivar
   integer :: iAtom, iTerm
   logical :: this_string_is_identified
   ! first expand the list of terms, and copy the old terms
   n_old = self%linear_vars_full(self%numberOf_linear_variables_full)%n_terms 
   self%linear_vars_full(self%numberOf_linear_variables_full)%n_terms = &
        self%linear_vars_full(self%numberOf_linear_variables_full)%n_terms + size(atom_linear_ops)
   allocate( buf_term( self%linear_vars_full(self%numberOf_linear_variables_full)%n_terms ))
   do iterm = 1, n_old
      buf_term(iterm) = self%linear_vars_full(self%numberOf_linear_variables_full)%term(iterm)
   end do
   ! now add the contributions encoded in atom_linear_ops
   do iAtom = 1, size(atom_linear_ops)
      buf_term(n_old+iAtom)%dx_exponent = atom_linear_ops(iAtom)%dx_exponent
      buf_term(n_old+iAtom)%dy_exponent = atom_linear_ops(iAtom)%dy_exponent
      buf_term(n_old+iAtom)%dz_exponent = atom_linear_ops(iAtom)%Iz_exponent
      buf_term(n_old+iAtom)%z_multiply  = atom_linear_ops(iAtom)%z_multiply 
      buf_term(n_old+iAtom)%dsca = atom_linear_ops(iAtom)%dsca
      buf_term(n_old+iAtom)%var_kind = 'null'
      do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
       this_string_is_identified = .False.
         ! perhaps the unsorted string is mean field?
         do ivar = 1, self%numberOf_Linear_variables_0
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%linear_vars_zero(ivar)%str) then 
              buf_term(n_old+iAtom)%var_index  = ivar
              buf_term(n_old+iAtom)%var_kind   = 'zero'
          this_string_is_identified = .True.
          exit
           end if
         end do
         ! or perhaps this is a fluctuation ? 
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
         do ivar = 1, self%numberOf_Linear_variables_k
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%linear_vars_kxky(ivar)%str) then 
             buf_term(n_old+iAtom)%var_index  = ivar
             buf_term(n_old+iAtom)%var_kind   = 'kxky'
          this_string_is_identified = .True.
          exit
           end if
         end do
 end if
         ! or perhaps this is a parameter ?   
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
         do ivar = 1, self%numberOf_parameters                  
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%list_parameters(ivar)%str) then 
             buf_term(n_old+iAtom)%dsca = buf_term(n_old+iAtom)%dsca * self%list_parameters(ivar)%val
          this_string_is_identified = .True.
          exit
           end if
         end do
          if (.not.this_string_is_identified) then
          print *, 'unknown string in equations:', atom_linear_ops(iAtom)%unsorted(i_str)%str
          error stop                                                                          
          end if
       end if

      end do
      if (buf_term(n_old+iAtom)%var_kind=='null') then
         STOP 'invalid full variable def.'
      end if

   end do
   call move_alloc(from=buf_term, to=self%linear_vars_full(self%numberOf_linear_variables_full)%term) 
 end subroutine

 subroutine build_zero_stif_matrix_in_full_recipe( self, atom_linear_ops)
   class(full_problem_recipe_T), intent(inOut) :: self
   type(atom_linearOp_T), dimension(:), allocatable, intent(in) :: atom_linear_ops
   integer :: n_old 
   integer :: n_new_pieces, counter_eqn
   integer :: ivar, i_str
   integer :: iAtom, this_recipe
   integer, allocatable :: buf_int(:)
   real(dp), allocatable :: buf_real(:)
   logical :: this_string_is_identified
   ! first expand the list of terms, and copy the old terms
   this_recipe = self%numberOf_coupled_zero_systems
   counter_eqn = self%zero_recipes(this_recipe)%n_equations
   n_new_pieces = size(atom_linear_ops,1)
   
   n_old = self%zero_recipes(this_recipe)%stif%n_pieces
   self%zero_recipes(this_recipe)%stif%n_pieces = &
   self%zero_recipes(this_recipe)%stif%n_pieces + n_new_pieces

   allocate( buf_real( n_old + n_new_pieces))
   buf_real(1:n_old) = self%zero_recipes(this_recipe)%stif%dsca
   call move_alloc(from=buf_real, to=self%zero_recipes(this_recipe)%stif%dsca)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%stif%kx_exponent
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%stif%kx_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%stif%ky_exponent
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%stif%ky_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%stif%Iz_exponent
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%stif%Iz_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%stif%z_multiply 
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%stif%z_multiply )

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%stif%jVar
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%stif%jVar)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%stif%iEqn
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%stif%iEqn)

   ! now add the contributions encoded in atom_linear_ops
   do iAtom = 1, size(atom_linear_ops)
     self%zero_recipes(this_recipe)%stif%kx_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%dx_exponent
     self%zero_recipes(this_recipe)%stif%ky_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%dy_exponent
     self%zero_recipes(this_recipe)%stif%Iz_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%Iz_exponent
     self%zero_recipes(this_recipe)%stif%z_multiply (n_old+iAtom) = atom_linear_ops(iAtom)%z_multiply 
     self%zero_recipes(this_recipe)%stif%dsca(n_old+iAtom) = atom_linear_ops(iAtom)%dsca         
     self%zero_recipes(this_recipe)%stif%jVar(n_old+iAtom) = 0                               
     self%zero_recipes(this_recipe)%stif%iEqn(n_old+iAtom) = counter_eqn
     do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
       this_string_is_identified = .False.
         ! perhaps the unsorted string is a field?
         do ivar = 1, self%zero_recipes(this_recipe)%n_coupled_vars
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%zero_recipes(this_recipe)%vars(ivar)%name) then 
            self%zero_recipes(this_recipe)%stif%jVar(n_old+iAtom) = ivar
          this_string_is_identified = .True.
          exit
           end if
         end do
         ! or perhaps this is a parameter ?   
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
         do ivar = 1, self%numberOf_parameters                  
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%list_parameters(ivar)%str) then 
             self%zero_recipes(this_recipe)%stif%dsca(n_old+iAtom) =  self%list_parameters(ivar)%val&
            *self%zero_recipes(this_recipe)%stif%dsca(n_old+iAtom)
            this_string_is_identified = .True.
            exit ! we identified the string. move on to the next
           end if
         end do
          if (.not.this_string_is_identified) then
          print *, 'unknown string in equations:', atom_linear_ops(iAtom)%unsorted(i_str)%str
          error stop                                                                          
          end if
       end if

      end do
      if (self%zero_recipes(this_recipe)%stif%jVar(n_old+iAtom)==0) then
         stop 'invalid zero stif matrix definition'
      end if

   end do
 end subroutine
   
 subroutine build_kxky_stif_matrix_in_full_recipe( self, atom_linear_ops)
   class(full_problem_recipe_T), intent(inOut) :: self
   type(atom_linearOp_T), dimension(:), allocatable, intent(in) :: atom_linear_ops
   integer :: n_old 
   integer :: n_new_pieces, counter_eqn
   integer :: ivar, i_str
   integer :: iAtom, this_recipe
   integer, allocatable :: buf_int(:)
   real(dp), allocatable :: buf_real(:)
   logical :: this_string_is_identified
   ! first expand the list of terms, and copy the old terms
   this_recipe = self%numberOf_coupled_kxky_systems
   counter_eqn = self%kxky_recipes(this_recipe)%n_equations
   n_new_pieces = size(atom_linear_ops,1)
   
   n_old = self%kxky_recipes(this_recipe)%stif%n_pieces
   self%kxky_recipes(this_recipe)%stif%n_pieces = &
   self%kxky_recipes(this_recipe)%stif%n_pieces + n_new_pieces


   allocate( buf_real( n_old + n_new_pieces))
   buf_real(1:n_old) = self%kxky_recipes(this_recipe)%stif%dsca
   call move_alloc(from=buf_real, to=self%kxky_recipes(this_recipe)%stif%dsca)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%stif%kx_exponent
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%stif%kx_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%stif%ky_exponent
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%stif%ky_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%stif%Iz_exponent
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%stif%Iz_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%stif%z_multiply
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%stif%z_multiply )

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%stif%jVar
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%stif%jVar)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%stif%iEqn
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%stif%iEqn)

   ! now add the contributions encoded in atom_linear_ops
   do iAtom = 1, size(atom_linear_ops)
     self%kxky_recipes(this_recipe)%stif%kx_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%dx_exponent
     self%kxky_recipes(this_recipe)%stif%ky_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%dy_exponent
     self%kxky_recipes(this_recipe)%stif%Iz_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%Iz_exponent
     self%kxky_recipes(this_recipe)%stif%z_multiply (n_old+iAtom) = atom_linear_ops(iAtom)%z_multiply 
     self%kxky_recipes(this_recipe)%stif%dsca(n_old+iAtom) = atom_linear_ops(iAtom)%dsca         
     self%kxky_recipes(this_recipe)%stif%jVar(n_old+iAtom) = 0                               
     self%kxky_recipes(this_recipe)%stif%iEqn(n_old+iAtom) = counter_eqn
     do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
       this_string_is_identified = .False.
         ! perhaps the unsorted string is a field?
         do ivar = 1, self%kxky_recipes(this_recipe)%n_coupled_vars
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%kxky_recipes(this_recipe)%vars(ivar)%name) then 
            self%kxky_recipes(this_recipe)%stif%jVar(n_old+iAtom) = ivar
          this_string_is_identified = .True.
          exit
           end if
         end do
         ! or perhaps this is a parameter ?   
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
         do ivar = 1, self%numberOf_parameters                  
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%list_parameters(ivar)%str) then 
             self%kxky_recipes(this_recipe)%stif%dsca(n_old+iAtom) =  self%list_parameters(ivar)%val&
            *self%kxky_recipes(this_recipe)%stif%dsca(n_old+iAtom)
            this_string_is_identified = .True.
            exit ! we identified the string. move on to the next
           end if
         end do
          if (.not.this_string_is_identified) then
          print *, 'unknown string in equations:', atom_linear_ops(iAtom)%unsorted(i_str)%str
          error stop                                                                          
          end if
       end if

      end do
      if (self%kxky_recipes(this_recipe)%stif%jVar(n_old+iAtom)==0) then
         do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
         print *, atom_linear_ops(iAtom)%unsorted(i_str)%str
         end do
         stop 'invalid kxky stif matrix definition'
      end if

   end do
 end subroutine
     
 subroutine build_zero_mass_matrix_in_full_recipe( self, atom_linear_ops)
   class(full_problem_recipe_T), intent(inOut) :: self
   type(atom_linearOp_T), dimension(:), allocatable, intent(in) :: atom_linear_ops
   integer :: n_old 
   integer :: n_new_pieces, counter_eqn
   integer :: ivar, i_str
   integer :: iAtom, this_recipe
   integer, allocatable :: buf_int(:)
   real(dp), allocatable :: buf_real(:)
   logical :: this_string_is_identified
   ! first expand the list of terms, and copy the old terms
   this_recipe = self%numberOf_coupled_zero_systems
   counter_eqn = self%zero_recipes(this_recipe)%n_equations
   n_new_pieces = size(atom_linear_ops,1)
   
   n_old = self%zero_recipes(this_recipe)%mass%n_pieces
   self%zero_recipes(this_recipe)%mass%n_pieces = &
   self%zero_recipes(this_recipe)%mass%n_pieces + n_new_pieces

   allocate( buf_real( n_old + n_new_pieces))
   buf_real(1:n_old) = self%zero_recipes(this_recipe)%mass%dsca
   call move_alloc(from=buf_real, to=self%zero_recipes(this_recipe)%mass%dsca)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%mass%kx_exponent
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%mass%kx_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%mass%ky_exponent
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%mass%ky_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%mass%Iz_exponent
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%mass%Iz_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%mass%z_multiply
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%mass%z_multiply )

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%mass%jVar
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%mass%jVar)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%zero_recipes(this_recipe)%mass%iEqn
   call move_alloc(from=buf_int, to=self%zero_recipes(this_recipe)%mass%iEqn)

   ! now add the contributions encoded in atom_linear_ops
   do iAtom = 1, size(atom_linear_ops)
     self%zero_recipes(this_recipe)%mass%kx_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%dx_exponent
     self%zero_recipes(this_recipe)%mass%ky_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%dy_exponent
     self%zero_recipes(this_recipe)%mass%Iz_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%Iz_exponent
     self%zero_recipes(this_recipe)%mass%z_multiply (n_old+iAtom) = atom_linear_ops(iAtom)%z_multiply 
     self%zero_recipes(this_recipe)%mass%dsca(n_old+iAtom) = atom_linear_ops(iAtom)%dsca         
     self%zero_recipes(this_recipe)%mass%jVar(n_old+iAtom) = 0                               
     self%zero_recipes(this_recipe)%mass%iEqn(n_old+iAtom) = counter_eqn
     do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
       this_string_is_identified = .False.
         ! perhaps the unsorted string is a field?
         do ivar = 1, self%zero_recipes(this_recipe)%n_coupled_vars
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%zero_recipes(this_recipe)%vars(ivar)%name) then 
            self%zero_recipes(this_recipe)%mass%jVar(n_old+iAtom) = ivar
          this_string_is_identified = .True.
          exit
           end if
         end do
         ! or perhaps this is a parameter ?   
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
         do ivar = 1, self%numberOf_parameters                  
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%list_parameters(ivar)%str) then 
             self%zero_recipes(this_recipe)%mass%dsca(n_old+iAtom) =  self%list_parameters(ivar)%val&
            *self%zero_recipes(this_recipe)%mass%dsca(n_old+iAtom)
            this_string_is_identified = .True.
            exit ! we identified the string. move on to the next
           end if
         end do
          if (.not.this_string_is_identified) then
          print *, 'unknown string in equations:', atom_linear_ops(iAtom)%unsorted(i_str)%str
          error stop                                                                          
          end if
       end if

      end do
      if (self%zero_recipes(this_recipe)%mass%jVar(n_old+iAtom)==0) then
         stop 'invalid zero mass matrix definition'
      end if

   end do
 end subroutine
     
   
 subroutine build_kxky_mass_matrix_in_full_recipe( self, atom_linear_ops)
   class(full_problem_recipe_T), intent(inOut) :: self
   type(atom_linearOp_T), dimension(:), allocatable, intent(in) :: atom_linear_ops
   integer :: n_old 
   integer :: n_new_pieces, counter_eqn
   integer :: ivar, i_str
   integer :: iAtom, this_recipe
   integer, allocatable :: buf_int(:)
   real(dp), allocatable :: buf_real(:)
   logical :: this_string_is_identified

   ! firstexpand the list of terms, and copy the old terms
   this_recipe = self%numberOf_coupled_kxky_systems
   counter_eqn = self%kxky_recipes(this_recipe)%n_equations
   n_new_pieces = size(atom_linear_ops,1)
   
   n_old = self%kxky_recipes(this_recipe)%mass%n_pieces
   self%kxky_recipes(this_recipe)%mass%n_pieces = &
   self%kxky_recipes(this_recipe)%mass%n_pieces + n_new_pieces

 
   allocate( buf_real( n_old + n_new_pieces))
   buf_real(1:n_old) = self%kxky_recipes(this_recipe)%mass%dsca
   call move_alloc(from=buf_real, to=self%kxky_recipes(this_recipe)%mass%dsca)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%mass%kx_exponent
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%mass%kx_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%mass%ky_exponent
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%mass%ky_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%mass%Iz_exponent
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%mass%Iz_exponent)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%mass%z_multiply
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%mass%z_multiply )

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%mass%jVar
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%mass%jVar)

   allocate( buf_int( n_old + n_new_pieces))
   buf_int(1:n_old) = self%kxky_recipes(this_recipe)%mass%iEqn
   call move_alloc(from=buf_int, to=self%kxky_recipes(this_recipe)%mass%iEqn)

   ! now add the contributions encoded in atom_linear_ops
   do iAtom = 1, size(atom_linear_ops)
     self%kxky_recipes(this_recipe)%mass%kx_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%dx_exponent
     self%kxky_recipes(this_recipe)%mass%ky_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%dy_exponent
     self%kxky_recipes(this_recipe)%mass%Iz_exponent(n_old+iAtom) = atom_linear_ops(iAtom)%Iz_exponent
     self%kxky_recipes(this_recipe)%mass%z_multiply (n_old+iAtom) = atom_linear_ops(iAtom)%z_multiply 
     self%kxky_recipes(this_recipe)%mass%dsca(n_old+iAtom) = atom_linear_ops(iAtom)%dsca         
     self%kxky_recipes(this_recipe)%mass%jVar(n_old+iAtom) = 0                               
     self%kxky_recipes(this_recipe)%mass%iEqn(n_old+iAtom) = counter_eqn
     do i_str = 1, atom_linear_ops(iAtom)%nOf_unsortedStrings
       this_string_is_identified = .False.
       ! perhaps the unsorted string is a field?
         ! perhaps the unsorted string is a field?
         do ivar = 1, self%kxky_recipes(this_recipe)%n_coupled_vars
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%kxky_recipes(this_recipe)%vars(ivar)%name) then 
            self%kxky_recipes(this_recipe)%mass%jVar(n_old+iAtom) = ivar
          this_string_is_identified = .True.
          exit
           end if
         end do
         ! or perhaps this is a parameter ?   
       ! or perhaps this is a parameter ?   
       if ( this_string_is_identified) then
          cycle ! we identified the string. move on to the next
       else
         do ivar = 1, self%numberOf_parameters                  
           if (atom_linear_ops(iAtom)%unsorted(i_str)%str==self%list_parameters(ivar)%str) then 
             self%kxky_recipes(this_recipe)%mass%dsca(n_old+iAtom) =  self%list_parameters(ivar)%val&
            *self%kxky_recipes(this_recipe)%mass%dsca(n_old+iAtom)
            this_string_is_identified = .True.
            exit ! we identified the string. move on to the next
           end if
         end do
          if (.not.this_string_is_identified) then
          print *, 'unknown string in equations:', atom_linear_ops(iAtom)%unsorted(i_str)%str
          error stop                                                                          
          end if
  end if

      end do
      if (self%kxky_recipes(this_recipe)%mass%jVar(n_old+iAtom)==0) then
         stop 'invalid kxky mass matrix definition'
      end if

   end do
 end subroutine
     


 subroutine determine_kind_of_information(text_line, bstep, eofSig, stopSig)
    character(len=1024),  intent(in) :: text_line
    integer, intent(inOut) :: bstep
    logical, intent(out) :: eofSig, stopSig
    eofSig  = .False.
    stopSig = .False.
    select case (bstep)
       ! from 'Param', accepted values are 'Param' or 'linear_variable_kxky' or 'sourceParam'
       case (1)
         if (text_line(3:7).eq.'Param') then 
            bstep = 1
         else if (text_line(3:13).eq.'sourceParam') then
            bstep = 11
         else if (text_line(3:22).eq.'linear_variable_kxky') then
            bstep = 2
         else
            stopSig = .True.
         end if
       ! from 'sourceParam', accepted values are 'Param' or 'linear_variable_kxky'
       case (11)
         if (text_line(3:7).eq.'Param') then 
            bstep = 1
         else if (text_line(3:13).eq.'sourceParam') then
            bstep = 11
         else if (text_line(3:22).eq.'linear_variable_kxky') then
            bstep = 2
         else
            stopSig = .True.
         end if
       ! from 'linear_variable_kxky', accepted values are 'linear_variable_{kxky,mean}'
       case (2)
         if (text_line(3:22).eq.'linear_variable_kxky') then
            bstep = 2
         else if (text_line(3:22).eq.'linear_variable_mean') then
            bstep = 3
         else
            stopSig = .True.
         end if
       ! from 'linear_variable_zero', accepted values are 'linear_variable_{mean,full}' 
       case (3)
         if (text_line(3:22).eq.'linear_variable_mean') then
            bstep = 3
         else if (text_line(3:22).eq.'linear_variable_full') then
            bstep = 4
         else
            stopSig = .True.
         end if
       ! from 'linear_variable_full', accepted values are 'linear_variable_full_build'
       case (4)
         if (text_line(3:28).eq.'linear_variable_full_build') then
            bstep = 5
         else
            stopSig = .True.
         end if
       ! from 'linear_variable_build', accepted values are 'linear_variable_{full_build/full}'
       !                                                   'quadratic_variable'
       case (5)
         if (text_line(3:28).eq.'linear_variable_full_build') then
            bstep = 5
         else if (text_line(3:22).eq.'linear_variable_full') then
            bstep = 4
         else if (text_line(3:20).eq.'quadratic_variable') then
            bstep = 6
         else
            stopSig = .True.
         end if
       ! from 'quadratic_variable', accepted values are 'quadratic_variable', 
       !                                                'add_set_of_coupled_kxky_equations'
       case (6)
         if (text_line(3:35).eq.'add_set_of_coupled_kxky_equations') then
            bstep = 7
         else if (text_line(3:20).eq.'quadratic_variable') then
            bstep = 6
         else
            stopSig = .True.
         end if
       ! from 'add_set_..kxky_eqns', accepted values are 'linearly_coupled_var'
       case (7)
         if (text_line(3:22).eq.'linearly_coupled_var') then
            bstep = 71
         else
            stopSig = .True.
         end if
       ! from 'linearly_coupled_var', accepted values are 'new_equation', 
       !                                              'add_BC_for_this_var'
       !                                              'linearly_coupled_var'
       case (71)
         if (text_line(3:22).eq.'linearly_coupled_var') then
            bstep = 71
         else if (text_line(3:14).eq.'new_equation') then
            bstep = 73
         else if (text_line(3:21).eq.'add_BC_for_this_var') then
            bstep = 72
         else
            stopSig = .True.
         end if
       ! from 'add_BC_for_this_var', accepted values are 'new_equation', 
       !                                              'linearly_coupled_var'
       case (72)
         if (text_line(3:22).eq.'linearly_coupled_var') then
            bstep = 71
         else if (text_line(3:14).eq.'new_equation') then
            bstep = 73
         else
            stopSig = .True.
         end if
       ! from 'new_equation', accepted values are 'add_d/dt_term'
       !                                          'add_rhs_linear'
       !                                          'add_rhs_NL'
       !                                          'new_equation'
       case (73)
         if (text_line(3:15).eq.'add_d/dt_term') then
            bstep = 731
         else if (text_line(3:16).eq.'add_rhs_linear') then
            bstep = 732
         else if (text_line(3:12).eq.'add_rhs_NL') then
            bstep = 733
         else
            stopSig = .True.
         end if
       ! from 'add_d/dt_term', accepted values are 'add_d/dt_term'
       !                                           'add_rhs_linear'
       !                                           'add_rhs_NL'
       !                                           'new_equation'
       case (731,732,733)
         if (text_line(3:15).eq.'add_d/dt_term') then
            bstep = 731
         else if (text_line(3:16).eq.'add_rhs_linear') then
            bstep = 732
         else if (text_line(3:12).eq.'add_rhs_NL') then
            bstep = 733
         else if (text_line(3:14).eq.'new_equation') then
            bstep = 73
         else if (text_line(3:35).eq.'add_set_of_coupled_kxky_equations') then
            bstep = 7
         else if (text_line(3:35).eq.'add_set_of_coupled_zero_equations') then
            bstep = 8 
         else
            stopSig = .True.
         end if
       ! from 'add_set_..zero_eqns', accepted values are 'linearly_coupled_var'
       case (8)
         if (text_line(3:22).eq.'linearly_coupled_var') then
            bstep = 81
         else
            stopSig = .True.
         end if
       ! from 'linearly_coupled_var', accepted values are 'new_equation', 
       !                                              'add_BC_for_this_var'
       !                                              'linearly_coupled_var'
       case (81)
         if (text_line(3:22).eq.'linearly_coupled_var') then
            bstep = 81
         else if (text_line(3:14).eq.'new_equation') then
            bstep = 83
         else if (text_line(3:21).eq.'add_BC_for_this_var') then
            bstep = 82
         else
            stopSig = .True.
         end if
       ! from 'add_BC_for_this_var', accepted values are 'new_equation', 
       !                                              'linearly_coupled_var'
       case (82)
         if (text_line(3:22).eq.'linearly_coupled_var') then
            bstep = 81
         else if (text_line(3:14).eq.'new_equation') then
            bstep = 83
         else
            stopSig = .True.
         end if
       ! from 'new_equation', accepted values are 'add_d/dt_term'
       !                                          'add_rhs_linear'
       !                                          'add_rhs_NL'
       !                                          'new_equation'
       case (83)
         if (text_line(3:15).eq.'add_d/dt_term') then
            bstep = 831
         else if (text_line(3:16).eq.'add_rhs_linear') then
            bstep = 832
         else if (text_line(3:12).eq.'add_rhs_NL') then
            bstep = 833
         else if (text_line(3:16).eq.'add_rhs_source') then
            bstep = 834
         else
            stopSig = .True.
         end if
       ! from 'add_d/dt_term', accepted values are 'add_d/dt_term'
       !                                           'add_rhs_linear'
       !                                           'add_rhs_NL'
       !                                           'new_equation'
       case (831,832,833,834)
         if (text_line(3:15).eq.'add_d/dt_term') then
            bstep = 831
         else if (text_line(3:16).eq.'add_rhs_linear') then
            bstep = 832
         else if (text_line(3:12).eq.'add_rhs_NL') then
            bstep = 833
         else if (text_line(3:16).eq.'add_rhs_source') then
            bstep = 834
         else if (text_line(3:14).eq.'new_equation') then
            bstep = 83
         else if (text_line(3:35).eq.'add_set_of_coupled_zero_equations') then
            bstep = 8 
         else if (text_line(3:5).eq.'EOF') then
            eofSig = .True.
         else
            stopSig = .True.
         end if
     end select
 
            
 end subroutine determine_kind_of_information

 End Module LP_equations
