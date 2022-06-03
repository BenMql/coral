module R_algebra_d_1D

 use fortran_kinds
 use sparse_formats
 use sparse_conversions
 use lapack95
 implicit none

 type :: dOperator_1d_1coupled_T
   type(ddia_matrix)  :: dia
   type(dcsr_matrix)  :: csr
   type(dLU_fact_T)   :: LU 
   logical :: has_dia, has_csr, has_lu
   integer :: nCoupled
  contains
   procedure :: set_size => set_size_dOp_DIA_1d_1coupled
   procedure :: init_CSR => initializeCSR_dOp_DIA_1d_1coupled
   procedure :: change_of_basis_CSR => change_of_basis_dCSR_dMat
   procedure :: copy_CSR_to_DIA => copy_dOp_csr_to_dia
   procedure :: build_dOp_DIA_1d_1coupled_from_csr
   procedure :: build_dOp_DIA_1d_1coupled_from_coo
   procedure :: factorize => factorize_dOp_DIA_1d_1coupled
   procedure :: dot_dOp_DIA_1d_1coupled_vec_oop
   procedure :: dot_dOp_DIA_1d_1coupled_vec_oop_normal
   generic   :: dot => dot_dOp_DIA_1d_1coupled_vec_oop_normal, &
                       dot_dOp_DIA_1d_1coupled_vec_oop
   procedure :: backsolve_dOp_DIA_1d_1coupled
   procedure :: backsolve_dOp_DIA_1d_1coupled_knownLdb
   procedure :: backsolve_dOp_DIA_1d_1coupled_unique
   procedure :: backsolve_dOp_DIA_1d_1coupled_PTR
   generic   :: backsolve => backsolve_dOp_DIA_1d_1coupled,&
                             backsolve_dOp_DIA_1d_1coupled_unique,&
                             backsolve_dOp_DIA_1d_1coupled_knownLdb, &
                             backsolve_dOp_DIA_1d_1coupled_PTR
   generic   :: constructor => build_dOp_DIA_1d_1coupled_from_csr, &
                               build_dOp_DIA_1d_1coupled_from_coo
 end type dOperator_1d_1coupled_T

 contains

 subroutine change_of_basis_dCSR_dMat(self, Qmatrix, Pmatrix)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   type(dcsr_matrix), intent(In) :: Qmatrix, Pmatrix
   type(dcsr_matrix) :: aux
   call self%csr%dot(Pmatrix, aux)
   call Qmatrix%dot(aux, self%csr)
 end subroutine

 subroutine set_size_dOp_DIA_1D_1coupled(self,ncoupled)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   integer, intent(in) :: ncoupled
   self%ncoupled = ncoupled
 end subroutine


 subroutine copy_dOp_csr_to_dia(self)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   call self%csr%to_dia( self%dia )
   self%has_csr = .True.
 end subroutine


 subroutine initializeCSR_dOp_DIA_1d_1coupled(self)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   call self%csr%empty(self%nCoupled)
   self%has_csr = .True.
 end subroutine



 subroutine build_dOp_DIA_1d_1coupled_from_coo(self, aCOO)
   type(dcoo_matrix), intent(in) :: aCOO
   class(dOperator_1d_1coupled_T) :: self
   call d_coo2csr(aCOO,self%CSR)            
   call self%csr%to_dia(self%dia)
   self%has_csr=.True.
   self%has_dia=.True.
   self%has_lu =.False.
 end subroutine

 subroutine build_dOp_DIA_1d_1coupled_from_csr(self, aCSR)
   type(dcsr_matrix), intent(in) :: aCSR
   class(dOperator_1d_1coupled_T) :: self
   call dcsr_copy(aCSR, self%csr)
   call self%CSR%to_dia(self%dia)
   self%has_csr=.True.
   self%has_dia=.True.
   self%has_lu =.False.
 end subroutine

 subroutine dot_dOp_DIA_1d_1coupled_vec_oop_normal(self,  x, y)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   real(dp), allocatable, intent(in)    :: x(:)
   real(dp), allocatable, intent(inout) :: y(:)
   real(dp) :: alpha = 1._dp
   real(dp) :: beta  = 0._dp
   
   if (.not.self%has_csr) STOP 'CSR matrix missing' 
   call mkl_dcsrmv('N', self%csr%nrow, self%csr%ncol,    alpha,&
                   'GxxFxx', self%csr%dat, &
                   self%csr%col, self%csr%row(1:self%csr%nrow), self%csr%row(2:(self%csr%nrow+1)), &
                   x, beta, y)
 end subroutine

 subroutine dot_dOp_DIA_1d_1coupled_vec_oop(self, transpose_A, x, y)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   character(len=1), intent(in) :: transpose_A
   real(dp), allocatable, intent(in)    :: x(:)
   real(dp), allocatable, intent(inout) :: y(:)
   real(dp) :: alpha = 1._dp
   real(dp) :: beta  = 0._dp
   
   if (.not.self%has_csr) STOP 'CSR matrix missing' 
   call mkl_dcsrmv(transpose_A, self%csr%nrow, self%csr%ncol,    alpha,&
                   'GxxFxx', self%csr%dat, &
                   self%csr%col, self%csr%row(1:self%csr%nrow), self%csr%row(2:(self%csr%nrow+1)), &
                   x, beta, y)
 end subroutine

 subroutine backsolve_dOp_DIA_1d_1coupled(self, right_hand_sides, nrhs, ldb)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   real(dp), allocatable :: right_hand_sides(:)
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldb
   integer :: info
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call dgbtrs('N', self%dia%ncol,  &   
                    self%dia%nl,    &   
                    self%dia%nu,    &   
                    nrhs,           &   
                    self%lu%factors,&  
                    2*self%dia%nl + self%dia%nu+1,&
                    self%lu %pivots,&
                    right_hand_sides, ldb, info)
    
   if (info.ne.0) print *, "Problem during LU routine dgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 end subroutine

 subroutine backsolve_dOp_DIA_1d_1coupled_knownLDB(self, right_hand_sides, nrhs)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   real(dp), allocatable :: right_hand_sides(:)
   integer, intent(in) :: nrhs
   integer :: info
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call dgbtrs('N', self%dia%ncol,  &   
                    self%dia%nl,    &   
                    self%dia%nu,    &   
                    nrhs,           &   
                    self%lu%factors,&  
                    2*self%dia%nl + self%dia%nu+1,&
                    self%lu %pivots,&
                    right_hand_sides, &
                    self%dia%ncol, info)
    
   if (info.ne.0) print *, "Problem during LU routine dgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 end subroutine

 subroutine backsolve_dOp_DIA_1d_1coupled_PTR(self, right_hand_sides, nrhs, ldb)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   real(dp), pointer :: right_hand_sides(:)
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldb
   integer :: info
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call dgbtrs('N', self%dia%ncol,  &   
                    self%dia%nl,    &   
                    self%dia%nu,    &   
                    nrhs,           &   
                    self%lu%factors,&  
                    2*self%dia%nl + self%dia%nu+1,&
                    self%lu %pivots,&
                    right_hand_sides, ldb, info)
    
   if (info.ne.0) print *, "Problem during LU routine dgbtrs"
   if (info.ne.0) print *, "in backsolve_zOp_DIA_1d_1coupled_PTR"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 end subroutine



 subroutine backsolve_dOp_DIA_1d_1coupled_unique(self, right_hand_sides)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   real(dp), allocatable, intent(inOut) :: right_hand_sides(:)
   integer :: nrhs = 1
   integer :: info
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call dgbtrs('N', self%dia%ncol,  &   
                    self%dia%nl,    &   
                    self%dia%nu,    &   
                    nrhs,           &   
                    self%lu%factors,&  
                    2*self%dia%nl + self%dia%nu+1,&
                    self%lu %pivots,&
                    right_hand_sides, &
                    self%dia%ncol, info)
   if (info.ne.0) print *, "Problem during LU routine dgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 end subroutine

 subroutine factorize_dOp_DIA_1d_1coupled(self)
   class(dOperator_1d_1coupled_T), intent(inOut) :: self
   integer :: info
   if (allocated(self%lu%factors)) deAllocate (self%lu%factors)
   allocate(self%lu%factors( self%dia%nl*2 + self%dia%nu + 1, self%dia%ncol))
   self%lu%factors(self%dia%nl+1:self%dia%nl*2 + self%dia%nu + 1, :) = self%dia%dat

   if (allocated(self%lu%pivots))  deAllocate (self%lu%pivots)
   allocate(self%lu%pivots( min(self%dia%nrow, self%dia%ncol)))
   call dgbtrf ( self%dia%nrow,&
                 self%dia%ncol,&
                 self%dia%nl,  &
                 self%dia%nu,  &
                 self%lu %factors,&
                 2*self%dia%nl + self%dia%nu+1,&
                 self%lu %pivots,&
                 info)
   if (info.ne.0) print *, "Problem during LU routine dgbtrf"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
   self%has_lu = .True.
 end subroutine

end module R_algebra_d_1D
