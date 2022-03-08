module PL_algebra_z_1D

 use fortran_kinds
 use sparse_formats
 use sparse_conversions
 use lapack95
 implicit none

 type :: zOperator_1d_1coupled_T
   type(zdia_matrix)  :: dia
   type(zcsr_matrix)  :: csr
   type(zLU_fact_T)   :: LU 
   logical :: has_dia, has_csr, has_lu
   integer :: nCoupled
  contains
   procedure :: set_size => set_size_zOp_DIA_1d_1coupled
   procedure :: init_CSR => initializeCSR_zOp_DIA_1d_1coupled
   procedure :: copy_CSR_to_DIA => copy_zOp_csr_to_dia
   procedure :: build_zOp_DIA_1d_1coupled_from_csr
   procedure :: build_zOp_DIA_1d_1coupled_from_coo
   procedure :: factorize => factorize_zOp_DIA_1d_1coupled
   procedure :: dot_zOp_DIA_1d_1coupled_vec_oop
   procedure :: dot_zOp_DIA_1d_1coupled_vec_oop_normal
   generic   :: dot => dot_zOp_DIA_1d_1coupled_vec_oop_normal, &
                       dot_zOp_DIA_1d_1coupled_vec_oop
   procedure :: backsolve_zOp_DIA_1d_1coupled
   procedure :: backsolve_zOp_DIA_1d_1coupled_PTR
   procedure :: backsolve_zOp_DIA_1d_1coupled_3dPTR
   procedure :: backsolve_zOp_DIA_1d_1coupled_3drhs
   procedure :: backsolve_zOp_DIA_1d_1coupled_knownLdb
   procedure :: backsolve_zOp_DIA_1d_1coupled_unique
   generic   :: backsolve => backsolve_zOp_DIA_1d_1coupled,&
                             backsolve_zOp_DIA_1d_1coupled_PTR,&
                             backsolve_zOp_DIA_1d_1coupled_3dPTR,&
                             backsolve_zOp_DIA_1d_1coupled_3drhs,&
                             backsolve_zOp_DIA_1d_1coupled_unique,&
                             backsolve_zOp_DIA_1d_1coupled_knownLdb
   generic   :: constructor => build_zOp_DIA_1d_1coupled_from_csr, &
                               build_zOp_DIA_1d_1coupled_from_coo
 end type zOperator_1d_1coupled_T

 contains


 subroutine set_size_zOp_DIA_1D_1coupled(self,ncoupled)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   integer, intent(in) :: ncoupled
   self%ncoupled = ncoupled
 end subroutine


 subroutine copy_zOp_csr_to_dia(self)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   call self%csr%to_dia( self%dia )
   self%has_csr = .True.
 end subroutine


 subroutine initializeCSR_zOp_DIA_1d_1coupled(self)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   call self%csr%empty(self%nCoupled)
   self%has_csr = .True.
 end subroutine



 subroutine build_zOp_DIA_1d_1coupled_from_coo(self, aCOO)
   type(zcoo_matrix), intent(in) :: aCOO
   class(zOperator_1d_1coupled_T) :: self
   call z_coo2csr(aCOO,self%CSR)            
   call self%csr%to_dia(self%dia)
   self%has_csr=.True.
   self%has_dia=.True.
   self%has_lu =.False.
 end subroutine

 subroutine build_zOp_DIA_1d_1coupled_from_csr(self, aCSR)
   type(zcsr_matrix), intent(in) :: aCSR
   class(zOperator_1d_1coupled_T) :: self
   call zcsr_copy(aCSR, self%csr)
   call self%CSR%to_dia(self%dia)
   self%has_csr=.True.
   self%has_dia=.True.
   self%has_lu =.False.
 end subroutine

 subroutine dot_zOp_DIA_1d_1coupled_vec_oop_normal(self,  x, y)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable, intent(in)    :: x(:)
   complex(dp), allocatable, intent(inout) :: y(:)
   complex(dp) :: alpha = cmplx(1._dp, 0._dp, kind=dp)
   complex(dp) :: beta  = cmplx(0._dp, 0._dp, kind=dp)
   
   if (.not.self%has_csr) STOP 'CSR matrix missing' 
   call mkl_zcsrmv('N', self%csr%nrow, self%csr%ncol,    alpha,&
                   'GxxFxx', self%csr%dat, &
                   self%csr%col, self%csr%row(1:self%csr%nrow), self%csr%row(2:(self%csr%nrow+1)), &
                   x, beta, y)
 end subroutine

 subroutine dot_zOp_DIA_1d_1coupled_vec_oop(self, transpose_A, x, y)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   character(len=1), intent(in) :: transpose_A
   complex(dp), allocatable, intent(in)    :: x(:)
   complex(dp), allocatable, intent(inout) :: y(:)
   complex(dp) :: alpha = cmplx(1._dp, 0._dp, kind=dp)
   complex(dp) :: beta  = cmplx(0._dp, 0._dp, kind=dp)
   
   if (.not.self%has_csr) STOP 'CSR matrix missing' 
   call mkl_zcsrmv(transpose_A, self%csr%nrow, self%csr%ncol,    alpha,&
                   'GxxFxx', self%csr%dat, &
                   self%csr%col, self%csr%row(1:self%csr%nrow), self%csr%row(2:(self%csr%nrow+1)), &
                   x, beta, y)
 end subroutine

 subroutine backsolve_zOp_DIA_1d_1coupled_3drhs(self, right_hand_sides_3d, nrhs, ldb)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable, intent(inOut), target :: right_hand_sides_3d(:,:,:)
   complex(dp), pointer :: right_hand_sides(:)
   type(c_ptr) :: dummptr
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldb
   integer :: info
   dummptr = c_loc(right_hand_sides_3d(1,1,1))
   call c_f_pointer(dummPtr, right_hand_sides, [nrhs*ldb])
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call zgbtrs('N', self%dia%ncol,  &   
                    self%dia%nl,    &   
                    self%dia%nu,    &   
                    nrhs,           &   
                    self%lu%factors,&  
                    2*self%dia%nl + self%dia%nu+1,&
                    self%lu %pivots,&
                    right_hand_sides, ldb, info)
    
   if (info.ne.0) print *, "Problem during LU routine zgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 end subroutine

 subroutine backsolve_zOp_DIA_1d_1coupled_3dPTR(self, right_hand_sides3d, nrhs, ldb)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   complex(dp), pointer :: right_hand_sides3d(:,:,:)
   complex(dp), pointer :: right_hand_sides  (:)
   type(C_ptr) :: dummptr                
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldb
   integer :: info
   dummptr = c_loc(right_hand_sides3d(1,1,1))
   call c_f_pointer(dummPtr, right_hand_sides, [nrhs*ldb])
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call zgbtrs('N', self%dia%ncol,  &   
                    self%dia%nl,    &   
                    self%dia%nu,    &   
                    nrhs,           &   
                    self%lu%factors,&  
                    2*self%dia%nl + self%dia%nu+1,&
                    self%lu %pivots,&
                    right_hand_sides, ldb, info)
    
   if (info.ne.0) print *, "Problem during LU routine zgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 end subroutine

 subroutine backsolve_zOp_DIA_1d_1coupled_PTR(self, right_hand_sides, nrhs, ldb)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   complex(dp), pointer :: right_hand_sides(:)
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldb
   integer :: info
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call zgbtrs('N', self%dia%ncol,  &   
                    self%dia%nl,    &   
                    self%dia%nu,    &   
                    nrhs,           &   
                    self%lu%factors,&  
                    2*self%dia%nl + self%dia%nu+1,&
                    self%lu %pivots,&
                    right_hand_sides, ldb, info)
    
   if (info.ne.0) print *, "Problem during LU routine zgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 end subroutine


 subroutine backsolve_zOp_DIA_1d_1coupled(self, right_hand_sides, nrhs, ldb)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable :: right_hand_sides(:)
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldb
   integer :: info
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call zgbtrs('N', self%dia%ncol,  &   
                    self%dia%nl,    &   
                    self%dia%nu,    &   
                    nrhs,           &   
                    self%lu%factors,&  
                    2*self%dia%nl + self%dia%nu+1,&
                    self%lu %pivots,&
                    right_hand_sides, ldb, info)
    
   if (info.ne.0) print *, "Problem during LU routine zgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 end subroutine

 subroutine backsolve_zOp_DIA_1d_1coupled_knownLDB(self, right_hand_sides, nrhs)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable :: right_hand_sides(:)
   integer, intent(in) :: nrhs
   integer :: info
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call zgbtrs('N', self%dia%ncol,  &   
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

 subroutine backsolve_zOp_DIA_1d_1coupled_unique(self, right_hand_sides)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable, intent(inOut) :: right_hand_sides(:)
   integer :: nrhs = 1
   integer :: info
   if (.not.self%has_lu) STOP 'LU factorization is missing'
   call zgbtrs('N', self%dia%ncol,  &   
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

 subroutine factorize_zOp_DIA_1d_1coupled(self)
   class(zOperator_1d_1coupled_T), intent(inOut) :: self
   integer :: info
   if (allocated(self%lu%factors)) deAllocate (self%lu%factors)
   allocate(self%lu%factors( self%dia%nl*2 + self%dia%nu + 1, self%dia%ncol))
   self%lu%factors(self%dia%nl+1:self%dia%nl*2 + self%dia%nu + 1, :) = self%dia%dat

   if (allocated(self%lu%pivots))  deAllocate (self%lu%pivots)
   allocate(self%lu%pivots( min(self%dia%nrow, self%dia%ncol)))
   call zgbtrf ( self%dia%nrow,&
                 self%dia%ncol,&
                 self%dia%nl,  &
                 self%dia%nu,  &
                 self%lu %factors,&
                 2*self%dia%nl + self%dia%nu+1,&
                 self%lu %pivots,&
                 info)
   if (info.ne.0) print *, "Problem during LU routine zgbtrf"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
   self%has_lu = .True.
 end subroutine

end module PL_algebra_z_1d
