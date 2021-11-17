module LP_algebra_d_3D

 use fortran_kinds
 use sparse_formats
 use sparse_conversions
 use lapack95
 implicit none


 type :: dOperator_3d_1coupled_T
   type(ddia_matrix), allocatable :: dia(:,:)
   type(dLU_fact_T), allocatable  :: LU (:,:) 
   integer :: n1, n2, ncoupled
   logical :: has_csr, has_dia, has_lu
  contains
   procedure :: factorize => factorize_dOp_DIA_3d_1coupled
   procedure :: backsolve_dOp_DIA_3d_1coupled
   procedure :: backsolve_dOp_DIA_3d_1coupled_knownLdb
   generic   :: backsolve => backsolve_dOp_DIA_3d_1coupled,&
                             backsolve_dOp_DIA_3d_1coupled_knownLDB
 end type dOperator_3d_1coupled_T

contains

 subroutine backsolve_dOp_DIA_3d_1coupled_knownLdb(self, right_hand_sides)
   class(dOperator_3d_1coupled_T), intent(inOut) :: self
   real(dp), allocatable, target :: right_hand_sides(:,:,:)
   integer :: nrhs = 1
   integer :: i1
   integer :: i2
   integer :: info
   do i2 = 1, self%n2
     do i1 = 1, self%n1
       call dgbtrs('N', self%dia(i1,i2)%ncol,  &   
                        self%dia(i1,i2)%nl,    &   
                        self%dia(i1,i2)%nu,    &   
                        nrhs,                  &   
                        self%lu(i1,i2)%factors,&  
                        2*self%dia(i1,i2)%nl + self%dia(i1,i2)%nu+1,&
                        self%lu (i1,i2)%pivots,&
                        right_hand_sides(:,i1,i2), &
                        self%dia(i1,i2)%ncol, info)
                        !b, ldb, info)
        if (info.ne.0) call xgbtrx_error_msg(info,.True.)
     end do
   end do
 end subroutine

 
 subroutine backsolve_dOp_DIA_3d_1coupled(self, right_hand_sides, ldb)
   class(dOperator_3d_1coupled_T), intent(inOut) :: self
   real(dp), allocatable, target :: right_hand_sides(:,:,:)
   integer :: nrhs = 1
   integer, intent(in) :: ldb
   integer :: i1
   integer :: i2
   integer :: info
   do i2 = 1, self%n2
     do i1 = 1, self%n1
       call dgbtrs('N', self%dia(i1,i2)%ncol,  &   
                        self%dia(i1,i2)%nl,    &   
                        self%dia(i1,i2)%nu,    &   
                        nrhs,                  &   
                        self%lu(i1,i2)%factors,&  
                        2*self%dia(i1,i2)%nl + self%dia(i1,i2)%nu+1,&
                        self%lu (i1,i2)%pivots,&
                        right_hand_sides(:,i1,i2), ldb, info)
        if (info.ne.0) call xgbtrx_error_msg(info,.True.)
     end do
   end do
 end subroutine

 subroutine factorize_dOp_DIA_3d_1coupled(self)
   class(dOperator_3d_1coupled_T), intent(inOut) :: self
   integer :: i1
   integer :: i2
   integer :: info
   do i2 = 1, self%n2
     do i1 = 1, self%n1
        if (allocated(self%lu(i1,i2)%pivots)) deAllocate (self%lu(i1,i2)%pivots)
        allocate(self%lu(i1,i2)%pivots( min(self%dia(i1,i2)%nrow, self%dia(i1,i2)%ncol)))
        call dgbtrf ( self%dia(i1,i2)%nrow,&
                      self%dia(i1,i2)%ncol,&
                      self%dia(i1,i2)%NL,  &
                      self%dia(i1,i2)%nu,  &
                      self%lu (i1,i2)%factors,&
                      2*self%dia(i1,i2)%nl + self%dia(i1,i2)%nu+1,&
                      self%lu (i1,i2)%pivots,&
                      info)
        if (info.ne.0) call xgbtrx_error_msg(info,.True.)
     end do
   end do
 end subroutine

end module LP_algebra_d_3D
