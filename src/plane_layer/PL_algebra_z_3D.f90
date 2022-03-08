module PL_algebra_z_3d

 use MPI_vars
 use fortran_kinds
 use sparse_formats
 use sparse_conversions
 use sparse_blas, only: dcsr_convert_zcsr
 use lapack95
 implicit none


 type :: zOperator_3d_1coupled_T
   type(zdia_matrix), allocatable :: dia(:,:)
   type(zLU_fact_T),  allocatable :: LU (:,:) 
   type(zcsr_matrix), allocatable :: csr(:,:)
   integer :: n1, n2, ncoupled
   logical :: has_csr, has_dia, has_lu
  contains
   procedure :: set_size => set_size_zOp_3d_1coupled
   procedure :: init_CSR => initializeCSR_zOp_3d_1coupled
   procedure :: alloc_DIA_and_LU => allocate_zOp_3d_and_LU
   procedure :: change_of_basis_CSR => change_of_basis_zCSR_dMat
   procedure :: copy_CSR_to_DIA => copy_zOp_csr_to_dia
   procedure :: factorize => factorize_zOp_3d_1coupled
   procedure :: backsolve_zOp_3d_1coupled
   procedure :: backsolve_zOp_3d_1coupled_knownLdb
   generic   :: backsolve => backsolve_zOp_3d_1coupled,&
                             backsolve_zOp_3d_1coupled_knownLDB
   procedure :: dot_zOp_3d_1coupled_vec_oop
   procedure :: dot_zOp_3d_1coupled_vec_oop_normal
   procedure :: dot_zOp_3d_1coupled_vec_oop_normal_to_pointer
   generic   :: dot => dot_zOp_3d_1coupled_vec_oop_normal, &
                       dot_zOp_3d_1coupled_vec_oop_normal_to_pointer, &
                       dot_zOp_3d_1coupled_vec_oop

 end type zOperator_3d_1coupled_T

contains

 subroutine change_of_basis_zCSR_dMat(self, Qmatrix, Pmatrix)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   type(dcsr_matrix), intent(In) :: Qmatrix, Pmatrix
   type(zcsr_matrix) :: zQmatrix, zPmatrix, aux
   integer :: i1, i2
   call dcsr_convert_zcsr(Qmatrix, zQmatrix)
   call dcsr_convert_zcsr(Pmatrix, zPmatrix)
   do i2 = 1,self%n2
   do i1 = 1,self%n1
      call self%csr(i1,i2)%dot(zPmatrix, aux)
      call zQmatrix%dot(aux, self%csr(i1,i2))
   end do
   end do
 end subroutine

 subroutine set_size_zOp_3D_1coupled(self, n1, n2, ncoupled)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   integer, intent(in) :: n1, n2, ncoupled
   self%n1 = n1
   self%n2 = n2
   self%ncoupled = ncoupled
 end subroutine

 subroutine copy_zOp_csr_to_dia(self)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   integer :: i1, i2
   allocate (self%dia(self%n1, self%n2))
   do i2 = 1, self%n2
   do i1 = 1, self%n1
     call self%csr( i1,i2 )%to_dia( self%dia(i1,i2) )
   end do
   end do
   self%has_dia = .True.
 end subroutine

 subroutine allocate_zOp_3d_and_LU(self)                   
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   allocate (self%dia(self%n1, self%n2))
   allocate (self%LU (self%n1, self%n2))
 end subroutine


 
 subroutine initializeCSR_zOp_3d_1coupled(self)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   integer :: i1, i2
   allocate (self%csr(self%n1, self%n2))
   do i2 = 1, self%n2                                                
   do i1 = 1, self%n1
      call self%csr(i1,i2)%empty( self%ncoupled )
   end do 
   end do 
   self%has_csr = .True.
 end subroutine

 subroutine dot_zOp_3d_1coupled_vec_oop_normal_to_pointer(self,  x, y)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable, intent(in)     :: x(:,:,:)
   complex(dp), pointer,     intent(inout)  :: y(:,:,:)
   complex(dp) :: alpha = cmplx(1._dp, 0._dp, kind=dp)
   complex(dp) :: beta  = cmplx(0._dp, 0._dp, kind=dp)
   complex(dp), allocatable :: columnVector_x(:)
   complex(dp), allocatable :: columnVector_y(:)
   integer :: i1, i2
   
   if (.not.self%has_csr) STOP 'CSR matrix missing' 
   allocate(columnVector_x(self%csr(1,1)%ncol))
   allocate(columnVector_y(self%csr(1,1)%ncol))
   do i2 = 1, self%n2
   do i1 = 1, self%n1
     columnVector_x = x(1:self%csr(i1,i2)%ncol, i1, i2)
     columnVector_y = cmplx(0._dp, 0._dp, kind=dp)
     call mkl_zcsrmv('N', &
                     self%csr(i1,i2)%nrow, &
                     self%csr(i1,i2)%ncol, alpha, 'GxxFxx', &
                     self%csr(i1,i2)%dat, &
                     self%csr(i1,i2)%col, &
                     self%csr(i1,i2)%row(1: self%csr(i1,i2)%nrow), &
                     self%csr(i1,i2)%row(2:(self%csr(i1,i2)%nrow+1)), &
                     columnVector_x, &
                     beta, &
                     columnVector_y)
     y(1:self%csr(i1,i2)%ncol, i1, i2) = columnVector_y 
   end do 
   end do 
 end subroutine

 subroutine dot_zOp_3d_1coupled_vec_oop_normal(self,  x, y)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable, intent(in)     :: x(:,:,:)
   complex(dp), allocatable, intent(inout)  :: y(:,:,:)
   complex(dp) :: alpha = cmplx(1._dp, 0._dp, kind=dp)
   complex(dp) :: beta  = cmplx(0._dp, 0._dp, kind=dp)
   complex(dp), allocatable :: columnVector_x(:)
   complex(dp), allocatable :: columnVector_y(:)
   integer :: i1, i2
   
   if (.not.self%has_csr) STOP 'CSR matrix missing' 
   allocate(columnVector_x(self%csr(1,1)%ncol))
   allocate(columnVector_y(self%csr(1,1)%ncol))
   do i2 = 1, self%n2
   do i1 = 1, self%n1
     columnVector_x = x(1:self%csr(i1,i2)%ncol, i1, i2)
     columnVector_y = cmplx(0._dp, 0._dp, kind=dp)
     call mkl_zcsrmv('N', &
                     self%csr(i1,i2)%nrow, &
                     self%csr(i1,i2)%ncol, alpha, 'GxxFxx', &
                     self%csr(i1,i2)%dat, &
                     self%csr(i1,i2)%col, &
                     self%csr(i1,i2)%row(1: self%csr(i1,i2)%nrow), &
                     self%csr(i1,i2)%row(2:(self%csr(i1,i2)%nrow+1)), &
                     columnVector_x, &
                     beta, &
                     columnVector_y)
     y(1:self%csr(i1,i2)%ncol, i1, i2) = columnVector_y 
   end do 
   end do 
 end subroutine

 subroutine dot_zOp_3d_1coupled_vec_oop(self, transpose_A, x, y)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   character(len=1), intent(in) :: transpose_A
   complex(dp), allocatable, intent(in), target     :: x(:,:,:)
   complex(dp), allocatable, intent(inout), target  :: y(:,:,:)
   complex(dp) :: alpha = cmplx(1._dp, 0._dp, kind=dp)
   complex(dp) :: beta  = cmplx(0._dp, 0._dp, kind=dp)
   complex(dp), pointer :: columnVector_x(:)=>null()
   complex(dp), pointer :: columnVector_y(:)=>null()
   integer :: i1, i2
   
   if (.not.self%has_csr) STOP 'CSR matrix missing' 
   do i2 = 1, self%n2
   do i1 = 1, self%n1
     columnVector_x => x(1:self%csr(i1,i2)%ncol, i1, i2)
     columnVector_y => y(1:self%csr(i1,i2)%ncol, i1, i2)
     call mkl_zcsrmv(transpose_A, &
                     self%csr(i1,i2)%nrow, &
                     self%csr(i1,i2)%ncol, alpha, 'GxxFxx', &
                     self%csr(i1,i2)%dat, &
                     self%csr(i1,i2)%col, &
                     self%csr(i1,i2)%row(1: self%csr(i1,i2)%nrow), &
                     self%csr(i1,i2)%row(2:(self%csr(i1,i2)%nrow+1)), &
                     columnVector_x, &
                     beta, &
                     columnVector_y)
   end do 
   end do 
 end subroutine

 subroutine backsolve_zOp_3d_1coupled_knownLdb(self, right_hand_sides)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable, intent(inOut) :: right_hand_sides(:,:,:)
   integer :: nrhs = 1
   integer :: i1
   integer :: i2
   integer :: info
   if (.not.self%has_lu) STOP 'missing LU factorization'
   do i2 = 1, self%n2
     do i1 = 1, self%n1
       !b => right_hand_sides(1:self%ncoupled, i1, i2) 
       call zgbtrs('N', self%dia(i1,i2)%ncol,  &   
                        self%dia(i1,i2)%nl,    &   
                        self%dia(i1,i2)%nu,    &   
                        nrhs,                  &   
                        self%lu(i1,i2)%factors,&  
                        2*self%dia(i1,i2)%nl + self%dia(i1,i2)%nu+1,&
                        self%lu (i1,i2)%pivots,&
                        right_hand_sides(:,i1,i2), &
                        self%dia(i1,i2)%ncol, info)
                        !b, ldb, info)
       if (info.ne.0) print *, "Problem during LU routine zgbtrs"
       if (info.ne.0) print *, "Error Diagnostic:", info
       if (info.ne.0) stop
     end do
   end do
 end subroutine


 subroutine backsolve_zOp_3d_1coupled(self, right_hand_sides, ldb)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   complex(dp), allocatable, target :: right_hand_sides(:,:,:)
   integer :: nrhs = 1
   integer, intent(in) :: ldb
   integer :: i1
   integer :: i2
   integer :: info
   if (.not.self%has_lu) STOP 'missing LU factorization'
   do i2 = 1, self%n2
     do i1 = 1, self%n1
       call zgbtrs('N', self%dia(i1,i2)%ncol,  &   
                        self%dia(i1,i2)%nl,    &   
                        self%dia(i1,i2)%nu,    &   
                        nrhs,                  &   
                        self%lu(i1,i2)%factors,&  
                        2*self%dia(i1,i2)%nl + self%dia(i1,i2)%nu+1,&
                        self%lu (i1,i2)%pivots,&
                        right_hand_sides(:,i1,i2), ldb, info)
                        !b, ldb, info)
       if (info.ne.0) print *, "Problem during LU routine zgbtrs, process", my_rank
       if (info.ne.0) print *, "Error Diagnostic:", info
       if (info.ne.0) stop
     end do
   end do
 end subroutine

 subroutine factorize_zOp_3d_1coupled(self)
   class(zOperator_3d_1coupled_T), intent(inOut) :: self
   integer :: i1
   integer :: i2
   integer :: info
   integer :: myNu, myNl, myNCol ! to improve legibility a little
   do i2 = 1, self%n2
     do i1 = 1, self%n1
        myNu   = self%dia(i1,i2)%nu
        myNl   = self%dia(i1,i2)%nl
        myNcol = self%dia(i1,i2)%nCol
        if (allocated(self%lu(i1,i2)%factors)) deAllocate (self%lu(i1,i2)%factors)
        allocate( self%lu(i1,i2)%factors( myNl*2 + myNu + 1 , myNcol ) ) 
        self%lu(i1,i2)%factors( myNl+1 : 2*myNl+myNu+1, :) = self%dia(i1,i2)%dat(:,:)

        if (allocated(self%lu(i1,i2)%pivots)) deAllocate (self%lu(i1,i2)%pivots)
        allocate(self%lu(i1,i2)%pivots( min(self%dia(i1,i2)%nrow, self%dia(i1,i2)%ncol)))
        call zgbtrf ( self%dia(i1,i2)%nrow,&
                      self%dia(i1,i2)%ncol,&
                      self%dia(i1,i2)%NL,  &
                      self%dia(i1,i2)%nu,  &
                      self%lu (i1,i2)%factors,&
                      2*self%dia(i1,i2)%nl + self%dia(i1,i2)%nu+1,&
                      self%lu (i1,i2)%pivots,&
                      info)
       if (info.ne.0) print *, "Problem during LU routine zgbtrf, (process, info, ix, iy)", my_rank, info, i1, i2
       if (info.ne.0) stop
     end do
   end do
   self%has_lu = .True.
 end subroutine

end module PL_algebra_z_3d
