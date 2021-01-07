module sparse_formats_d
 use fortran_kinds
 use chdir_mod
 implicit none

            !======================!
            !                      !
            !     CSR Matrices     !
            !                      !
            !======================!

 !> @brief
 !> CSR Matrix, double real.         
 !> @details
 !> We implement the 3-array variation of the complex sparse-row format.
 !! More specifically, the entries of the i-th row of the matrix are stored
 !! in \p dat between the indices \p row(i) and \p row(i+1)-1 . Their column
 !! is contained in the array \p col at the same positions. For efficiency,
 !! it is preferable that entries of each row are sorted by increasing 
 !! column. If this is the case, \p is_sorted is True.
 type :: dcsr_matrix
    real(dp), dimension(:), allocatable :: dat !< Double Precision Entries.
    integer, dimension(:), allocatable :: row  !< Row index.
    integer, dimension(:), allocatable :: col  !< Column.
    integer :: nrow !< Number of rows.
    integer :: ncol !< Number of columns.
    integer :: nelems !< number of non zero entries.
    logical :: is_sorted = .False. !< Are columns sorted by increasing order?
  contains 
    procedure :: write2disk => output_dcsr_matrix
    procedure :: dot_dcsr => wrap_dcsrmultcsr
    generic   :: dot => dot_dcsr, &
                 dot_dcsr_1d_singleVec_and_insert_in_3d_111, &
                 dot_dcsr_1d_1coupled_vec_oop_normal_general,&
                 dot_dcsr_1d_singleVec
    procedure :: copy => dcsr_copy
    procedure :: sort => dcsr_sort              
    procedure :: truncate => dcsr_truncate
    procedure :: truncate_rows => dcsr_truncate_rows
    procedure :: to_dia => d_csr2dia
    procedure :: empty => dcsr_empty
    procedure :: dot_dcsr_1d_singleVec
    procedure :: dot_dcsr_1d_singleVec_and_insert_in_3d_111
    procedure :: dot_dcsr_1d_1coupled_vec_oop_normal_general
    procedure :: identity => dcsr_create_identity
 end type dcsr_matrix

 !> @brief
 !> COO Matrix, double real.    
 !> @details
 !> Coordinate sparse matrix format: each non zero entry is defined by its
 !! position \p row(i), \p col(i) and its value \p dat(i).
 Type :: dcoo_matrix
    Real(dp), Dimension(:), Allocatable :: dat !< Double Precision Entries.
    Integer, Dimension(:), Allocatable :: row  !< Row index.
    Integer, Dimension(:), Allocatable :: col  !< Column.
    Integer :: nrow !< Number of rows.
    Integer :: ncol !< Number of columns.
    Integer :: nelems !< number of non zero entries.
  Contains
    Procedure :: copy => dcoo_copy
    Procedure :: transpose_out_of_place => dcoo_transpose_out_of_place
    Procedure :: transpose_in_place     => dcoo_transpose_in_place
    Generic   :: transpose => transpose_out_of_place, transpose_in_place
 End type dcoo_matrix

            !======================!
            !                      !
            !     DIA Matrices     !
            !                      !
            !======================!
 
 !> @brief
 !> DIA Matrix, double real.    
 !> @details
 !> See Intel MKL Sparse BLAS webpage.
 Type :: ddia_matrix
    Real(dp), Dimension(:,:), Allocatable :: dat!< Double real entries.
    Integer :: nrow !< Number of rows.
    Integer :: ncol !< Number of columns.
    Integer :: nl   !< Number of lower diagonals
    Integer :: nu   !< Number of upper diagonals
 End type ddia_matrix

 !> @brief
 !> DIA Matrix and its LU factorization, double precision.
 Type :: ddia_and_lu 
    Type(ddia_matrix) :: dia  !< Matrix in DIA format
    Real(dp), Dimension(:,:), Allocatable :: lu !< LU factorization
    Integer,  Dimension(:),   Allocatable :: piv!< pivots
 End Type ddia_and_lu

 type :: dLU_fact_T
   real(dp), allocatable :: factors(:,:)
   integer, allocatable :: pivots(:)
 end type dLU_fact_T

 contains

 subroutine dot_dcsr_1d_singleVec_and_insert_in_3d_111(self, x, C, o_or_c)
   class(dcsr_matrix), intent(in) :: self
   real(dp), allocatable, intent(in) :: x(:)
   real(dp), allocatable             :: y(:)
   complex(dp), allocatable, intent(inOut) :: C(:,:,:)
   character(len=5) :: o_or_c
   allocate( y(size(C, dim=1) ) )
   call self%dot(x, y, 1._dp, 0._dp)
   
   select case (o_or_c)
          case('overW')
              C(:,1,1) = cmplx(0._dp, 0._dp, kind=dp)
              C(:,1,1) = cmplx(y, 0._dp, kind=dp)
          case('cumul')
              C(:,1,1) = C(:,1,1) + cmplx(y, 0._dp, kind=dp)
          case default
              print *, 'invalid overwrite_or_cumul var'
              print *, o_or_c, ' was passed to subroutine dot_dcsr_1d_singleVec_and_insert_in_3d_111'
              stop
   end select
 end subroutine
   
 subroutine dot_dcsr_1d_singleVec(self, x, y, o_or_c)
   class(dcsr_matrix), intent(in) :: self
   real(dp), allocatable, intent(in) :: x(:)
   real(dp), allocatable, intent(inOut) :: y(:)
   character(len=5) :: o_or_c
   select case (o_or_c)
          case('overW')
              call self%dot(x, y, 1._dp, 0._dp)
          case('cumul')
              call self%dot(x, y, 1._dp, 1._dp)
          case default
              print *, 'invalid overwrite_or_cumul var'
              print *, o_or_c, ' was passed to subroutine dot_dcsr_1d_singleVec_and_insert_in_3d_111'
              stop
   end select
 end subroutine

 subroutine dot_dcsr_1d_1coupled_vec_oop_normal_general(self,  x, y, alpha, beta)
   class(dcsr_matrix), intent(in) :: self
   real(dp), allocatable, intent(in)    :: x(:)
   real(dp), allocatable, intent(inout) :: y(:)
   real(dp), intent(in) :: alpha 
   real(dp), intent(in) :: beta  
   
   call mkl_dcsrmv('N', self%nrow, self%ncol,    alpha,&
                   'GxxFxx', self%dat, &
                   self%col, self%row(1:self%nrow), self%row(2:(self%nrow+1)), &
                   x, beta, y)
 end subroutine

 subroutine dcsr_expExp(A, expA, B, expB, C)

 !
 ! computes C = A^(expA) + B^(expB)
 ! All matrices have to be *square*, in CSR format
 !
 ! Input:
 !   . A: Matrix A, CSR format
 !   . B: Matrix B, CSR format
 !   . expA, expB : exponents for matrices A and B, respectively. 
 !     It is assumed that expA > 0 (otherwise, use SparMat_EXPO)
 !     expB = 0 is not a problem.
 !
 ! Output:
 !   . C: Matrix C = A^(expA) * B^(expB)
 !
 
   type(dcsr_matrix), intent(in)  :: A, B
   type(dcsr_matrix), intent(out) :: C
   type(dcsr_matrix) :: aux1, aux2
   integer, intent(In) :: expA, expB
   integer :: i 

   ! some checks
   if (A%nRow .ne. A%nCol) print *, 'error in dcsr_expExp (shape of A)'
   if (A%nRow .ne. A%nCol) error stop
   if (B%nRow .ne. B%nCol) print *, 'error in dcsr_expExp (shape of B)'
   if (B%nRow .ne. B%nCol) error stop
   if (B%nRow .ne. A%nCol) print *, 'error in dcsr_expExp (shapes of A and B)'
   if (B%nRow .ne. A%nCol) error stop

   if (expA.eq.0) then
   call dcsr_create_identity(aux1, A%nrow)
   else
   
   call A%copy(aux1)
   do i = 1, ((expA-1)/2)
      call aux1%dot(A, aux2)
      call aux2%dot(A, aux1)
   end do

   if (expA.gt.1) then
   if (mod(expA,2).eq.0) then
      call aux1%dot(A, aux2)
      call aux2%copy(aux1)
   end if
   end if
   end if
   !................................................................
   ! now do the same for B
   !................................................................
   if (expB.ge.1) then
      call aux1%dot(B, aux2)
      call aux2%copy(aux1)
   end if
   do i = 1,((expB-1)/2)
      call aux1%dot(B, aux2)
      call aux2%dot(B, aux1)
   end do
   if ((expB.gt.1).and.(mod(expB,2).eq.0)) then
      call aux1%dot(B, aux2)
      call aux2%copy( aux1)
   end if

   call aux1%copy(C)

 end subroutine dcsr_ExpExp
 
 
 subroutine dcsr_empty(A,N)
   class(dcsr_matrix), Intent(Out) :: A
   Integer, Intent(In) :: N
   Allocate( A%dat(0))
   Allocate( A%col(0))
   Allocate( A%row(N+1))
   A%row = 1
   A%nrow = N
   A%ncol = N
   A%nelems = 0
   call A%sort()
 end subroutine dcsr_empty


 Subroutine wrap_dcsrmultcsr(self, other, CCsr)

   Class(dcsr_matrix), Intent(In) :: self 
   Type(dcsr_matrix), Intent(In) :: other
   Type(dcsr_matrix), Intent(Out):: CCsr 
   Integer :: info
   
   If (.not.(self%is_sorted))  print *, 'wrap_dcsrmultcsr unsrt. self mat.'         
   If (.not.(other%is_sorted)) print *, 'wrap_dcsrmultcsr unsrt. other mat.'          
   If (.not.(self%is_sorted.and.other%is_sorted)) error stop

   If (self%ncol.eq.other%nrow) then

   CCsr%nrow = self%nrow
   CCsr%ncol = other%ncol

   allocate( CCsr%row (CCsr%nrow+1))
   Allocate( CCsr%col(1), CCsr%dat(1)) ! Keeps the compiler happy.
   Call mkl_dcsrmultcsr('n', 1, 8, self%nrow, self%ncol, other%ncol,&
                                   self%dat,  self%col,  self%row,&
                                   other%dat,  other%col,  other%row,&
                                   CCsr%dat,  CCsr%col,  CCsr%row,&
                                   0, info)
   CCSR%nelems = CCsr%row (CCSr%Nrow+1) -1
   DeAllocate( CCsr%col, CCsr%dat)
   allocate( CCsr%col (CCsr%nelems))
   allocate( CCsr%dat (CCsr%nelems))
   Call mkl_dcsrmultcsr('n', 2, 8, self%nrow, self%ncol, other%ncol,&
                                   self%dat,  self%col,  self%row,&
                                   other%dat,  other%col,  other%row,&
                                   CCsr%dat,  CCsr%col,  CCsr%row,&
                                   0, info)
   Call CCSR%sort()
   else
      print *, 'error in wrap_dcsrmultcsr'
      print *, 'self%nCol :', self %nCol
      print *, 'other%nRow:', other%nRow
      print *, 'dCsr matrices dimensions not compatible'
      Stop
   end if
 End Subroutine wrap_dcsrmultcsr

 Subroutine dcsr_sort(self)

   class(dcsr_matrix), Intent(InOut) :: self
   Integer :: irow
   Integer :: imax, iel
   Logical :: are_we_finished
   Integer  :: left_col, right_col
   Real(dp) :: left_dat, right_dat
   
   Do irow = 1,self%Nrow
      imax = self%row(irow+1)-self%row(irow)
      if (imax.gt.1) then
      are_we_finished = .False.
      Do While (.not.are_we_finished)
         imax = imax - 1 
         are_we_finished = .true.
         Do iel = 1,imax       
            left_col  = self%col(self%row(irow)+iel-1)
            right_col = self%col(self%row(irow)+iel)
            if (left_col.gt.right_col) then
               self%col(self%row(irow)+iel-1) = right_col
               self%col(self%row(irow)+iel)   = left_col
               left_dat  = self%dat(self%row(irow)+iel-1)
               right_dat = self%dat(self%row(irow)+iel)
               self%dat(self%row(irow)+iel-1) = right_dat
               self%dat(self%row(irow)+iel)   = left_dat
               are_we_finished = .false.
            end if 
          end do
       end do
       end if
     end do
    
   self%is_sorted = .True.

 End Subroutine dcsr_sort


 Subroutine dcsr_copy (self,other)

   class(dcsr_matrix), Intent(In)  :: self  
   Type(dcsr_matrix), Intent(Out) :: other
 
   other%nelems = self%nelems
   other%nrow   = self%nrow
   other%ncol   = self%ncol
   Allocate(other%dat(other%nelems))
   Allocate(other%col(other%nelems))
   Allocate(other%row(other%nrow+1))   
   other%dat = self%dat
   other%col = self%col
   other%row = self%row
   Call other%sort()           

 end subroutine dcsr_copy
  
 subroutine dcsr_truncate_rows(self, other, iMin, iMax)
   class(dcsr_matrix), Intent(In)  :: self
   type(dcsr_matrix), Intent(Out) :: other
   integer, intent(in) :: iMin, iMax
   type(dcsr_matrix) :: aux
   integer :: i, ni
   !elementary checks:
   if (iMin.lt.1) then
     print *, 'error in dcsr_truncate_rows'
     print *, 'iMin has to be greater than 0'
     error stop
   end if
   if (iMax.gt.self%nRow) then
     print *, 'error in dcsr_truncate_rows'
     print *, 'iMax cannot be greater than self%nRow'
     error stop
   end if
   ni = iMax+1-iMin
   aux%nelems = ni
   aux%ncol = self%nrow
   aux%nrow = ni
   allocate (aux%dat(ni))
   allocate (aux%col(ni))
   allocate (aux%row(ni+1))
   do i = 1, ni
      aux%dat(i) = 1._dp
      aux%col(i) = iMin + i - 1
      aux%row(i) = i
   end do
   aux%row(ni+1) = ni+1
   ! multiply matrix A on the left to extract the lines
   call aux%sort()
   call aux%dot(self, other)
 end subroutine

 subroutine dcsr_truncate(self, other, imin, imax, jmin, jmax)
   
   class(dcsr_matrix), Intent(In)  :: self
   Type(dcsr_matrix), Intent(Out) :: other
   Type(dcsr_matrix) :: B, aux
   Integer, Intent(In) ::  imin, imax, jmin, jmax
   Integer :: nj
   Integer :: i

   ! >>> elementary checks: ---------------------------------------
   if (iMin.lt.1) then
     print *, 'error in dcsr_truncate_rows'
     print *, 'iMin has to be greater than 0'
     error stop
   end if
   if (iMax.gt.self%nRow) then
     print *, 'error in dcsr_truncate_rows'
     print *, 'iMax cannot be greater than self%nRow'
     error stop
   end if
   !elementary checks:
   if (jMin.lt.1) then
     print *, 'error in dcsr_truncate_rows'
     print *, 'jMin has to be greater than 0'
     error stop
   end if
   if (jMax.gt.self%nCol) then
     print *, 'error in dcsr_truncate_rows'
     print *, 'jMax cannot be greater than self%nCol'
     error stop
   end if
   ! /// elementary checks: ---------------------------------------

   if ((jMin.eq.1).and.(jMax.eq.self%nCol)) then
      call self%truncate_rows(other, iMin, iMax)
   else
      call self%truncate_rows(B,     iMin, iMax)
      ! now create a 'tall' auxilliary matrix
      nj = jmax+1-jmin
      aux%ncol = nj
      aux%nrow = B%ncol
      aux%nelems = nj
      allocate( aux%dat(nj))
      allocate( aux%col(nj))
      allocate( aux%row(B%Ncol+1))
      do i = 1,jmin
          aux%row(i) = 1
      end do
      do i = 1,nj
          aux%dat(i) = 1._dp
          aux%col(i) = i
          aux%row(jmin+i) = i+1
      end do
      do i= (jmax+2),(self%ncol+1)
          aux%row(i) = nj+1
      end do
      ! multiply matrix B on the right to extract the columns
      call aux%sort()         
      call B%dot(aux,other)
   end if
 End Subroutine dcsr_Truncate
  

 !> @brief
 !! Export a real CSR matrix on the disk.                       
 !> @details
 !! Matrices are exported in folder \p dir_str. 
 !! The entries of the matrix are saved as double precision 
 !! real (for numpy: \p np.float_ or \p np.float64) in the file 
 !! \p prefix_str.coreXXXX.dat, where \p XXXX is the number of
 !! the process, in a 4 digit format (left-padded with zeros if
 !! necessary). The file \p prefix_str.coreXXXX.indices contains
 !! 4-byte integers (\p np.int32) that represent:
 !! the number of elemens \p nelems [0]; the number of rows \p nrow [1];
 !! the number of columns \p ncol [2]; the row pointer (see CSR format def.) [3:3+nrow+1];
 !! the list of columns (see CSR format def.) [3+nrow+1:3+nrow+1+nelems] 
 Subroutine output_dcsr_matrix(self_CSR, prefix_str, prefix_len, my_id, &
                                         dir_str, dir_len) 
   Class(dcsr_matrix), Intent(In) :: self_CSR  !< the matrix we output
   Integer, Intent(In) :: dir_len                   !< length of variable dir_str
   Character(len=dir_len),    Intent(In) :: dir_str !< output directory
   Integer, Intent(In) :: prefix_len !< length of variable prefix_str
   Character(len=prefix_len), Intent(In) :: prefix_str  !< prefix for filename
   Integer, Intent(In) :: my_id   !< MPI rank
   Character(len=4) :: corestr    !< internal var.: MPI rank  string
   Character(len=100) :: file_str !< internal var.: full filename string
   Call chdir(dir_str)
   Call chdir('./Matrices')
   write(corestr,"(i4.4)") my_id 
   ! >>>> write the entries first. <<<<
   file_str = prefix_str//".core"//corestr//".dat"
   Open (Unit=9, File=file_str, Status='replace', Access='stream')
   Write(9) self_CSR%dat
   Close (Unit=9)
   ! >>>> write nelems, nrow, ncol, row, col <<<<<
   file_str = prefix_str//".core"//corestr//".indices"
   Open (Unit=9, File=file_str, Status='replace', Access='stream')
   Write(9) self_CSR%nelems
   Write(9) self_CSR%nrow  
   Write(9) self_CSR%ncol  
   Write(9) self_CSR%row  
   Write(9) self_CSR%col  
   Close (Unit=9)
 End Subroutine output_dcsr_matrix
 
 Subroutine d_csr2dia(aCSR, aDIA)
   class(dcsr_matrix), Intent(In) :: aCSR
   Type(ddia_matrix), Intent(Out):: aDIA
   Integer :: nl, nu
   Integer :: my_dia, my_col, my_row
   Integer :: jcol,  irow

   if (aCSR%ncol.ne.aCSR%nrow) print *, "d_csr2dia not meant for non-square matrices"
   if (aCSR%ncol.ne.aCSR%nrow) print *, "Stopping now."                              
   if (aCSR%ncol.ne.aCSR%nrow) STOP

   !determin nl and nu, the number of lower and upper diagonals
   nl = 0
   nu = 0
   Do irow = 1,aCSR%nrow
      nl = max(nl, irow- aCSR%col(aCSR%row(irow))          )
      nu = max(nu,       aCSR%col(aCSR%row(irow+1)-1) -irow)
   End Do
   aDIA%nl = nl
   aDIA%nu = nu
   aDIA%ncol = aCSR%ncol
   aDIA%nrow = aCSR%nrow
   Allocate(aDIA%dat(nl+nu+1,aCSR%ncol))
   aDIA%dat = 0._dp

   !read the CSR mat and fill-in
   Do irow = 1,aCSR%nrow
   Do jcol = aCSR%row(irow),aCSR%row(irow+1)-1
      my_row = irow
      my_col = aCSR%col(jcol)
      my_dia = my_col-my_row !value between nu and -nl
      aDIA%dat(nu-my_dia+1, my_dia+irow) = aCSR%dat(jcol)
   end do
   end do

 End Subroutine d_csr2dia


 Subroutine dcoo_copy (A,B)
   Class(dcoo_matrix), Intent(In)  :: A
   Type(dcoo_matrix), Intent(Out) :: B
   B%nelems = A%nelems
   B%nrow   = A%nrow
   B%ncol   = A%ncol
   Allocate(B%dat(B%nelems))
   Allocate(B%col(B%nelems))
   Allocate(B%row(B%nelems))   
   B%dat = A%dat
   B%col = A%col
   B%row = A%row
 End Subroutine dcoo_copy

 Subroutine dcoo_transpose_out_of_place (A,B)
   Class(dcoo_matrix), Intent(In)  :: A
   Type(dcoo_matrix), Intent(Out) :: B
   Call dcoo_copy(A, B)
   B%row = A%col
   B%col = A%row
   B%nCol = A%nRow
   B%nRow = A%nCol
 End Subroutine dcoo_transpose_out_of_place

 Subroutine dcoo_transpose_in_place (A)
   Class(dcoo_matrix), Intent(InOut)  :: A
   Type(dcoo_matrix) :: B
   integer :: bufInt
   Call dcoo_copy(A, B)
   B%row = A%col
   B%col = A%row
   A%col = B%col
   A%row = B%row
   bufInt = A%nCol
   A%nCol = A%nRow
   A%nRow = bufInt
 End Subroutine dcoo_transpose_in_place

 Subroutine dcsr_create_identity(A,N)
   class(dcsr_matrix), intent(out) :: A
   integer, intent(in) :: N
   integer :: i 
   allocate (A%dat(N))
   allocate (A%col(N))
   allocate (A%row(N+1))
   Do i=1,N
      A%dat(i) = 1._dp
      A%col(i) = i
      A%row(i) = i
   End do
   A%row(N+1) = N+1
   A%nelems = N
   A%ncol = N
   A%nrow = N
   Call dcsr_sort(A)
 end subroutine dcsr_create_identity
 
end module sparse_formats_d
