 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: Sparse_formats
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> Definition of derived type for storing sparse matrices in the common
 !! CSR, COO, and DIA format for both double real and double complex precision.
 !
 !=============================================================================

 module sparse_formats_z
 use fortran_kinds, Only: dp
 use chdir_mod
 use, Intrinsic :: ISO_C_BINDING
 implicit none

            !======================!
            !                      !
            !     CSR Matrices     !
            !                      !
            !======================!
 
 !> @brief
 !> CSR Matrix, double complex.
 !> @details
 !> We implement the 3-array variation of the complex sparse-row format.
 !! More specifically, the entries of the i-th row of the matrix are stored
 !! in \p dat between the indices \p row(i) and \p row(i+1)-1 . Their column
 !! is contained in the array \p col at the same positions. For efficiency,
 !! it is preferable that entries of each row are sorted by increasing 
 !! column. If this is the case, \p is_sorted is True.
 type :: zcsr_matrix
    complex(dp), dimension(:), allocatable :: dat !< Double complex entries.
    integer, dimension(:), allocatable :: row     !< Row index.
    integer, dimension(:), allocatable :: col     !< Column.
    integer :: nrow !< Number of rows.
    integer :: ncol !< Number of columns.
    integer :: nelems !< number of non zero entries.
    logical :: is_sorted = .False. !< Are columns sorted by increasing order?
  contains 
    procedure :: write2disk => output_zcsr_matrix
    procedure :: truncate => zcsr_truncate
    procedure :: to_dia => z_csr2dia
    procedure :: dot_zcsr => wrap_zcsrmultcsr
    procedure :: Tdot_zcsr => wrap_zcsrmultcsr_T
    procedure :: empty => zcsr_empty
    procedure :: identity => zcsr_create_identity
    procedure :: sort => zcsr_sort              
    procedure :: dot_zCSR_1d_manyCoupled_vec_rescale
    procedure :: dot_zCSR_1d_manyCoupled_vec
    procedure :: dot_zCSR_1d_manyCoupled_vec_general
    procedure :: dot_zCSR_1d_manyCoupled_vec_rescale_T
    procedure :: dot_zCSR_1d_manyCoupled_vec_T
    !procedure :: dot_zCSR_1d_manyCoupled_vec_transpose_general
    generic :: dot => dot_zCSR_1d_manyCoupled_vec_rescale, &
                      dot_zCSR_1d_manyCoupled_vec, &
                      dot_zCSR_1d_manyCoupled_vec_general, &
                      dot_zCSR_1d_manyCoupled_vec_rescale_T, &
                      dot_zCSR_1d_manyCoupled_vec_T, &
                      dot_zcsr
    generic :: T_dot => Tdot_zcsr
 end type zcsr_matrix

            !======================!
            !                      !
            !     COO Matrices     !
            !                      !
            !======================!


 !> @brief
 !> COO Matrix, double complex.
 !> @details
 !> Coordinate sparse matrix format: each non zero entry is defined by its
 !! position \p row(i), \p col(i) and its value \p dat(i).
 Type :: zcoo_matrix
    Complex(dp), Dimension(:), Allocatable :: dat !< Double complex entries.
    Integer, Dimension(:), Allocatable :: row     !< Row index.
    Integer, Dimension(:), Allocatable :: col     !< Column.
    Integer :: nrow !< Number of rows.
    Integer :: ncol !< Number of columns.
    Integer :: nelems !< number of non zero entries.
  Contains
    Procedure :: copy => zcoo_copy
    Procedure :: transpose_out_of_place => zcoo_transpose_out_of_place
    Procedure :: transpose_in_place     => zcoo_transpose_in_place
    Generic   :: transpose => transpose_out_of_place, transpose_in_place
 End type zcoo_matrix
            !======================!
            !                      !
            !     DIA Matrices     !
            !                      !
            !======================!
 
 !> @brief
 !> DIA Matrix, double complex. 
 !> @details
 !> See Intel MKL Sparse BLAS webpage.
 Type :: zdia_matrix
    Complex(dp), Dimension(:,:), Allocatable :: dat !< Double complex entries.
    Integer :: nrow !< Number of rows.
    Integer :: ncol !< Number of columns.
    Integer :: nl   !< Number of lower diagonals
    Integer :: nu   !< Number of upper diagonals
 End type zdia_matrix

 !> @brief
 !> DIA Matrix and its LU factorization, double complex.
 Type :: zdia_and_lu 
    Type(zdia_matrix) :: dia  !< Matrix in DIA format
    Complex(dp), Dimension(:,:), Allocatable :: lu !< LU factorization
    Integer,     Dimension(:),   Allocatable :: piv!< pivots
 End Type zdia_and_lu


 type :: zLU_fact_T
   complex(dp), allocatable :: factors(:,:)
   integer, allocatable :: pivots(:)
 end type zLU_fact_T

 Contains 

 subroutine dot_zCSR_1d_manyCoupled_vec_rescale_T(self, transpose_flag, B, C, o_or_c, zAlpha)
   class(zcsr_matrix), intent(in) :: self
   character(len=1), intent(in) :: transpose_flag
   Complex(Kind=dp), allocatable,  intent(inOut) :: B(:,:,:)                                        
   Complex(Kind=dp), allocatable,  intent(inOut) :: C(:,:,:)                                        
   character (len=5), intent(in) :: o_or_c
   complex(kind=dp), intent(in) :: zAlpha
   complex(kind=dp) :: zBeta
   if (transpose_flag.ne.'T') print *, 'bad transpose flag'
   if (transpose_flag.ne.'T') print *, transpose_flag, 'was passed to dot_zCSR_1d_manyCoupled_vec_rescale_T'
   if (transpose_flag.ne.'T') stop 
 
   
   select case (o_or_c)
          case ('overW')
               zbeta = cmplx(0._dp, 0._dp, kind=dp)
          case ('cumul')
               zbeta = cmplx(1._dp, 0._dp, kind=dp)
          case default
               print *, 'bad overwrite_or_cumul flag in dot_zOp_1d_manyCoupled_vec'
               print *, 'passed value is: ', o_or_c
               stop
   end select
   call self%dot_zcsr_1d_manyCoupled_vec_general(&
                      'T', &
                      zAlpha,&                          
                      B,    & 
                      zbeta,&
                      C,&
                      size(B, dim=2),  &
                      size(B, dim=3), 1)
                      
 end subroutine
 

 subroutine dot_zCSR_1d_manyCoupled_vec_T(self, transpose_flag, B, C, o_or_c)
   class(zcsr_matrix), intent(in) :: self
   character(len=1), intent(in) :: transpose_flag
   Complex(Kind=dp), allocatable,  intent(inOut) :: B(:,:,:)                                        
   Complex(Kind=dp), allocatable,  intent(inOut) :: C(:,:,:)                                        
   character (len=5) :: o_or_c
   complex(kind=dp) :: zBeta
   if (transpose_flag.ne.'T') print *, 'bad transpose flag'
   if (transpose_flag.ne.'T') print *, transpose_flag, 'was passed to dot_zCSR_1d_manyCoupled_vec_rescale_T'
   if (transpose_flag.ne.'T') stop 
   select case (o_or_c)
          case ('overW')
               zbeta = cmplx(0._dp, 0._dp, kind=dp)
          case ('cumul')
               zbeta = cmplx(1._dp, 0._dp, kind=dp)
          case default
               print *, 'bad overwrite_or_cumul flag in dot_zOp_1d_manyCoupled_vec'
               print *, 'passed value is: ', o_or_c
               stop
   end select
   call self%dot_zcsr_1d_manyCoupled_vec_general(&
                      'T', &
                      cmplx(1._dp, 0._dp, kind=dp),&
                      B,& 
                      zbeta,&
                      C,&
                      size(B, dim=2),              &
                      size(B, dim=3), 1)
                      
 end subroutine
 
 !> @brief
 !> \callgraph Apply a sparse operator \p A on multiple right-hand sides \p B and
 !> store the result in \p C
 subroutine dot_zcsr_1d_manyCoupled_vec_general(self_csr, transpose_A, &
                  zalpha, B, zbeta, C, N2, N3, sta3)
   !................................................................
   ! Intended use:
   ! Apply an operator (discretized by csr_matrix A) along the first
   ! dimension of a 3D array B. Update C with the result. Visually:
   ! Do i2=1,N2 and i3=sta3,N3
   ! C(:,i2,i3) = zalpha * A(:,j) * B(j,i2,i3) + zbeta * C(:,i2,i3)
   !................................................................
   class(zcsr_matrix), intent(in) :: self_csr
   Complex(Kind=dp), allocatable, intent(inOut), target :: B(:,:,:) 
   Complex(Kind=dp), allocatable, intent(inOut), target :: C(:,:,:) 
   Complex(Kind=dp), Intent(in) :: zalpha, zbeta
   Integer, Intent(In) :: N2, N3
   Integer, Intent(In) :: sta3
   Character(len=1), intent(in) :: transpose_A 
   type(C_Ptr) :: dumPtr
   Complex(Kind=dp), Pointer :: ptr2B(:,:)
   Complex(Kind=dp), Pointer :: ptr2C(:,:)
   Character :: matdescra(4)
   matdescra = ['G','X', 'N', 'F']
   if (sta3.le.N3) then
   dumPtr = C_Loc(B(1,1,sta3))
   Call C_F_Pointer (dumPtr, ptr2B, [self_csr%nCol, N2*(N3+1-sta3)] )
   dumPtr = C_Loc(C(1,1,sta3))
   Call C_F_Pointer (dumPtr, ptr2C, [self_csr%nCol, N2*(N3+1-sta3)] )
   Call mkl_zcsrmm(transpose_A,&
                   self_csr%nrow, N2*(N3+1-sta3), self_csr%ncol, zalpha, matdescra,&
                   self_csr%dat, self_csr%col, &
                   self_csr%row(1: self_csr%nrow   ), &
                   self_csr%row(2:(self_csr%nrow+1)), &
                   ptr2B, self_csr%ncol, zbeta, &
                   ptr2C, self_csr%nrow)
   end if
 end subroutine 


 subroutine dot_zCSR_1d_manyCoupled_vec_rescale(self, B, C, o_or_c, zAlpha)
   class(zcsr_matrix), intent(in) :: self
   Complex(Kind=dp), allocatable,  intent(inOut) :: B(:,:,:)                                        
   Complex(Kind=dp), allocatable,  intent(inOut) :: C(:,:,:)                                        
   character (len=5), intent(in) :: o_or_c
   complex(kind=dp), intent(in) :: zAlpha
   complex(kind=dp) :: zBeta
   select case (o_or_c)
          case ('overW')
               zbeta = cmplx(0._dp, 0._dp, kind=dp)
          case ('cumul')
               zbeta = cmplx(1._dp, 0._dp, kind=dp)
          case default
               print *, 'bad overwrite_or_cumul flag in dot_zOp_1d_manyCoupled_vec'
               print *, 'passed value is: ', o_or_c
               stop
   end select
   call self%dot_zcsr_1d_manyCoupled_vec_general(&
                      'N', &
                      zAlpha,&                          
                      B,    & 
                      zbeta,&
                      C,&
                      size(B, dim=2),  &
                      size(B, dim=3), 1)
                      
 end subroutine
 

 subroutine dot_zCSR_1d_manyCoupled_vec(self, B, C, o_or_c)
   class(zcsr_matrix), intent(in) :: self
   Complex(Kind=dp), allocatable,  intent(inOut) :: B(:,:,:)                                        
   Complex(Kind=dp), allocatable,  intent(inOut) :: C(:,:,:)                                        
   character (len=5) :: o_or_c
   complex(kind=dp) :: zBeta
   select case (o_or_c)
          case ('overW')
               zbeta = cmplx(0._dp, 0._dp, kind=dp)
          case ('cumul')
               zbeta = cmplx(1._dp, 0._dp, kind=dp)
          case default
               print *, 'bad overwrite_or_cumul flag in dot_zOp_1d_manyCoupled_vec'
               print *, 'passed value is: ', o_or_c
               stop
   end select
   if (size(B, dim=2) .ne. size(C, dim=2)) then
      print *, 'dot_zCSR_1d_manyCoupled_vec, uncompatible sizes along dim.2'
      print *, 'C:', size(C, dim=2),'; B:', size(B, dim=2)
      stop
   end if
   if (size(B, dim=3) .ne. size(C, dim=3)) then
      print *, 'dot_zCSR_1d_manyCoupled_vec, uncompatible sizes along dim.3'
      print *, 'C:', size(C, dim=3),'; B:', size(B, dim=3)
      stop
   end if
   call self%dot_zcsr_1d_manyCoupled_vec_general(&
                      'N', &
                      cmplx(1._dp, 0._dp, kind=dp),&
                      B,& 
                      zbeta,&
                      C,&
                      size(B, dim=2),              &
                      size(B, dim=3), 1)
                      
 end subroutine
 
 !!> @brief
 !!> \callgraph Apply a sparse operator \p A on multiple right-hand sides \p B and
 !!> store the result in \p C
 !subroutine dot_zcsr_1d_manyCoupled_vec_normal_general(self_csr,&
                  !zalpha, B, zbeta, C, N2, N3, sta3)
   !!................................................................
   !! Intended use:
   !! Apply an operator (discretized by csr_matrix A) along the first
   !! dimension of a 3D array B. Update C with the result. Visually:
   !! Do i2=1,N2 and i3=sta3,N3
   !! C(:,i2,i3) = zalpha * A(:,j) * B(j,i2,i3) + zbeta * C(:,i2,i3)
   !!................................................................
   !class(zcsr_matrix), intent(in) :: self_csr
   !Complex(Kind=dp), allocatable, intent(inOut), target :: B(:,:,:) 
   !Complex(Kind=dp), allocatable, intent(inOut), target :: C(:,:,:) 
   !Complex(Kind=dp), Intent(in) :: zalpha, zbeta
   !Integer, Intent(In) :: N2, N3
   !Integer, Intent(In) :: sta3
   !Character(len=1) :: transpose_A = 'N'
   !type(C_Ptr) :: dumPtr
   !Complex(Kind=dp), Pointer :: ptr2B(:,:)
   !Complex(Kind=dp), Pointer :: ptr2C(:,:)
   !Character :: matdescra(4)
   !matdescra = ['G','X', 'N', 'F']
   !if (sta3.le.N3) then
   !dumPtr = C_Loc(B(1,1,sta3))
   !Call C_F_Pointer (dumPtr, ptr2B, [self_csr%nCol, N2*(N3+1-sta3)] )
   !dumPtr = C_Loc(C(1,1,sta3))
   !Call C_F_Pointer (dumPtr, ptr2C, [self_csr%nCol, N2*(N3+1-sta3)] )
   !Call mkl_zcsrmm(transpose_A,&
                   !self_csr%nrow, N2*(N3+1-sta3), self_csr%ncol, zalpha, matdescra,&
                   !self_csr%dat, self_csr%col, &
                   !self_csr%row(1: self_csr%nrow   ), &
                   !self_csr%row(2:(self_csr%nrow+1)), &
                   !ptr2B, self_csr%ncol, zbeta, &
                   !ptr2C, self_csr%nrow)
   !end if
 !end subroutine 

 subroutine zcsr_empty(A,N)
   class(zcsr_matrix), Intent(Out) :: A
   Integer, Intent(In) :: N
   Allocate( A%dat(0))
   Allocate( A%col(0))
   Allocate( A%row(N+1))
   A%row = 1
   A%nrow = N
   A%ncol = N
   A%nelems = 0
   call A%sort()
 end subroutine zcsr_empty

 Subroutine zcsr_create_identity(A,N)
   class(zcsr_matrix), intent(out) :: A
   integer, intent(in) :: N
   integer :: i 
   allocate (A%dat(N))
   allocate (A%col(N))
   allocate (A%row(N+1))
   Do i=1,N
      A%dat(i) = cmplx(1._dp, 0._dp, kind=dp)
      A%col(i) = i
      A%row(i) = i
   End do
   A%row(N+1) = N+1
   A%nelems = N
   A%ncol = N
   A%nrow = N
   Call zcsr_sort(A)
 end subroutine zcsr_create_identity
 
 !> @brief
 !! Export a complex CSR matrix on the disk.                       
 !> @details
 !! Matrices are exported in folder \p dir_str. 
 !! The entries of the matrix are saved as double precision 
 !! complex (for numpy: \p np.complex_ or \p np.complex128) in the file 
 !! \p prefix_str.coreXXXX.dat, where \p XXXX is the number of
 !! the process, in a 4 digit format (left-padded with zeros if
 !! necessary). The file \p prefix_str.coreXXXX.indices contains
 !! 4-byte integers (\p np.int32) that represent:
 !! the number of elemens \p nelems [0]; the number of rows \p nrow [1];
 !! the number of columns \p ncol [2]; the row pointer (see CSR format def.) [3:3+nrow+1];
 !! the list of columns (see CSR format def.) [3+nrow+1:3+nrow+1+nelems] 
 Subroutine output_zcsr_matrix(self_CSR, prefix_str, prefix_len, my_id, &
                                         dir_str, dir_len) 
   Class(zcsr_matrix), Intent(In) :: self_CSR  !< the matrix we output
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
 End Subroutine output_zcsr_matrix

 Subroutine zcoo_copy (A,B)

   Class(zcoo_matrix), Intent(In)  :: A
   Type (zcoo_matrix), Intent(Out) :: B
   B%nelems = A%nelems
   B%nrow   = A%nrow
   B%ncol   = A%ncol
   Allocate(B%dat(B%nelems))
   Allocate(B%col(B%nelems))
   Allocate(B%row(B%nelems))   
   B%dat = A%dat
   B%col = A%col
   B%row = A%row

 End Subroutine zcoo_copy

 Subroutine zcoo_transpose_out_of_place (A,B)

   Class(zcoo_matrix), Intent(In)  :: A
   Type (zcoo_matrix), Intent(Out) :: B
   Call zcoo_copy(A, B)
   B%row = A%col
   B%col = A%row
 End Subroutine zcoo_transpose_out_of_place

 Subroutine zcoo_transpose_in_place (A)

   Class(zcoo_matrix), Intent(InOut)  :: A
   Type (zcoo_matrix) :: B
   Call zcoo_copy(A, B)
   B%row = A%col
   B%col = A%row
   A%col = B%col
   A%row = B%row
 End Subroutine zcoo_transpose_in_place

 Subroutine wrap_zcsrmultcsr_T(self, other, CCsr)

   class(zcsr_matrix), Intent(In) :: self 
   Type(zcsr_matrix), Intent(In) :: other 
   Type(zcsr_matrix), Intent(Out):: CCsr 
   Integer :: info
   If (.not.(self%is_sorted.and.other%is_sorted)) print *, 'wrap_zcsrmultcsr unsrt. matrices'
   If (.not.(self%is_sorted.and.other%is_sorted)) STOP

   CCsr%nrow = self%nrow
   CCsr%ncol = other%ncol

   allocate( CCsr%row (CCsr%nrow+1))
   Allocate( CCsr%col(1), CCsr%dat(1)) ! Keeps the compiler happy.
   Call mkl_zcsrmultcsr('T', 1, 8, self%nrow, self%ncol, other%ncol,&
                                   self%dat,  self%col,  self%row,&
                                   other%dat,  other%col,  other%row,&
                                   CCsr%dat,  CCsr%col,  CCsr%row,&
                                   0, info)
   CCSR%nelems = CCsr%row (CCSr%Nrow+1) -1
   DeAllocate( CCsr%col, CCsr%dat)
   allocate( CCsr%col (CCsr%nelems))
   allocate( CCsr%dat (CCsr%nelems))
   Call mkl_zcsrmultcsr('T', 2, 8, self%nrow, self%ncol, other%ncol,&
                                   self%dat,  self%col,  self%row,&
                                   other%dat,  other%col,  other%row,&
                                   CCsr%dat,  CCsr%col,  CCsr%row,&
                                   0, info)
   Call zcsr_sort(CCSR)
 end subroutine wrap_zcsrmultcsr_T

 subroutine wrap_zcsrmultcsr(self, other, CCsr)

   class(zcsr_matrix), intent(in) :: self 
   type(zcsr_matrix), intent(in) :: other 
   type(zcsr_matrix), intent(out):: CCsr 
   integer :: info
   if (.not.(self%is_sorted.and.other%is_sorted)) print *, 'wrap_zcsrmultcsr unsrt. matrices'
   if (.not.(self%is_sorted.and.other%is_sorted)) STOP

   CCsr%nrow = self%nrow
   CCsr%ncol = other%ncol

   allocate( CCsr%row (CCsr%nrow+1))
   allocate( CCsr%col(1), CCsr%dat(1)) ! Keeps the compiler happy.
   call mkl_zcsrmultcsr('n', 1, 8, self%nrow, self%ncol, other%ncol,&
                                   self%dat,  self%col,  self%row,&
                                   other%dat,  other%col,  other%row,&
                                   CCsr%dat,  CCsr%col,  CCsr%row,&
                                   0, info)
   CCSR%nelems = CCsr%row (CCSr%Nrow+1) -1
   deAllocate( CCsr%col, CCsr%dat)
   allocate( CCsr%col (CCsr%nelems))
   allocate( CCsr%dat (CCsr%nelems))
   call mkl_zcsrmultcsr('n', 2, 8, self%nrow, self%ncol, other%ncol,&
                                   self%dat,  self%col,  self%row,&
                                   other%dat,  other%col,  other%row,&
                                   CCsr%dat,  CCsr%col,  CCsr%row,&
                                   0, info)
   call zcsr_sort(CCSR)
 end subroutine wrap_zcsrmultcsr


 subroutine zcsr_sort(self)

   class(zcsr_matrix), Intent(InOut) :: self
   Integer :: irow
   Integer :: imax, iel
   Logical :: are_we_finished
   Integer  :: left_col, right_col
   Complex(dp) :: left_dat, right_dat
   
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

 End Subroutine zcsr_sort



 Subroutine zcsr_truncate(self, other, imin, imax, jmin, jmax)
   
   class(zcsr_matrix), Intent(In)  :: self
   Type(zcsr_matrix), Intent(Out) :: other
   Type(zcsr_matrix) :: B, aux
   Integer, Intent(In) ::  imin, imax, jmin, jmax
   Integer :: ni, nj
   Integer :: i
   
   ! create a 'flat' auxilliary matrix to extract the relevant lines 
   ni = imax+1-imin
   allocate( aux%dat(ni))
   allocate( aux%col(ni))
   allocate( aux%row(ni+1))
   do i = 1,ni
       aux%dat(i) = Cmplx(1._dp, 0._dp, kind=dp)
       aux%col(i) = imin + i - 1
       aux%row(i) = i
   end do
   aux%row(ni+1) = ni+1
   aux%nrow = ni
   aux%ncol = self%nrow
   aux%nelems = ni
   ! multiply matrix A on the left to extract the lines
   Call wrap_zcsrmultcsr(aux, self, B)
   ! now create a 'tall' auxilliary matrix
   nj = jmax+0-jmin
   deallocate(aux%dat, aux%col, aux%row)
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
       aux%dat(i) = Cmplx(1._dp, 0._dp, kind=dp)
       aux%col(i) = i
       aux%row(jmin+i) = i+1
   end do
   do i= (jmax+2),(self%ncol+1)
       aux%row(i) = nj+1
   end do
   ! multiply matrix B on the right to extract the columns
   Call zcsr_sort(aux)
   Call wrap_zcsrmultcsr(B,aux, other)
 End Subroutine zcsr_Truncate
  
 !

 !
 !
 !
 !
 !
 !
 Subroutine z_csr2dia(aCSR, aDIA)
   class(zcsr_matrix), Intent(In) :: aCSR
   Type(zdia_matrix), Intent(Out):: aDIA
   Integer :: nl, nu
   Integer :: my_dia, my_col, my_row
   Integer :: jcol,  irow

   if (aCSR%ncol.ne.aCSR%nrow) print *, "z_csr2dia not meant for non-square matrices"
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
   aDIA%dat = Cmplx(0._dp, 0._dp, Kind=dp)

   !read the CSR mat and fill-in
   Do irow = 1,aCSR%nrow
   Do jcol = aCSR%row(irow),aCSR%row(irow+1)-1
      my_row = irow
      my_col = aCSR%col(jcol)
      my_dia = my_col-my_row !value between nu and -nl
      aDIA%dat(nu-my_dia+1, my_dia+irow) = aCSR%dat(jcol)
   end do
   end do

 End Subroutine z_csr2dia

 End Module sparse_formats_z

