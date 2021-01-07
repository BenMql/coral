 Module sparse_blas
   Use fortran_kinds
   Use sparse_formats
   Use sparse_conversions

 Implicit None

  interface dcsr_convert_zcsr
     module procedure dcsr_convert_zcsr_rescale
     module procedure dcsr_convert_zcsr_noRescaling
  end interface

 Contains
 
 ! ================================================ !
 !                                                  !
 !         DOUBLE PRECISION ROUTINES                !
 !                                                  !
 ! ================================================ !

 Subroutine wrap_dcsrmultcsr_T(ACsr, BCsr, CCsr)

   Type(dcsr_matrix), Intent(In) :: ACsr 
   Type(dcsr_matrix), Intent(In) :: BCsr 
   Type(dcsr_matrix), Intent(Out):: CCsr 
   Integer :: info
   If (.not.(ACsr%is_sorted.and.BCSR%is_sorted)) print *, 'wrap_dcsrmultcsr unsrt. matrices'
   If (.not.(ACsr%is_sorted.and.BCSR%is_sorted)) STOP

   CCsr%nrow = ACsr%nrow
   CCsr%ncol = BCsr%ncol

   allocate( CCsr%row (CCsr%nrow+1))
   Allocate( CCsr%col(1), CCsr%dat(1)) ! Keeps the compiler happy.
   Call mkl_dcsrmultcsr('T', 1, 8, ACsr%nrow, ACsr%ncol, BCSR%ncol,&
                                   ACsr%dat,  ACsr%col,  ACsr%row,&
                                   BCsr%dat,  BCsr%col,  BCsr%row,&
                                   CCsr%dat,  CCsr%col,  CCsr%row,&
                                   0, info)
   CCSR%nelems = CCsr%row (CCSr%Nrow+1) -1
   DeAllocate( CCsr%col, CCsr%dat)
   allocate( CCsr%col (CCsr%nelems))
   allocate( CCsr%dat (CCsr%nelems))
   Call mkl_dcsrmultcsr('T', 2, 8, ACsr%nrow, ACsr%ncol, BCSR%ncol,&
                                   ACsr%dat,  ACsr%col,  ACsr%row,&
                                   BCsr%dat,  BCsr%col,  BCsr%row,&
                                   CCsr%dat,  CCsr%col,  CCsr%row,&
                                   0, info)
   Call dcsr_sort(CCSR)
 End Subroutine wrap_dcsrmultcsr_T
 
 Subroutine wrap_dcsradd(ACsr, BCsr, CCsr, sca)

   Type(dcsr_matrix), Intent(In)  :: ACsr 
   Type(dcsr_matrix), Intent(In)  :: BCsr 
   Type(dcsr_matrix), Intent(Out) :: CCsr 
   Real(dp), Intent(In) :: sca
   Integer :: info

   If (.not.(ACsr%is_sorted.and.BCSR%is_sorted)) print *, 'wrap_dcsradd unsrt. matrices'
   If (.not.(ACsr%is_sorted.and.BCSR%is_sorted)) STOP

   CCsr%nrow = ACsr%nrow
   CCsr%ncol = ACsr%ncol

   If ((ACsr%nrow.ne.BCsr%nrow).or.(Acsr%ncol.ne.BCsr%Ncol)) print *, 'inconsistent dims.'
   If ((ACsr%nrow.ne.BCsr%nrow).or.(Acsr%ncol.ne.BCsr%Ncol)) STOP

   allocate( CCsr%row (CCsr%nrow+1))
   Allocate( CCsr%col(1), CCsr%dat(1)) ! Keeps the compiler happy.
   Call mkl_dcsradd('N', 1, 3, ACsr%nrow, ACsr%ncol,&
                               ACsr%dat,  ACsr%col,  ACsr%row, sca,&
                               BCsr%dat,  BCsr%col,  BCsr%row, &
                               CCsr%dat,  CCsr%col,  CCsr%row, 0, info)
   CCsr%nelems = CCsr%row(CCsr%nrow+1)-1
   DeAllocate(CCsr%dat, CCsr%col)
   Allocate (CCsr%dat(CCsr%nelems))
   Allocate (CCsr%col(CCsr%nelems))
   Call mkl_dcsradd('N', 2, 3, ACsr%nrow, ACsr%ncol,&
                               ACsr%dat,  ACsr%col,  ACsr%row, sca,&
                               BCsr%dat,  BCsr%col,  BCsr%row, &
                               CCsr%dat,  CCsr%col,  CCsr%row, 0, info)
   Call dcsr_sort(CCSR)
 End Subroutine wrap_dcsradd

  
   
  
 Subroutine dcsr_convert_zcsr_noRescaling (A, B)
   Type(dcsr_matrix), Intent(In)  :: A
   Type(zcsr_matrix), Intent(Out) :: B
   Integer :: j
   B%nelems = A%nelems
   B%nrow   = A%nrow
   B%ncol   = A%ncol
   Allocate(B%dat(B%nelems))
   Allocate(B%col(B%nelems))
   Allocate(B%row(B%nrow+1))   
   do j=1,B%nelems
   B%dat(j) = Cmplx(A%dat(j), 0._dp, Kind=dp)
   B%col(j) = A%col(j)
   end do
   do j=1,(B%nrow +1)
   B%row(j) = A%row(j)
   end do
   Call zcsr_sort(B)
 End Subroutine dcsr_convert_zcsr_noRescaling
   
  
  
 Subroutine dcsr_convert_zcsr_rescale (A, B, sca)
   Type(dcsr_matrix), Intent(In)  :: A
   Type(zcsr_matrix), Intent(Out) :: B
   Complex(kind=dp), Intent(In)   :: sca
   Integer :: j
   B%nelems = A%nelems
   B%nrow   = A%nrow
   B%ncol   = A%ncol
   Allocate(B%dat(B%nelems))
   Allocate(B%col(B%nelems))
   Allocate(B%row(B%nrow+1))   
   do j=1,B%nelems
   B%dat(j) = Cmplx(A%dat(j), 0._dp, Kind=dp) * sca
   B%col(j) = A%col(j)
   end do
   do j=1,(B%nrow +1)
   B%row(j) = A%row(j)
   end do
   Call zcsr_sort(B)
 End Subroutine dcsr_convert_zcsr_rescale
   
   
 
 ! ================================================ !
 !                                                  !
 !           DOUBLE COMPLEX ROUTINES                !
 !                                                  !
 ! ================================================ !



 Subroutine wrap_zcsradd(ACsr, BCsr, CCsr, sca)

   Type(zcsr_matrix), Intent(In) :: ACsr 
   Type(zcsr_matrix), Intent(In) :: BCsr 
   Type(zcsr_matrix), Intent(Out):: CCsr 
   Complex(dp), Intent(In) :: sca
   Integer :: info
   If (.not.(ACsr%is_sorted)) print *, 'wrap_zcsradd unsrt. matrix A'
   If (.not.(BCsr%is_sorted)) print *, 'wrap_zcsradd unsrt. matrix B'
   If (.not.(ACsr%is_sorted.and.BCSR%is_sorted)) STOP

   CCsr%nrow = ACsr%nrow
   CCsr%ncol = ACsr%ncol

   If ((ACsr%nrow.ne.BCsr%nrow).or.(Acsr%ncol.ne.BCsr%Ncol)) print *, 'inconsistent dims.'
   If ((ACsr%nrow.ne.BCsr%nrow).or.(Acsr%ncol.ne.BCsr%Ncol)) STOP

   allocate( CCsr%row (CCsr%nrow+1))
   Allocate( CCsr%col(1), CCsr%dat(1)) ! Keeps the compiler happy.
   Call mkl_zcsradd('N', 1, 3, ACsr%nrow, ACsr%ncol,&
                               ACsr%dat,  ACsr%col,  ACsr%row, sca,&
                               BCsr%dat,  BCsr%col,  BCsr%row, &
                               CCsr%dat,  CCsr%col,  CCsr%row, 0, info)
   CCsr%nelems = CCsr%row(CCsr%nrow+1)-1
   DeAllocate(CCsr%dat, CCsr%col)
   Allocate (CCsr%dat(CCsr%nelems))
   Allocate (CCsr%col(CCsr%nelems))
   Call mkl_zcsradd('N', 2, 3, ACsr%nrow, ACsr%ncol,&
                               ACsr%dat,  ACsr%col,  ACsr%row, sca,&
                               BCsr%dat,  BCsr%col,  BCsr%row, &
                               CCsr%dat,  CCsr%col,  CCsr%row, 0, info)
   Call zcsr_sort(CCsr)
 End Subroutine wrap_zcsradd

 End Module sparse_blas

