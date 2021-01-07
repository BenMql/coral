 Module sparse_conversions
 
 Use fortran_kinds
 Use sparse_formats

 Implicit None
 
 ! WARNING: ============================================================
 ! WARNING: Most of these routines *demand* that matrices passed as 
 ! WARNING: arguments in the CSR format are sorted (i.e. for each
 ! WARNING: row, the column indices are increasing). To be safe, always
 ! WARNING: sort CSR matrices *right after* their creation, by invoking
 ! WARNING: dcsr_sort(..) or zcsr_sort(..).
 ! WARNING: ============================================================

 Contains

 ! ================================================ !
 !                                                  !
 !         DOUBLE PRECISION ROUTINES                !
 !                                                  !
 ! ================================================ !


 Subroutine d_csr2coo(aCSR, aCOO)

   Type(dcsr_matrix), Intent(In) :: aCSR
   Type(dcoo_matrix), Intent(Out):: aCOO
   Integer :: irow

   ! WARNING: the CSR matrix *has to* be sorted
   
   aCOO%nelems = aCSR%nelems
   aCOO%nrow   = aCSR%nrow
   aCOO%ncol   = aCSR%ncol
   Allocate(aCOO%dat(aCSR%nelems))
   Allocate(aCOO%col(aCSR%nelems))
   Allocate(aCOO%row(aCSR%nelems))
   aCOO%dat = aCSR%dat      
   aCOO%col = aCSR%col      

   Do irow = 1, aCOO%nrow
      aCOO%row(aCSR%row(irow):(aCSR%row(irow+1)-1)) = irow
   End do

 End Subroutine d_csr2coo


 Subroutine d_coo2csr(aCOO, aCSR)

   Type(dcoo_matrix), Intent(In) :: aCOO
   Type(dcsr_matrix), Intent(Out):: aCSR
   Integer :: irow
   Integer :: fill_in, iel, this_row_length
   
   Allocate(aCSR%row(aCOO%nrow+1))
   Allocate(aCSR%col(aCOO%nelems))
   Allocate(aCSR%dat(aCOO%nelems))
   aCSR%nrow = aCOO%nrow
   aCSR%ncol = aCOO%ncol
   aCSR%nelems = aCOO%nelems

   aCSR%row(1) = 1
   fill_in = 0
   Do irow = 1, aCOO%nrow
      this_row_length = 0 
      Do iel = 1, aCOO%nelems
         if (aCOO%row(iel).eq.irow) then
            this_row_length = this_row_length + 1 
            fill_in = fill_in + 1 
            aCSR%dat(fill_in) = aCOO%dat(iel)
            aCSR%col(fill_in) = aCOO%col(iel)
         end if
         aCSR%row(irow+1) = aCSR%row(irow) + this_row_length 
      end Do
   End Do  

   Call dcsr_sort(aCSR)

   ! a routine check (better be safe...)
   if (aCSR%nelems.ne.(aCSR%row(aCSR%nrow+1)-1)) print *, "Pb conversion coo2csr"
   if (aCSR%nelems.ne.(aCSR%row(aCSR%nrow+1)-1)) STOP


 End Subroutine d_coo2csr

 Subroutine d_coo2dense(aCOO, aFull)
   Type(dcoo_matrix), Intent(In) :: aCOO
   Real(dp), Dimension(:,:), Allocatable, intent(Out) :: aFull
   Integer :: iel

   Allocate(aFull(aCOO%nrow, aCOO%ncol))
   aFull = 0._dp
   Do iel = 1, aCOO%nelems
      aFull(aCOO%row(iel), aCOO%col(iel)) = aCOO%dat(iel)
   End Do
 End Subroutine d_coo2dense 

 Subroutine d_csr2dense(aCSR, aFull)

   Type(dcsr_matrix), Intent(In) :: aCsr
   Real(dp), Dimension(:,:), Allocatable, Intent(Out) :: aFull
   Type(dcsr_matrix) :: aux
   Type(dcoo_matrix) :: aCOO

   Call dcsr_copy(aCSR, Aux)
   Call dcsr_sort(aux)
   Call d_csr2coo  (aux,  aCOO)
   Call d_coo2dense(aCOO, aFull)

 End Subroutine d_csr2dense 

 Subroutine d_old_csr2dia(aCSR, aDIA)
   Type(dcsr_matrix), Intent(In) :: aCSR
   Type(ddia_matrix), Intent(Out):: aDIA
   Real(dp), Dimension(:,:), Allocatable :: aFull
   Integer :: nl, nu
   Integer :: idiag, irow

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
   print *, 'nl=', nl
   print *, 'nu=', nu
   aDIA%nl = nl
   aDIA%nu = nu
   aDIA%ncol = aCSR%ncol
   aDIA%nrow = aCSR%nrow
   Allocate(aDIA%dat(nl+nu+1,aCSR%ncol))
   
   Call d_csr2dense(aCSR, aFull)

   !fill-in the upper-diagonals
   Do idiag = 1,nu
      Do irow = 1,aCSR%nrow-idiag
         aDIA%dat(nu-idiag+1,idiag+irow) = aFull(irow,irow+idiag) 
      End Do
   End Do
   !fill-in diagonal
   Do irow = 1,aCSR%nrow
   aDIA%dat(nu+1,irow) = aFull(irow,irow) 
   End Do
   !fill-in the upper-diagonals
   Do idiag = 1,nl
      Do irow = idiag,aCSR%nrow
         aDIA%dat(nu+1+idiag,irow-idiag+1) = aFull(irow,irow-idiag+1) 
      End Do
   End Do
   
 End Subroutine d_old_csr2dia

 Subroutine d_dia4gbtrf(aDIA, luDIA)
   type(ddia_matrix), Intent(In) :: aDIA
   Real(dp), Dimension(:,:), Allocatable, Intent(Out) :: luDIA
   Allocate(luDIA(aDIA%nl*2 + aDIA%nu + 1, aDIA%ncol))
   luDIA(aDIA%nl+1:aDIA%nl*2 + aDIA%nu + 1, :) = aDIA%dat
 End Subroutine d_dia4gbtrf

   

 ! ================================================ !
 !                                                  !
 !           DOUBLE COMPLEX ROUTINES                !
 !                                                  !
 ! ================================================ !
  
 Subroutine zcsr_copy (A,B)
   Type(zcsr_matrix), Intent(In)  :: A
   Type(zcsr_matrix), Intent(Out) :: B
   B%nelems = A%nelems
   B%nrow   = A%nrow
   B%ncol   = A%ncol
   Allocate(B%dat(B%nelems))
   Allocate(B%col(B%nelems))
   Allocate(B%row(B%nrow+1))   
   B%dat = A%dat
   B%col = A%col
   B%row = A%row
   Call B%sort()            
 End Subroutine zcsr_copy

 Subroutine z_csr2coo(aCSR, aCOO)

   Type(zcsr_matrix), Intent(In) :: aCSR
   Type(zcoo_matrix), Intent(Out):: aCOO
   Integer :: irow

   aCOO%nelems = aCSR%nelems
   aCOO%nrow   = aCSR%nrow
   aCOO%ncol   = aCSR%ncol
   Allocate(aCOO%dat(aCSR%nelems))
   Allocate(aCOO%col(aCSR%nelems))
   Allocate(aCOO%row(aCSR%nelems))
   aCOO%dat = aCSR%dat      
   aCOO%col = aCSR%col      

   Do irow = 1, aCOO%nrow
      aCOO%row(aCSR%row(irow):(aCSR%row(irow+1)-1)) = irow
   End do

 End Subroutine z_csr2coo


 Subroutine z_coo2csr(aCOO, aCSR)

   Type(zcoo_matrix), Intent(In) :: aCOO
   Type(zcsr_matrix), Intent(Out):: aCSR
   Integer :: irow
   Integer :: fill_in, iel, this_row_length
   
   Allocate(aCSR%row(aCOO%nrow+1))
   Allocate(aCSR%col(aCOO%nelems))
   Allocate(aCSR%dat(aCOO%nelems))
   aCSR%nrow = aCOO%nrow
   aCSR%ncol = aCOO%ncol
   aCSR%nelems = aCOO%nelems

   aCSR%row(1) = 1
   fill_in = 0
   Do irow = 1, aCOO%nrow
      this_row_length = 0 
      Do iel = 1, aCOO%nelems
         if (aCOO%row(iel).eq.irow) then
            this_row_length = this_row_length + 1 
            fill_in = fill_in + 1 
            aCSR%dat(fill_in) = aCOO%dat(iel)
            aCSR%col(fill_in) = aCOO%col(iel)
         end if
         aCSR%row(irow+1) = aCSR%row(irow) + this_row_length 
      end Do
   End Do  

   Call zcsr_sort(aCSR)
   aCSR%is_sorted = .True.

 End Subroutine z_coo2csr

 Subroutine z_coo2dense(aCOO, aFull)
   Type(zcoo_matrix), Intent(In) :: aCOO
   Complex(dp), Dimension(:,:), Allocatable, intent(Out) :: aFull
   Integer :: iel

   Allocate(aFull(aCOO%nrow, aCOO%ncol))
   aFull = Cmplx( 0._dp, 0._dp, kind=dp)
   Do iel = 1, aCOO%nelems
      aFull(aCOO%row(iel), aCOO%col(iel)) = aCOO%dat(iel)
   End Do
 End Subroutine z_coo2dense 

 Subroutine z_csr2dense(aCSR, aFull)

   Type(zcsr_matrix), Intent(In) :: aCsr
   Type(zcsr_matrix) :: aux 
   Complex(dp), Dimension(:,:), Allocatable, Intent(Out) :: aFull
   Type(zcoo_matrix) :: aCOO

   Call zcsr_copy(aCSR, aux)
   Call zcsr_sort(aux)
   Call z_csr2coo  (aux, aCOO)
   Call z_coo2dense(aCOO, aFull)

 End Subroutine z_csr2dense 


 Subroutine z_dia4gbtrf(aDIA, luDIA)
   type(zdia_matrix), Intent(In) :: aDIA
   Complex(dp), Dimension(:,:), Allocatable, Intent(Out) :: luDIA
   Allocate(luDIA(aDIA%nl*2 + aDIA%nu + 1, aDIA%ncol))
   luDIA(aDIA%nl+1:aDIA%nl*2 + aDIA%nu + 1, :) = aDIA%dat
 End Subroutine z_dia4gbtrf

 End Module Sparse_conversions
