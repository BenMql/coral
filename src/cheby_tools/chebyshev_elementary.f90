 module Chebyshev_elementary
 use Sparse_manipulations
 use MPI_Vars, only: my_rank
 use fortran_kinds, only: dp
 implicit None
   
 Contains

 subroutine chebyshev_first_N_primitives(N_op, ChebyOp, N, CEM_I)
   type(dcsr_matrix), intent(in) :: CEM_I
   integer, intent(in) :: N_op
   integer, intent(in) :: N
   type(dcsr_matrix), dimension(:), allocatable, intent(out) :: ChebyOp
   integer :: i 
   type(dcsr_matrix) :: aux
   allocate( ChebyOp(N_op+1))
   call dcsr_create_identity(ChebyOp(1),N) 
   if (N_op.ge.1) then
   call dcsr_truncate(CEM_I, ChebyOp(2), 1, N, 1, N)
   do i = 2,N_op 
      call dcsr_exp(CEM_I,i, aux )  
      call dcsr_truncate(aux, ChebyOp(i+1),   1, N, 1, N)
   end do
   end if
 end subroutine chebyshev_first_N_primitives


 subroutine chebyshev_elementary_matrices(N, center, gap, CEM_I, CEM_M, padding)
   type(dcsr_matrix), intent(out) :: CEM_I
   type(dcsr_matrix), intent(out) :: CEM_M
   integer, intent(in) :: padding
   integer, intent(in) :: N
   real(dp), intent(in) :: center, gap
   ! auxilliary arrays to build the matrices in COO format:
   type(dcoo_matrix) :: aux
   ! misc. stuff:
   integer :: Npad   ! size of the padded matrix N+padding
   integer :: NM, NI ! number of elements in M and I, resp. 
   integer :: i
   integer :: fill_in_index

   Npad = N + padding
   if (my_rank.eq.0) then
   write (*,*)  '================================================================='
   write (*,*)  '...'
   write (*,*)  '... Creating Chebyshev Elementary Matrices for integration and'
   110 format (' ... product, with total size ',I4,' after padding with ',I2,' elmts.')
   write (*,110) Npad, padding 
   write (*,*)  '...'
   write (*,*)  '================================================================='
   write (*,*) New_Line('s')
   end if

   NM = 3 * Npad - 2   ! -1, 0, and +1 diagonals (main, sub- and super-)
   NI = 2 * Npad - 3   ! -1 and +1 diagonals, first row = 0 
   
   fill_in_index = 1   ! use that to keep trace of where we are
   
   ! ###################################################################################
   ! Building CEM_M
   ! ...................................................................................
   Allocate(aux%dat(NM), CEM_M%dat(NM))
   Allocate(aux%col(NM), CEM_M%col(NM))
   Allocate(aux%row(NM), CEM_M%row(Npad+1))      ! CSR: size Npad+1 (*not* a typo)
   aux%nrow = Npad
   aux%ncol = Npad
   aux%nelems = NM                     
   ! CORRECTING SOME STUFF HERE?
   ! super diagonal:
      ! the first element is gap/2
      ! correction: this is actually the *SUB*diagonal
      i=1
      aux%dat(fill_in_index) = gap/2.0_dp
      aux%row(fill_in_index) = i+1
      aux%col(fill_in_index) = i
      fill_in_index = fill_in_index + 1
      ! and the rest of the super-diagonal is gap/4 
      Do i=1,Npad-1
         aux%dat(fill_in_index) = gap/4._dp
         aux%row(fill_in_index) = i
         aux%col(fill_in_index) = i+1
         fill_in_index = fill_in_index + 1
      End Do
   ! ...................................................................................
   ! sub diagonal is gap/4 (no exception):
      Do i=2,Npad-1
         aux%dat(fill_in_index) = gap/4._dp
         aux%row(fill_in_index) = i+1
         aux%col(fill_in_index) = i
         fill_in_index = fill_in_index + 1
      End Do
   ! ...................................................................................
   ! main diagonal is center:
      Do i=1,Npad
         aux%dat(fill_in_index) = center
         aux%row(fill_in_index) = i
         aux%col(fill_in_index) = i
         fill_in_index = fill_in_index + 1
      End Do
   ! convert and clean behind us:
   CALL d_coo2csr(aux,CEM_M)
   Deallocate(aux%dat, aux%row, aux%col)
   ! ...................................................................................
   ! End of building CEM_M
   ! ###################################################################################

   
   fill_in_index = 1   ! use that to keep trace of where we are

   ! ###################################################################################
   ! Building CEM_I
   ! ...................................................................................
   Allocate(aux%dat(NI), CEM_I%dat(NI))
   Allocate(aux%col(NI), CEM_I%col(NI))
   Allocate(aux%row(NI), CEM_I%row(Npad+1))      ! CSR: size Npad+1 (*not* a typo)
   aux%nrow = Npad
   aux%ncol = Npad
   aux%nelems = NI                     
   ! sub diagonal: 
      ! the first element is 2
      i=1
      aux%dat (fill_in_index) = 2.0_dp
      aux%row (fill_in_index) = i+1
      aux%col (fill_in_index) = i
      fill_in_index = fill_in_index + 1
      ! other elements are 1/col (where col is 1-based, fortran-like index)
      Do i = 2,Npad-1
         aux%col (fill_in_index) = i
         aux%row (fill_in_index) = i+1
         aux%dat (fill_in_index) = 1._dp/real(i,kind=dp)
         fill_in_index = fill_in_index + 1
      End Do
   ! super diagonal: nothing on row 1. Other entries are -1/(row-1) [Fortran]
      Do i = 2,Npad-1
         aux%row (fill_in_index) = i
         aux%col (fill_in_index) = i+1
         aux%dat (fill_in_index) = -1._dp/(real(i-1,kind=dp))
         fill_in_index = fill_in_index + 1
      End Do
   ! renormalize everything with gap/4
   aux%dat = aux%dat*gap / 4.0_dp
   ! convert and clean behind us:
   CALL d_coo2csr(aux,CEM_I)
   DeAllocate(aux%dat, aux%row, aux%col)
   ! ...................................................................................
   ! End of building CEM_I
   ! ###################################################################################


   
 end subroutine chebyshev_elementary_matrices   

 end module chebyshev_elementary
      
