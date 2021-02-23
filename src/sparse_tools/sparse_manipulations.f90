 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: Sparse_manipulations
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> To do later
 !
 !=============================================================================

 module sparse_manipulations
 use fortran_kinds
 use sparse_formats
 use sparse_conversions
 use sparse_blas
 use, Intrinsic :: ISO_C_BINDING
 implicit none


 contains

 Subroutine dcsr_mul_ddense_3D(dalpha, A, transpose_A, B, dbeta, C, N2, N3, sta3)
   !................................................................
   ! Intended use:
   ! Apply an operator (discretized by csr_matrix A) along the first
   ! dimension of a 3D array B. Update C with the result. Visually:
   ! Do i2=1,N2 and i3=sta3,N3
   ! C(:,i2,i3) = dalpha * A(:,j) * B(j,i2,i3) + dbeta * C(:,i2,i3)
   !................................................................
   Type(dcsr_matrix), Intent(in) :: A
   Real(Kind=dp), Pointer :: B(:,:,:)                                        
   Real(Kind=dp), Pointer :: C(:,:,:)                                        
   Real(Kind=dp), Intent(in) :: dalpha, dbeta
   Integer, Intent(In) :: N2, N3
   Integer, Intent(In) :: sta3
   Character(len=1), Intent(In) :: transpose_A
   type(C_Ptr) :: dumptr
   Real(Kind=dp), Pointer :: ptr2B(:,:)
   Real(Kind=dp), Pointer :: ptr2C(:,:)
   Character :: matdescra(4)
   matdescra = ['G','X', 'N', 'F']
   if (sta3.le.N3) Then
   dumptr = C_Loc(B(1,1,sta3))
   Call C_F_Pointer (dumptr, ptr2B, [A%Ncol, N2*(N3+1-sta3)] )
   dumptr = C_Loc(C(1,1,sta3))
   Call C_F_Pointer (dumptr, ptr2C, [A%Ncol, N2*(N3+1-sta3)] )
   Call mkl_dcsrmm(transpose_A,&
                   A%nrow, N2*(N3+1-sta3), A%ncol, dalpha, matdescra,&
                   A%dat, A%col, A%row(1:A%nrow), A%row(2:(A%nrow+1)), &
                   ptr2B, A%ncol, dbeta, ptr2C, A%nrow)
   End if
 End Subroutine dcsr_mul_ddense_3D
  


 !> @brief
 !> \callgraph Apply a sparse operator \p A on multiple right-hand sides \p B and
 !> store the result in \p C
 Subroutine zcsr_mul_zdense_3D(zalpha, A, transpose_A, B, zbeta, C, N2, N3, sta3)
   !................................................................
   ! Intended use:
   ! Apply an operator (discretized by csr_matrix A) along the first
   ! dimension of a 3D array B. Update C with the result. Visually:
   ! Do i2=1,N2 and i3=sta3,N3
   ! C(:,i2,i3) = zalpha * A(:,j) * B(j,i2,i3) + zbeta * C(:,i2,i3)
   !................................................................
   Type(zcsr_matrix), Intent(in) :: A
   Complex(Kind=dp), Pointer :: B(:,:,:)                                        
   Complex(Kind=dp), Pointer :: C(:,:,:)                                        
   Complex(Kind=dp), Intent(in) :: zalpha, zbeta
   Integer, Intent(In) :: N2, N3
   Integer, Intent(In) :: sta3
   Character(len=1), Intent(In) :: transpose_A
   type(C_Ptr) :: dumptr
   Complex(Kind=dp), Pointer :: ptr2B(:,:)
   Complex(Kind=dp), Pointer :: ptr2C(:,:)
   Character :: matdescra(4)
   matdescra = ['G','X', 'N', 'F']
   if (sta3.le.N3) Then
   dumptr = C_Loc(B(1,1,sta3))
   Call C_F_Pointer (dumptr, ptr2B, [A%Ncol, N2*(N3+1-sta3)] )
   dumptr = C_Loc(C(1,1,sta3))
   Call C_F_Pointer (dumptr, ptr2C, [A%Ncol, N2*(N3+1-sta3)] )
   Call mkl_zcsrmm(transpose_A,&
                   A%nrow, N2*(N3+1-sta3), A%ncol, zalpha, matdescra,&
                   A%dat, A%col, A%row(1:A%nrow), A%row(2:(A%nrow+1)), &
                   ptr2B, A%ncol, zbeta, ptr2C, A%nrow)
   End if
 End Subroutine zcsr_mul_zdense_3D
  
 Subroutine dcsr_mul_dvec_1D(dalpha, A, transpose_A, B, dbeta, C)
   !................................................................
   ! Intended use:
   ! Apply an operator (discretized by csr_matrix A) along the first
   ! dimension of a 3D array B. Update C with the result. Visually:
   ! Do i2=1,N2 and i3=sta3,N3
   ! C(:) = dalpha * A(:,j) * B(j) + dbeta * C(:)
   !................................................................
   Type(dcsr_matrix), Intent(in) :: A
   real(Kind=dp), Pointer :: B(:)
   real(Kind=dp), Pointer :: C(:)
   Real(Kind=dp), Intent(in) :: dalpha, dbeta
   Integer :: N2, N3, sta3
   Character(len=1), Intent(In) :: transpose_A
   type(C_Ptr) :: dumptr
   Real(Kind=dp), Pointer :: ptr2B(:,:)
   Real(Kind=dp), Pointer :: ptr2C(:,:)
   Character :: matdescra(4)
   matdescra = ['G','X', 'N', 'F']
   sta3 = 1
   N3 = 1
   N2 = 1
   dumptr = C_Loc(B(1))
   Call C_F_Pointer (dumptr, ptr2B, [A%Ncol, N2*(N3+1-sta3)] )
   dumptr = C_Loc(C(1))
   Call C_F_Pointer (dumptr, ptr2C, [A%Ncol, N2*(N3+1-sta3)] )
   Call mkl_dcsrmm(transpose_A,&
                   A%nrow, N2*(N3+1-sta3), A%ncol, dalpha, matdescra,&
                   A%dat, A%col, A%row(1:A%nrow), A%row(2:(A%nrow+1)), &
                   ptr2B, A%ncol, dbeta, ptr2C, A%nrow)
 End Subroutine dcsr_mul_dvec_1D

 Subroutine zcsr_mul_zdense_2D(zalpha, A, transpose_A, B, zbeta, C, N3, sta3)

   !................................................................
   ! Intended use:
   ! Apply an operator (discretized by csr_matrix A) along the first
   ! dimension of a 3D array B. Update C with the result. Visually:
   ! Do i2=1,N2 and i3=sta3,N3
   ! C(:,i2,i3) = zalpha * A(:,j) * B(j,i2,i3) + zbeta * C(:,i2,i3)
   !................................................................
   Type(zcsr_matrix), Intent(in) :: A
   Complex(Kind=dp), Pointer :: B(:,:)                                        
   Complex(Kind=dp), Pointer :: C(:,:)                                        
   Complex(Kind=dp), Intent(in) :: zalpha, zbeta
   Integer, Intent(In) :: N3
   Integer, Intent(In) :: sta3
   Character(len=1), Intent(In) :: transpose_A
   type(C_Ptr) :: dumptr
   Complex(Kind=dp), Pointer :: ptr2B(:,:)
   Complex(Kind=dp), Pointer :: ptr2C(:,:)
   Character :: matdescra(4)
   matdescra = ['G','X', 'N', 'F']
   if (sta3.le.N3) Then
   dumptr = C_Loc(B(1,sta3))
   Call C_F_Pointer (dumptr, ptr2B, [A%Ncol, (N3+1-sta3)] )
   dumptr = C_Loc(C(1,sta3))
   Call C_F_Pointer (dumptr, ptr2C, [A%Ncol, (N3+1-sta3)] )
   Call mkl_zcsrmm(transpose_A,&
                   A%nrow, (N3+1-sta3), A%ncol, zalpha, matdescra,&
                   A%dat, A%col, A%row(1:A%nrow), A%row(2:(A%nrow+1)), &
                   ptr2B, A%ncol, zbeta, ptr2C, A%nrow)
   End if
 End Subroutine zcsr_mul_zdense_2D
  
 Subroutine dcsr_mul_ddense_2D(dalpha, A, transpose_A, B, dbeta, C, N3, sta3)

   !................................................................
   ! Intended use:
   ! Apply an operator (discretized by csr_matrix A) along the first
   ! dimension of a 3D array B. Update C with the result. Visually:
   ! Do i2=1,N2 and i3=sta3,N3
   ! C(:,i2,i3) = zalpha * A(:,j) * B(j,i2,i3) + zbeta * C(:,i2,i3)
   !................................................................
   Type(dcsr_matrix), Intent(in) :: A
   Real(Kind=dp), Pointer :: B(:,:)                                        
   Real(Kind=dp), Pointer :: C(:,:)                                        
   Real(Kind=dp), Intent(in) :: dalpha, dbeta
   Integer, Intent(In) :: N3
   Integer, Intent(In) :: sta3
   Character(len=1), Intent(In) :: transpose_A
   type(C_Ptr) :: dumptr
   Real(Kind=dp), Pointer :: ptr2B(:,:)
   Real(Kind=dp), Pointer :: ptr2C(:,:)
   Character :: matdescra(4)
   matdescra = ['G','X', 'N', 'F']
   if (sta3.le.N3) Then
   dumptr = C_Loc(B(1,sta3))
   Call C_F_Pointer (dumptr, ptr2B, [A%Ncol, (N3+1-sta3)] )
   dumptr = C_Loc(C(1,sta3))
   Call C_F_Pointer (dumptr, ptr2C, [A%Ncol, (N3+1-sta3)] )
   Call mkl_dcsrmm(transpose_A,&
                   A%nrow, (N3+1-sta3), A%ncol, dalpha, matdescra,&
                   A%dat, A%col, A%row(1:A%nrow), A%row(2:(A%nrow+1)), &
                   ptr2B, A%ncol, dbeta, ptr2C, A%nrow)
   End if
 End Subroutine dcsr_mul_ddense_2D
  

 Subroutine zcsr_add_block_dcsr_z(A, B, sca, rowpos, colpos)
   Type(zcsr_matrix), Intent(InOut) :: A
   Type(dcsr_matrix), Intent(In)    :: B
   Complex(Kind=dp), Intent(In) :: sca
   Integer, Intent(In) :: rowpos, colpos
   Type(zcsr_matrix) :: zcopy
   Type(zcsr_matrix) :: buf
   Type(dcsr_matrix) :: aux
   Call dcsr_pad(B,A%nrow, A%ncol, rowpos, colpos, aux)
   Call dcsr_convert_zcsr(aux,zcopy, sca)
   Call wrap_zcsradd (A, zcopy, buf, Cmplx(1.0_dp, 0.0_dp, kind=dp))
   Call zcsr_copy(buf, A)
 End Subroutine zcsr_add_block_dcsr_z
 !
 !
 !
 !
 !
 Subroutine dcsr_add_block_dcsr_d(A, B, sca, rowpos, colpos)
   type(dcsr_matrix), Intent(InOut) :: A
   Type(dcsr_matrix), Intent(In) :: B
   Real(Kind=dp), Intent(In) :: sca
   Integer, Intent(In) :: rowpos, colpos
   Type(dcsr_matrix) :: buf
   Type(dcsr_matrix) :: aux
   Call dcsr_pad(B,A%nrow, A%ncol, rowpos, colpos, aux)
   Call wrap_dcsradd(A, aux, buf, sca)
   Call dcsr_copy(buf, A)
 End Subroutine dcsr_add_block_dcsr_d
 !
 !
 !
 !
 !
 !
 !

 SUBROUTINE dcsr_zero_out_n_lines(A,nzeros, C)
 
   Type(dcsr_matrix), Intent(In)  :: A 
   Type(dcsr_matrix), Intent(Out) :: C 
   Type(dcsr_matrix) :: aux
   Integer, Intent(In) :: nzeros
   Integer :: i
   ! mkl stuff
   Integer :: n 

   n = A%nrow   

   aux%nrow = N
   aux%ncol = N
   aux%nelems = N-nzeros

   ! create a diagonal matrix with nzeros zeros
   allocate( aux%dat(n-nzeros))
   allocate( aux%col(n-nzeros))
   allocate( aux%row(n+1))
   do i=1,nzeros
       aux%row(i) = 1
   end do
   do i = 1,(n-nzeros)
       aux%dat(i) = 1._dp
       aux%col(i) = nzeros + i 
       aux%row(nzeros+i) = i
   end do
   aux%row(n+1) = n-nzeros+1

   
   !Multiply A on the left
   Call dcsr_sort(aux)
   Call wrap_dcsrmultcsr(aux, A, C)

 End Subroutine dcsr_zero_out_n_lines
 !
 !
 !
 !
 !


 Subroutine dcsr_zero_out_last_n_cols(A,nzeros, C)
 
   Type(dcsr_matrix), Intent(IN)  :: A 
   Type(dcsr_matrix), Intent(out) :: C 
   Type(dcsr_matrix) :: aux
   Integer, Intent(In) :: nzeros
   Integer :: i
   ! mkl stuff
   Integer :: N 

   N = A%ncol   
   aux%nrow = N
   aux%ncol = N
   aux%nelems = N-nzeros


   Allocate( aux%dat(N-nzeros))
   Allocate( aux%col(N-nzeros))
   Allocate( aux%row(N+1))
   Do i = 1,(N-nzeros)
      aux%dat(i) = 1._dp
      aux%col(i) = i
      aux%row(i) = i
   End Do
   Do i = (N-nzeros+1), (N+1)
      aux%row(i) = (N-nzeros+1)
   End Do
   ! multiply matrix A on the right to zero out the nzeros last columns
   Call dcsr_sort(aux)
   Call wrap_dcsrmultcsr(A, aux, C)

 End Subroutine dcsr_zero_out_last_n_cols

 !
 !
 !
 !
 !

 !subroutine dcsr_mul_dcsr(self, other, third)
   !class(dcsr_matrix), intent(in) :: self
   !type (dcsr_matrix), intent(in) :: other
   !type (dcsr_matrix), intent(out):: third
   !real(kind=dp), dimension(:), allocatable :: Aux1, Aux2
   !integer,       dimension(:), allocatable :: IAux1, IAux2, JAux1, JAux2
   !integer :: Nmat
   !integer :: info
   !integer :: i
!
   !Nmat = A%nrow
   !allocate( Aux1(A%row(Nmat+1)-1))
   !allocate(JAux1(A%row(Nmat+1)-1))
   !allocate(IAux1(   Nmat+1   ))
   !allocate(IAux2(   Nmat+1   ))
   !allocate( Aux2(1), Jaux2(1)) ! just in case ...
   !Aux1  = A%dat
   !IAux1 = A%row
   !JAux1 = A%col
!
      !call mkl_dcsrmultcsr('n',1,8,self%nrow, self%nrow, selfNmat, Nmat, Nmat, Aux1,JAux1,IAux1,&
                                             !A%dat,A%col,A%row,Aux2,JAux2,IAux2,0,info)
      !deAllocate(aux2, jaux2)
      !allocate (Jaux2(Iaux2(Nmat+1)-1))
      !allocate ( aux2(Iaux2(Nmat+1)-1))
      !call mkl_dcsrmultcsr('n',2,8,Nmat, Nmat, Nmat, Aux1,JAux1,IAux1,&
                                             !A%dat,A%col,A%row,Aux2,JAux2,IAux2,0,info)
      !! copy everything back in aux1
      !deAllocate (aux1,jaux1)
      !iaux1 = iaux2
      !allocate (Jaux1(Iaux1(Nmat+1)-1))
      !allocate ( aux1(Iaux1(Nmat+1)-1))
      !Jaux1 = Jaux2
      !Aux1  = Aux2




 !
 !
 !
 !


 !
 !
 !
 !
 !

 subroutine dcsr_exp(A,expA, C)
   Type(dcsr_matrix), Intent(In)  :: A
   Type(dcsr_matrix), Intent(Out) :: C
   Type(dcsr_matrix) :: aux
   Integer, Intent(in) :: expA
   Call dcsr_create_identity(aux, A%nrow)
   Call dcsr_expexp(A, expA, aux, 0, C)
 end subroutine dcsr_exp
 
 
 !
 !
 !
 !
 !


 !
 !
 !
 !
 !

 Subroutine dcsr_pad(A,rowmax, colmax, rowpos, colpos,C)
   ! insert matrix A into a bigger matrix C of size (rowmax x colmax)
   ! at the position (rowpos, colpos). Other entries are zero.
   type(dcsr_matrix), Intent(in)  :: A ! the 'small' matrix
   type(dcsr_matrix), Intent(out) :: C ! the 'big' host matrix
   Integer, Intent(in) :: rowmax, colmax, rowpos, colpos
   integer :: i 
   C%nelems = A%nelems
   C%nrow = rowmax
   C%ncol = colmax
   allocate( C%col(C%nelems))
   allocate( C%dat(C%nelems))
   allocate( C%row(C%nrow+1))
   C%col  = A%col + colpos - 1
   C%dat = A%dat
   do i = 1,rowpos
      C%row(i) = 1
   end do
   do i = (rowpos+1), (rowpos+A%nrow)
      C%row(i) = A%row(i-rowpos+1)
   end do
   do i = (rowpos+A%nrow+1),(C%nrow+1)
      C%row(i) = A%row(A%nrow+1)
   end do
   Call dcsr_sort(C)
 end subroutine dcsr_pad
 
   
 !
 !
 !
 !
 !

 end module sparse_manipulations
