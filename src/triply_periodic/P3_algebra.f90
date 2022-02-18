 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: P3_algebra            
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel, benjamin.miquel@tutanota.com
 !
 ! DESCRIPTION
 !> Linear algebra operations for handling (relatively small) sets of equations
 !! decoupled along (kx,ky,kz).
 !
 !=============================================================================
 Module P3_algebra

 Use Fortran_kinds
 use P3_lapack_wrappers

 Implicit None

 Type :: operator_3d_T
    Complex(kind=dp), dimension(:,:,:,:,:), Allocatable :: mat
    Complex(kind=dp), dimension(:,:,:,:,:), Allocatable :: inv
    Integer :: n1, n2, n3, nvar
  Contains
    Procedure :: initialise => initialise_zOp5
    Procedure :: invert => invert_zOp5_lapack
    Procedure :: action_zOp5_upon_zOp5_derivedTypesOnly_outOfPlace
    Procedure :: action_zOp5_upon_zOp5_TypeArrayArray_outOfPlace
    Procedure :: action_zOp5_upon_zVec4_outOfPlace
    Procedure :: action_zOp5_upon_zVec4_inPlace
    Procedure :: action_inv_zOp5_upon_zVec4_outOfPlace
    Procedure :: action_inv_zOp5_upon_zVec4_inPlace
    Procedure :: check  => check_inverse_zOp5
    Procedure :: zOp5_sum_zOp5_zarrays
    Procedure :: zOp5_sum_zOp5_dtypes
    generic   :: equals_sum_of => zOp5_sum_zOp5_dtypes, zOp5_sum_zOp5_zarrays
    Procedure :: fill_in => fill_zOp5
    generic   :: dot  => action_zOp5_upon_zOp5_derivedTypesOnly_outOfPlace,&
                         action_zOp5_upon_zOp5_typeArrayArray_outOfPlace, &
                         action_zOp5_upon_zVec4_outOfPlace, &
                         action_zOp5_upon_zVec4_inPlace
    generic   :: backsolve => action_inv_zOp5_upon_zVec4_outOfPlace, &
                              action_inv_zOp5_upon_zVec4_inPlace
 End type operator_3d_T


 Contains
 
 !> @brief
 !! creates an empty operator
 subroutine initialise_zOp5(self, nnvar, nn1, nn2, nn3)

   class(operator_3d_T) :: self
   integer, intent(in) :: nnvar !< number of coupled variables
   integer, intent(in) :: nn1   !< fast index 
   integer, intent(in) :: nn2   !< intermediate index
   integer, intent(in) :: nn3   !< slow index 

   if (Allocated(self%mat)) DeAllocate (self%mat) 
   if (Allocated(self%inv)) DeAllocate (self%inv)

   self%n1 = nn1
   self%n2 = nn2
   self%n3 = nn3
   self%nVar = nnvar

   allocate(self%mat(nn1, nn2, nn3, nnvar, nnvar), source=cmplx(0._dp, 0._dp, kind=dp))

 end subroutine initialise_zOp5

 !> @brief
 !! compute the inverse for all operators
 !> @details
 !! for all k1, k2, k3, we compute self%inv(k1,k2,k3,:,:)
 !! so that matmul( self%inv(k1,k2,k3,:,:), self%mat(k1,k2,k3,:,:) ) = identity 
 subroutine invert_zOp5_lapack (self)
   class(operator_3d_T) :: self
   complex(kind=dp), dimension(:,:), allocatable :: buf2d !< a temporary buffer
                                         !! for storing each matrix contiguously

   integer :: n1, n2, n3             !< shorthand
   integer :: i1, i2, i3, ivar, jvar !< dummy indices
   N1 = self%n1       
   N2 = self%n2       
   N3 = self%n3       
   if (allocated(self%inv)) DeAllocate(self%inv)
   Allocate(self%inv(n1, n2, n3, self%nvar, self%nvar), &
            source = Cmplx(0._dp, 0._dp, kind=dp))
   
   allocate(buf2d(self%nvar, self%nvar))
   do i3 = 1, n3
   do i2 = 1, n2
   do i1 = 1, n1
      do ivar = 1, self%nVar
      do jvar = 1, self%nVar
      buf2d(ivar, jvar) = self%mat(i1,i2,i3,ivar,jvar)
      end do
      end do
      call wrap_LUinvert_z(buf2d, .False.)
      do ivar = 1, self%nVar
      do jvar = 1, self%nVar
      self%inv(i1,i2,i3,ivar,jvar) = buf2d(ivar, jvar)
      end do
      end do
   end do
   end do
   end do

 end subroutine invert_zOp5_lapack

 subroutine action_zOp5_upon_zOp5_typeArrayArray_outOfPlace (self, other, self_dot_other)
   class(operator_3d_T) :: self
   complex(kind=dp), dimension(:,:,:,:,:), allocatable, intent(in)  :: other
   complex(kind=dp), dimension(:,:,:,:,:), allocatable, intent(out) :: self_dot_other
   integer :: ivar, kvar, jvar

   allocate(self_dot_other(self%N1,self%N2,self%N3,self%Nvar,self%Nvar),&
            source = Cmplx(0._dp, 0._dp, kind=dp))

   do ivar = 1,self%Nvar
   do jvar = 1,self%Nvar
   do kvar = 1,self%Nvar
   self_dot_other(:,:,:,ivar, jvar) = self_dot_other(:,:,:,ivar, jvar) +&
                           self%mat(:,:,:,ivar,kvar) * other(:,:,:,kvar, jvar)
   end do
   end do
   end do
 end subroutine action_zOp5_upon_zOp5_typeArrayArray_outOfPlace

 subroutine action_zOp5_upon_zOp5_derivedTypesOnly_outOfPlace(self, other, self_dot_other)
   class(operator_3d_T) :: self
   type(operator_3d_T), intent(in) :: other
   type(operator_3d_T), intent(out) :: self_dot_other
   Integer :: ivar, kvar, jvar
   integer :: nref

   nref = self%nVar
   if (other%nVar .ne. nref) print *, 'inconsistent sizes, matmul Op5 x Op5'
   if (other%nVar .ne. nref) stop
   
   call self_dot_other%initialise( self%nVar, self%n1, self%n2, self%n3)

   do ivar = 1,self%Nvar
   do jvar = 1,self%Nvar
   do kvar = 1,self%Nvar
   self_dot_other%mat(:,:,:,ivar, jvar) = self_dot_other%mat(:,:,:,ivar, jvar) +&
                           self%mat(:,:,:,ivar,kvar) * other%mat(:,:,:,kvar, jvar)
   end do
   end do
   end do
 end subroutine action_zOp5_upon_zOp5_derivedTypesOnly_outOfPlace

 subroutine action_zOp5_upon_zVec4_outOfPlace (self, other, self_dot_other)
   class(operator_3d_T) :: self
   complex(kind=dp), dimension(:,:,:,:), allocatable, intent(in)  :: other
   complex(kind=dp), dimension(:,:,:,:), allocatable, intent(out) :: self_dot_other
   integer :: ivar, kvar

   allocate(self_dot_other(self%N1,self%N2,self%N3,self%nVar))
   self_dot_other = Cmplx(0._dp, 0._dp, Kind=dp)

   do ivar = 1,self%nVar
   do kvar = 1,self%nVar
   self_dot_other(:,:,:,ivar) = self_dot_other(:,:,:,ivar) +&
                           self%mat(:,:,:,ivar,kvar) * other(:,:,:,kvar)
   end do
   end do
 end subroutine action_zOp5_upon_zVec4_outOfPlace

 subroutine action_zOp5_upon_zVec4_inPlace (self, other)
   class(operator_3d_T) :: self
   complex(kind=dp), dimension(:,:,:,:), allocatable, intent(inOut)  :: other
   complex(kind=dp), dimension(:,:,:,:), allocatable :: self_dot_other
   integer :: ivar, kvar

   allocate(self_dot_other(self%N1,self%N2,self%N3,self%nVar))
   self_dot_other = cmplx(0._dp, 0._dp, Kind=dp)

   do ivar = 1,self%nVar
   do kvar = 1,self%nVar
   self_dot_other(:,:,:,ivar) = self_dot_other(:,:,:,ivar) +&
                           self%mat(:,:,:,ivar,kvar) * other(:,:,:,kvar)
   end do
   end do
   other = self_dot_other
   deAllocate (self_dot_other)
 end subroutine action_zOp5_upon_zVec4_inPlace

 subroutine action_inv_zOp5_upon_zVec4_outOfPlace (self, other, self_dot_other)
   class(operator_3d_T) :: self
   complex(kind=dp), dimension(:,:,:,:), allocatable, intent(in)  :: other
   complex(kind=dp), dimension(:,:,:,:), allocatable, intent(out) :: self_dot_other
   integer :: ivar, kvar

   allocate(self_dot_other(self%N1,self%N2,self%N3,self%nVar))
   self_dot_other = Cmplx(0._dp, 0._dp, Kind=dp)

   do ivar = 1,self%nVar
   do kvar = 1,self%nVar
   self_dot_other(:,:,:,ivar) = self_dot_other(:,:,:,ivar) +&
                           self%inv(:,:,:,ivar,kvar) * other(:,:,:,kvar)
   end do
   end do
 end subroutine action_inv_zOp5_upon_zVec4_outOfPlace

 subroutine action_inv_zOp5_upon_zVec4_inPlace (self, other)
   class(operator_3d_T) :: self
   complex(kind=dp), dimension(:,:,:,:), allocatable, intent(inOut)  :: other
   complex(kind=dp), dimension(:,:,:,:), allocatable :: self_dot_other
   integer :: ivar, kvar

   allocate(self_dot_other(self%N1,self%N2,self%N3,self%nVar))
   self_dot_other = cmplx(0._dp, 0._dp, Kind=dp)

   do ivar = 1,self%nVar
   do kvar = 1,self%nVar
   self_dot_other(:,:,:,ivar) = self_dot_other(:,:,:,ivar) +&
                           self%inv(:,:,:,ivar,kvar) * other(:,:,:,kvar)
   end do
   end do
   other = self_dot_other
   deAllocate (self_dot_other)
 end subroutine action_inv_zOp5_upon_zVec4_inPlace

 Subroutine zOp5_sum_zOp5_zArrays(self, Amat, Bmat, zSca)
   class(operator_3d_T) :: self
   complex(kind=dp), dimension(:,:,:,:,:), allocatable, intent(in) :: Amat
   complex(kind=dp), dimension(:,:,:,:,:), allocatable, intent(in) :: Bmat
   complex(kind=dp), intent(in) :: zSca
   
   if (allocated(self%mat)) DeAllocate (self%mat) 
   call self%initialise( size(Amat,4), &
                         size(Amat,1), &
                         size(Amat,2), &
                         size(Amat,3)  ) 

   if (self%nVar.ne. size(Amat,5)) print *, "Error in adding operators: dim 4 must equal dim 5"
   if (self%nVar.ne. size(Amat,5)) stop

   self%mat = Amat + zSca * Bmat

 End Subroutine zOp5_sum_zOp5_zArrays

 Subroutine zOp5_sum_zOp5_dtypes(self, typeA, typeB, zSca)
   class(operator_3d_T) :: self
   class(operator_3d_T), intent(in) :: typeA
   class(operator_3d_T), intent(in) :: typeB
   complex(kind=dp), intent(in) :: zSca
   
   call self%initialise( typeA%nVar, typeA%n1, typeA%n2, typeA%n3)

   if ((self%n1  .ne. typeB%n1  )  .or. &
       (self%n2  .ne. typeB%n2  )  .or. &
       (self%n3  .ne. typeB%n3  )  .or. &
       (self%nVar.ne. typeB%nVar)) then
      print *, "Error in adding operators: non conformant"
      stop
   end if                                    

   self%mat = typeA%mat + zSca * typeB%mat

 End Subroutine zOp5_sum_zOp5_dtypes

 Subroutine check_inverse_zOp5(self)
   Class(operator_3d_T) :: self
   Complex(kind=dp), Dimension(:,:,:,:,:), Allocatable :: auxMat 
   Integer :: i1, i2, i3, ivar
   Allocate(auxMat(self%N1,self%N2,self%N3,self%Nvar,self%Nvar))
   auxMat = Cmplx(0._dp, 0._dp, Kind=dp)
   Call self%dot( self%inv, auxMat)
   !Print *, auxMat 
   Do ivar = 1,self%Nvar
   Do i3 = 1,self%n3
   Do i2 = 1,self%n2
   Do i1 = 1,self%n1
   auxMat(i1,i2,i3,ivar,ivar) = auxMat(i1,i2,i3,ivar,ivar) &
                               - Cmplx(1._dp, 0._dp, kind=dp)
   end do
   end do
   end do
   end do
   Print *, "Inverse correct to accuracy: ", sqrt(sum(abs(auxMat)**2))/&
                       Real(self%N1*self%N2*self%N3*self%nvar**2, kind=dp)
 End Subroutine check_inverse_zOp5
 
 !> @brief
 !! In the coupling matrix self, add a coupling term for variable ivar 
 !! equation ieq.
 Subroutine fill_zOp5(self, arr_3d, ivar, ieq)
   Class(operator_3d_T) :: self
   Complex(kind=dp), Dimension(:,:,:), Allocatable, Intent(In) :: arr_3d
   Integer, intent(in) :: ivar, ieq
   integer :: i1,i2,i3

   Do i3 = 1,self%n3
   Do i2 = 1,self%n2
   Do i1 = 1,self%n1
      self%mat(i1,i2,i3,ieq,ivar) = self%mat(i1,i2,i3,ieq,ivar) + arr_3d(i1,i2,i3)
   End Do
   End Do
   End Do

 End Subroutine fill_zOp5

 End Module P3_algebra



 !!@brief 
 !!< a naive (disastrously inefficient?) implementation
 !Subroutine invert_zOp5 (self)
   !Class(operator_3d_T) :: self
   !Complex(kind=dp), Dimension(:,:,:,:,:), Allocatable :: Comat
   !Integer :: n1, n2, n3
   !Integer :: i1, i2, i3, ivar, jvar
   !N1 = self%n1       
   !N2 = self%n2       
   !N3 = self%n3       
   !if (allocated(self%inv)) DeAllocate(self%inv)
   !Allocate(self%inv(n1, n2, n3, self%nvar, self%nvar), &
            !source = Cmplx(0._dp, 0._dp, kind=dp))
   !Select Case (self%Nvar)
   !Case(1)
      !self%inv(:,:,:,1,1) = 1._dp/self%mat(:,:,:,1,1)   
   !Case(2)
      !if (allocated(self%det)) DeAllocate(self%det)
      !Allocate(self%det(n1, n2, n3), source=cmplx(0._dp, 0._dp, kind=dp))
      !self%det = (self%mat(:,:,:,1,1)*self%mat(:,:,:,2,2)&
                !!- self%mat(:,:,:,1,2)*self%mat(:,:,:,2,1))
      !self%inv(:,:,:,1,1) = self%mat(:,:,:,2,2) / self%Det
      !self%inv(:,:,:,2,2) = self%mat(:,:,:,1,1) / self%Det
      !self%inv(:,:,:,1,2) =-self%mat(:,:,:,1,2) / self%Det
      !self%inv(:,:,:,2,1) =-self%mat(:,:,:,2,1) / self%Det
   !Case DEFAULT
      !Call self%compute_determinant()
      !Call self%compute_comatrix   (Comat)
      !Do ivar = 1,self%Nvar
      !Do jvar = 1,self%Nvar
         !self%inv(:,:,:,ivar,jvar) = Comat(:,:,:,jvar,ivar)/self%det 
      !End Do
      !End Do
   !End Select
 !End Subroutine invert_zOp5

 !Subroutine Compute_the_comatrix_zOp5(self, Comat)
   !Class(operator_3d_T) :: self
   !Complex(kind=dp), Dimension(:,:,:,:,:), Allocatable, Intent(out) :: Comat
   !Complex(kind=dp), Dimension(:,:), Allocatable:: My_slice
   !Complex(kind=dp), Dimension(:,:), Allocatable:: this_wn 
   !integer :: i1, i2, i3, ivar, jvar
   !Real(dp) :: dsca
   !Complex(dp) :: this_det
   !Allocate(Comat(self%N1,  & 
                  !self%N2,  &
                  !self%N3,  &
                  !self%Nvar,&
                  !self%Nvar))
   !allocate(this_wn(self%nvar, self%nvar))
   !Do i3=1,self%n3
   !Do i2=1,self%n2
   !Do i1=1,self%n1
   !
   !Deallocate(this_wn)
   !allocate(this_wn(self%nvar, self%nvar))
   !Do jvar = 1,self%Nvar
   !Do ivar = 1,self%Nvar
   !this_wn(ivar, jvar) = self%mat(i1,i2,i3,ivar, jvar)
   !End Do 
   !End Do 
!
   !Do jvar = 1,self%Nvar
   !Do ivar = 1,self%Nvar
     !Call remove_row_col( this_wn, ivar, jvar, my_slice)
     !dsca = (-1.0_dp)**(ivar+jvar)
     !this_det=recursive_determinant(my_slice, self%nvar-1)
     !Comat(i1,i2,i3,ivar,jvar) = dsca * this_det
   !End Do
   !End Do
   !End Do
   !End Do
   !End Do
!
   !
 !End Subroutine Compute_the_comatrix_zOp5
!
!
 !Subroutine Compute_the_determinant_zOp5(self)
   !Class(operator_3d_T) :: self
   !Complex(kind=dp), Dimension(:,:), Allocatable:: My_slice
   !integer :: i1, i2, i3
   !if (allocated(self%det)) DeAllocate(self%det)
   !Allocate(self%det(self%N1, self%N2, self%N3))
   !Allocate(my_slice(self%nvar,self%nvar))
   !Do i1=1,self%n1
   !Do i2=1,self%n2
   !Do i3=1,self%n3
      !my_slice = self%mat(i1,i2,i3,:,:)
      !self%Det(i1,i2,i3) = recursive_determinant(my_slice,self%nvar)
   !End Do
   !End Do
   !End Do
 !End Subroutine Compute_the_determinant_zop5
   !
 !Recursive Function recursive_determinant(a2,nv) result(my_det)
 !Complex(kind=dp), Dimension(:,:), Allocatable, Intent(IN) :: a2
 !Integer, Intent(In) :: nv 
 !Integer :: iv
 !Complex(kind=dp) :: my_det
 !Complex(kind=dp) :: aux
 !Complex(kind=dp), Dimension(:,:), Allocatable :: a_slice
 !Real   (kind=dp) :: mult_factor
 !If (nv.eq.1) then
   !my_det = a2(1,1)
 !Else
   !! we expand on the first line:
   !aux = Cmplx(0._dp, 0._dp, kind=dp)
   !do iv = 1,nv
   !Call remove_row_col(a2,1,iv,a_slice) 
   !if (modulo(iv+1,2).eq.0) then
      !mult_factor = 1._dp
   !else 
      !mult_factor =-1._dp
   !end if
   !aux = aux + mult_factor * a2(1,iv)* recursive_determinant(a_slice,nv-1)
   !end do
   !my_det = aux
 !End If
 !End Function recursive_determinant
!
 !Subroutine remove_row_col(a, irow, icol, b)
 !Complex(kind=dp), Dimension(:,:), Allocatable, Intent(IN)  :: a
 !Complex(kind=dp), Dimension(:,:), Allocatable, Intent(out) :: b 
 !Integer, Intent(IN) :: irow, icol
 !Integer :: ix, iy, ixb, iyb
 !Allocate(b(size(a,1)-1, size(a,2)-1))
 !iyb = 0
 !Do iy = 1,size(a,1)
 !if (iy.ne.icol) then
 !iyb = iyb+1
 !ixb = 0
 !Do ix = 1,size(a,2)
   !if (ix.ne.irow) then
   !ixb = ixb +1
   !b(ixb,iyb) = a(ix,iy)
   !end if
 !end do
 !End if
 !end do
 !End Subroutine remove_row_col
 
