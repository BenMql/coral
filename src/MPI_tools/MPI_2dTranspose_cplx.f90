Module MPI_2dTranspose_cplx
Use Double
Use MPI_VARS
Implicit None

Type :: mpi2dT_buffer_cplx
   Complex(dp), Pointer :: AN(:,:)
   Complex(dp), Pointer :: AT(:,:)
   Complex(dp), Allocatable, Dimension(:,:) :: BN 
   Complex(dp), Allocatable, Dimension(:,:) :: BT 
   Integer :: NZ, NX, locZ, locX, block_size
   Character(len=2) :: state
   Contains
   Procedure :: init =>mpi2Dt_Check_and_init_cplx
   Procedure :: pack_T2N =>mpi2Dt_Pack_T2N_cplx
   Procedure :: unpack_T2N =>mpi2Dt_unPack_T2N_cplx
   Procedure :: pack_N2T =>mpi2Dt_Pack_N2T_cplx
   Procedure :: unpack_N2T =>mpi2Dt_unPack_N2T_cplx
   Procedure :: Forward_transpose =>mpi2Dt_transpose_N2T_cplx
   Procedure :: inverse_transpose =>mpi2Dt_transpose_T2N_cplx
End Type mpi2dT_buffer_cplx
   
Contains

Subroutine mpi2Dt_Check_and_init_cplx(self, N_Cheb, N_four)
                                
  Class(mpi2dT_buffer_cplx), intent(InOut) :: self
  Integer, Intent(In) :: N_Cheb, N_Four

  ! >>> first check that the dimensions are OK
  if (modulo(N_cheb, world_size).ne.0) Then
     Print*, '>>> ERROR'
     Print*, 'The number of Chebyshev Pols. (including the 3/2'
     Print*, 'dealiasing factor) is not a multiple of the num-'
     Print*, 'ber of processes. Stopping now.'                      
     Stop 
  End If
  if (modulo(N_four, 2*world_size).ne.0) Then
     Print*, '>>> ERROR'
     Print*, 'The number of  Fourier  modes (including the 3/2'
     Print*, 'dealiasing factor) is not a multiple of the num-'
     Print*, 'ber of processes. Stopping now.'                      
     Stop 
  End If
  self%NX = N_four/2
  self%NZ = N_cheb
  self%locX = N_four/2/world_size
  self%locZ = N_cheb/world_size
  self%block_size = self%locX * self%locZ

  Allocate(self%BN (self%block_size, world_size))
  Allocate(self%BT (self%block_size, world_size))
  Allocate(self%AN (self%NZ, self%locX))
  Allocate(self%AT (self%NX, self%locZ))

End Subroutine mpi2Dt_Check_and_init_cplx

Subroutine mpi2Dt_Pack_N2T_Cplx(self)
  Class(mpi2dT_buffer_cplx), intent(InOut) :: self
  ! --- local variables ---
  Integer :: iblock, index_a, index_b
  Integer :: ix, iz

   self%BN(:,:) = Cmplx(0._dp, 0._dp, kind=dp)            
   Do iblock = 1,world_size
      Do ix = 1, self%locX
         Do iz = 1, self%locZ
            index_b = (ix-1)*self%locZ + iz
            index_a = (iblock-1) * self%locZ + iz
            self%BN (index_b, iblock) = self%AN(index_a, ix)
         End Do
      End Do
   End Do
End Subroutine mpi2Dt_Pack_N2T_Cplx


Subroutine mpi2Dt_UnPack_T2N_Cplx(self)
  Class(mpi2dT_buffer_cplx), intent(InOut) :: self
  ! --- local variables ---
  Integer :: iblock, index_a, index_b
  Integer :: ix, iz

   self%AN(:,:) = Cmplx(0._dp, 0._dp, kind=dp)            
   Do iblock = 1,world_size
      Do ix = 1, self%locX
         Do iz = 1, self%locZ
            index_b = (ix-1)*self%locZ + iz
            index_a = (iblock-1) * self%locZ + iz
            self%AN(index_a, ix) = self%BN (index_b, iblock) 
         End Do
      End Do
   End Do
End Subroutine mpi2Dt_UnPack_T2N_Cplx

Subroutine mpi2dt_fill_in_dummy(self)
  Class(mpi2dT_buffer_cplx), intent(InOut) :: self
  Integer :: ix, iz
  Do ix = 1, self%locX
  Do iz = 1, self%NZ  
    self%AN(iz,ix) = Cmplx((my_rank+1)*1000+ix+self%NX*iz, &
                           (my_rank+1)*2000+ix+self%NX*iz, kind=dp)
  End do
  End do
End Subroutine mpi2dt_fill_in_dummy

Subroutine mpi2Dt_unPack_N2T_Cplx(self)
  Class(mpi2dT_buffer_cplx), intent(InOut) :: self
  ! --- local variables ---
  Integer :: iblock, index_a, index_b
  Integer :: ix, iz

   self%AT(:,:) = Cmplx(0._dp, 0._dp, kind=dp)            
   Do iblock = 1,world_size
      Do iz = 1, self%locz
         Do ix = 1, self%locx
            index_b = (ix-1)*self%locZ + iz
            index_a = (iblock-1) * self%locX + ix
            self%AT(index_a, iz) = self%BT (index_b, iblock)
         End Do
      End Do
   End Do
End Subroutine mpi2Dt_unPack_N2T_Cplx

Subroutine mpi2Dt_Pack_T2N_Cplx(self)
  Class(mpi2dT_buffer_cplx), intent(InOut) :: self
  ! --- local variables ---
  Integer :: iblock, index_a, index_b
  Integer :: ix, iz

   self%BT(:,:) = Cmplx(0._dp, 0._dp, kind=dp)            
   Do iblock = 1,world_size
      Do iz = 1, self%locz
         Do ix = 1, self%locx
            index_b = (ix-1)*self%locZ + iz
            index_a = (iblock-1) * self%locX + ix
            self%BT (index_b, iblock) = self%AT(index_a, iz) 
         End Do
      End Do
   End Do
End Subroutine mpi2Dt_Pack_T2N_Cplx


Subroutine mpi2Dt_transpose_N2T_Cplx(self)
  Class(mpi2dT_buffer_cplx), intent(InOut) :: self
   
   Call self%pack_N2T()
   Call MPI_AllToAll(self%BN, 2*self%block_size, MPI_DOUBLE, &
                     self%BT, 2*self%block_size, MPI_DOUBLE, MPI_COmm_world, ierr)
   Call self%unpack_N2T()
   
End Subroutine mpi2Dt_transpose_N2T_Cplx

Subroutine mpi2Dt_transpose_T2N_Cplx(self)
  Class(mpi2dT_buffer_cplx), intent(InOut) :: self

   Call self%pack_T2N()
   !the do the transpose
   Call MPI_AllToAll(self%BT, 2*self%block_size, MPI_DOUBLE, &
                     self%BN, 2*self%block_size, MPI_DOUBLE, MPI_COmm_world, ierr)
   Call self%unpack_T2N()
   
End Subroutine mpi2Dt_transpose_T2N_Cplx
End Module MPI_2dTranspose_cplx
