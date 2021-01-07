Module MPI_2dTranspose_real
Use Fortran_kinds
Use MPI_VARS
Implicit None

Type :: mpi2dT_buffer_real
   Real(dp), Pointer :: AN(:,:)
   Real(dp), Pointer :: AT(:,:)
   Real(dp), Allocatable, Dimension(:,:) :: BN 
   Real(dp), Allocatable, Dimension(:,:) :: BT 
   Integer :: NZ, NX, locZ, locX, block_size
   Character(len=2) :: state
   Contains
   Procedure :: init =>mpi2Dt_Check_and_init_Real
   Procedure :: pack_T2N =>mpi2Dt_Pack_T2N_Real
   Procedure :: unpack_T2N =>mpi2Dt_unPack_T2N_Real
   Procedure :: pack_N2T =>mpi2Dt_Pack_N2T_Real
   Procedure :: unpack_N2T =>mpi2Dt_unPack_N2T_Real
   Procedure :: Forward_transpose =>mpi2Dt_transpose_N2T_Real
   Procedure :: inverse_transpose =>mpi2Dt_transpose_T2N_Real
End Type mpi2dT_buffer_real
   
Contains

Subroutine mpi2Dt_Check_and_init_real(self, N_Cheb, N_four)
                                
  Class(mpi2dT_buffer_real), intent(InOut) :: self
  Integer, Intent(In) :: N_Cheb, N_Four

  ! >>> first check that the dimensions are OK
  if (modulo(N_cheb, world_size).ne.0) Then
     Print*, '>>> ERROR'
     Print*, 'The number of Chebyshev Pols. (including the 3/2'
     Print*, 'dealiasing factor) is not a multiple of the num-'
     Print*, 'ber of processes. Stopping now.'                      
     Stop 
  End If
  if (modulo(N_four, world_size).ne.0) Then
     Print*, '>>> ERROR'
     Print*, 'The number of  Fourier  modes (including the 3/2'
     Print*, 'dealiasing factor) is not a multiple of the num-'
     Print*, 'ber of processes. Stopping now.'                      
     Stop 
  End If
  self%NX = N_four
  self%NZ = N_cheb
  self%locX = N_four/world_size
  self%locZ = N_cheb/world_size
  self%block_size = self%locX * self%locZ

  Allocate(self%BN (self%block_size, world_size))
  Allocate(self%BT (self%block_size, world_size))
  Allocate(self%AN (self%NZ, self%locX))
  Allocate(self%AT (self%NX, self%locZ))

End Subroutine mpi2Dt_Check_and_init_real 

Subroutine mpi2Dt_Pack_N2T_Real(self)
  Class(mpi2dT_buffer_real), intent(InOut) :: self
  ! --- local variables ---
  Integer :: iblock, index_a, index_b
  Integer :: ix, iz

   self%BN(:,:) = 0._dp
   Do iblock = 1,world_size
      Do ix = 1, self%locX
         Do iz = 1, self%locZ
            index_b = (ix-1)*self%locZ + iz
            index_a = (iblock-1) * self%locZ + iz
            self%BN (index_b, iblock) = self%AN(index_a, ix)
         End Do
      End Do
   End Do
End Subroutine mpi2Dt_Pack_N2T_real


Subroutine mpi2Dt_UnPack_T2N_Real(self)
  Class(mpi2dT_buffer_real), intent(InOut) :: self
  ! --- local variables ---
  Integer :: iblock, index_a, index_b
  Integer :: ix, iz

   self%AN(:,:) = 0._dp
   Do iblock = 1,world_size
      Do ix = 1, self%locX
         Do iz = 1, self%locZ
            index_b = (ix-1)*self%locZ + iz
            index_a = (iblock-1) * self%locZ + iz
            self%AN(index_a, ix) = self%BN (index_b, iblock) 
         End Do
      End Do
   End Do
End Subroutine mpi2Dt_UnPack_T2N_real

Subroutine mpi2Dt_unPack_N2T_Real(self)
  Class(mpi2dT_buffer_real), intent(InOut) :: self
  ! --- local variables ---
  Integer :: iblock, index_a, index_b
  Integer :: ix, iz

   self%AT(:,:) = 0._dp
   Do iblock = 1,world_size
      Do iz = 1, self%locz
         Do ix = 1, self%locx
            index_b = (ix-1)*self%locZ + iz
            index_a = (iblock-1) * self%locX + ix
            self%AT(index_a, iz) = self%BT (index_b, iblock)
         End Do
      End Do
   End Do
End Subroutine mpi2Dt_unPack_N2T_Real

Subroutine mpi2Dt_Pack_T2N_Real(self)
  Class(mpi2dT_buffer_real), intent(InOut) :: self
  ! --- local variables ---
  Integer :: iblock, index_a, index_b
  Integer :: ix, iz

   self%BT(:,:) = 0._dp
   Do iblock = 1,world_size
      Do iz = 1, self%locz
         Do ix = 1, self%locx
            index_b = (ix-1)*self%locZ + iz
            index_a = (iblock-1) * self%locX + ix
            self%BT (index_b, iblock) = self%AT(index_a, iz) 
         End Do
      End Do
   End Do
End Subroutine mpi2Dt_Pack_T2N_Real


Subroutine mpi2Dt_transpose_N2T_Real(self)
  Class(mpi2dT_buffer_real), intent(InOut) :: self
   
   Call self%pack_N2T()
   Call MPI_AllToAll(self%BN, self%block_size, MPI_DOUBLE, &
                     self%BT, self%block_size, MPI_DOUBLE, MPI_COmm_world, ierr)
   Call self%unpack_N2T()
   
End Subroutine mpi2Dt_transpose_N2T_real

Subroutine mpi2Dt_transpose_T2N_real(self)
  Class(mpi2dT_buffer_real), intent(InOut) :: self

   Call self%pack_T2N()
   Call MPI_AllToAll(self%BT, self%block_size, MPI_DOUBLE, &
                     self%BN, self%block_size, MPI_DOUBLE, MPI_COmm_world, ierr)
   Call self%unpack_T2N()
   
End Subroutine mpi2Dt_transpose_T2N_real
End Module MPI_2dTranspose_real
