program mpi_2ddriver
Use mpi_2dtranspose
!Type :: mpi2dT_buffer_cplx
   !Complex(dp), Allocatable, Dimension(:,:) :: AN
   !Complex(dp), Allocatable, Dimension(:,:) :: AT
   !Complex(dp), Allocatable, Dimension(:,:) :: BN 
   !Complex(dp), Allocatable, Dimension(:,:) :: BT 
   !Integer :: NZ, NX, locZ, locX, block_size
   !Character(len=2) :: state
   !Contains
   !Procedure :: init =>mpi2Dt_Check_and_init_cplx
   !Procedure :: pack_dat =>mpi2Dt_Pack_cplx
   !Procedure :: unpack_dat =>mpi2Dt_unPack_cplx
   !Procedure :: Forward_transpose =>mpi2Dt_transpose_cplx
!End Type mpi2dT_buffer_cplx


Implicit NOne
Type(mpi2dT_buffer_cplx) :: mpiBuf
Integer, Parameter :: NZ = 512
Integer, Parameter :: NX = 768
Complex(dp), Allocatable, Dimension(:,:) :: backup_a
Real(dP) :: relerror

 Call MPI_init(ierr)
 Call MPI_Comm_Size(MPI_Comm_world, world_size, ierr) 
 Call MPI_Comm_Rank(MPI_Comm_world, My_rank,    ierr) 

 Call mpiBuf%init(NZ, NX)
 Call mpi2dt_fill_in_dummy(mpiBuf)
 Allocate (Backup_a(NZ, mpiBuf%locX))
 Backup_a=mpiBuf%AN
 Print *, 'for core', my_rank,', local NZ =',mpibuf%LocZ,', local NX =',mpibuf%locX
! Print *, 'for core', my_rank,', self%AN=',abs(mpibuf%AN)                            

 Call mpiBuf%forward_transpose()
! Print *, 'for core', my_rank,', self%AT=',abs(mpibuf%AT)                            
 mpiBUF%AN = cmplx(0._dp, 0._dp, kind=dp)
 mpiBUF%BT = cmplx(0._dp, 0._dp, kind=dp)
 mpiBUF%BN = cmplx(0._dp, 0._dp, kind=dp)
 Call mpiBuf%inverse_transpose()
! Print *, 'for core', my_rank,', self%AN=',abs(mpibuf%AN)                            
 Call MPI_barrier(MPI_Comm_world, ierr)
 relerror = sum(abs(backup_a-mpibuf%AN))/sum(abs(backup_a))
 Print *, 'for core', my_rank,', error =', relerror
 

 Call MPI_finalize(ierr)

End program mpi_2ddriver
