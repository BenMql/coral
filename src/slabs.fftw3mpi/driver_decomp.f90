program driver_decomp
use fortran_kinds
use mpi_vars
use domain_decomposition
implicit none

 call MPI_init(ierr)
 call MPI_Comm_Size(MPI_Comm_world, world_size, ierr)
 call MPI_Comm_Rank(MPI_Comm_world, My_rank,    ierr) 
 call domain_decomp%init( 3, 4, 12)
print *, 'coucou'
End program driver_decomp
