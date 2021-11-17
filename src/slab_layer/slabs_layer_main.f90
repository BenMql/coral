program slab_layer_main
use SL_geometry
use SL_algebra
implicit none


call mpi_init(ierr)
call MPI_Comm_Size(MPI_Comm_world, world_size, ierr)
call MPI_Comm_Rank(MPI_Comm_world, My_rank,    ierr)
call fftw_MPI_init()
call domain_decomp%decomp_init(128, 128, 100)
if (my_rank.eq.0) print *, 'my_rank, phys_offset, spec_offset, localSize_phys, localSize_spec'
print *, my_rank,&
         domain_decomp%phys_offset, & 
         domain_decomp%spec_offset, &
         domain_decomp%localSize_phys,&
         domain_decomp%localSize_spec

call geom%init()

call MPI_Finalize( ierr)


end program slab_layer_main
