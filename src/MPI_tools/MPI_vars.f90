 module MPI_vars
 use fortran_kinds
 use MPI
 implicit none

 integer :: my_rank
 integer :: world_size
 integer :: ierr
 integer :: tag
 integer :: status(MPI_STATUS_SIZE)
 integer :: irank


 end module MPI_Vars
