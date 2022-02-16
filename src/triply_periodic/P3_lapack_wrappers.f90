 Module P3_lapack_wrappers
   use Fortran_kinds
   use lapack95     
   Implicit None

   interface wrap_LUfactor_square
     module procedure wrap_LUfactor_square_z
     module procedure wrap_LUfactor_square_d
   end interface wrap_LUfactor_square

   contains

   Subroutine wrap_LUinvert_z(aFull, verb)
     Complex(dp), dimension(:,:), allocatable, intent(InOut) :: aFull
     logical, intent(In) :: verb
     Integer, dimension(:), allocatable :: aPiv
     integer :: aSize
     Complex(dp), dimension(:), allocatable :: work
     integer :: info

     if (size(aFull,1).ne.size(aFull,2)) STOP "nonSquare matrix!"
     aSize = size(aFull,1)
     allocate(work(2*aSize))
     call wrap_LUfactor_square_z(aFull, aPiv, verb)
     call zgetri(aSize, aFull, aSize, aPiv, work, 2*aSize, info)
     if (info.ne.0) STOP "INV problem"
     
   End Subroutine wrap_LUinvert_z

   Subroutine wrap_LUfactor_square_z(aFull, aPiv, verb)
     Complex(dp), dimension(:,:), allocatable, intent(InOut) :: aFull
     Logical, intent(In) :: verb
     Integer, dimension(:), allocatable, intent(out) :: aPiv
     integer :: aSize
     integer :: info

     if (size(aFull,1).ne.size(aFull,2)) STOP "nonSquare matrix!"

     aSize = size(aFull,1)
     allocate (aPiv(aSize))             

     call zgetrf(aSize, aSize, aFull, aSize, aPiv, info)

     if (info.ne.0) STOP "LU problem"

   End Subroutine wrap_LUfactor_square_z
     
   Subroutine wrap_LUfactor_square_d(aFull, aPiv, verb)
     Real(dp), dimension(:,:), allocatable, intent(InOut) :: aFull
     Integer, dimension(:), allocatable, intent(out) :: aPiv
     integer :: aSize
     Logical, intent(In) :: verb
     integer :: info

     if (size(aFull,1).ne.size(aFull,2)) STOP "nonSquare matrix!"

     aSize = size(aFull,1)
     allocate (aPiv(aSize))             

     call dgetrf(aSize, aSize, aFull, aSize, aPiv, info)

     if (info.ne.0) STOP "LU problem"

   End Subroutine wrap_LUfactor_square_d
     
   

 End Module P3_lapack_wrappers
