 Module R_Cheby_misc   
 Use fortran_kinds
 Implicit None

 Contains


 Subroutine Gauss_Chebyshev_grid_cos(GC_cos, N1)
   Integer, Intent(In) :: N1
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_cos
   Real(kind=dp), Parameter :: pi = 4._dp*atan(1._dp)
   Integer :: i
   Allocate(GC_cos(N1))
   Do i=1,N1
      GC_cos(i) = Cos(pi*Real(2*i-1, kind=dp)/Real(2*N1, kind=dp))
   End Do
 End Subroutine Gauss_Chebyshev_grid_cos

 Subroutine Gauss_Chebyshev_grid_1d(GC_grid1, N1, zcenter, zgap)
   Integer, Intent(In) :: N1
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_grid1
   Real(Kind=dp), Dimension(:), Allocatable :: GC_cos
   Real(Kind=dp), Intent(In) :: zcenter, zgap
   Allocate(GC_grid1(N1))
   Call Gauss_Chebyshev_grid_cos(GC_cos, N1)
   GC_grid1 = zcenter + GC_cos*0.5_dp*zgap 
 End Subroutine Gauss_Chebyshev_grid_1d

 Subroutine Gauss_Chebyshev_weight_1d(GC_weight1, N1, zgap)
   Integer, Intent(In) :: N1
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_weight1
   Real(Kind=dp), Dimension(:), Allocatable :: GC_grid1
   Real(Kind=dp), Dimension(:), Allocatable :: augmented_midpoints
   Real(Kind=dp), intent(in) :: zgap
   Integer :: i
   Allocate(GC_weight1(N1))
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !compute a Gauss-Chebyshev grid between -1 and 1
   Call Gauss_Chebyshev_grid_1d(GC_grid1,N1, 0._dp, 2._dp)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !find the midpoints, including the endpoints -1 and 1
   Allocate(augmented_midpoints(N1+1))
   augmented_midpoints(1)    = 1._dp
   augmented_midpoints(N1+1) =-1._dp
   Do i = 1,(N1-1)
      augmented_midpoints(i+1) = 0.5_dp*(GC_grid1(i+1)+GC_grid1(i))
   End Do
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !find the distance between midpoints and normalize:
   Do i = 1,N1
      GC_weight1(i) = augmented_midpoints(i)-augmented_midpoints(i+1)
   End Do
   GC_weight1 = GC_weight1 /2._dp*zgap
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 End Subroutine Gauss_Chebyshev_weight_1d


 Subroutine Gauss_Chebyshev_weight_2d(GC_weight2, GC_weight1, N1, NX,zgap)
   Integer, Intent(In) :: N1, NX
   Real(Kind=dp), Dimension(:,:), Allocatable, Intent(Out) :: GC_weight2
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_weight1
   Real(Kind=dp), Intent(In) ::  zgap
   Integer :: i,j
   Call Gauss_Chebyshev_weight_1d(GC_weight1, N1, zgap)
   !!! FIX ME FIX ME FIX ME Should be transposed?????
   !!! FIXME FIXME FIXME    Should be transposed?????
   !!! FIX ME FIX ME FIX ME Should be transposed?????
   !!! FIXME FIXME FIXME    Should be transposed?????
   Allocate(GC_weight2(N1, NX))
   Do j=1,NX
   Do i=1,N1
      GC_weight2(i,j) = GC_weight1(i)
   End do
   End do
 End Subroutine Gauss_Chebyshev_weight_2d



 Subroutine Gauss_Chebyshev_grid_2d(GC_grid2, GC_grid1, N1, NX, zcenter, zgap)
   Integer, Intent(In) :: N1, NX
   Real(Kind=dp), Dimension(:,:), Allocatable, Intent(Out) :: GC_grid2
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_grid1
   Real(Kind=dp), Intent(In) :: zcenter, zgap
   Integer :: i,j
   Call Gauss_Chebyshev_grid_1d(GC_grid1,N1, zcenter, zgap)
   !!! FIX ME FIX ME FIX ME Should be transposed?????
   !!! FIXME FIXME FIXME    Should be transposed?????
   !!! FIX ME FIX ME FIX ME Should be transposed?????
   !!! FIXME FIXME FIXME    Should be transposed?????
   Allocate(GC_grid2(N1, NX))
   Do j=1,NX
   Do i=1,N1
      GC_grid2(i,j) = GC_grid1(i)
   End do
   End do
 End Subroutine Gauss_Chebyshev_grid_2d
  

  
 End Module R_Cheby_misc   
