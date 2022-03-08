 Module PL_Cheby_misc   
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


 Subroutine Gauss_Chebyshev_weight_3d(GC_weight3, GC_weight1, N1, NY, NX,zgap)
   Integer, Intent(In) :: N1, NY, NX
   Real(Kind=dp), Dimension(:,:,:), Allocatable, Intent(Out) :: GC_weight3
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_weight1
   Real(Kind=dp), Intent(In) ::  zgap
   Integer :: i,j,k
   Call Gauss_Chebyshev_weight_1d(GC_weight1, N1, zgap)
   Allocate(GC_weight3(N1, NX, NY))
   Do k=1,NY
   Do j=1,NX
   Do i=1,N1
      GC_weight3(i,j,k) = GC_weight1(i)
   End do
   End do
   End do
 End Subroutine Gauss_Chebyshev_weight_3d



 Subroutine Gauss_Chebyshev_grid_3d(GC_grid3, GC_grid1, N1, NY, NX, zcenter, zgap)
   Integer, Intent(In) :: N1, NY, NX
   Real(Kind=dp), Dimension(:,:,:), Allocatable, Intent(Out) :: GC_grid3
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_grid1
   Real(Kind=dp), Intent(In) :: zcenter, zgap
   Integer :: i,j,k
   Call Gauss_Chebyshev_grid_1d(GC_grid1,N1, zcenter, zgap)
   Allocate(GC_grid3(N1, NX, NY))
   Do k=1,NY
   Do j=1,NX
   Do i=1,N1
      GC_grid3(i,j,k) = GC_grid1(i)
   End do
   End do
   End do
 End Subroutine Gauss_Chebyshev_grid_3d
  

  
 End Module PL_Cheby_misc   
