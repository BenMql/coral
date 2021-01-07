Program initial_temperature_background

   Use, Intrinsic :: ISO_C_Binding
   Implicit None

   Integer, Parameter :: dp = C_double
   Real(C_double), Dimension(:), Allocatable :: z_grid
   Real(C_double), Dimension(:), Allocatable :: T_back
 
   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! SOME HARD CODED PARAMETERS
   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Real(C_double), Parameter :: layer_gap   =1.0_dp
   Real(C_double), Parameter :: layer_center=0.5_dp
   Real(C_double), Parameter :: eta = 10._dp         
   Character(len=44), Parameter :: output_dir='./'
   Integer, Parameter :: NZ=128
   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! END OF HARD CODED PARAMETERS
   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>



   Call Gauss_Lobatto_grid_1d(z_grid, NZ, layer_center, layer_gap)

   Allocate (T_back(NZ))

   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !>> ENTER BELOW THE FORMULA FOR THE BACKGROUND >>>>>
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !T_back =- z_grid**2/2._dp &
   !        - z_grid / (exp(-eta)-1._dp) &
   !        - exp(-eta*z_grid) / eta / (exp(-eta)-1._dp) 
   T_back =- z_grid**2/2._dp &
           + z_grid 
   !T_back = z_grid &
          !+ exp(-z_grid*15._dp)/15._dp & 
          !+ exp((z_grid-1.0_dp)*15._dp)/15._dp  
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !<< END OF THE BACKGROUND TEMPERATURE FORMULA <<<<<<
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   Call chdir(output_dir)
   Open (Unit=9, File="NZ.lut", Status='replace', Access='stream')
   Write(9) NZ  
   Close (Unit=9)
   Open (Unit=9, File="temperature_background.dat", Status='replace', Access='stream')
   Write(9) z_grid
   Write(9) T_back
   Close (Unit=9)

Contains

 Subroutine Gauss_Lobatto_grid_cos(GC_cos, N1)
   Integer, Intent(In) :: N1
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_cos
   Real(kind=dp), Parameter :: pi = 4._dp*atan(1._dp)
   Integer :: i
   Allocate(GC_cos(N1))
   Do i=1,N1
      GC_cos(i) = Cos(pi*Real(i-1, kind=dp)/Real(N1-1, kind=dp))
   End Do
 End Subroutine Gauss_Lobatto_grid_cos

 Subroutine Gauss_Lobatto_grid_1d(GC_grid1, N1, zcenter, zgap)
   Integer, Intent(In) :: N1
   Real(Kind=dp), Dimension(:), Allocatable, Intent(Out) :: GC_grid1
   Real(Kind=dp), Dimension(:), Allocatable :: GC_cos
   Real(Kind=dp), Intent(In) :: zcenter, zgap
   Allocate(GC_grid1(N1))
   Call Gauss_Lobatto_grid_cos(GC_cos, N1)
   GC_grid1 = zcenter + GC_cos*0.5_dp*zgap 
 End Subroutine Gauss_Lobatto_grid_1d

 End Program Initial_temperature_Background
