module chebyshev_galerkin_2
 use sparse_formats
 implicit none

contains

 Subroutine Chebyshev_Galerkin_unique_stencil( N, bc_type, S)
   
   Integer, Intent(In) :: N 
   Integer, intent(In) :: bc_type
   Type(dcsr_matrix), Intent(Out) :: S
   ! the stencil matrix is prepared in COO format, stored in:
   Real(dp), Dimension(:), Allocatable :: aux_dat
   Integer, Dimension(:), Allocatable :: aux_col, aux_row
   Integer :: j 
   Integer :: fill_in_index
   Real(dp) :: naux, maux
   Integer :: job(8)
   Integer :: info
   
   S%nrow = N 
   fill_in_index = 1

   !##########################################################
   ! bc_type = 00: No boundary condition
   !..........................................................
   Select Case (bc_type)
      Case(0)
         S%ncol  = N
         S%nelems = N
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
         ! build the identity
         Do j=1,N
            aux_dat(j) = 1._dp
            aux_row(j) = j
            aux_col(j) = j
         End do
   !..........................................................
   ! bc_type = 10: Dirichlet at +1, f=0
   !..........................................................
      Case(10)
         S%ncol = N-1
         S%nelems = 2*(N-1)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-1
               aux_dat(fill_in_index) =-1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               aux_dat(fill_in_index) = 1.0_dp
               aux_row(fill_in_index) = j+1
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   !..........................................................
   ! bc_type = 11: Dirichlet at -1, f=0
   !..........................................................
      Case(11)
         S%ncol  = N-1
         S%nelems = 2*(N-1)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-1
               aux_dat(fill_in_index) = 1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               aux_dat(fill_in_index) = 1.0_dp
               aux_row(fill_in_index) = j+1
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
         S%nelems = fill_in_index - 1  
   !..........................................................
   ! bc_type = 12: Neuman at +1, Df=0
   !..........................................................
      Case(12)
         S%ncol  = N-1
         S%nelems = 2*(N-1)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-1
               aux_dat(fill_in_index) = -(j+1._dp)**2
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               aux_dat(fill_in_index) = real(j, kind=dp)**2
               aux_row(fill_in_index) = j+1
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
         S%nelems = fill_in_index - 1  
   !..........................................................
   ! bc_type = 13: Neuman at -1, Df=0
   !..........................................................
      Case(13)
         S%ncol  = N-1
         S%nelems = 2*(N-1)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-1
               aux_dat(fill_in_index) =  (j+1._dp)**2
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               aux_dat(fill_in_index) = real(j, kind=dp)**2
               aux_row(fill_in_index) = j+1
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
         S%nelems = fill_in_index - 1  
   !..........................................................
   ! bc_type = 20: Both Dirichlet, f=0 
   !..........................................................
      Case(20)
         S%ncol  = N-2
         S%nelems = 2*(N-2)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-2
               aux_dat(fill_in_index) =-1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               aux_dat(fill_in_index) = 1.0_dp
               aux_row(fill_in_index) = j+2
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   !..........................................................
   ! bc_type = 21: Both Neuman, Df=0
   !..........................................................
      Case(21)
         S%ncol  = N-2
         S%nelems = 2*(N-2)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-2
               aux_dat(fill_in_index) =-1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               naux = Real(j+1, kind=dp)
               aux_dat(fill_in_index) = (naux-2.0_dp)**2/naux**2
               aux_row(fill_in_index) = j+2
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   !..........................................................
   ! bc_type = 22:  Top No-Slip, Bottom Stress-Free 
   !                Top Dirichlet, Bottom Neumann   
   !..........................................................
      Case(22)
         S%ncol  = N-2
         S%nelems = 3*(N-2)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-2
               naux = Real(j-1, kind=dp)
               ! diagonal of ones
               aux_dat(fill_in_index) = 1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal 
               aux_dat(fill_in_index) = - (4._dp*naux + 4._dp)/&
                           (2._dp*naux**2 + 6._dp*naux + 5._dp)
               aux_row(fill_in_index) = j+1
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subsubdiagonal 
               aux_dat(fill_in_index) = -(naux**2 + (naux+1._dp)**2)/&
                         ((naux+1._dp)**2 + (naux+2._dp)**2)
               aux_row(fill_in_index) = j+2
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   !..........................................................
   ! bc_type = 23:  Top Stress-Free, Bottom No-Slip 
   !                Top Neumann, Bottom Dirichlet   
   !..........................................................
      Case(23)
         S%ncol  = N-2
         S%nelems = 3*(N-2)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-2
               naux = Real(j-1, kind=dp)
               ! diagonal of ones
               aux_dat(fill_in_index) = 1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal 
               aux_dat(fill_in_index) =   (4._dp*naux + 4._dp)/&
                           (2._dp*naux**2 + 6._dp*naux + 5._dp)
               aux_row(fill_in_index) = j+1
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subsubdiagonal 
               aux_dat(fill_in_index) = -(naux**2 + (naux+1._dp)**2)/&
                         ((naux+1._dp)**2 + (naux+2._dp)**2)
               aux_row(fill_in_index) = j+2
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   !..........................................................
   ! bc_type = 40: Both no-slip (f=0 and Df=0)
   !..........................................................
      Case(40)
         S%ncol  = N-4
         S%nelems = 3*(N-4)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-4
               aux_dat(fill_in_index) =-1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               naux = real(j+3, kind=dp)
               aux_dat(fill_in_index) = (2._dp*naux-4._dp)/(naux-1._dp)
               aux_row(fill_in_index) = j+2
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               aux_dat(fill_in_index) =-(naux-3._dp)/(naux-1._dp)
               aux_row(fill_in_index) = j+4
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   !..........................................................
   ! bc_type = 41: Both Stress-free (f=0 and DDf=0)
   !..........................................................
      Case(41)
         S%ncol  = N-4
         S%nelems = 3*(N-4)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-4
               aux_dat(fill_in_index) =-1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               naux = Real(j+3, kind=dp)
               maux = 24._dp*naux/(2._dp*naux**2-4._dp*naux+3._dp)
               maux = maux - 18._dp/(naux-1._dp)
               aux_dat(fill_in_index) = maux +2._dp
               aux_row(fill_in_index) = j+2
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               aux_dat(fill_in_index) =-maux -1._dp
               aux_row(fill_in_index) = j+4
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   !..........................................................
   ! bc_type = 42:  Top No-Slip, Bottom Stress-Free 
   !                Top Dirichlet, Bottom Neumann   
   !..........................................................
      Case(42)
         S%ncol  = N-4
         S%nelems = 5*(N-4)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-4
               naux = Real(j-1, kind=dp)
               ! diagonal of ones
               aux_dat(fill_in_index) = 1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal j+1
               aux_dat(fill_in_index) = - (2._dp*naux + 2._dp) /&
                        (naux**2 + 5._dp*naux + 7._dp)
               aux_row(fill_in_index) = j+1
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal j+2
               aux_dat(fill_in_index) = - ( 2._dp*naux**3 +&
                                           12._dp*naux**2 +&
                                           28._dp*naux    +&
                                           24._dp)/&
                                          (       naux**3 +&
                                            8._dp*naux**2 +&
                                           22._dp*naux    +&
                                           21._dp)
               aux_row(fill_in_index) = j+2
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal j+3
               aux_dat(fill_in_index) =   (2._dp*naux + 2._dp) /&
                        (naux**2 + 5._dp*naux + 7._dp)
               aux_row(fill_in_index) = j+3
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal j+4
               aux_dat(fill_in_index) =   (       naux**3 +&
                                            4._dp*naux**2 +&
                                            6._dp*naux    +&
                                            3._dp)/&
                                          (       naux**3 +&
                                            8._dp*naux**2 +&
                                           22._dp*naux    +&
                                           21._dp)
               aux_row(fill_in_index) = j+4
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   !..........................................................
   ! bc_type = 43:  Top Stress-Free, Bottom No-Slip
   !..........................................................
      Case(43)
         S%ncol  = N-4
         S%nelems = 5*(N-4)
         Allocate(aux_dat(S%nelems))
         Allocate(aux_row(S%nelems))
         Allocate(aux_col(S%nelems))
            Do j = 1,N-4
               naux = Real(j-1, kind=dp)
               ! diagonal of ones
               aux_dat(fill_in_index) = 1.0_dp
               aux_row(fill_in_index) =  j
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal j+1
               aux_dat(fill_in_index) =   (2._dp*naux + 2._dp) /&
                        (naux**2 + 5._dp*naux + 7._dp)
               aux_row(fill_in_index) = j+1
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal j+2
               aux_dat(fill_in_index) = - ( 2._dp*naux**3 +&
                                           12._dp*naux**2 +&
                                           28._dp*naux    +&
                                           24._dp)/&
                                          (       naux**3 +&
                                            8._dp*naux**2 +&
                                           22._dp*naux    +&
                                           21._dp)
               aux_row(fill_in_index) = j+2
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal j+3
               aux_dat(fill_in_index) = - (2._dp*naux + 2._dp) /&
                        (naux**2 + 5._dp*naux + 7._dp)
               aux_row(fill_in_index) = j+3
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
               ! subdiagonal j+4
               aux_dat(fill_in_index) =   (       naux**3 +&
                                            4._dp*naux**2 +&
                                            6._dp*naux    +&
                                            3._dp)/&
                                          (       naux**3 +&
                                            8._dp*naux**2 +&
                                           22._dp*naux    +&
                                           21._dp)
               aux_row(fill_in_index) = j+4
               aux_col(fill_in_index) =  j
               fill_in_index = fill_in_index + 1 
            End do
   End Select
   !..........................................................
   ! Done with boundary condition            
   !##########################################################
      
   !
   !

   !##########################################################
   ! Convert aux to CSR format, and output
   !..........................................................
   ! set-up the mkl COO<->CSR conversion routine
   job(1) = 2 ! COO to CSR, sorted rows (increasing column)
   job(2) = 1 ! Fortran indexing, CSR
   job(3) = 1 ! Fortran indexing, COO
   job(4) = 0 ! ? unused here ?
   job(5) = 0 ! ? allegedly unused if job(1) <> 0
   job(6) = 0 ! fill-in the three arrays
   job(7) = 0 ! ? unused here ?
   job(8) = 0 ! ? unused here ?
   Allocate (S%dat( S%nelems ))
   Allocate (S%col( S%nelems ))
   Allocate (S%row(    N + 1    ))
   ! convert and clean behind us:
   CALL mkl_dcsrcoo(job, N, S%dat, S%col, S%row,&
              S%nelems,   aux_dat,  aux_row,  aux_col, info) 
   DEALLOCATE(aux_dat, aux_row, aux_col)
   !..........................................................
   ! Done with conversion
   !##########################################################

 end subroutine Chebyshev_Galerkin_unique_stencil

end module chebyshev_galerkin_2

