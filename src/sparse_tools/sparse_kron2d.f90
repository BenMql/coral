 Module Sparse_kron2d
 Use Double, Only: dp
 Use CSR_matrices
 implicit none

 Contains

 ! ================================================ !
 !                                                  !
 !         DOUBLE PRECISION ROUTINES                !
 !                                                  !
 ! ================================================ !

 Subroutine d_CSRkronCSR (Acsr, Bcsr, Ccsr)
   Type(dcsr_matrix), Intent(InOut) :: ACsr
   Type(dcsr_matrix), Intent(InOut) :: BCsr
   Type(dcsr_matrix), Intent(InOut) :: CCsr
   Type(dcoo_matrix), Intent(InOut) :: ACoo
   Type(dcoo_matrix), Intent(InOut) :: BCoo
   Type(dcoo_matrix), Intent(InOut) :: CCoo

   Call d_ordercsr(Acsr)
   Call d_ordercsr(Bcsr)
   Call d_csr2coo (Acsr, Acoo)
   Call d_csr2coo (Bcsr, Bcoo)
   Call d_cookroncoo(Acoo, Bcoo, Ccoo)
   Call d_coo2csr (CCoo, CCsr)

 End Subroutine d_CSRkronCSR

 Subroutine d_COOkronCOO (M1, M2, M3)
   Type(dcoo_matrix), Intent(In) :: M1
   Type(dcoo_matrix), Intent(In) :: M2
   Type(dcoo_matrix), Intent(Out):: M3
   Integer :: compteur, ila, ilb
   Integer :: Nnz1, NNZ2, neq2
   NNZ1 = M1%nelems
   NNZ2 = M2%nelems
   Neq2 = M2%nrow 
   Allocate (M3%dat(NNZ2*NNZ1))
   Allocate (M3%row(NNZ2*NNZ1))
   Allocate (M3%col(NNZ2*NNZ1))
   compteur = 0
   Do ila = 1, NNZ1
   Do ilb = 1, NNZ2
       compteur = compteur + 1
       M3%dat (compteur) = M1%dat (ila) * M2%dat (ilb)           
       M3%row (compteur) =(M1%row (ila)-1) * Neq2 + M2%row(ilb)
       M3%col (compteur) =(M1%col (ila)-1) * Neq2 + M2%col(ilb)
  End Do 
  End Do

 End Subroutine d_COOkronCOO



 ! ================================================ !
 !                                                  !
 !           DOUBLE COMPLEX ROUTINES                !
 !                                                  !
 ! ================================================ !


 Subroutine z_CSRkronCSR (Acsr, Bcsr, Ccsr)
   Type(zcsr_matrix), Intent(InOut) :: ACsr
   Type(zcsr_matrix), Intent(InOut) :: BCsr
   Type(zcsr_matrix), Intent(InOut) :: CCsr
   Type(zcoo_matrix), Intent(InOut) :: ACoo
   Type(zcoo_matrix), Intent(InOut) :: BCoo
   Type(zcoo_matrix), Intent(InOut) :: CCoo

   Call z_ordercsr(Acsr)
   Call z_ordercsr(Bcsr)
   Call z_csr2coo (Acsr, Acoo)
   Call z_csr2coo (Bcsr, Bcoo)
   Call z_cookroncoo(Acoo, Bcoo, Ccoo)
   Call z_coo2csr (CCoo, CCsr)

 End Subroutine z_CSRkronCSR

 Subroutine z_COOkronCOO (M1, M2, M3)
   Type(zcoo_matrix), Intent(In) :: M1
   Type(zcoo_matrix), Intent(In) :: M2
   Type(zcoo_matrix), Intent(Out):: M3
   Integer :: compteur, ila, ilb
   Integer :: Nnz1, NNZ2, neq2
   NNZ1 = M1%nelems
   NNZ2 = M2%nelems
   Neq2 = M2%nrow 
   Allocate (M3%dat(NNZ2*NNZ1))
   Allocate (M3%row(NNZ2*NNZ1))
   Allocate (M3%col(NNZ2*NNZ1))
   compteur = 0
   Do ila = 1, NNZ1
   Do ilb = 1, NNZ2
       compteur = compteur + 1
       M3%dat (compteur) = M1%dat (ila) * M2%dat (ilb)           
       M3%row (compteur) =(M1%row (ila)-1) * Neq2 + M2%row(ilb)
       M3%col (compteur) =(M1%col (ila)-1) * Neq2 + M2%col(ilb)
  End Do 
  End Do

 End Subroutine z_COOkronCOO



 End Module Sparse_kron2d
