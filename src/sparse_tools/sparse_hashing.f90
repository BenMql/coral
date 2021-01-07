 Module sparse_hashing
 Use double
 Use sparse_format
 Implicit None

 Contains
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! some hashing functions for debugging sparse matrices.
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Subroutine hash_z(A,val)
    Type(zcsr_matrix), Intent(In) :: A
    Real(dp), Intent(Out) :: val
    Integer :: k, ik
    val = 0._dp
    Do k = 1,A%nrow
      Do ik = A%row(k),(A%row(k+1)-1)
      val = val + cos(real(k,kind=dp)) *real(A%dat(ik)) * sin(real(A%col(ik), kind=dp))*abs(A%dat(ik))
      End Do
    End Do
 End Subroutine

 Subroutine hash_d(A,val)
    Type(dcsr_matrix), Intent(In) :: A
    Real(dp), Intent(Out) :: val
    Integer :: k, ik
    val = 0._dp
    Do k = 1,A%nrow
      Do ik = A%row(k),(A%row(k+1)-1)
      val = val + cos(real(k,kind=dp)) *real(A%dat(ik)) * sin(real(A%col(ik), kind=dp))*abs(A%dat(ik))
      End Do
    End Do
 End Subroutine
 End Module sparse_hashing
