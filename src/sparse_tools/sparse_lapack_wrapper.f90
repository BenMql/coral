 Module Sparse_LAPACK_wrapper
 Use Sparse_formats
 Use Lapack95
 Implicit None

 Contains

 Subroutine wrap_LUfact_band_z(A, verb)
   Type(zdia_and_lu), Intent(InOut) :: A
   Logical, Intent(In) :: verb
   Integer :: info
   if (allocated(A%piv)) DeAllocate(A%piv)
   Allocate(A%piv(min(A%dia%nrow, A%dia%ncol)))
   Call zgbtrf(A%dia%nrow, A%dia%ncol, A%dia%nl, A%dia%nu, A%lu, 2*A%dia%nl+A%dia%nu+1, A%piv, info)
   if (info.ne.0) print *, "Problem during LU routine zgbtrf"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 End Subroutine wrap_LUfact_band_z


 Subroutine wrap_LUfact_band_d(A, verb)
   Type(ddia_and_lu), Intent(InOut) :: A
   Logical, Intent(In) :: verb
   Integer :: info
   if (allocated(A%piv)) DeAllocate(A%piv)
   Allocate(A%piv(min(A%dia%nrow, A%dia%ncol)))
   Call dgbtrf(A%dia%nrow, A%dia%ncol, A%dia%nl, A%dia%nu, A%lu, 2*A%dia%nl+A%dia%nu+1, A%piv, info)
   if (info.ne.0) print *, "Problem during LU routine dgbtrf"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 End Subroutine wrap_LUfact_band_d

 Subroutine wrap_LUsolve_band_z(A, b, nrhs, ldb, verb)
   Type(zdia_and_lu), Intent(In) :: A
   Complex(dp), Pointer :: b(:)
   Integer, Intent(In) :: nrhs
   Integer, Intent(In) :: ldb 
   Logical, Intent(In) :: verb
   Integer :: info
   Call zgbtrs('N', A%dia%ncol, A%dia%nl, A%dia%nu, nrhs, A%LU, &
                  2*A%dia%nl+A%dia%nu+1, A%piv, b, ldb, info)
   if (info.ne.0) print *, "Problem during LU routine zgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 End Subroutine wrap_LUsolve_band_z
 
 Subroutine wrap_LUsolve_band_d(A, b, nrhs, ldb, verb)
   Type(ddia_and_lu), Intent(In) :: A
   Real(dp), Pointer :: b(:)
   Integer, Intent(In) :: nrhs
   Integer, Intent(In) :: ldb 
   Logical, Intent(In) :: verb
   Integer :: info
   Call dgbtrs('N', A%dia%ncol, A%dia%nl, A%dia%nu, nrhs, A%LU, &
                  2*A%dia%nl+A%dia%nu+1, A%piv, b, ldb, info)
   if (info.ne.0) print *, "Problem during LU routine dgbtrs"
   if (info.ne.0) print *, "Error Diagnostic:", info
   if (info.ne.0) stop
 End Subroutine wrap_LUsolve_band_d
 
 End Module Sparse_LAPACK_wrapper
