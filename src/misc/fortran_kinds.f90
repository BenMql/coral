 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: Fortran_kinds
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel
 !
 ! DESCRIPTION
 !> Definition of some fortran kinds based upon the ISO C standard.
 !
 !=============================================================================

 Module Fortran_kinds
 Use, Intrinsic :: ISO_C_Binding
 Integer, Parameter :: dp = C_Double !< Double precision number (8 bytes)
 Integer, Parameter :: ci32 = C_INT32_T  !< 4-byte integer
 Integer, Parameter :: ci64 = C_INTPTR_T !< 8-byte integer

 contains

  subroutine benAbort()
    integer, allocatable :: a(:)
    allocate( a(1) ) 
    a(2) = 3
  end subroutine

 end Module Fortran_kinds
