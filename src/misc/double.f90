 !=============================================================================
 !                            C O R A L
 !=============================================================================
 !
 ! MODULE: Double           
 ! 
 !> @author
 !> Ben Miquel, www.normalesup.org/~benmiquel
 !
 ! DESCRIPTION
 !> Definition of some fortran kinds based upon the ISO C standard.
 !
 !=============================================================================

 Module double          
 Use, Intrinsic :: ISO_C_Binding
 Integer, Parameter :: dp = C_Double !< Double precision number (8 bytes)
 Integer, Parameter :: ci32 = C_INT32_T !< 4-bytes integer
 end Module double
