Program test_idiot
  Use, Intrinsic :: Iso_C_binding
  character(len=10) :: val_str
  Real(C_Double) :: val

  Integer, Pointer :: i => null()

  
  val_str = trim('7.d0     ')
  Read(val_str,*) val
  write(*,*) val_str
  write(*,*) val
  write(*,*) associated(i)
End Program test_idiot

