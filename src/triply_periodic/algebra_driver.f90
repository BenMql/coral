program algebra_driver
  
  use P3_algebra
  implicit none

  integer :: N1 = 6
  integer :: NV = 2
  integer :: ivar, jeq
  complex(kind=dp), allocatable :: buf3d(:,:,:)

  type(operator_3d) :: a 

  call a%initialise(NV,N1, N1, N1)
  !////////////////////////////
  allocate (buf3d(n1, n1, n1))
  do ivar = 1,NV
  buf3d = ivar**2 
  call a%fill_in (buf3d, ivar, ivar)
  end do
  !////////////////////////////
  call a%invert()
  call a%check() 


  
end program algebra_driver
