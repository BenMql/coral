
Program timesteppers_test

   use ts_test_classes

   Type(timestepping_operators) :: to
   type(buffers_for_solution) :: b4sol
   real(dp), dimension(:), allocatable :: f
   Integer, parameter :: NT = 10000
   real(dp), parameter :: dt = 0.01
   Type(scheme_handle) :: shandle
   Integer :: it


   call shandle%init('ARS_111', .True.)

   allocate(f(NT))

   to%mass%val = 10._dp
   to%stif%val = 2._dp
   call to%factorize(dt, shandle)
   call b4sol%alloc(to, shandle)

   b4sol%field = 0.2_dp

   do it = 1, NT
     call b4sol%march_forward(to, shandle, dt, .False.)
     f(it) = b4sol%field
   end do

   print *, f
   
   Open (Unit=9, File='out.dat', Status='replace', Access='stream')
   Write(9) f                              
   Close (Unit=9)
   


end Program timesteppers_test

