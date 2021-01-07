 Module IMEX_schemes
 use fortran_kinds, only: dp
 Implicit None

 type :: scheme_handle
   Integer :: imp_stages
   Integer :: exp_stages
   Integer :: accuracy
   Real(dp), Dimension(:,:), Allocatable :: A_arr
   Real(dp), Dimension(:),   Allocatable :: b_arr
   Real(dp), Dimension(:),   Allocatable :: c_arr
   Real(dp), Dimension(:,:), Allocatable :: A_hat
   Real(dp), Dimension(:),   Allocatable :: b_hat
   Real(dp), Dimension(:),   Allocatable :: c_hat
  Contains
   Procedure :: init => initialise_IMEX_scheme
 end type

 Contains

 Subroutine initialise_IMEX_scheme(self, scheme_name, verb)
   class(scheme_handle), intent(out) :: self
   character(len=7), intent(in) :: scheme_name
   logical, intent(in) :: verb
   
   select case(scheme_name)
    Case ('ARS_111')
    Call Define_Forward_Backward_Euler_111(self, verb)
    Case ('ARS_122')
    Call Define_Implicit_Explicit_Midpoint_122(self, verb)
    Case ('ARS_233')
    Call Define_233(self, verb)
    Case ('ARS_232')
    Call Define_Lstable_232(self, verb)
    Case ('ARS_222')
    Call Define_Lstable_222(self, verb)
    Case ('ARS_343')
    Call Define_Lstable_343(self, verb)
    Case ('ARS_443')
    Call Define_Lstable_443(self, verb)
    Case Default
    if (verb) print *, "Invalid time-stepping scheme code."
    if (verb) print *, "Please modify the corresponding line in TFP_parameters.in"
    if (verb) print *, "Valid 7-character strings of the form ''ARS_???'':"
    if (verb) print *, "    ??? = 111, 122, 233, 232, 222, 343, 443"
    if (verb) print *, "..."
    if (verb) print *, "Stopping now."
    Stop
 End Select
 end subroutine

   
 Subroutine Define_Forward_Backward_Euler_111(shandle, verb)
   Type(Scheme_handle), Intent(Out) :: shandle 
   logical, intent(In) :: verb
   if (verb) then
   write (*,*) '======================================================'
   write (*,*) '>>> IMEX SCHEME:'
   write (*,*) '... First order Backward-Euler.'         
   write (*,*) '======================================================'
   write (*,*) New_line('a')
   end if
   shandle%imp_stages = 1
   shandle%exp_stages = 1
   shandle%accuracy   = 1
   Allocate( shandle%A_arr(1,1) )
   Allocate( shandle%B_arr(1)   )
   Allocate( shandle%C_arr(1)   )
   shandle%A_arr(1,1) = 1._dp
   shandle%B_arr(1)   = 1._dp
   shandle%C_arr(1)   = 1._dp
   Allocate( shandle%A_hat(2,2) )
   Allocate( shandle%B_hat(2)   )
   Allocate( shandle%C_hat(2)   )
   shandle%A_hat = 0._dp
   shandle%B_hat = 0._dp
   shandle%C_hat = 0._dp
   shandle%A_hat(2,1) = 1._dp
   shandle%B_hat(1)   = 1._dp
   shandle%C_hat(1)   = 1._dp
 End Subroutine Define_Forward_Backward_Euler_111

 Subroutine Define_Implicit_Explicit_Midpoint_122(shandle,verb)
   Logical, intent(in) :: verb
   Type(Scheme_handle), Intent(Out) :: shandle 
   if (verb) then
   write (*,*) '======================================================'
   write (*,*) '>>> IMEX SCHEME:'
   write (*,*) '... Ascher-Ruuth-Spiteri, Runge-kutta 2'
   Write (*,*) '... Stages: 1 implicits, 2 explicits.'
   write (*,*) '======================================================'
   write (*,*) New_line('a')
   end if
   shandle%imp_stages = 1
   shandle%exp_stages = 2
   shandle%accuracy   = 2
   Allocate( shandle%A_arr(1,1) )
   Allocate( shandle%B_arr(1)   )
   Allocate( shandle%C_arr(1)   )
   shandle%A_arr(1,1) = 0.5_dp
   shandle%B_arr(1)   = 1.0_dp
   shandle%C_arr(1)   = 0.5_dp
   Allocate( shandle%A_hat(2,2) )
   Allocate( shandle%B_hat(2)   )
   Allocate( shandle%C_hat(2)   )
   shandle%A_hat = 0._dp
   shandle%B_hat = 0._dp
   shandle%C_hat = 0._dp
   shandle%A_hat(2,1) = 0.5_dp
   shandle%C_hat(2)   = 0.5_dp
   shandle%B_hat(2)   = 1.0_dp
 End Subroutine
 
 Subroutine Define_233(shandle,verb)
   Logical, intent(in) :: verb
   Type(Scheme_handle), Intent(Out) :: shandle 
   Real(dp) :: gam
   if (verb) then
   write (*,*) '======================================================'
   write (*,*) '>>> IMEX SCHEME:'
   write (*,*) '... Ascher-Ruuth-Spiteri, Runge-kutta 3'
   Write (*,*) '... Stages: 2 implicits, 3 explicits.'
   write (*,*) '======================================================'
   write (*,*) New_line('a')
   end if
   gam = (3._dp+sqrt(3._dp))/6._dp
   shandle%imp_stages = 2
   shandle%exp_stages = 3
   shandle%accuracy   = 3
   Allocate( shandle%A_arr(2,2) )
   Allocate( shandle%B_arr(2)   )
   Allocate( shandle%C_arr(2)   )
   shandle%A_arr = 0._dp
   shandle%B_arr = 0._dp
   shandle%C_arr = 0._dp
   shandle%A_arr(1,1) = gam  
   shandle%A_arr(2,2) = gam  
   shandle%A_arr(2,1) = 1.0_dp-2._dp*gam 
   shandle%B_arr(1)   = 0.5_dp
   shandle%B_arr(2)   = 0.5_dp
   shandle%C_arr(1)   = gam
   shandle%C_arr(2)   = 1.0_dp-gam
   Allocate( shandle%A_hat(3,3) )
   Allocate( shandle%B_hat(3)   )
   Allocate( shandle%C_hat(3)   )
   shandle%A_hat = 0._dp
   shandle%B_hat = 0._dp
   shandle%C_hat = 0._dp
   shandle%A_hat(2,1) = gam  
   shandle%A_hat(3,1) = gam  -1._dp
   shandle%A_hat(3,2) = 2._dp*(1.0_dp-gam)
   shandle%B_hat(1)   = 0.0_dp
   shandle%B_hat(2)   = 0.5_dp
   shandle%B_hat(3)   = 0.5_dp
   shandle%C_hat(1)   = 0._dp
   shandle%C_hat(2)   = gam
   shandle%C_hat(3)   = 1.0_dp-gam
 End Subroutine Define_233


 Subroutine Define_Lstable_232(shandle,verb)

   Logical, intent(in) :: verb
   Type(Scheme_handle), Intent(Out) :: shandle 
   Real(dp) :: gam
   Real(dp) :: delta
   if (verb) then
   write (*,*) '======================================================'
   write (*,*) '>>> IMEX SCHEME:'
   write (*,*) '... Ascher-Ruuth-Spiteri, Runge-kutta 2, L-Stable'
   Write (*,*) '... Stages: 2 implicits, 3 explicits.'
   write (*,*) '======================================================'
   write (*,*) New_line('a')
   end if
   delta = -2._dp*sqrt(2.0_dp)/3._dp
   gam = 1._dp-1._dp/sqrt(2.d0)   
   shandle%imp_stages = 2
   shandle%exp_stages = 3
   shandle%accuracy   = 2
   Allocate( shandle%A_arr(2,2) )
   Allocate( shandle%B_arr(2)   )
   Allocate( shandle%C_arr(2)   )
   shandle%A_arr = 0._dp
   shandle%B_arr = 0._dp
   shandle%C_arr = 0._dp
   shandle%A_arr(1,1) = gam  
   shandle%A_arr(2,2) = gam  
   shandle%A_arr(2,1) = 1.0_dp - gam 
   shandle%B_arr(1)   = 1.0_dp - gam
   shandle%B_arr(2)   = gam   
   shandle%C_arr(1)   = gam
   shandle%C_arr(2)   = 1.0_dp
   Allocate( shandle%A_hat(3,3) )
   Allocate( shandle%B_hat(3)   )
   Allocate( shandle%C_hat(3)   )
   shandle%A_hat = 0._dp
   shandle%B_hat = 0._dp
   shandle%C_hat = 0._dp
   shandle%A_hat(2,1) = gam  
   shandle%A_hat(3,1) = delta
   shandle%A_hat(3,2) = 1._dp-delta
   shandle%B_hat(1)   = 0.0_dp 
   shandle%B_hat(2)   = 1._dp-gam
   shandle%B_hat(3)   = gam
   shandle%C_hat(1)   = 0._dp
   shandle%C_hat(2)   = gam
   shandle%C_hat(3)   = 1.0_dp
 End Subroutine Define_Lstable_232

 Subroutine Define_Lstable_222(shandle,verb)
   Logical, intent(in) :: verb
   Type(Scheme_handle), Intent(Out) :: shandle 
   Real(8) :: gam
   Real(8) :: delta
   if (verb) then
   write (*,*) '======================================================'
   write (*,*) '>>> IMEX SCHEME:'
   write (*,*) '... Ascher-Ruuth-Spiteri, Runge-kutta 2, L-Stable'
   Write (*,*) '... Stages: 2 implicits, 2 explicits.'
   write (*,*) '======================================================'
   write (*,*) New_line('a')
   end if
   gam = 1._dp-1._dp/sqrt(2._dp)   
   delta = 1.0_dp-1._dp/(2._dp*gam)
   shandle%imp_stages = 2
   shandle%exp_stages = 2
   shandle%accuracy   = 2
   Allocate( shandle%A_arr(2,2) )
   Allocate( shandle%B_arr(2)   )
   Allocate( shandle%C_arr(2)   )
   shandle%A_arr = 0._dp
   shandle%B_arr = 0._dp
   shandle%C_arr = 0._dp
   shandle%A_arr(1,1) = gam
   shandle%A_arr(2,2) = gam
   shandle%A_arr(2,1) = 1._dp-gam
   shandle%c_arr(1)   = gam
   shandle%c_arr(2)   = 1._dp
   shandle%b_arr(1)   = 1._dp-gam
   shandle%b_arr(2)   = gam
   Allocate( shandle%A_hat(3,3) )
   Allocate( shandle%B_hat(3)   )
   Allocate( shandle%C_hat(3)   )
   shandle%A_hat = 0._dp
   shandle%B_hat = 0._dp
   shandle%C_hat = 0._dp
   shandle%A_hat(2,1) = gam
   shandle%A_hat(3,1) = delta
   shandle%A_hat(3,2) = 1._dp-delta
   shandle%b_hat(1) = delta
   shandle%b_hat(2) = 1._dp-delta
   shandle%c_hat(2) = gam
   shandle%c_hat(3) = 1._dp
 End Subroutine Define_Lstable_222

 Subroutine Define_Lstable_343(shandle, verb)
   Logical, intent(in) :: verb
   Type(Scheme_handle), Intent(Out) :: shandle 
   Real(8) :: gam
   Real(8) :: upg
   Real(8) :: umg
   Real(8) :: b1
   Real(8) :: b2
   Real(8) :: a42hat
   if (verb) then
   write (*,*) '======================================================'
   write (*,*) '>>> IMEX SCHEME:'
   write (*,*) '... Ascher-Ruuth-Spiteri, Runge-kutta 3, L-Stable'
   Write (*,*) '... Stages: 3 implicits, 4 explicits.'
   Write (*,*) '###  ~~~ WARNING: ...................... '
   Write (*,*) '###  ~~~ WARNING: The selected scheme is '
   Write (*,*) '###  ~~~ WARNING: single precision only! '
   Write (*,*) '###  ~~~ WARNING: ...................... '
   write (*,*) '======================================================'
   write (*,*) New_line('a')
   end if
   gam = 0.43586652150845895
   upg = 0.71793326075422947
   umg = 0.28206673924577053
   b1  = 1.2084966491760099
   b2  = -0.64436317068446880 
   a42hat = 0.5529291479d0
   shandle%imp_stages = 3
   shandle%exp_stages = 4
   shandle%accuracy   = 3
   Allocate( shandle%A_arr(3,3) )
   Allocate( shandle%B_arr(3)   )
   Allocate( shandle%C_arr(3)   )
   shandle%A_arr = 0.d0
   shandle%B_arr = 0.d0
   shandle%C_arr = 0.d0
   shandle%A_arr(1,1) = gam  
   shandle%A_arr(2,2) = gam  
   shandle%A_arr(3,3) = gam  
   shandle%A_arr(2,1) = umg          
   shandle%A_arr(3,1) = b1          
   shandle%A_arr(3,2) = b2          
   shandle%B_arr(1)   = b1
   shandle%B_arr(2)   = b2
   shandle%B_arr(3)   = gam
   shandle%C_arr(1)   = gam
   shandle%C_arr(2)   = upg
   shandle%C_arr(3)   = 1.d0
   Allocate( shandle%A_hat(4,4) )
   Allocate( shandle%B_hat(4)   )
   Allocate( shandle%C_hat(4)   )
   shandle%A_hat = 0.d0
   shandle%B_hat = 0.d0
   shandle%C_hat = 0.d0
   shandle%A_hat(2,1) = gam  
   shandle%A_hat(3,1) = 0.3212788860
   shandle%A_hat(3,1) = (15.d0/4.d0 - 15.d0*gam +21.d0/4.d0*gam**2)*a42hat &
                        - 7.d0/2.d0+13.d0*gam - 4.5d0*gam**2
   print *, 'a31 :', shandle%A_hat(3,1) 
   shandle%A_hat(3,2) = 0.3966543747
   shandle%A_hat(3,2) =-(15.d0/4.d0 - 15.d0*gam +21.d0/4.d0*gam**2)*a42hat &
                        + 4.d0 - 12.5d0  *gam + 4.5d0*gam**2
   print *, 'a32 :', shandle%A_hat(3,2) 
   shandle%A_hat(4,1) =-0.105858296 
   shandle%A_hat(4,1) =1.0d0 - a42hat - a42hat                         
   print *, 'a41 :', shandle%A_hat(4,1) 
   shandle%A_hat(4,2) = 0.5529291479
   shandle%A_hat(4,2) = a42hat        
   print *, 'a42 :', shandle%A_hat(4,2) 
   shandle%A_hat(4,3) = 0.5529291479
   shandle%A_hat(4,3) = a42hat        
   print *, 'a43 :', shandle%A_hat(4,3) 
   shandle%B_hat(1)   = 0.0d0 
   shandle%B_hat(2)   = b1
   shandle%B_hat(3)   = b2
   shandle%B_hat(4)   = gam
   shandle%C_hat(1)   = 0.d0
   shandle%C_hat(2)   = gam
   shandle%C_hat(3)   = upg  
   shandle%C_hat(4)   = 1.0d0
 End Subroutine Define_Lstable_343


 Subroutine Define_Lstable_443(shandle, verb)
   Logical, intent(in) :: verb
   Type(Scheme_handle), Intent(Out) :: shandle 
   if (verb) then
   write (*,*) '======================================================'
   write (*,*) '>>> IMEX SCHEME:'
   write (*,*) '... Ascher-Ruuth-Spiteri, Runge-kutta 3, L-Stable'
   Write (*,*) '... Stages: 4 implicits, 4 explicits.'
   write (*,*) '======================================================'
   write (*,*) New_line('a')
   end if
   shandle%imp_stages = 4
   shandle%exp_stages = 4
   shandle%accuracy   = 3
   Allocate( shandle%A_arr(4,4) )
   Allocate( shandle%B_arr(4)   )
   Allocate( shandle%C_arr(4)   )
   shandle%A_arr = 0._dp
   shandle%B_arr = 0._dp
   shandle%C_arr = 0._dp
   shandle%A_arr(1,1) = 0.5_dp
   shandle%A_arr(2,2) = 0.5_dp
   shandle%A_arr(3,3) = 0.5_dp
   shandle%A_arr(4,4) = 0.5_dp
   shandle%A_arr(2,1) = 1._dp/6._dp
   shandle%A_arr(3,1) =-0.5_dp
   shandle%A_arr(3,2) = 0.5_dp
   shandle%A_arr(4,1) = 1.5_dp
   shandle%A_arr(4,2) =-1.5_dp
   shandle%A_arr(4,3) = 0.5_dp
   shandle%B_arr(1)   = 1.5_dp
   shandle%B_arr(2)   =-1.5_dp
   shandle%B_arr(3)   = 0.5_dp
   shandle%B_arr(4)   = 0.5_dp
   shandle%C_arr(1)   = 0.5_dp
   shandle%C_arr(2)   = 2.0_dp/3._dp
   shandle%C_arr(3)   = 0.5_dp
   shandle%C_arr(4)   = 1.0_dp
   Allocate( shandle%A_hat(5,5) )
   Allocate( shandle%B_hat(5)   )
   Allocate( shandle%C_hat(5)   )
   shandle%A_hat = 0._dp
   shandle%B_hat = 0._dp
   shandle%C_hat = 0._dp
   shandle%A_hat(2,1) = 0.5_dp
   shandle%A_hat(3,1) = 11._dp/18._dp
   shandle%A_hat(3,2) = 1._dp/18._dp
   shandle%A_hat(4,1) = 5._dp/ 6._dp
   shandle%A_hat(4,2) =-5._dp/ 6._dp
   shandle%A_hat(4,3) = 0.5_dp
   shandle%A_hat(5,1) = 0.25_dp
   shandle%A_hat(5,2) = 7._dp/4._dp
   shandle%A_hat(5,3) = 0.75_dp    
   shandle%A_hat(5,4) =-7._dp/4._dp 
   shandle%B_hat(1)   = 0.25_dp
   shandle%B_hat(2)   = 7._dp/4._dp
   shandle%B_hat(3)   = 0.75_dp
   shandle%B_hat(4)   =-7._dp/4._dp
   shandle%C_hat(1)   = 0._dp
   shandle%C_hat(2)   = 0.5_dp
   shandle%C_hat(3)   = 2.0_dp/3.0_dp
   shandle%C_hat(4)   = 0.5_dp
   shandle%C_hat(5)   = 1.0_dp
 End Subroutine Define_Lstable_443


   

 !Subroutine scalar_eq_onetimestep( L, u0, shandle, dt)
   !Real(dp), Intent(In) :: L
   !Real(dp), Intent(InOut) :: u0
   !Type(scheme_handle), Intent(In) :: shandle
   !Real(dp), Intent(In) :: dt
   !! INTERNAL VARIABLES:
   !Integer :: i, j
   !Real(dp), Dimension(:), Allocatable :: Kvar
   !Real(dp), Dimension(:), Allocatable :: Khat
   !Real(dp) :: ui
   !Real(dp) :: aux
!
   !Allocate( Kvar(shandle%imp_stages) )
   !Allocate( Khat(shandle%exp_stages) )
!
   !Do i = 1, shandle%imp_stages
      !If (i.eq.1) Then
        !Call scalar_nonlinearity(u0,Khat(i))
      !Else 
        !ui = aux + dt * shandle%A_arr((i-1),(i-1)) * Kvar(i-1)
        !Call scalar_nonlinearity(ui,Khat(i))
      !End If
      !aux = u0
      !Do j=1,i-1
        !aux = aux + dt * shandle%a_arr(i,  j) * Kvar(j) 
        !aux = aux + dt * shandle%a_hat(i+1,j) * Khat(j)
      !End Do
      !aux = aux + dt * shandle%a_hat(i+1,i) * Khat(i)
      !!solve for kvar(i):
      !kvar(i) = l*aux / (1.0d0- l*dt * shandle%a_arr(i,i) )
   !End Do
   !If (shandle%Exp_stages.eq.(shandle%imp_stages+1)) then
        !ui = aux + dt * shandle%A_arr(shandle%imp_stages,shandle%imp_stages) &
                                               !* Kvar(shandle%imp_stages)
        !Call scalar_nonlinearity(ui,Khat(shandle%exp_stages))
   !End If
!
   !Do j=1,shandle%imp_stages
      !u0 = u0 + dt*shandle%b_arr(j)*Kvar(j)
   !End Do
   !Do j=1,shandle%exp_stages
      !u0 = u0 + dt*shandle%b_hat(j)*Khat(j)
   !End Do
 !End Subroutine scalar_eq_onetimestep
!
 !Subroutine scalar_nonlinearity(u,q)
   !Real(dp), Intent(In)  :: u
   !Real(dp), Intent(Out) :: q
   !q =-u**3
 !End Subroutine scalar_nonlinearity


 End Module IMEX_schemes
