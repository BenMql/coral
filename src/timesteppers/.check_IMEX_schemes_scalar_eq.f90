 Program test_timestep_scheme_scalar
 Use IMEX_schemes
 Implicit None

 Type(Scheme_handle) :: my_sche
 Integer :: NT
 Real(8) :: dt 
 Real(8) :: cur_val
 Real(8), Dimension(:), Allocatable :: arr_val
 Integer :: i 
 Real(8) :: tmax = 2.75d0
 Character(len=6) :: suffix2sav

 Call chdir('/home/ben/dev/sandbox/')

!########################################################
!########################################################
 suffix2sav = 'Sch222'
 Call Define_Lstable_222(my_sche)
!........................................................
 NT = 100 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 1000 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 100000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000000
 print *, NT
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)












!########################################################
!########################################################
 suffix2sav = 'Sch111'
 Call Define_Forward_Backward_Euler_111(my_sche)
!........................................................
 NT = 100 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 1000 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 100000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000000
 print *, NT
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!########################################################
!########################################################
 suffix2sav = 'Sch122'
 Call Define_Implicit_Explicit_Midpoint_122(my_sche)
!........................................................
 NT = 100 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 1000 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 100000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000000
 print *, NT
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!########################################################
!########################################################
 suffix2sav = 'Sch233'
 Call Define_233(my_sche)
!........................................................
 NT = 100 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 1000 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 100000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000000
 print *, NT
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!########################################################
!########################################################
 suffix2sav = 'Sch232'
 Call Define_Lstable_232(my_sche)
!........................................................
 NT = 100 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 1000 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 100000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000000
 print *, NT
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)

!########################################################
!########################################################
 suffix2sav = 'Sch343'
 Call Define_Lstable_343(my_sche)
!........................................................
 NT = 10 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 25 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 60  
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 100 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 200 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 400 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 1000 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 100000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000000
 print *, NT
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)






!########################################################
!########################################################
 suffix2sav = 'Sch443'
 Call Define_Lstable_443(my_sche)
!........................................................
 NT = 100 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 1000 
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 100000
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)
!........................................................
 NT = 10000000
 print *, NT
 Allocate( arr_val(NT) )
 cur_val = 0.01d0 
 dt = tmax/dble(NT)
 Do i = 1,NT
   Call scalar_eq_onetimestep( +2.d0, cur_val, my_sche, dt)
   arr_val(i) = cur_val
 End Do
 Call save_1d_arr(arr_val, suffix2sav,NT)
 Deallocate (arr_val)













 
 Contains 
 Subroutine save_1d_arr(arr, radical_save, numero)
   Character(len=6), Intent(In)  :: radical_save
   Integer, Intent(in) :: numero
   Real(8), Dimension(:), Allocatable, Intent(In) :: arr
   Character(len=100) :: fichiercom
   Character(len=8)   :: nom
   WRITE(nom, "(i8.8)") numero
   fichiercom = TRIM(radical_save)//"_"//nom//".dat"
   open(unit=8,file=fichiercom,form='unformatted',access='stream', Status = 'replace')
   write(unit=8) arr(:)
   close(unit=8)
 End Subroutine save_1d_arr
 End Program test_timestep_scheme_scalar
