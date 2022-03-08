program driver_algebra
 use LP_algebra
 implicit none

 type(dcoo_matrix) :: coo
 type(dOperator_1d_1coupled_T) :: Op

 integer :: N, iterm

 real(dp), allocatable :: vec(:)
 real(dp), allocatable :: rhs(:)

 N = 5
 ! COO matrix definition
 coo%nelems = 15
 allocate(coo%row(coo%nelems)) 
 allocate(coo%col(coo%nelems)) 
 allocate(coo%dat(coo%nelems)) 
 coo%nrow = N
 coo%ncol = N
 
 iterm = 1

 call add_term(1, 1, 0._dp, coo, iterm)    
 call add_term(1, 2, 1._dp, coo, iterm)    
 call add_term(1, 3, 2._dp, coo, iterm)    
 call add_term(2, 1, 9._dp, coo, iterm)    
 call add_term(2, 2, 4._dp, coo, iterm)    
 call add_term(2, 3, 6._dp, coo, iterm)    
 call add_term(2, 4, 8._dp, coo, iterm)    
 call add_term(3, 2, 1._dp, coo, iterm)    
 call add_term(3, 3,-3._dp, coo, iterm)    
 call add_term(3, 4, 5._dp, coo, iterm)    
 call add_term(3, 5, 7._dp, coo, iterm)    
 call add_term(4, 3,-2._dp, coo, iterm)    
 call add_term(4, 5, 4._dp, coo, iterm)    
 call add_term(5, 4, 1._dp, coo, iterm)    
 call add_term(5, 5, 3._dp, coo, iterm)    

 call Op%constructor(coo)

 allocate(vec(5))
 allocate(rhs(5))
 rhs = [-4._dp,13._dp, -4._dp, -14._dp, -11._dp]

 call Op%factorize()
 call Op%backsolve(rhs, 1, N)
 print *, rhs
 call Op%dot('N',rhs, vec)
 print *, (vec(1) - (-4._dp))**2 + &
          (vec(2) - (13._dp))**2 + &
          (vec(3) - (-4._dp))**2 + &
          (vec(4) -(-14._dp))**2 + &
          (vec(5) -(-11._dp))**2 
 call Op%dot('N',rhs, vec)
 print *, (vec(1) - (-4._dp))**2 + &
          (vec(2) - (13._dp))**2 + &
          (vec(3) - (-4._dp))**2 + &
          (vec(4) -(-14._dp))**2 + &
          (vec(5) -(-11._dp))**2 

 contains

 subroutine add_term(irow, icol, ddat, mycoo, myterm)
 integer, intent(in) :: irow, icol
 integer, intent(inOut) :: myTerm
 type(dcoo_matrix) :: myCoo
 real(dp), intent(in) :: ddat
 myCoo%dat(myterm) = ddat
 myCoo%row(myterm) = irow
 myCoo%col(myterm) = icol
 myterm = myterm + 1
 end subroutine 

end program driver_algebra

