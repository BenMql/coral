module chebyshev_tau
  use fortran_kinds
  implicit none

  contains

  subroutine cheby_tau_line_coefficients( N, pos, derivative_order, tau_line)
     integer, intent(in) :: N
     integer, intent(in) :: pos ! should be -1 (bottom) or +1 (top)
     integer, intent(in) :: derivative_order
     real(dp), allocatable, intent(inOut) :: tau_line(:)
     integer, allocatable :: pDegree(:)
     integer :: i 

     !check the input
     if (pos**2.ne.1) then
        error stop 'cheby_tau_line_coefficients: pos should be -1 (bottom) or +1 (top)'
     end if

     allocate (pDegree(N))
     do i = 1,N
       pDegree(i) = i-1
     end do

     allocate (tau_line(N))

     select case (derivative_order)
        case (0)
           if (pos.eq.+1) tau_line = 1._dp
           if (pos.eq.-1) tau_line = (-1)**pDegree ! polynomials parity 
        case (1)
           tau_line = pDegree**2
           if (pos.eq.-1) tau_line = tau_line * (-1)**(pDegree+1) ! derivative have flipped parity
        case (2)
           tau_line = pDegree**2 * (pDegree**2 - 1) / 3._dp
           if (pos.eq.-1) tau_line = tau_line * (-1)**(pDegree) ! (d/dz)**2 have same parity as T(z)'s
        case (3:)
           tau_line = 1._dp 
           do i=0,(derivative_order-1)
           tau_line = tau_line * (pDegree**2 - i **2) / (2._dp*i + 1._dp)
           end do
           if (pos.eq.-1) tau_line = tau_line * (-1)**(pDegree+derivative_order) ! 
        case (:-1)
           error stop 'cheby_tau_line_coefficients: derivative_order should be non-negative'
     end select

     deAllocate (pDegree)

  end subroutine

end module chebyshev_tau
