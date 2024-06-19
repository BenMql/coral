module chebyshev_tau
  use fortran_kinds
  implicit none

  contains

  subroutine cheby_tau_line_coefficients( N, pos, derivative_order)
     integer, intent(in) :: N
     integer, intent(in) :: pos ! should be -1 (bottom) or +1 (top)
     integer, intent(in) :: derivative_order
  end subroutine

end module chebyshev_tau
