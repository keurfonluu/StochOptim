program test

  use stochoptim, only: Evolutionary

  implicit none

  real(kind = 8), parameter :: pi = 3.141592653589793238460d0
  integer(kind = 4), parameter :: n_dim = 5, popsize = 10, max_iter = 500
  real(kind = 8), dimension(:), allocatable :: lower, upper
  type(Evolutionary) :: ea

  ! Define search boundaries
  allocate(lower(n_dim), upper(n_dim))
  lower = -5.12
  upper = 5.12

  ! Initialize evolutionary optimizer
  ea = Evolutionary(rastrigin, lower, upper, &
                    popsize = popsize, max_iter = max_iter)

  ! Optimize
  call ea % optimize(solver = "cpso")

  ! Display results
  call ea % print()

contains

  function rastrigin(x) result(f)
    real(kind = 8), dimension(:), intent(in) :: x
    integer(kind = 4) :: n_dim
    real(kind = 8) :: f, sum1

    n_dim = size(x)
    sum1 = sum( x**2 - 10.0d0 * cos( 2.0d0 * pi * x ) )
    f = 10.0d0 * n_dim + sum1
    return
  end function rastrigin

end program test
