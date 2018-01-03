!=======================================================================
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!=======================================================================

module benchmark_functions

  use forlab, only : IPRE, RPRE, pi, ones, linspace, randu

  implicit none

  private
  public :: BenchmarkFunction

  abstract interface
    real(kind = RPRE) function optifunc(x)
      use forlab, only : RPRE
      import
      real(kind = RPRE), dimension(:), intent(in) :: x
    end function optifunc
  end interface

  type BenchmarkFunction
    procedure(optifunc), pointer, nopass :: func
    real(kind = RPRE), dimension(:), allocatable :: lower, upper
    integer(kind = IPRE) :: n_dim = 2
    real(kind = RPRE) :: min = 0.
  end type BenchmarkFunction

  interface BenchmarkFunction
    module procedure init_BenchmarkFunction
  end interface BenchmarkFunction

contains

  function init_BenchmarkFunction(func, n_dim) result(bf)
    type(BenchmarkFunction) :: bf
    character(len = *), intent(in) :: func
    integer(kind = IPRE), intent(in), optional :: n_dim

    if ( present(n_dim) ) then
      if ( n_dim .lt. 1 ) then
        print *, "Error: n_dim must be an integer > 0"
        stop
      else
        bf % n_dim = n_dim
      end if
    end if

    allocate(bf % lower(bf % n_dim), bf % upper(bf % n_dim))
    select case(func)
    case("ackley")
      bf % func => ackley
      bf % lower = -32.768
      bf % upper = 32.768
    case("griewank")
      bf % func => griewank
      bf % lower = -600.
      bf % upper = 600.
    case("quartic")
      bf % func => quartic
      bf % lower = -1.28
      bf % upper = 1.28
    case("quartic_noise")
      bf % func => quartic_noise
      bf % lower = -1.28
      bf % upper = 1.28
    case("rastrigin")
      bf % func => rastrigin
      bf % lower = -5.12
      bf % upper = 5.12
    case("rosenbrock")
      bf % func => rosenbrock
      bf % lower = -5.12
      bf % upper = 5.12
    case("sphere")
      bf % func => sphere
      bf % lower = -5.12
      bf % upper = 5.12
    case("styblinski-tang")
      bf % func => styblinski_tang
      bf % lower = -5.12
      bf % upper = 5.12
    case default
      print *, "Error: unknown benchmark function '" // trim(func) // "'"
      stop
    end select
    return
  end function init_BenchmarkFunction

  real(kind = RPRE) function ackley(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1, sum2, e = 2.7182818284590451d0

    nd = size(x)
    sum1 = sqrt( 1.0d0 / nd * sum( x**2 ) )
    sum2 = 1.0d0 / nd * sum( cos( 2.0d0 * pi * x ) )
    ackley = 20.0d0 + e - 20.0d0 * exp( -0.2d0 * sum1 ) - exp(sum2)
    return
  end function ackley

  real(kind = RPRE) function griewank(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1, prod1

    nd = size(x)
    sum1 = sum( x**2 ) / 4000.0d0
    prod1 = product( cos( x / sqrt( linspace(1, nd, nd) ) ) )
    griewank = 1.0d0 + sum1 - prod1
    return
  end function griewank

  real(kind = RPRE) function quartic(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd

    nd = size(x)
    quartic = sum( linspace(1, nd, nd) * x**4 )
    return
  end function quartic

  real(kind = RPRE) function quartic_noise(x)
    real(kind = RPRE), dimension(:), intent(in) :: x

    quartic_noise = quartic(x) + randu()
    return
  end function quartic_noise

  real(kind = RPRE) function rastrigin(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1

    nd = size(x)
    sum1 = sum( x**2 - 10.0d0 * cos( 2.0d0 * pi * x ) )
    rastrigin = 10.0d0 * nd + sum1
    return
  end function rastrigin

  real(kind = RPRE) function rosenbrock(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1, sum2

    nd = size(x)
    sum1 = sum( ( x(2:) - x(:nd-1)**2 )**2 )
    sum2 = sum( ( 1.0d0 - x(:nd-1) )**2 )
    rosenbrock = 100.0d0 * sum1 + sum2
    return
  end function rosenbrock

  real(kind = RPRE) function schwefel(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1

    nd = size(x)
    sum1 = sum( x * sin( sqrt( abs(x) ) ) )
    schwefel = 418.9829d0 * nd - sum1
    return
  end function schwefel

  real(kind = RPRE) function sphere(x)
    real(kind = RPRE), dimension(:), intent(in) :: x

    sphere = sum( x**2 )
    return
  end function sphere

  real(kind = RPRE) function styblinski_tang(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    real(kind = RPRE) :: sum1

    sum1 = sum( x**4 - 16.0d0 * x**2 + 5.0d0 * x )
    styblinski_tang = sum1 / 2.0d0 + 39.16599 * size(x)
    return
  end function styblinski_tang

end module benchmark_functions
