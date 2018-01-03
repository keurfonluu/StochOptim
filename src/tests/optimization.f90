!=======================================================================
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!=======================================================================

program optimization

  use forlab, only: IPRE, RPRE, CLEN, split_argument, num2str, rng
  use benchmark_functions, only: BenchmarkFunction
  use stochoptim, only: Evolutionary
#ifdef do_mpi
  use mpi
#endif

  implicit none

  ! Inputs
  integer(kind = IPRE) :: n_dim, popsize, max_iter, seed, verbose
  real(kind = RPRE) :: w, c1, c2, gamma, F, CR, sigma, mu_perc
  character(len = :), allocatable :: solver, strategy, func

  ! Local variables
  character(len = 18) :: attributes
  type(Evolutionary) :: ea
  type(BenchmarkFunction) :: bf

  ! MPI variables
  integer(kind = IPRE) :: master = 0, rank, ierr, iprv

  ! Initialize MPI
#ifdef do_mpi
  call mpi_init_thread(mpi_thread_funneled, iprv, ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
#else
  rank = 0
#endif

  ! Read command line arguments
  call command_arguments()

  ! Initialize random number generator
  if ( seed .ge. 0 ) then
    call rng(seed)
  else
    call rng()
  end if

  ! Optimize
  bf = BenchmarkFunction(func, n_dim)
  ea = Evolutionary(bf % func, bf % lower, bf % upper, &
                    popsize = popsize, max_iter = max_iter, &
                    constrain = .true., snap = .false., verbose = verbose)

  call ea % optimize(&
    solver = solver, &
    w = w, &
    c1 = c1, &
    c2 = c2, &
    gamma = gamma, &
    F = F, &
    CR = CR, &
    strategy = strategy, &
    sigma = sigma, &
    mu_perc = mu_perc)

  ! Display results
  if ( rank .eq. master ) then
    print *
    attributes = "function:"
    print *, adjustr(attributes) // " '" // trim(func) // "'"
    call ea % print()
    print *
  end if

#ifdef do_mpi
  call mpi_finalize(ierr)
#endif

contains

  subroutine command_arguments()
    integer(kind = IPRE) :: i, n_argin
    character(len = CLEN) :: argin
    character(len = :), allocatable :: argname, argval

    solver = "cpso"
    func = "sphere"
    n_dim = 2
    popsize = 16
    w = 0.7298
    c1 = 1.49618
    c2 = 1.49618
    gamma = 1.
    F = 0.5
    CR = 0.1
    strategy = "rand1"
    sigma = 1.
    mu_perc = 0.5
    max_iter = 200
    seed = -1
    verbose = 0

    n_argin = command_argument_count()
    do i = 1, n_argin
      call get_command_argument(i, argin)
      call split_argument(argin, argname, argval)

      select case(argname)
      ! Number of dimensions
      case("n_dim")
        read(argval, *) n_dim

      ! Solver
      case("solver")
        solver = argval

      ! Function
      case("func")
        func = argval

      ! Population size
      case("popsize")
        read(argval, *) popsize

      ! Inertia weight
      case("w")
        read(argval, *) w

      ! Cognition parameter
      case("c1")
        read(argval, *) c1

      ! Sociability parameter
      case("c2")
        read(argval, *) c2

      ! Competitivity parameter
      case("gamma")
        read(argval, *) gamma

      ! Differential weight
      case("F")
        read(argval, *) F

      ! Crossover probability
      case("CR")
        read(argval, *) CR

      ! Strategy
      case("strategy")
        strategy = argval

      ! Step size
      case("sigma")
        read(argval, *) sigma

      ! Percentage of offsprings
      case("mu_perc")
        read(argval, *) mu_perc

      ! Maximum number of iterations
      case("max_iter")
        read(argval, *) max_iter

      ! Random number seed
      case("seed")
        read(argval, *) seed

      ! Verbosity
      case("verbose")
        read(argval, *) verbose

      case default
        print *, "Error: unknown argument '" // trim(argname) // "'"
        stop
      end select
    end do
    return
  end subroutine command_arguments

end program optimization
