!=======================================================================
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!=======================================================================

program benchmark

  use forlab, only: IPRE, RPRE, CLEN, progress_perc, split_argument, &
    num2str, rng, mean, median, std, tic, toc
  use benchmark_functions, only: BenchmarkFunction
  use stochoptim, only: Evolutionary
#ifdef do_mpi
  use mpi
#endif

  implicit none

  ! Inputs
  integer(kind = IPRE) :: n_dim, popsize, max_iter, max_run, seed
  real(kind = RPRE) :: w, c1, c2, gamma, F, CR, sigma, mu_perc
  character(len = :), allocatable :: solver, func

  ! Local variables
  integer(kind = IPRE) :: n_run
  real(kind = RPRE) :: min_fit, max_fit, mean_fit, median_fit, std_fit
  real(kind = RPRE), dimension(:), allocatable :: fit
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

  ! Optimize max_run times
  bf = BenchmarkFunction(func, n_dim)
  ea = Evolutionary(bf % func, bf % lower, bf % upper, &
                    popsize = popsize, max_iter = max_iter, &
                    eps2 = 1d-30, constrain = .false., snap = .false.)

  allocate(fit(max_run))
  if ( rank .eq. master ) then
    print *
    call tic()
    call progress_perc(0, max_run, " Processing: ")
  end if
  do n_run = 1, max_run
    call ea % optimize(&
      solver = solver, &
      w = w, &
      c1 = c1, &
      c2 = c2, &
      gamma = gamma, &
      F = F, &
      CR = CR, &
      sigma = sigma, &
      mu_perc = mu_perc)
    fit(n_run) = ea % gfit
    if ( rank .eq. master ) call progress_perc(n_run, max_run, " Processing: ")
  end do

  ! Display results
  if ( rank .eq. master ) then
    min_fit = minval(fit)
    max_fit = maxval(fit)
    mean_fit = mean(fit)
    median_fit = median(fit)
    std_fit = std(fit)

    print *
    print *

    attributes = "function:"
    print *, adjustr(attributes) // " '" // trim(func) // "'"
    call ea % print_parameters()
    print *
    attributes = "max_run:"
    print *, adjustr(attributes) // " " // num2str(max_run)
    attributes = "min:"
    print *, adjustr(attributes) // " " // num2str(min_fit)
    attributes = "max:"
    print *, adjustr(attributes) // " " // num2str(max_fit)
    attributes = "mean:"
    print *, adjustr(attributes) // " " // num2str(mean_fit)
    attributes = "median:"
    print *, adjustr(attributes) // " " // num2str(median_fit)
    attributes = "std:"
    print *, adjustr(attributes) // " " // num2str(std_fit)

    print *
    call toc()
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
    sigma = 0.5
    mu_perc = 0.5
    max_iter = 200
    max_run = 100
    seed = -1

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

      ! Step size
      case("sigma")
        read(argval, *) sigma

      ! Percentage of offsprings
      case("mu_perc")
        read(argval, *) mu_perc

      ! Maximum number of iterations
      case("max_iter")
        read(argval, *) max_iter

      ! Maximum number of iterations
      case("max_run")
        read(argval, *) max_run

      ! Random number seed
      case("seed")
        read(argval, *) seed

      case default
        print *, "Error: unknown argument '" // trim(argname) // "'"
        stop
      end select
    end do
    return
  end subroutine command_arguments

end program benchmark
