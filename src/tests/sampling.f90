!=======================================================================
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!=======================================================================

program sampling

  use forlab, only: IPRE, RPRE, CLEN, split_argument, ones, num2str, rng
  use benchmark_functions, only: BenchmarkFunction
  use stochoptim, only: MonteCarlo
#ifdef do_mpi
  use mpi
#endif

  implicit none

  ! Inputs
  integer(kind = IPRE) :: n_dim, max_iter, seed, verbose
  real(kind = RPRE) :: stepsize, perc
  character(len = :), allocatable :: sampler, func

  ! Local variables
  character(len = 18) :: attributes
  type(MonteCarlo) :: mc
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

  bf = BenchmarkFunction(func, n_dim)
  mc = MonteCarlo(bf % func, bf % lower, bf % upper, max_iter = max_iter, &
                  perc = perc, constrain = .true., verbose = verbose)

  call mc % sample(&
    sampler = sampler, &
    stepsize = stepsize * ones(n_dim))

  ! Display results
  if ( rank .eq. master ) then
    print *
    attributes = "function:"
    print *, adjustr(attributes) // " '" // trim(func) // "'"
    call mc % print()
    print *
  end if

  ! Save results
  if ( rank .eq. master ) then
    call system(&
      "rm -rf output;" // &
      "mkdir -p output")
    call mc % save_models("output/")
  end if

#ifdef do_mpi
  call mpi_finalize(ierr)
#endif

contains

  subroutine command_arguments()
    integer(kind = IPRE) :: i, n_argin
    character(len = CLEN) :: argin
    character(len = :), allocatable :: argname, argval

    sampler = "mcmc"
    func = "sphere"
    n_dim = 2
    max_iter = 1000
    stepsize = 0.1
    perc = 1.
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

      ! Sampler
      case("sampler")
        sampler = argval

      ! Function
      case("func")
        func = argval

      ! Maximum number of iterations
      case("max_iter")
        read(argval, *) max_iter

      ! Stepsize
      case("stepsize")
        read(argval, *) stepsize

      ! Number of dimensions to perturbe as a percentage
      case("perc")
        read(argval, *) perc

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

end program sampling
