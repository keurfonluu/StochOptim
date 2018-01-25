!=======================================================================
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!=======================================================================

program sensitivity

  use forlab, only: IPRE, RPRE, CLEN, progress_perc, split_argument, &
    num2str, rng, linspace, mean, savebin, tic, toc
  use benchmark_functions, only: BenchmarkFunction
  use stochoptim, only: Evolutionary
#ifdef do_mpi
  use mpi
#endif

  implicit none

  ! Inputs
  integer(kind = IPRE) :: n_dim, popfactor, max_iter, max_run, seed, verbose, nx, ny
  real(kind = RPRE) :: w, c, gamma, F, CR, sigma, mu_perc, xmin, xmax, ymin, ymax
  character(len = :), allocatable :: solver, func, xaxis_name, yaxis_name

  ! Local variables
  integer(kind = IPRE) :: popsize, n_run, i, j, counter
  real(kind = RPRE), dimension(:), allocatable :: xaxis, yaxis
  real(kind = RPRE), dimension(:,:,:), allocatable :: Z
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

  ! Analyze sensitivity to parameters
  allocate(Z(nx,ny,max_run))
  xaxis = linspace(xmin, xmax, nx)
  yaxis = linspace(ymin, ymax, ny)
  counter = 0

  if ( rank .eq. master ) then
    print *
    call tic()
    call progress_perc(0, nx*ny, " Processing: ")
  end if
  select case(solver)
  case("pso", "cpso")
    do i = 1, nx
      do j = 1, ny
        counter = counter + 1
        !$omp parallel default(shared) private(popsize, bf, ea)
        !$omp do schedule(runtime)
        do n_run = 1, max_run
          if ( xaxis_name .eq. "c" .and. yaxis_name .eq. "w" ) then
            popsize = n_dim * popfactor
            bf = BenchmarkFunction(func, n_dim)
            ea = Evolutionary(bf % func, bf % lower, bf % upper, &
                              popsize = popsize, max_iter = max_iter, &
                              constrain = .false., snap = .false., verbose = verbose)
            call ea % optimize(solver = solver, w = yaxis(j), c1 = xaxis(i), &
                               c2 = xaxis(i), gamma = gamma)
          else if ( xaxis_name .eq. "c" .and. yaxis_name .eq. "n_dim" ) then
            popsize = yaxis(j) * popfactor
            bf = BenchmarkFunction(func, int(yaxis(j)))
            ea = Evolutionary(bf % func, bf % lower, bf % upper, &
                              popsize = popsize, max_iter = max_iter, &
                              constrain = .false., snap = .false., verbose = verbose)
            call ea % optimize(solver = solver, w = w, c1 = xaxis(i), &
                               c2 = xaxis(i), gamma = gamma)
          else if ( xaxis_name .eq. "w" .and. yaxis_name .eq. "n_dim" ) then
            popsize = yaxis(j) * popfactor
            bf = BenchmarkFunction(func, int(yaxis(j)))
            ea = Evolutionary(bf % func, bf % lower, bf % upper, &
                              popsize = popsize, max_iter = max_iter, &
                              constrain = .false., snap = .false., verbose = verbose)
            call ea % optimize(solver = solver, w = xaxis(i), c1 = c, &
                               c2 = c, gamma = gamma)
          end if
          Z(i,j,n_run) = ea % gfit
        end do
        !$omp end parallel
        if ( rank .eq. master ) call progress_perc(counter, nx*ny, " Processing: ")
      end do
    end do
  end select

  ! Save results
  if ( rank .eq. master ) then
    call system("mkdir -p output")
    call system("rm -f output/" // solver // "_" // func // "_" // num2str(n_dim) // "_Z.bin")
    call system("rm -f output/" // solver // "_" // func // "_" // num2str(n_dim) // "_xaxis.bin")
    call system("rm -f output/" // solver // "_" // func // "_" // num2str(n_dim) // "_yaxis.bin")
    call savebin("output/" // solver // "_" // func // "_" // num2str(n_dim) // "_Z.bin", Z)
    call savebin("output/" // solver // "_" // func // "_" // num2str(n_dim) // "_xaxis.bin", xaxis)
    call savebin("output/" // solver // "_" // func // "_" // num2str(n_dim) // "_yaxis.bin", yaxis)
  end if

  ! Display results
  if ( rank .eq. master ) then
    print *
    print *
    attributes = "function:"
    print *, adjustr(attributes) // " '" // trim(func) // "'"
    call ea % print()

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
    popfactor = 5
    w = 0.7298
    c = 1.49618
    gamma = 1.
    F = 0.5
    CR = 0.1
    sigma = 0.5
    mu_perc = 0.5
    max_iter = 200
    seed = -1
    xaxis_name = "c"
    xmin = 0.
    xmax = 3.
    nx = 10
    yaxis_name = "w"
    ymin = 0.
    ymax = 1.
    ny = 10
    max_run = 100
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

      ! Population size factor (popsize = popfactor * n_dim)
      case("popfactor")
        read(argval, *) popfactor

      ! Inertia weight
      case("w")
        read(argval, *) w

      ! Acceleration parameters
      case("c")
        read(argval, *) c

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

      ! Random number seed
      case("seed")
        read(argval, *) seed

      ! X axis
      case("xaxis")
        xaxis_name = argval

      ! X min
      case("xmin")
        read(argval, *) xmin

      ! X max
      case("xmax")
        read(argval, *) xmax

      ! NX
      case("nx")
        read(argval, *) nx

      ! Y axis
      case("yaxis")
        yaxis_name = argval

      ! Y min
      case("ymin")
        read(argval, *) ymin

      ! Y max
      case("ymax")
        read(argval, *) ymax

      ! NY
      case("ny")
        read(argval, *) ny

      ! Maximum number of iterations
      case("max_run")
        read(argval, *) max_run

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

  subroutine select_axis(axis, axis_name)
    real(kind = RPRE), dimension(:), intent(inout) :: axis
    character(len = *), intent(in) :: axis_name

    return
  end subroutine select_axis

end program sensitivity
