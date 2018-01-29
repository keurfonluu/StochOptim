# StochOptim
StochOptim provides user friendly functions to solve optimization problems using stochastic algorithms in a parallel environment (**MPI**).

| | |  
|:-:|---|
| **Version:** | 1.0.3 |
| **Author:** | Keurfon Luu |
| **Web site:** | https://github.com/keurfonluu/stochoptim |
| **Copyright:** | This document has been placed in the public domain. |
| **License:** | StochOptim is released under the MIT License. |

**NOTE:** StochOptim has been implemented in the frame of my Ph. D. thesis. If you find any error or bug, or if you have any suggestion, please don't hesitate to contact me.


## Features

StochOptim provides routines written in oriented-object Fortran to sample a parameter model space and optimize objective functions:

* Metropolis-Hastings algorithm
* Differential Evolution
* Particle Swarm Optimization
* Competitive Particle Swarm Optimization
* Covariance Matrix Adaptation - Evolution Strategy


## Usage

Place files in *src/lib* in your source directory and import *stochoptim* module:

```fortran
  use stochoptim, only: Evolutionary, MonteCarlo
```

Methods available:

* optimize / sample
* print_parameters
* print_results
* print
* save_models

For evolutionary optimizers, population individuals can be parallelized using **MPI** simply by adding the preprocessor option *-Ddo_mpi*.


### Example

Examples codes can be found in *src/tests*, or you can try the following code:

```fortran
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
```


### Results

```bash
            solver: 'cpso'
        parameters: w = 0.72, c1 = 1.49, c2 = 1.49, gamma = 1.00
           popsize: 10
          max_iter: 500
             n_dim: 5
          solution:
                           2.4281118235742959E-006
                          -1.8335474048010150E-006
                           1.0076754625181591E-006
                          -1.8331826109239015E-009
                           9.6263231204725598E-007
           fitness: 2.2219381889954093E-009
            n_iter: 358
            n_eval: 3580
         n_restart: 4
              flag: fitness is lower than threshold eps2 (9.9999999392252903E-009)
```
