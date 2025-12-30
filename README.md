# Mesh Adaptive Direct Search (MADS)

A pure Julia implementation of (Robust)LtMADS and (Robust)OrthoMADS for derivative-free blackbox optimization.

See [NOMAD.jl](https://github.com/ppascal97/NOMAD.jl) for a Julia wrapper of [NOMAD](https://www.gerad.ca/nomad/).

## Installation

Type `]` in the Julia REPL to enter the package REPL, then
```julia
add https://github.com/jbrea/MeshAdaptiveDirectSearch.jl
```
and backspace or ^C to leave it again.

---

## Quick Start

### Simple API

```julia
using MeshAdaptiveDirectSearch

f(x) = sum(abs2, x)  # Sphere function

# Basic usage - finds minimum near [0, 0]
result = solve(f, [5.0, -3.0]; lb=-10.0, ub=10.0, max_evals=3000)

println("Best f(x) = $(result.f)")
println("Best x    = $(result.x)")
```

### Problem Builder API

```julia
using MeshAdaptiveDirectSearch

rosenbrock(x) = (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2

# Create and configure problem
p = Problem(rosenbrock, 2)
set_initial!(p, [-1.0, 1.0])
set_bounds!(p, -5.0, 5.0)
add_constraint!(p, x -> x[1]^2 + x[2]^2 < 4.0)

# Solve
result = solve(p; max_evals=5000, verbosity=Final)
```

---

## Algorithms

| Method | Description |
|--------|-------------|
| `LtMADS(n)` | LT-MADS with random lower triangular directions |
| `OrthoMADS(n)` | Deterministic MADS with orthogonal Halton directions |
| `RobustLtMADS(n)` | LT-MADS with kernel smoothing for noisy objectives |
| `RobustOrthoMADS(n)` | OrthoMADS with kernel smoothing for noisy objectives |

### Using with `solve()`

```julia
# Default: OrthoMADS
result = solve(f, x0; lb=-10.0, ub=10.0)

# Specify method
result = solve(f, x0; method=:LtMADS, lb=-10.0, ub=10.0)

# For noisy functions
result = solve(noisy_f, x0; method=:RobustOrthoMADS, lb=-10.0, ub=10.0)
```

### Using with `minimize()` (Advanced)

```julia
result = minimize(LtMADS(2), f, x0; lowerbound=[-10, -10], upperbound=[10, 10])
result = minimize(OrthoMADS(2), f, x0; lowerbound=[-10, -10], upperbound=[10, 10])
result = minimize(RobustLtMADS(2), noisy_f, x0; lowerbound=[-10, -10], upperbound=[10, 10])
```

---

## Verbosity Levels

Control output with the `verbosity` parameter:

| Level | Description |
|-------|-------------|
| `Silent` | No output (default) |
| `Final` | Print summary at termination |
| `Iter` | One line per iteration |
| `Step` | Detailed search/poll stage info |
| `Debug` | Diagnostics (direction count, cache size) |

```julia
result = solve(f, x0; verbosity=Iter, log_interval=10)  # Log every 10 iterations
```

---

## Options

All stopping criteria and logging options available via keyword arguments:

```julia
result = solve(f, x0;
    # Bounds
    lb = -10.0,                    # Lower bound (scalar or vector)
    ub = 10.0,                     # Upper bound (scalar or vector)
    
    # Stopping Criteria
    max_iters = 10_000,            # Maximum iterations
    max_evals = 10_000,            # Maximum function evaluations
    max_time = 60.0,               # Maximum time (seconds)
    min_mesh_size = eps()/2,       # Minimum mesh size (convergence)
    f_tol = 0.0,                   # Stop if improvement < f_tol
    x_tol = 0.0,                   # Stop if step size < x_tol
    
    # Output
    verbosity = Silent,            # Logging level
    log_interval = 1,              # Log every N iterations
    store_trace = false,           # Store optimization history
    
    # Reproducibility
    seed = 12345,                  # Random seed
)
```

---

## Constraints

Define constraints as boolean functions that return `true` when satisfied:

```julia
# With solve()
result = solve(f, x0; 
    lb=-10.0, ub=10.0,
    constraints=[x -> x[1] + x[2] < 1.0, x -> x[1]^2 < 2.0]
)

# With Problem API
p = Problem(f, 2)
set_bounds!(p, -10.0, 10.0)
add_constraint!(p, x -> sum(x) < 0.5)
add_constraint!(p, x -> x[1] > -3.0)
result = solve(p)

# With minimize()
result = minimize(OrthoMADS(2), f, x0;
    lowerbound=[-10, -10], 
    upperbound=[10, 10],
    constraints=[x -> sum(x) < 0.5]
)
```

---

## Result

The result is a named tuple containing:

```julia
result.f                # Best objective value found
result.x                # Best point found
result.stopping_reason  # Why optimization stopped
result.iterations       # Number of iterations
result.evaluations      # Number of function evaluations
result.trace            # History (if store_trace=true)
```

Stopping reasons: `MinMeshSize`, `MaxIterations`, `MaxEvaluations`, `MaxTime`, `FTolReached`, `XTolReached`

---

## Complete Example

```julia
using MeshAdaptiveDirectSearch, Random

# Rosenbrock function (global minimum at [1, 1])
rosenbrock(x) = (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2

# Noisy version
noisy_rosenbrock(x) = rosenbrock(x) + 0.1*randn()

Random.seed!(42)

# Standard optimization
result = solve(rosenbrock, [-2.0, 2.0];
    lb = -5.0,
    ub = 5.0,
    max_evals = 5000,
    verbosity = Final
)

# For noisy objective
result = solve(noisy_rosenbrock, [-2.0, 2.0];
    method = :RobustOrthoMADS,
    lb = -5.0,
    ub = 5.0,
    max_evals = 5000,
    verbosity = Final
)
```

---

## References

- Audet, Charles and Dennis, J. E., "Mesh Adaptive Direct Search Algorithms for Constrained Optimization", 2006, [doi](http://dx.doi.org/10.1137/040603371)

- Abramson, Mark A. and Audet, Charles and Dennis, J. E. and Le Digabel, Sébastien, "OrthoMADS: A Deterministic MADS Instance with Orthogonal Directions", 2009, [doi](http://dx.doi.org/10.1137/080716980)

- Audet, Charles and Ianni, Andrea and Le Digabel, Sébastien and Tribes, Christophe, "Reducing the Number of Function Evaluations in Mesh Adaptive Direct Search Algorithms", 2014, [doi](http://dx.doi.org/10.1137/120895056)

- Audet, Charles and Ihaddadene, Amina and Le Digabel, Sébastien and Tribes, Christophe, "Robust optimization of noisy blackbox problems using the Mesh Adaptive Direct Search algorithm", 2018, [doi](http://dx.doi.org/10.1007/s11590-017-1226-6)
