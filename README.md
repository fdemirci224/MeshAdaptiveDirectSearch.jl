# Mesh Adaptive Direct Search (MADS)

A pure Julia implementation of (Robust)LtMADS and (Robust)OrthoMADS for derivative-free blackbox optimization.

## Installation

```julia
] add https://github.com/fdemirci224/MeshAdaptiveDirectSearch.jl.git
```

---

## Quick Start

```julia
using MeshAdaptiveDirectSearch

f(x) = (x[1] - 1)^2 + 100*(x[2] - x[1]^2)^2  # Rosenbrock

result = solve(f, [0.0, 0.0]; lb=-5.0, ub=5.0, max_evals=5000, verbosity=Final)
```

---

## Two API Styles

### 1. Functional API (`solve`)

```julia
f(x) = sum(abs2, x)

result = solve(f, [5.0, -3.0]; 
    lb = -10.0, 
    ub = 10.0, 
    max_evals = 3000,
    method = :OrthoMADS,
    verbosity = Iter
)
```

### 2. Problem Builder API

```julia
rosenbrock(x) = (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2

p = Problem(rosenbrock, 2)
set_initial!(p, [-1.0, 1.0])
set_bounds!(p, -5.0, 5.0)
add_constraint!(p, x -> x[1] + x[2] < 1.5)

result = solve(p; max_evals=5000, verbosity=Final)
```

### 3. Using Options Struct

```julia
# Create reusable options
opts = Options(
    max_iterations = 10_000,
    max_evaluations = 5000,
    f_target = 1e-6,        # Stop when f(x) < target
    min_mesh_size = 1e-12,
    verbosity = Final,
    seed = 42               # For reproducibility
)

# Use with minimize()
result = minimize(OrthoMADS(2), f, [0.0, 0.0];
    lowerbound = [-5.0, -5.0],
    upperbound = [5.0, 5.0],
    options = opts
)
```

### 4. Low-Level API (`minimize`)

```julia
# Direct method construction
result = minimize(LtMADS(2), f, [0.0, 0.0];
    lowerbound = [-10.0, -10.0],
    upperbound = [10.0, 10.0],
    max_iterations = 5000,
    verbosity = Silent
)
```

---

## Available Methods

| Method | Symbol | Description |
|--------|--------|-------------|
| `LtMADS(n)` | `:LtMADS` | Random lower triangular directions |
| `OrthoMADS(n)` | `:OrthoMADS` | Deterministic orthogonal Halton directions |
| `RobustLtMADS(n)` | `:RobustLtMADS` | LT-MADS with kernel smoothing for noisy objectives |
| `RobustOrthoMADS(n)` | `:RobustOrthoMADS` | OrthoMADS with kernel smoothing for noisy objectives |

---

## Stopping Criteria

| Parameter | Description |
|-----------|-------------|
| `max_iters` | Maximum iterations |
| `max_evals` | Maximum function evaluations |
| `max_time` | Wall-clock time limit (seconds) |
| `min_mesh_size` | Stop when mesh size < threshold (convergence) |
| `f_target` | Stop when f(x) < target value |
| `f_tol` | Stop after 5 consecutive stagnating iterations (f improvement < threshold) |
| `x_tol` | Stop after 5 consecutive stagnating iterations (step size < threshold) |

---

## Verbosity Levels

| Level | Description |
|-------|-------------|
| `Silent` | No output |
| `Final` | Summary at termination (iterations, evaluations, time, result) |
| `Iter` | One line per iteration with eval count |
| `Step` | Detailed search/poll stage info |
| `Debug` | Diagnostics per iteration |

---

## Result

```julia
result.f                # Best objective value
result.x                # Best point found
result.stopping_reason  # Why optimization stopped
result.iterations       # Total iterations
result.evaluations      # Total function evaluations
result.time             # Total time (seconds)
result.trace            # History (if store_trace=true)
```

**Stopping reasons:** `MinMeshSize`, `MaxIterations`, `MaxEvaluations`, `MaxTime`, `FTargetReached`, `FTolReached`, `XTolReached`

---

## Constraints

```julia
result = solve(f, x0; 
    lb = -10.0, ub = 10.0,
    constraints = [
        x -> x[1] + x[2] < 1.0,   # Linear
        x -> x[1]^2 + x[2]^2 < 4  # Nonlinear
    ]
)
```

---

## References

- Audet & Dennis (2006), "Mesh Adaptive Direct Search Algorithms for Constrained Optimization", [doi](http://dx.doi.org/10.1137/040603371)
- Abramson et al. (2009), "OrthoMADS: A Deterministic MADS Instance with Orthogonal Directions", [doi](http://dx.doi.org/10.1137/080716980)
- Audet et al. (2014), "Reducing the Number of Function Evaluations in MADS", [doi](http://dx.doi.org/10.1137/120895056)
- Audet et al. (2018), "Robust optimization of noisy blackbox problems using MADS", [doi](http://dx.doi.org/10.1007/s11590-017-1226-6)
