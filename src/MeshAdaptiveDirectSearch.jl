"""
    MeshAdaptiveDirectSearch

Julia implementation of the Mesh Adaptive Direct Search (MADS) algorithm
for derivative-free optimization.

# Algorithms
- `LtMADS`: LT-MADS (Audet & Dennis, 2006)
- `OrthoMADS`: OrthoMADS (Abramson et al., 2009)
- `RobustLtMADS`: Robust LT-MADS for noisy functions (Audet et al., 2018)
- `RobustOrthoMADS`: Robust OrthoMADS for noisy functions

# Basic Usage
```julia
using MeshAdaptiveDirectSearch

# Simple API
f(x) = sum(abs2, x)
result = solve(f, [1.0, 2.0]; lb=-10.0, ub=10.0, max_evals=1000)

# Problem API
p = Problem(rosenbrock, 2)
set_initial!(p, [0.0, 0.0])
set_bounds!(p, -5.0, 5.0)
result = solve(p; max_evals=3000)

# Advanced: direct method construction
result = minimize(OrthoMADS(2), f, [0.0, 0.0]; 
                  lowerbound=[-10.0, -10.0], 
                  upperbound=[10.0, 10.0])
```
"""
module MeshAdaptiveDirectSearch

using StaticArrays, Random, ElasticArrays, LinearAlgebra, Primes
using Printf: @sprintf

# Core algorithm components (in dependency order)
include("meshes.jl")
include("directions_lt.jl")
include("directions_ortho.jl")
include("kernels.jl")
include("cache.jl")
include("constraints.jl")
include("transformations.jl")
include("search.jl")
include("mads_types.jl")
include("poll.jl")

# Logging and options
include("logging.jl")
include("options.jl")

# Optimization loop
include("optimize.jl")

# High-level APIs
include("problem.jl")
include("solve.jl")

# Exports - Algorithms
export MADS, LtMADS, OrthoMADS
export RobustMADS, RobustLtMADS, RobustOrthoMADS

# Exports - Main functions
export minimize, solve

# Exports - Problem API
export Problem, set_initial!, set_bounds!, add_constraint!, clear_constraints!

# Exports - Options and Verbosity
export Options, Verbosity
export Silent, Final, Iter, Step, Debug

# Exports - Stopping reasons
export StoppingReason
export MaxIterations, MinMeshSize, MaxEvaluations, MaxTime, FTolReached, XTolReached, FTargetReached

# Exports - Reduction strategies
export NegReduction, NoReduction

# Exports - Mesh types (for advanced users)
export Mesh, LogMesh

end # module
