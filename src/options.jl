# Options struct for MADS
# Strictly non-algorithmic: stopping criteria and output only

"""
    Options(; kwargs...)

Non-algorithmic options for MADS optimization.
Controls stopping criteria, logging, and trace storage.

Algorithm behavior (like poll strategy) is controlled by the MADS method constructor,
not by Options.

# Stopping Criteria
- `max_iterations::Int = 10_000`: Maximum number of main loop iterations
- `max_evaluations::Int = typemax(Int)`: Maximum function evaluations
- `max_time::Float64 = Inf`: Maximum wall-clock time in seconds
- `min_mesh_size::Float64 = eps(Float64)/2`: Stop when mesh size falls below this
- `f_target::Float64 = -Inf`: Stop when f(x) < f_target (target value reached)
- `f_tol::Float64 = 0.0`: Stop if improvement in f(x) is below this threshold (after 5 consecutive stagnating iters)
- `x_tol::Float64 = 0.0`: Stop if step size ||x_new - x_old|| is below this threshold (after 5 consecutive stagnating iters)

# Output & Logging
- `verbosity::Verbosity = Silent`: Logging verbosity level
- `log_interval::Int = 1`: Log every N iterations (for Iter level)
- `store_trace::Bool = false`: Store history of all visited points

# Reproducibility
- `seed::Union{Nothing,Int} = nothing`: Random seed for direction generation
"""
Base.@kwdef struct Options
    # Stopping Criteria
    max_iterations::Int = 10_000
    max_evaluations::Int = typemax(Int)
    max_time::Float64 = Inf
    min_mesh_size::Float64 = eps(Float64) / 2
    f_target::Float64 = -Inf
    f_tol::Float64 = 0.0
    x_tol::Float64 = 0.0

    # Output & Logging
    verbosity::Verbosity = Silent
    log_interval::Int = 1
    store_trace::Bool = false

    # Reproducibility
    seed::Union{Nothing,Int} = nothing
end

"""
    OptimizationTrace

Stores the history of optimization for later analysis.
"""
mutable struct OptimizationTrace
    x::Vector{Vector{Float64}}
    f::Vector{Float64}
    iterations::Vector{Int}
end

OptimizationTrace() = OptimizationTrace(Vector{Float64}[], Float64[], Int[])

function push_trace!(trace::OptimizationTrace, x::Vector{Float64}, fx::Float64, k::Int)
    push!(trace.x, copy(x))
    push!(trace.f, fx)
    push!(trace.iterations, k)
end
