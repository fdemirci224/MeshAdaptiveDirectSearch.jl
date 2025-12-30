# High-level solve API for MADS
# Simple functional interface with centralized method dispatch

"""
    METHOD_CONSTRUCTORS

Centralized dispatch table mapping method symbols to constructors.
This is the single place where symbol → constructor mapping is defined.
"""
const METHOD_CONSTRUCTORS = Dict{Symbol,Any}(
    :LtMADS => LtMADS,
    :OrthoMADS => OrthoMADS,
    :RobustLtMADS => RobustLtMADS,
    :RobustOrthoMADS => RobustOrthoMADS,
    :MADS => MADS,
    :RobustMADS => RobustMADS,
)

"""
    solve(f, x0; kwargs...)

High-level interface to solve an optimization problem.

# Arguments
- `f`: Objective function to minimize
- `x0`: Initial point (vector)

# Keyword Arguments
- `lb = -Inf`: Lower bound (scalar or vector)
- `ub = Inf`: Upper bound (scalar or vector)
- `max_evals = 10_000`: Maximum function evaluations
- `max_iters = 10_000`: Maximum iterations
- `method = :OrthoMADS`: Algorithm (:LtMADS, :OrthoMADS, :RobustLtMADS, :RobustOrthoMADS)
- `verbosity = Silent`: Logging verbosity
- `constraints = []`: Constraint functions
- `seed = nothing`: Random seed for reproducibility
- Additional options can be passed as keyword arguments

# Example
```julia
f(x) = sum(abs2, x)
result = solve(f, [1.0, 2.0]; lb=-10.0, ub=10.0, max_evals=1000)
```
"""
function solve(f::Function, x0::AbstractVector{<:Real};
    lb::Union{Real,AbstractVector{<:Real},Nothing}=nothing,
    ub::Union{Real,AbstractVector{<:Real},Nothing}=nothing,
    max_evals::Int=10_000,
    max_iters::Int=10_000,
    method::Symbol=:OrthoMADS,
    verbosity::Verbosity=Silent,
    constraints::Vector=[],
    seed::Union{Nothing,Int}=nothing,
    min_mesh_size::Float64=eps(Float64) / 2,
    f_target::Float64=-Inf,
    f_tol::Float64=0.0,
    x_tol::Float64=0.0,
    max_time::Float64=Inf,
    log_interval::Int=1,
    store_trace::Bool=false)

    # Validate bounds are provided
    if lb === nothing || ub === nothing
        error("Bounds are required. Please provide both `lb` and `ub`.")
    end

    n = length(x0)
    x0_vec = collect(Float64, x0)

    # Process bounds
    lb_vec = lb isa Real ? fill(Float64(lb), n) : collect(Float64, lb)
    ub_vec = ub isa Real ? fill(Float64(ub), n) : collect(Float64, ub)

    # Validate bounds
    length(lb_vec) == n || throw(DimensionMismatch("lb must have length $n"))
    length(ub_vec) == n || throw(DimensionMismatch("ub must have length $n"))

    # Clamp initial point to bounds
    x0_vec = clamp.(x0_vec, lb_vec, ub_vec)

    # Get method constructor from centralized dispatch table
    haskey(METHOD_CONSTRUCTORS, method) ||
        throw(ArgumentError("Unknown method :$method. Valid methods: $(keys(METHOD_CONSTRUCTORS))"))
    constructor = METHOD_CONSTRUCTORS[method]

    # Create MADS instance
    mads = constructor(n)

    # Create options
    opts = Options(
        max_iterations=max_iters,
        max_evaluations=max_evals,
        max_time=max_time,
        min_mesh_size=min_mesh_size,
        f_target=f_target,
        f_tol=f_tol,
        x_tol=x_tol,
        verbosity=verbosity,
        log_interval=log_interval,
        store_trace=store_trace,
        seed=seed
    )

    # Run optimization
    return minimize(mads, f, x0_vec;
        lowerbound=lb_vec,
        upperbound=ub_vec,
        constraints=constraints,
        options=opts)
end

"""
    solve(p::Problem; kwargs...)

Solve an optimization problem defined by a Problem object.

# Example
```julia
p = Problem(rosenbrock, 2)
set_initial!(p, [0.0, 0.0])
set_bounds!(p, -5.0, 5.0)
result = solve(p; max_evals=3000)
```
"""
function solve(p::Problem;
    max_evals::Int=10_000,
    max_iters::Int=10_000,
    method::Symbol=:OrthoMADS,
    verbosity::Verbosity=Silent,
    seed::Union{Nothing,Int}=nothing,
    min_mesh_size::Float64=eps(Float64) / 2,
    f_target::Float64=-Inf,
    f_tol::Float64=0.0,
    x_tol::Float64=0.0,
    max_time::Float64=Inf,
    log_interval::Int=1,
    store_trace::Bool=false)

    return solve(p.f, p.x0;
        lb=p.lb,
        ub=p.ub,
        max_evals=max_evals,
        max_iters=max_iters,
        method=method,
        verbosity=verbosity,
        constraints=p.constraints,
        seed=seed,
        min_mesh_size=min_mesh_size,
        f_target=f_target,
        f_tol=f_tol,
        x_tol=x_tol,
        max_time=max_time,
        log_interval=log_interval,
        store_trace=store_trace)
end
