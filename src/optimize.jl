# Optimization loop for MADS
# Refactored minimize with Options integration and professional logging

"""
    EvalCounter

Wrapper to count function evaluations.
"""
mutable struct EvalCounter{F}
    f::F
    count::Int
end

EvalCounter(f) = EvalCounter(f, 0)

function (ec::EvalCounter)(x)
    ec.count += 1
    return ec.f(x)
end

"""
    minimize(m::AbstractMADS, f, x0 = zeros(length(x0));
             lowerbound = nothing,
             upperbound = nothing,
             max_iterations = 10^4,
             min_mesh_size = eps(Float64)/2,
             verbosity = Silent,
             constraints = [],
             options = nothing)

Minimize function `f` with method `m`.

For possible methods `m` see [`LtMADS`](@ref), [`OrthoMADS`](@ref), [`MADS`](@ref),
[`RobustLtMADS`](@ref), [`RobustOrthoMADS`](@ref), [`RobustMADS`](@ref).

Constraints can be defined by boolean functions, e.g.
`constraints = [x -> sum(x) > 1, x -> x[1]^2 < 3]`.

# Arguments
- `m`: MADS method (LtMADS, OrthoMADS, etc.)
- `f`: Objective function to minimize
- `x0`: Initial point

# Keyword Arguments
- `lowerbound`: Lower bounds (required - scalar or vector)
- `upperbound`: Upper bounds (required - scalar or vector)
- `max_iterations`: Maximum iterations (default: 10000)
- `min_mesh_size`: Minimum mesh size for convergence (default: eps()/2)
- `verbosity`: Logging level (Silent, Final, Iter, Step, Debug)
- `constraints`: Vector of constraint functions
- `options`: Options struct (overrides individual keyword arguments)
"""
function minimize(m::AbstractMADS, f, x0=zeros(length(x0));
    lowerbound=nothing,
    upperbound=nothing,
    max_iterations=10^4,
    min_mesh_size=eps(Float64) / 2,
    verbosity=Silent,
    constraints=[],
    options::Union{Nothing,Options}=nothing)

    # Validate bounds are provided
    if lowerbound === nothing || upperbound === nothing
        error("Bounds are required. Please provide both `lowerbound` and `upperbound`.")
    end

    # If options provided, use those values
    if options !== nothing
        max_iterations = options.max_iterations
        min_mesh_size = options.min_mesh_size
        verbosity = options.verbosity

        # Set random seed if specified
        if options.seed !== nothing
            Random.seed!(options.seed)
        end
    end

    # Create logger
    log_interval = options !== nothing ? options.log_interval : 1
    logger = MADSLogger(verbosity=verbosity, log_interval=log_interval)

    # Validate initial point
    isvalid(constraints, x0) || error("x0 = $x0 doesn't satisfy all constraints.")

    # Setup transformations
    to, from = standardtransformation(lowerbound, upperbound)

    # Wrap function with evaluation counter
    eval_counter = EvalCounter(x -> f(to(x)))
    cinternal = [x -> c(to(x)) for c in constraints]

    # Initialize
    incumbent = from(x0)
    init!(m.poll, incumbent)
    fincumbent = eval_counter(incumbent)

    # Get method name for logging
    method_name = string(typeof(m).name.wrapper)
    log_start!(logger, length(x0), method_name)

    # Trace storage
    trace = options !== nothing && options.store_trace ? OptimizationTrace() : nothing

    # Additional stopping criteria from options
    max_evals = options !== nothing ? options.max_evaluations : typemax(Int)
    max_time = options !== nothing ? options.max_time : Inf
    f_target = options !== nothing ? options.f_target : -Inf
    f_tol = options !== nothing ? options.f_tol : 0.0
    x_tol = options !== nothing ? options.x_tol : 0.0

    # Tracking
    start_time = time()
    prev_incumbent = copy(incumbent)
    prev_fincumbent = fincumbent

    # Stagnation tracking - require consecutive stagnation before stopping
    stagnation_count = 0
    min_stagnation_iters = 5  # Require 5 consecutive stagnating iterations

    # Main optimization loop
    for k in 1:max_iterations
        # Check max_evals before doing more work
        if eval_counter.count >= max_evals
            elapsed_time = time() - start_time
            result = (f=fincumbent, x=to(incumbent), stopping_reason=MaxEvaluations,
                iterations=k - 1, evaluations=eval_counter.count, time=elapsed_time, trace=trace)
            log_final!(logger, result)
            return result
        end

        # Search stage
        search_result = search(m, eval_counter, cinternal, incumbent, fincumbent)
        search_incumbent, search_fincumbent, search_improved = search_result.incumbent, search_result.fx, search_result.hasimproved

        # Poll stage (if search didn't find full success)
        if search_improved != 1
            poll_result = poll(m, eval_counter, cinternal, incumbent, fincumbent)
            incumbent, fincumbent, i = poll_result.incumbent, poll_result.fx, poll_result.hasimproved
        else
            incumbent, fincumbent, i = search_incumbent, search_fincumbent, search_improved
        end

        # Update mesh
        update!(m.mesh, i)

        # Store trace
        if trace !== nothing
            push_trace!(trace, to(incumbent), fincumbent, k)
        end

        # Logging
        current_Δ = Δ(m.mesh)
        log_iteration!(logger, k, fincumbent, current_Δ, i, eval_counter.count)
        log_step!(logger, k, i == search_improved ? :search : :poll, i,
            search_fincumbent, fincumbent, current_Δ)

        # Debug logging (gated, per-iteration summary only)
        if should_log(logger, Debug)
            cache_size = m isa RobustMADS ? length(m.cache.y) : 0
            log_debug!(logger, k,
                directions_count=length(m.poll) + 1,
                cache_size=cache_size,
                mesh_level=ℓ(m.mesh))
        end

        # Check stopping criteria

        # Min mesh size
        if current_Δ < min_mesh_size
            elapsed_time = time() - start_time
            result = (f=fincumbent, x=to(incumbent), stopping_reason=MinMeshSize,
                iterations=k, evaluations=eval_counter.count, time=elapsed_time, trace=trace)
            log_final!(logger, result)
            return result
        end

        # F target reached
        if fincumbent < f_target
            elapsed_time = time() - start_time
            result = (f=fincumbent, x=to(incumbent), stopping_reason=FTargetReached,
                iterations=k, evaluations=eval_counter.count, time=elapsed_time, trace=trace)
            log_final!(logger, result)
            return result
        end

        # Max evaluations
        if eval_counter.count >= max_evals
            elapsed_time = time() - start_time
            result = (f=fincumbent, x=to(incumbent), stopping_reason=MaxEvaluations,
                iterations=k, evaluations=eval_counter.count, time=elapsed_time, trace=trace)
            log_final!(logger, result)
            return result
        end

        # Max time
        if time() - start_time > max_time
            elapsed_time = time() - start_time
            result = (f=fincumbent, x=to(incumbent), stopping_reason=MaxTime,
                iterations=k, evaluations=eval_counter.count, time=elapsed_time, trace=trace)
            log_final!(logger, result)
            return result
        end

        # F tolerance (stagnation) - require consecutive stagnation
        f_stagnating = f_tol > 0 && abs(fincumbent - prev_fincumbent) < f_tol && i <= 0
        x_stagnating = x_tol > 0 && norm(incumbent - prev_incumbent) < x_tol && i <= 0

        if f_stagnating || x_stagnating
            stagnation_count += 1
        else
            stagnation_count = 0
        end

        if stagnation_count >= min_stagnation_iters
            elapsed_time = time() - start_time
            stopping_reason = f_stagnating ? FTolReached : XTolReached
            result = (f=fincumbent, x=to(incumbent), stopping_reason=stopping_reason,
                iterations=k, evaluations=eval_counter.count, time=elapsed_time, trace=trace)
            log_final!(logger, result)
            return result
        end

        # Update previous values
        prev_incumbent = copy(incumbent)
        prev_fincumbent = fincumbent
    end

    elapsed_time = time() - start_time
    result = (f=fincumbent, x=to(incumbent), stopping_reason=MaxIterations,
        iterations=max_iterations, evaluations=eval_counter.count, time=elapsed_time, trace=trace)
    log_final!(logger, result)
    return result
end
