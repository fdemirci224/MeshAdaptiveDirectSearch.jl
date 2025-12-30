# Optimization loop for MADS
# Refactored minimize with Options integration and professional logging

"""
    minimize(m::AbstractMADS, f, x0 = zeros(length(x0));
             lowerbound = -ones(length(x0)),
             upperbound = ones(length(x0)),
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
- `lowerbound`: Lower bounds (scalar or vector)
- `upperbound`: Upper bounds (scalar or vector)
- `max_iterations`: Maximum iterations (default: 10000)
- `min_mesh_size`: Minimum mesh size for convergence (default: eps()/2)
- `verbosity`: Logging level (Silent, Final, Iter, Step, Debug)
- `constraints`: Vector of constraint functions
- `options`: Options struct (overrides individual keyword arguments)
"""
function minimize(m::AbstractMADS, f, x0=zeros(length(x0));
    lowerbound=-ones(length(x0)),
    upperbound=ones(length(x0)),
    max_iterations=10^4,
    min_mesh_size=eps(Float64) / 2,
    verbosity=Silent,
    constraints=[],
    options::Union{Nothing,Options}=nothing)

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
    finternal = x -> f(to(x))
    cinternal = [x -> c(to(x)) for c in constraints]

    # Initialize
    incumbent = from(x0)
    init!(m.poll, incumbent)
    fincumbent = finternal(incumbent)

    # Get method name for logging
    method_name = string(typeof(m).name.wrapper)
    log_start!(logger, length(x0), method_name)

    # Trace storage
    trace = options !== nothing && options.store_trace ? OptimizationTrace() : nothing

    # Additional stopping criteria from options
    max_evals = options !== nothing ? options.max_evaluations : typemax(Int)
    max_time = options !== nothing ? options.max_time : Inf
    f_tol = options !== nothing ? options.f_tol : 0.0
    x_tol = options !== nothing ? options.x_tol : 0.0

    # Tracking
    eval_count = 1  # Initial evaluation
    start_time = time()
    prev_incumbent = copy(incumbent)
    prev_fincumbent = fincumbent

    # Main optimization loop
    for k in 1:max_iterations
        # Search stage
        search_result = search(m, finternal, cinternal, incumbent, fincumbent)
        search_incumbent, search_fincumbent, search_improved = search_result.incumbent, search_result.fx, search_result.hasimproved

        # Poll stage (if search didn't find full success)
        if search_improved != 1
            poll_result = poll(m, finternal, cinternal, incumbent, fincumbent)
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
        log_iteration!(logger, k, fincumbent, current_Δ, i)
        log_step!(logger, k, i == search_improved ? :search : :poll, i,
            search_fincumbent, fincumbent, current_Δ)

        # Debug logging (gated, per-iteration summary only)
        if should_log(logger, Debug)
            cache_size = m isa RobustMADS ? length(m.cache.y) : 0
            log_debug!(logger, k,
                directions_count=length(m.poll) + 1,  # Approximate
                cache_size=cache_size,
                mesh_level=ℓ(m.mesh))
        end

        # Check stopping criteria

        # Min mesh size
        if current_Δ < min_mesh_size
            result = (f=fincumbent, x=to(incumbent), stopping_reason=MinMeshSize,
                iterations=k, evaluations=eval_count, trace=trace)
            log_final!(logger, result)
            return result
        end

        # Max time
        if time() - start_time > max_time
            result = (f=fincumbent, x=to(incumbent), stopping_reason=MaxTime,
                iterations=k, evaluations=eval_count, trace=trace)
            log_final!(logger, result)
            return result
        end

        # F tolerance (stagnation)
        if f_tol > 0 && abs(fincumbent - prev_fincumbent) < f_tol && i <= 0
            result = (f=fincumbent, x=to(incumbent), stopping_reason=FTolReached,
                iterations=k, evaluations=eval_count, trace=trace)
            log_final!(logger, result)
            return result
        end

        # X tolerance
        if x_tol > 0 && norm(incumbent - prev_incumbent) < x_tol && i <= 0
            result = (f=fincumbent, x=to(incumbent), stopping_reason=XTolReached,
                iterations=k, evaluations=eval_count, trace=trace)
            log_final!(logger, result)
            return result
        end

        # Update previous values
        prev_incumbent = copy(incumbent)
        prev_fincumbent = fincumbent
    end

    result = (f=fincumbent, x=to(incumbent), stopping_reason=MaxIterations,
        iterations=max_iterations, evaluations=eval_count, trace=trace)
    log_final!(logger, result)
    return result
end
