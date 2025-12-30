# Professional Logging System for MADS
# Provides 5 verbosity levels with clean API (no println in main code)

"""
    Verbosity

Verbosity levels for MADS logging.
- `Silent`: No output
- `Final`: Print only termination summary (best f, x, stopping reason, iterations, evaluations, time)
- `Iter`: One line per iteration (iteration index, incumbent f, mesh size Δ, success flag, eval count)
- `Step`: Detailed stage-level info (search vs poll outcome, mesh update, success type)
- `Debug`: Diagnostics per iteration (direction count, cache size) - NOT per evaluation
"""
@enum Verbosity Silent Final Iter Step Debug

"""
    StoppingReason

Enum for optimization termination reasons.
"""
@enum StoppingReason MaxIterations MinMeshSize MaxEvaluations MaxTime FTolReached XTolReached FTargetReached

"""
    MADSLogger

Logger for MADS optimization with configurable verbosity.
"""
struct MADSLogger
    verbosity::Verbosity
    io::IO
    log_interval::Int
end

MADSLogger(; verbosity::Verbosity=Silent, io::IO=stdout, log_interval::Int=1) =
    MADSLogger(verbosity, io, log_interval)

# Check if logging at a given level is enabled
should_log(logger::MADSLogger, level::Verbosity) = Int(logger.verbosity) >= Int(level)

"""
    log_iteration!(logger, k, fx, Δ, success, eval_count)

Log iteration summary (Iter level and above).
Only logs every `log_interval` iterations.
"""
function log_iteration!(logger::MADSLogger, k::Int, fx::Float64, Δ::Float64, success::Int, eval_count::Int)
    should_log(logger, Iter) || return
    k % logger.log_interval != 0 && return

    success_str = success > 0 ? "✓" : success == 0 ? "○" : "✗"
    println(logger.io, "iter=$k  evals=$eval_count  f=$(@sprintf("%.6e", fx))  Δ=$(@sprintf("%.2e", Δ))  [$success_str]")
end

"""
    log_step!(logger, stage, outcome, details...)

Log detailed stage information (Step level and above).
"""
function log_step!(logger::MADSLogger, k::Int, stage::Symbol, improved::Int,
    search_fx::Float64, poll_fx::Float64, Δ::Float64)
    should_log(logger, Step) || return

    stage_str = stage == :search ? "SEARCH" : "POLL  "
    outcome = improved > 0 ? "success" : improved == 0 ? "cache" : "fail"

    println(logger.io, "  [$k] $stage_str → $outcome  search_f=$(@sprintf("%.4e", search_fx))  poll_f=$(@sprintf("%.4e", poll_fx))  Δ=$(@sprintf("%.2e", Δ))")
end

"""
    log_debug!(logger, k, directions_count, cache_size, mesh_level)

Log debug diagnostics (Debug level only).
Logs summary per iteration, NOT per evaluation.
"""
function log_debug!(logger::MADSLogger, k::Int;
    directions_count::Int=0,
    cache_size::Int=0,
    mesh_level::Int=0)
    should_log(logger, Debug) || return

    println(logger.io, "    [DEBUG] iter=$k  dirs=$directions_count  cache=$cache_size  ℓ=$mesh_level")
end

"""
    log_final!(logger, result)

Log final optimization result (Final level and above).
Includes evaluations count and total time.
"""
function log_final!(logger::MADSLogger, result::NamedTuple)
    should_log(logger, Final) || return

    # Format time nicely
    time_sec = haskey(result, :time) ? result.time : 0.0
    if time_sec < 60
        time_str = @sprintf("%.2f seconds", time_sec)
    elseif time_sec < 3600
        time_str = @sprintf("%.2f minutes", time_sec / 60)
    else
        time_str = @sprintf("%.2f hours", time_sec / 3600)
    end

    println(logger.io, "")
    println(logger.io, "═══════════════════════════════════════════════════════")
    println(logger.io, " MADS Optimization Complete")
    println(logger.io, "───────────────────────────────────────────────────────")
    println(logger.io, " Stopping reason : $(result.stopping_reason)")
    println(logger.io, " Iterations      : $(result.iterations)")
    println(logger.io, " Evaluations     : $(result.evaluations)")
    println(logger.io, " Total time      : $time_str")
    println(logger.io, " Best f(x)       : $(@sprintf("%.10e", result.f))")
    println(logger.io, " Best x          : $(result.x)")
    println(logger.io, "═══════════════════════════════════════════════════════")
end

"""
    log_start!(logger, n, method_name)

Log optimization start (Iter level and above).
"""
function log_start!(logger::MADSLogger, n::Int, method_name::String)
    should_log(logger, Iter) || return

    println(logger.io, "")
    println(logger.io, "Starting MADS optimization (n=$n, method=$method_name)")
    println(logger.io, "───────────────────────────────────────────────────────")
end
