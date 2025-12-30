# Problem struct for builder-style API

"""
    Problem(f, n::Int)

A mutable configuration object for defining optimization problems.

This is a **builder** — mutate via `set_initial!`, `set_bounds!`, `add_constraint!`.
**NOT** an immutable mathematical object.

# Example
```julia
p = Problem(rosenbrock, 2)
set_initial!(p, [0.0, 0.0])
set_bounds!(p, -5.0, 5.0)
add_constraint!(p, x -> x[1] + x[2] < 1.0)
result = solve(p)
```
"""
mutable struct Problem
    f::Function
    n::Int
    x0::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    constraints::Vector{Function}
end

"""
    Problem(f, n::Int)

Create a new optimization problem for function `f` with `n` dimensions.
Defaults: x0 = zeros(n), lb = -ones(n), ub = ones(n), no constraints.
"""
function Problem(f::Function, n::Int)
    Problem(f, n, zeros(n), -ones(n), ones(n), Function[])
end

"""
    set_initial!(p::Problem, x0)

Set the initial point for optimization.
"""
function set_initial!(p::Problem, x0::AbstractVector{<:Real})
    length(x0) == p.n || throw(DimensionMismatch("x0 must have length $(p.n), got $(length(x0))"))
    p.x0 = collect(Float64, x0)
    return p
end

"""
    set_bounds!(p::Problem, lb, ub)

Set bounds for all dimensions.
- If `lb` is a scalar, applies to all dimensions.
- If `lb` is a vector, must have length `n`.
"""
function set_bounds!(p::Problem, lb::Real, ub::Real)
    p.lb = fill(Float64(lb), p.n)
    p.ub = fill(Float64(ub), p.n)
    return p
end

function set_bounds!(p::Problem, lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real})
    length(lb) == p.n || throw(DimensionMismatch("lb must have length $(p.n)"))
    length(ub) == p.n || throw(DimensionMismatch("ub must have length $(p.n)"))
    p.lb = collect(Float64, lb)
    p.ub = collect(Float64, ub)
    return p
end

"""
    set_bounds!(p::Problem, i::Int, lb, ub)

Set bounds for a single dimension `i`.
"""
function set_bounds!(p::Problem, i::Int, lb::Real, ub::Real)
    1 <= i <= p.n || throw(BoundsError("Dimension $i out of range 1:$(p.n)"))
    p.lb[i] = Float64(lb)
    p.ub[i] = Float64(ub)
    return p
end

"""
    add_constraint!(p::Problem, c::Function)

Add a constraint function. Constraint must return `true` when satisfied.
"""
function add_constraint!(p::Problem, c::Function)
    push!(p.constraints, c)
    return p
end

"""
    clear_constraints!(p::Problem)

Remove all constraints from the problem.
"""
function clear_constraints!(p::Problem)
    empty!(p.constraints)
    return p
end
