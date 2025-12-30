# MADS Type Definitions
# Implements general MADS, p. 193 from Audet & Dennis, 2006

"""
    AbstractMADS

Abstract type for all MADS variants.
"""
abstract type AbstractMADS end

"""
    struct MADS{Tmesh,Tsearch,Tpoll} <: AbstractMADS

Standard MADS algorithm structure.

# Fields
- `mesh`: Mesh size manager (Mesh or LogMesh)
- `search`: Search strategy (e.g., NoSearch)
- `poll`: Poll direction generator (e.g., LTDirectionGenerator, OrthoDirectionGenerator)
"""
struct MADS{Tmesh,Tsearch,Tpoll} <: AbstractMADS
    mesh::Tmesh
    search::Tsearch
    poll::Tpoll
end

"""
    MADS(N; search=NoSearch(), poll=LTDirectionGenerator(N), mesh=LogMesh())

Construct a MADS optimizer for N-dimensional problems.
"""
function MADS(N; search=NoSearch(),
    poll=LTDirectionGenerator(N),
    mesh=LogMesh())
    MADS(mesh, search, poll)
end

"""
     LtMADS(N; search = NoSearch(), mesh = LogMesh())

Returns a `MADS` object with `poll = LTDirectionGenerator(N)` where `N` is the dimensionality of the problem.
See Audet & Dennis (2006), section 4, LTMADS.
"""
function LtMADS(N; search=NoSearch(), mesh=LogMesh())
    MADS(mesh, search, LTDirectionGenerator(N))
end

"""
    OrthoMADS(N; search = NoSearch(), mesh = LogMesh(), reduction = NegReduction(N))

Returns a `MADS` object with `poll = OrthoDirectionGenerator(N)` where `N` is
the dimensionality of the problem and `reduction` can be
[`NegReduction(N)`](@ref) or [`NoReduction()`](@ref).
See Abramson et al. (2009), ORTHOMADS and Audet et al. 2014 for NegReduction.
"""
function OrthoMADS(N; search=NoSearch(),
    mesh=LogMesh(),
    reduction=NegReduction(N))
    MADS(mesh, search, OrthoDirectionGenerator(N, reduction=reduction))
end

"""
    struct RobustMADS{Tmesh,Tsearch,Tpoll,Tkernel} <: AbstractMADS

Robust MADS for noisy objective functions.

# Fields
- `mesh`: Mesh size manager
- `search`: Search strategy
- `poll`: Poll direction generator
- `kernel`: Gaussian kernel for smoothing
- `cache`: Cache of evaluated points
- `f`: Smoothed function values
- `P`: Kernel density estimates
"""
struct RobustMADS{Tmesh,Tsearch,Tpoll,Tkernel} <: AbstractMADS
    mesh::Tmesh
    search::Tsearch
    poll::Tpoll
    kernel::Tkernel
    cache::Cache
    f::Vector{Float64}
    P::Vector{Float64}
end

"""
    RobustMADS(N; search = NoSearch(), poll = LTDirectionGenerator(N),
                  mesh = LogMesh(), kernel = GaussKernel(1, 1), cache = Cache(N))

Returns a `RobustMADS` object where `N` is the dimensionality of the problem.
See Audet et al. 2018.
"""
function RobustMADS(N; search=NoSearch(), poll=LTDirectionGenerator(N),
    mesh=LogMesh(), kernel=GaussKernel(1, 1), cache=Cache(N))
    RobustMADS(mesh, search, poll, kernel, cache, Float64[], Float64[])
end

"""
    RobustLtMADS(N; kwargs...)

Returns a `RobustMADS` object with `poll = LTDirectionGenerator(N)` where `N` is the dimensionality of the problem.
"""
function RobustLtMADS(N; kwargs...)
    RobustMADS(N; poll=LTDirectionGenerator(N), kwargs...)
end

"""
    RobustOrthoMADS(N; kwargs...)

Returns a `RobustMADS` object with `poll = OrthoDirectionGenerator(N)` where `N` is the dimensionality of the problem.
"""
function RobustOrthoMADS(N; kwargs...)
    RobustMADS(N; poll=OrthoDirectionGenerator(N), kwargs...)
end
