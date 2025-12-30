# OrthoMADS Direction Generator
# Implements OrthoMADS from Abramson et al. 2009

"""
    haltonnumber(base::Integer, index::Integer)::Float64

Compute the Halton sequence value at given index for given base.
"""
@inline function haltonnumber(base::Integer, index::Integer)::Float64
    res = 0.
    f = 1 / base
    i = index
    while i > 0
        res += f * (i % base)
        i = div(i, base)
        f = f / base
    end
    res
end

"""
    HaltonIterator{N}

Iterator for N-dimensional Halton sequence.
"""
struct HaltonIterator{N} end

function iterate(it::HaltonIterator{N}, i=1) where N
    [haltonnumber(prime(k), i) for k in 1:N], i + 1
end

"""
    normalized_halton_direction(u, l)

Compute normalized Halton direction for mesh level l.
"""
function normalized_halton_direction(u, l)
    α = 2^(l / 2) / sqrt(length(u)) - 1 / 2
    q = normalize(2 * u .- ones(length(u)))
    while norm(round.(α * q)) < 2^(l / 2)
        α += 0.1
    end
    round.(Int, (α - 0.1) * q)
end

"""
    scaledhouseholder(q)

Compute scaled Householder matrix from direction q.
"""
scaledhouseholder(q) = sum(q .^ 2) * I - 2 * q * q'

"""
    NoReduction()

No reduction strategy for OrthoMADS.
See Abramson et al. 2009.
"""
struct NoReduction end

"""
    NegReduction

(suc, neg) reduction strategy for OrthoMADS.
See Audet et al. 2014.
"""
struct NegReduction
    oldincumbent::Vector{Float64}
    w::Vector{Float64}
end
NegReduction(N) = NegReduction(zeros(N), fill(NaN, N))

"""
    OrthoDirectionGenerator{N,R}

OrthoMADS direction generator with reduction strategy R.
"""
mutable struct OrthoDirectionGenerator{N,R}
    reduction::R
    t₀::Int
    ℓmax::Int
    tmax::Int
end

function OrthoDirectionGenerator(N; t0=2 * N, reduction=NegReduction(N))
    OrthoDirectionGenerator{N,typeof(reduction)}(reduction, t0, 0, 0)
end

iterator(g::OrthoDirectionGenerator, l) = OrthoDirectionIterator(g, l)

init!(g::OrthoDirectionGenerator{N,NegReduction}, x) where N = g.reduction.oldincumbent .= x

function update!(g::OrthoDirectionGenerator{N,NegReduction}, x) where N
    @. g.reduction.w = x - g.reduction.oldincumbent
    @. g.reduction.oldincumbent = x
    return nothing
end

"""
    OrthoDirectionIterator{N,R}

Iterator over OrthoMADS poll directions.
"""
struct OrthoDirectionIterator{N,R}
    reduction::R
    H::Matrix{Float64}
end

function determine_t!(g::OrthoDirectionGenerator, ℓ)
    if ℓ >= g.ℓmax
        g.ℓmax = ℓ
        ℓ + g.t₀ > g.tmax && (g.tmax = ℓ + g.t₀)
        return ℓ + g.t₀
    else
        g.tmax += 1
        return g.tmax
    end
end

function OrthoDirectionIterator(g::OrthoDirectionGenerator{N,R}, ℓ) where {N,R}
    t = determine_t!(g, ℓ)
    u = first(iterate(HaltonIterator{N}(), t))
    q = normalized_halton_direction(u, ℓ)
    H = scaledhouseholder(q)
    if R === NegReduction && isnan(g.reduction.w[1])
        reduction = NoReduction()
    else
        reduction = g.reduction
    end
    OrthoDirectionIterator{N,typeof(reduction)}(reduction, H)
end

length(::OrthoDirectionIterator{N}) where N = 2N
length(::OrthoDirectionIterator{N,NegReduction}) where N = N + 1

function iterate(it::OrthoDirectionIterator{N}, i=1) where N
    i > 2N && return nothing
    i <= N && return @view(it.H[:, i]), i + 1
    return -@view(it.H[:, i-N]), i + 1
end

function iterate(it::OrthoDirectionIterator{N,NegReduction}, state=(1, zeros(N))) where N
    i, sumd = state
    i > N + 1 && return nothing
    if i <= N
        d = @view(it.H[:, i])
        sign = dot(d, it.reduction.w) >= 0 ? 1 : -1
        d .*= sign
        sumd .-= d
        return d, (i + 1, sumd)
    else
        return sumd, (i + 1, sumd)
    end
end
