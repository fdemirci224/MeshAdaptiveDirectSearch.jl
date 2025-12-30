# Mesh implementations for MADS
# Implements Eq. 2.1 from Audet & Dennis 2006

"""
    Mesh(; τ=4., Δᵐ=1., w⁺=1, w⁻=-1)

Standard mesh structure with random mesh size updates.
"""
mutable struct Mesh
    τ::Float64
    w⁺::Int
    w⁻::Int
    Δᵐ::Float64
end
Mesh(; τ=4., Δᵐ=1., w⁺=1, w⁻=-1) = Mesh(τ, w⁺, w⁻, Δᵐ)

"""
    update!(m::Mesh, i)

Update mesh size based on poll result `i`:
- i > 0: success, potentially increase mesh size
- i < 0: failure, decrease mesh size
- i == 0: no change
"""
function update!(m::Mesh, i)
    i == 0 && return m
    if i > 0
        m.Δᵐ == 1. && return m
        w = rand(0:m.w⁺)
    else
        w = rand(m.w⁻:-1)
    end
    m.Δᵐ *= m.τ^w
    m
end

"""
    Δ(m::Mesh)

Return current mesh size parameter.
"""
Δ(m::Mesh) = m.Δᵐ

"""
    ℓ(m::Mesh)

Return mesh index (log-scale level).
"""
ℓ(m::Mesh) = round(Int, -log(m.τ, m.Δᵐ))

"""
    LogMesh(; τ=4, Δᵐ=1)

Efficient log-space mesh for w⁺ = -w⁻ = 1.
Stores mesh level directly instead of computing logarithms.
"""
mutable struct LogMesh
    τ::Int
    neglogΔᵐ::Int
end
LogMesh(; τ=4, Δᵐ=1) = LogMesh(τ, Int(-log(τ, Δᵐ)))

"""
    update!(m::LogMesh, i)

Update log-mesh based on poll result.
"""
function update!(m::LogMesh, i)
    i == 0 && return m
    if i > 0
        m.neglogΔᵐ -= 1
    else
        m.neglogΔᵐ += 1
    end
    m
end

Δ(m::LogMesh) = min(1, (1 / m.τ)^m.neglogΔᵐ)
ℓ(m::LogMesh) = m.neglogΔᵐ
