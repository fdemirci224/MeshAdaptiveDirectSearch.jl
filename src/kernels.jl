# Gaussian Kernel for RobustMADS
# Implements kernel smoothing from Audet et al. 2018

"""
    GaussKernel(β, σ²)

Gaussian kernel for robust MADS smoothing.
- `β`: bandwidth scaling factor
- `σ²`: variance parameter (updated based on mesh size)
"""
mutable struct GaussKernel
    β::Float64
    σ²::Float64
end

"""
    (g::GaussKernel)(x, y)

Evaluate kernel between points x and y.
"""
(g::GaussKernel)(x, y) = g() * reshape(exp.(-sum((x .- y) .^ 2, dims=1) / (2 * g.σ²)), :)

"""
    (g::GaussKernel)()

Return kernel normalization constant.
"""
(g::GaussKernel)() = 1 / sqrt(2 * π * g.σ²)

"""
    update!kernel!(g::GaussKernel, mesh)

Update kernel variance based on current mesh size.
"""
update!kernel!(g::GaussKernel, mesh) = g.σ² = (g.β * Δ(mesh))^2
