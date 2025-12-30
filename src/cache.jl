# Cache for RobustMADS
# Stores evaluated points and incumbents

"""
    Cache(N)

Cache for storing evaluated points in RobustMADS.
- `x`: matrix of evaluated points (N × num_points)
- `y`: vector of function values
- `incumbents`: indices of incumbent points
"""
struct Cache
    x::ElasticArray{Float64,2,1}
    y::Vector{Float64}
    incumbents::Vector{Int}
end

Cache(N) = Cache(ElasticArray{Float64}(undef, N, 0), Float64[], Int[])

import Base.push!

"""
    push!(c::Cache, x, y)

Add a new point x with value y to the cache.
"""
push!(c::Cache, x, y) = begin
    append!(c.x, x)
    push!(c.y, y)
end
