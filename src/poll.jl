# Poll stage implementations for MADS

"""
    update!(::Any, ::Any)

Default no-op update for poll direction generators.
"""
update!(::Any, ::Any) = nothing

"""
    init!(::Any, ::Any)

Default no-op initialization for poll direction generators.
"""
init!(::Any, ::Any) = nothing

"""
    poll(m, f, constraints, x, fx)

Execute poll stage of MADS algorithm using current mesh settings.
"""
poll(m, f, constraints, x, fx) = poll(m, iterator(m.poll, ℓ(m.mesh)), Δ(m.mesh), f, constraints, x, fx)

"""
    isnewincumbent(m::MADS, x, fx, oldfx)

Determine if x is a new incumbent for standard MADS.
Returns (success_code, incumbent_point, incumbent_value).
"""
isnewincumbent(m::MADS, x, fx, oldfx) = fx < oldfx ? 1 : -1, x, fx

"""
    isnewincumbent(m::RobustMADS, x, fx, oldfx)

Determine if x is a new incumbent for RobustMADS using kernel smoothing.
Returns (success_code, incumbent_point, smoothed_value).

Success codes:
- 1: new incumbent found (current point is best)
- 0: cache success (different point in cache became best)
- -1: no improvement
"""
function isnewincumbent(m::RobustMADS, x, fx, oldfx)
    update!kernel!(m.kernel, m.mesh)
    if length(m.cache.x) > 0
        ψ = m.kernel(m.cache.x, x)
        m.f .*= m.P
        m.f .+= ψ * fx
        m.P .+= ψ
        m.f ./= m.P
        push!(m.P, sum(ψ) + m.kernel())
        push!(m.f, (dot(m.cache.y, ψ) + m.kernel() * fx) / m.P[end])
    else
        push!(m.P, m.kernel())
        push!(m.f, fx)
    end
    append!(m.cache.x, x)
    push!(m.cache.y, fx)
    i = argmin(m.f)
    success = i == length(m.f)
    cachesuccess = length(m.cache.incumbents) == 0 || i != m.cache.incumbents[end]
    success || cachesuccess && push!(m.cache.incumbents, i)
    success ? 1 : cachesuccess ? 0 : -1, m.cache.x[:, i], m.f[i]
end

"""
    poll(m, it, Δᵐ, f, constraints, x, fx)

Poll over directions from iterator `it` at mesh size `Δᵐ`.
Returns named tuple with incumbent, function value, and improvement flag.
"""
@inline function poll(m, it, Δᵐ, f, constraints, x, fx)
    for v in it
        newx = clamp!(x .+ Δᵐ * v, -1., 1.)
        isvalid(constraints, newx) || continue
        fnewx = f(newx)
        success, newincumbent, fnewx = isnewincumbent(m, newx, fnewx, fx)
        if success >= 0
            update!(m.poll, newincumbent)
            return (incumbent=newincumbent, fx=fnewx, hasimproved=success)
        end
    end
    (incumbent=x, fx=fx, hasimproved=-1)
end
