# Constraint handling for MADS

"""
    isvalid(constraints, x)

Check if point x satisfies all constraints.
Returns true if all constraints are satisfied, false otherwise.
"""
function isvalid(cs, x)
    for c in cs
        c(x) || return false
    end
    return true
end
