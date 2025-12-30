# Search stage implementations for MADS

"""
    NoSearch

Default search strategy that performs no search step.
"""
struct NoSearch end

"""
    search(m, f, constraints, x, fx)

Execute search stage of MADS algorithm.
"""
search(m, f, constraints, x, fx) = search(m.search, f, constraints, x, fx)

"""
    search(::NoSearch, f, constraints, x, fx)

No-op search that returns current incumbent unchanged.
"""
search(::NoSearch, f, constraints, x, fx) = (incumbent=x, fx=fx, hasimproved=0)
