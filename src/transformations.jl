# Coordinate transformations for bound handling

"""
    standardtransformation(lowerbound, upperbound)

Create transformation functions to/from the standard [-1, 1] domain.

Returns named tuple with:
- `to`: transforms from [-1, 1] to [lowerbound, upperbound]
- `from`: transforms from [lowerbound, upperbound] to [-1, 1]
"""
function standardtransformation(lowerbound, upperbound)
    d = upperbound - lowerbound
    (to=x -> @.((x + 1) / 2 * (upperbound - lowerbound) + lowerbound),
        from=x -> @.(2 * (x - lowerbound) / (upperbound - lowerbound) - 1))
end
