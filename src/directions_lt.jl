# LT Direction Generator
# Implements generation of directions, box on p. 203-204 from Audet & Dennis 2006

"""
    LTDirectionGenerator{N}

Generates LT (Lower Triangular) directions for MADS poll step.
See Audet & Dennis (2006), section 4, LTMADS.
"""
struct LTDirectionGenerator{N}
    b::Dict{Int,NamedTuple{(:i, :v),Tuple{Int,SVector{N,Int}}}}
end
LTDirectionGenerator(N) = LTDirectionGenerator{N}(Dict{Int,NamedTuple{(:i, :v),Tuple{Int,SVector{N,Int}}}}())

"""
    b(d::LTDirectionGenerator{N}, l::Int)

Generate direction b(l) for mesh level l.
"""
function b(d::LTDirectionGenerator{N}, l::Int) where N
    haskey(d.b, l) && return d.b[l]
    i = rand(1:N)
    d.b[l] = (i=i,
        v=SVector{N}([rand((-1, 1)) * (k == i ? 2^l : (2^l - 1)) for k in 1:N]))
end

"""
    iterator(g::LTDirectionGenerator, l)

Create an iterator over poll directions for mesh level l.
"""
iterator(g::LTDirectionGenerator, l) = LTDirectionIterator(g, max(0, l))

"""
    LTDirectionIterator{N,Np}

Iterator that generates positive spanning set using LT directions.
Implements generation of the positive basis, box on p. 204 from Audet & Dennis 2006.
"""
struct LTDirectionIterator{N,Np}
    l::Int
    i::Int
    b::SVector{N,Int}
    rowperm::SVector{Np,Int}
    colperm::SVector{Np,Int}
    sum::MVector{N,Int}
end

function LTDirectionIterator(generator::LTDirectionGenerator{N}, l) where N
    i, v = b(generator, l)
    LTDirectionIterator{N,N - 1}(l, i, v,
        SVector{N - 1}(randperm(N - 1)),
        SVector{N - 1}(randperm(N - 1)),
        zeros(Int, N))
end

"""
    L(i, j, l)

Lower triangular matrix element generator.
"""
function L(i, j, l)
    j > i && return 0
    j == i && return rand((-1, 1)) * 2^l
    return rand(-2^l+1:2^l-1)
end

import Base: iterate, length

length(::LTDirectionIterator{N,Np}) where {N,Np} = N + 1

function iterate(it::LTDirectionIterator{N,Np}, j=1) where {N,Np}
    j == N + 2 && return nothing
    j == N + 1 && return -SVector(it.sum), j + 1
    j == 1 && (it.sum .*= 0)
    if j == it.i
        v = it.b
    else
        v = SVector{N}([i == it.i ? 0 : L(it.rowperm[i-(i>it.i)],
            it.colperm[j-(j>it.i)], it.l)
                        for i in 1:N])
    end
    it.sum .+= v
    v, j + 1
end
