struct TopologyEdge{T<:Integer, U<:Real} <: AbstractSimpleWeightedEdge{T}
    src::T
    dst::T
    weight::U
end

#TopologyEdge(t::NTuple{2}) = TopologyEdge(t[1], t[2], one(Float64))

#function TopologyEdge{T,U}(t::NTuple{2}) where {T<:Integer,U<:Real}
#    return TopologyEdge(T(t[1]), T(t[2]), one(U))
#end

#TopologyEdge(p::Pair) = TopologyEdge(p.first, p.second, one(Float64))

function TopologyEdge{T, U}(p::Pair) where {T<:Integer, U<:Real}
    return TopologyEdge(T(p.first), T(p.second), one(U))
end


TopologyEdge(t::NTuple{3}) = TopologyEdge(t[1], t[2], t[3])

function TopologyEdge{T, U}(t::NTuple{3}) where {T<:Integer, U<:Real}
    return TopologyEdge(T(t[1]), T(t[2]), U(t[3]))
end

function TopologyEdge{T, U}(p::Pair{T1, T2}) where {T<:Integer, U<:Real, T1<:Integer, T2<:Integer}
    return TopologyEdge(p[1], p[2], one(U))
end

function TopologyEdge{T, U}(u::T1, v::T2) where {T<:Integer, T1<:Integer, T2<:Integer, U<:Real}
    return TopologyEdge(u, T1(v), one(U))
end

Base.eltype(::TopologyEdge{T}) where {T} = T

# Accessors
SimpleWeightedGraphs.src(e::TopologyEdge) = e.src
SimpleWeightedGraphs.dst(e::TopologyEdge) = e.dst
SimpleWeightedGraphs.weight(e::TopologyEdge) = e.weight

# I/O
function Base.show(io::IO, e::TopologyEdge)
    return print(io, "Edge $(e.src) => $(e.dst) with weight $(e.weight)")
end

# Conversions
Base.Tuple(e::TopologyEdge) = (src(e), dst(e), weight(e))

function (::Type{TopologyEdge{T, U}})(
    e::TopologyEdge
) where {T<:Integer, U<:Real}
    return TopologyEdge{T,U}(T(e.src), T(e.dst), U(e.weight))
end

# Convenience functions - note that these do not use weight.
Base.reverse(e::TopologyEdge{T, U}) where {T, U} = TopologyEdge{T, U}(dst(e), src(e), weight(e))

function Base.:(==)(e1::TopologyEdge, e2::TopologyEdge)
    return (src(e1) == src(e2) && dst(e1) == dst(e2))
end
