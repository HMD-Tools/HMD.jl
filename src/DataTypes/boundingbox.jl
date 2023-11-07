#####
##### Type `BoundingBox` definition
#####

struct BoundingBox{D, F <: AbstractFloat, L} <: AbstractBbox{D, F, L}
    origin::SVector{D, F}
    axis::SMatrix{D, D, F, L}
end

function BoundingBox{D, F}(origin::AbstractVector{T}, axis::AbstractMatrix{T}) where {D, F<:AbstractFloat, T<:Real}
    d = det(axis)
    if d < 0
        error("left-handed system not allowed")
    elseif d == 0
        error("box is degenerated")
    end
    return BoundingBox{D, F, D*D}(SVector{D, F}(origin), SMatrix{D, D, F, D*D}(axis))
end

function BoundingBox{D, F}() where {D, F<:AbstractFloat}
    BoundingBox{D, F, D*D}(zeros(F, 3), Matrix{F}(I, D, D))
end
