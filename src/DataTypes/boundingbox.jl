#####
##### Type `BoundingBox` definition
#####

struct BoundingBox{D, F <: AbstractFloat} <: AbstractBbox{D, F}
    origin::SVector{D, F}
    axis::SMatrix{D, D, F}
end

function BoundingBox{D, F}(origin::Vector{F}, axis::Matrix{F}) where {D, F<:AbstractFloat}
    d = det(axis)
    if d < 0
        error("left-handed system not allowed")
    elseif d == 0
        error("box is degenerated")
    end
    return BoundingBox{D, F}(SVector{D, F}(origin), SMatrix{D, D, F}(axis))
end

function BoundingBox{D, F}() where {D, F<:AbstractFloat}
    BoundingBox{D, F}(zeros(F, 3), Matrix{F}(I, D, D))
end
