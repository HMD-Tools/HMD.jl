#####
##### Type `Position` definition
#####

const Position{D, F} = Vector{SVector{D, F}} where {D, F <: AbstractFloat}

function Position{D, F}(natom::Integer) where {D, F}
    return Vector(SVector{D, F})(undef, natom)
end

function Position{D, F}() where {D, F}
    return SVector{D, F}[]
end

function Position{D, F}(x::Matrix{F}) where {D, F}
    if size(x, 1) != D
        error("Dimension mismatch")
    end
    #return [SVector{D, F}(x[:, i]) for i in 1:size(x, 2)]
    return [SVector{D, F}(v...) for v in 1:eachcol(x)]
end

function Position()
    Position{3, Float64}()
end

function Position(n::Integer)
    Position{3, Float64}(n)
end
