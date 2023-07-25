function prop_names(s::System)
    s.props |> keys |> collect
end

function prop(s::System, pname::AbstractString)
    s.props[pname]
end

function set_prop!(s::System, pname::AbstractString, p::AbstractArray{<:Integer})
    s.props[pname] = p
end

function set_prop!(s::System, pname::AbstractString, p::AbstractArray{Float32})
    s.props[pname] = p
end

function set_prop!(s::System, pname::AbstractString, p::AbstractArray{Float64})
    s.props[pname] = p
end

function set_prop!(s::System, pname::AbstractString, p::Union{Float64, Float32})
    s.props[pname] = [p]
end
