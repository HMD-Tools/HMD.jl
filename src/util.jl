#function sort(eiter::AbstractEdgeIter)
#    edge_array = collect(eiter)
#    sort!(edge_array; by = x -> src(x))
#
#    current = src(edge_array[1])
#    index = 1
#    for (i, edge) in enumerate(edge_array)
#        if src(edge) != current
#            sub = view(edge_array, index:i-1)
#            sort!(sub; by = x -> dst(x))
#            index = i
#            current = src(edge)
#        end
#    end
#
#    return edge_array
#end

function oblique_coord(x::AbstractVector, origin, axis)
    # x = α * axis[:, 1] + β * axis[:, 2] + γ * axis[:, 3]
    #x_relative = similar(x)
    x_relative = x .- origin

    #return ([axis[:, 1], axis[:, 2], axis[:, 3]] \ x_oblique) .+ origin
    return (axis \ x_relative) .+ origin
end

# serialization for data storage

"""
    serialize(strings)

convert Vector{String} to chars::Vector{UInt8} and bounds::Vector{Int64}.
The first vector contains all the characters in the strings,
and the second vector contains the starting index of each string.

"""
function serialize(strings::Vector{String})
    if isempty(strings)
        return UInt8[], Int64[]
    end


    chars = Vector{UInt8}(undef, sum(length(str) for str in strings))
    bounds = Vector{Int64}(undef, length(strings))
    nchar = 1
    for i in 1:length(strings)
        bounds[i] = nchar
        for c in codeunits(strings[i])
            chars[nchar] = c
            nchar += 1
        end
    end
    push!(bounds, nchar)

    return chars, bounds
end

"""
    deserialize(chars, bounds)

convert chars::Vector{UInt8} and bounds::Vector{Int64} to Vector{string}.
`char` contains all the characters in the strings,
and `bounds` contains the starting index of each string in `chars`.
"""
function deserialize(chars::Vector{UInt8}, bounds::Vector{Int64})
    return map(1:length(bounds)-1) do i
        String(chars[bounds[i]:bounds[i+1]-1])
    end
end
