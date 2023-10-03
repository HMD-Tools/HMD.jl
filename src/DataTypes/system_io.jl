#####
##### HDF5 IO
#####

function serialize(vec::Vector{SVector{D, T}}) where {D, T<:Real}
    #return [vec[atom_id][dim] for dim in 1:D, atom_id in eachindex(vec)]
    spos = Matrix{T}(undef, D, length(vec))
    for atom_id in eachindex(vec)
        spos[:, atom_id] = vec[atom_id]
    end
    return spos
end

function deserialize(D::Integer, mat::Matrix{T}) where {T<:Real}
    pos = SVector{D, T}[]
    resize!(pos, size(mat, 2))
    for atom_id in eachindex(pos)
        pos[atom_id] = mat[:, atom_id]
    end
    return pos
end

struct SerializedTopology
    num_node::Int64
    edges_org::Vector{Int64}
    edges_dst::Vector{Int64}
    denominator::Vector{BO_Precision}
    numerator::Vector{BO_Precision}
end

function serialize(topo::SimpleWeightedGraph)
    num_node = nv(topo)
    edges_org = Vector{Int64}(undef, ne(topo))
    edges_dst = Vector{Int64}(undef, ne(topo))
    denominator = Vector{Int16}(undef, ne(topo))
    numerator = Vector{Int16}(undef, ne(topo))

    for (i, edge) in enumerate(edges(topo))
        edges_org[i], edges_dst[i] = src(edge), dst(edge)
        weight = get_weight(topo, edges_org[i], edges_dst[i])
        denominator[i], numerator[i] = weight.den, weight.num
    end

    return SerializedTopology(num_node, edges_org, edges_dst, denominator, numerator)
end

function deserialize(ser_topo::SerializedTopology)
    topo = SimpleWeightedGraph{Int64, Rational{BO_Precision}}()
    add_vertices!(topo, ser_topo.num_node)
    for i in 1:length(ser_topo.edges_org)
        add_edge!(topo, ser_topo.edges_org[i], ser_topo.edges_dst[i], Rational{BO_Precision}(ser_topo.numerator[i], ser_topo.denominator[i]))
        #add_edge!(topo, ser_topo.edges_org[i], ser_topo.edges_dst[i], ser_topo.numerator[i]//ser_topo.denominator[i])
    end

    return topo
end

mutable struct H5system <: AbstractFileFormat
    file::Union{HDF5.File, HDF5.Group}
end

function h5system(name::AbstractString, mode::AbstractString)
    file_handler = H5system(h5open(name, mode))
    file = get_file(file_handler)
    if mode != "w" && read(file, "infotype") != "System"
        close(file_handler)
        error("file $(name) is not a System file. ")
    end

    return file_handler
end

function close(file_handler::H5system)
    close(get_file(file_handler))
end

function get_file(file_handler::H5system)
    return file_handler.file
end

function hmdsave(
    file_handler::H5system,
    s::System{D, F, SysType, L};
    compress = false
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if compress
        @warn "warning: compression is not supported yet. ignoring..."
    end

    file = get_file(file_handler)

    # metadata
    file["infotype"] = "System"
    file["dimension"] = dimension(s)
    file["precision"] = precision(s) |> string
    file["system_type"] = system_type(s) |> string

    #data
    file["time"] = time(s)
    file["position"] = serialize(all_positions(s))
    file["travel"] = serialize(all_travels(s))
    file["wrapped"] = wrapped(s)

    file["box/origin"] = Vector(box(s).origin)
    file["box/axis"] = Matrix(box(s).axis)

    #chars, bounds = serialize(all_elements(s))
    #file["element/chars"]  = chars
    #file["element/bounds"] = bounds
    file["element"] = all_elements(s)

    # topology
    stopo = serialize(topology(s))
    for fname in fieldnames(SerializedTopology)
        file["topology/$(fname)"] = getfield(stopo, fname)
    end

    # label hierarchy
    file["hierarchy_names"] = hierarchy_names(s)
    for hname in hierarchy_names(s)
        ser_hierarchy = serialize(hierarchy(s, hname))
        for fname in fieldnames(PackedHierarchy)
            file["hierarchy/$hname/$(fname)"] = getfield(ser_hierarchy, fname)
        end
    end

    # properties
    file["property_names"] = prop_names(s)
    for pname in prop_names(s)
        file["props/$pname"] = get_prop(s, pname)
    end

    return nothing
end

function read_system(system_file::H5system, template::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    # metadata check
    D_file, F_file, SysType_file = get_metadata(system_file)
    if (D_file, F_file) != (D, F)
        close(system_file)
        error("System type mismatch. (file: $(D_file), $(F_file), $(SysType_file), system: $(D), $(F), $(SysType)")
    end

    # data
    s = similar(template)
    import_dynamic!(s, system_file)
    import_static!(s, system_file)
    return s
end

function import_dynamic!(s::System{D, F, SysType, L}, system_file::H5system) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    file = get_file(system_file)
    set_time!(s, read(file, "time"))
    set_box!(s, BoundingBox{D, F}(read(file, "box/origin"), read(file, "box/axis")))
    s.position = deserialize(D, read(file, "position"))
    s.travel = deserialize(D, read(file, "travel"))
    s.wrapped = read(file, "wrapped")

    for pname in read(file, "property_names")
        set_prop!(s, pname, read(file, "props/$pname"))
    end

    return nothing
end

function import_static!(s::System{D, F, SysType, L}, system_file::H5system) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    file = get_file(system_file)
    #s.element = deserialize(read(file, "element/chars"), read(file, "element/bounds"))
    s.element = read(file, "element")
    s.topology = SerializedTopology(read(file, "topology/num_node"),
                                    read(file, "topology/edges_org"),
                                    read(file, "topology/edges_dst"),
                                    read(file, "topology/denominator"),
                                    read(file, "topology/numerator")) |> deserialize
    for hname in read(file, "hierarchy_names")
        if hname ∉ hierarchy_names(s)
            add_hierarchy!(s, hname)
        end
        ph = PackedHierarchy(read(file, "hierarchy/$hname/num_node"),
                            read(file, "hierarchy/$hname/edges_org"),
                            read(file, "hierarchy/$hname/edges_dst"),
                            read(file, "hierarchy/$hname/label_ids"),
                            read(file, "hierarchy/$hname/chars"),
                            read(file, "hierarchy/$hname/bounds"))
        s.hierarchy[hname] = deserialize(ph)
    end

    return nothing
end

function get_metadata(system_file::H5system)
    file = get_file(system_file)
    D = read(file, "dimension")
    F = read(file, "precision") |> Symbol |> eval
    SysType = read(file, "system_type")

    return D, F, SysType
end
