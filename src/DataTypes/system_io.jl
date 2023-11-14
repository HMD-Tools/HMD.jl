#####
##### HDF5 IO
#####

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
    s::System{D, F, SysType, L},
    precision = F;
    compress = false,
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if compress
        @warn "warning: compression is not supported yet. ignoring..."
    end

    file = get_file(file_handler)

    # metadata
    file["infotype"] = "System"
    file["dimension"] = dimension(s)
    file["precision"] = precision |> string
    file["system_type"] = system_type(s) |> string

    #data
    file["time"] = time(s)
    file["position"] = serialize(all_positions(s), precision)
    file["travel"] = serialize(all_travels(s))
    file["wrapped"] = wrapped(s)

    file["box/origin"] = Vector{F}(box(s).origin)
    file["box/axis"] = Matrix{F}(box(s).axis)

    file["element"] = all_elements(s)

    if haskey(file, "velocity")
        file["velocity"] = serialize(all_velocities(s), precision)
    end
    if haskey(file, "force")
        file["force"] = serialize(all_forces(s), precision)
    end

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

    return nothing
end

function read_system(system_file::H5system, template::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    # metadata check
    D_file, F_file, SysType_file = get_metadata(system_file)
    if D_file != D
        close(system_file)
        error("System type mismatch. (file: $(D_file), $(SysType_file), system: $(D), $(SysType)")
    elseif F_file > F
        @info "warning: system precision $(F) is lower than the file precision $(F_file). \n" *
            "This cause information loss."
    end

    # data
    s = similar(template)
    import_dynamic!(s, system_file)
    import_static!(s, system_file)
    return s
end

function import_dynamic!(
    s::System{D, F, SysType, L},
    system_file::H5system
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    file = get_file(system_file)
    set_time!(s, read(file, "time"))
    set_box!(s, BoundingBox{D, F}(read(file, "box/origin"), read(file, "box/axis")))
    s.position = let
        mat = read(file, "position")
        reinterpret(SVector{D, F}, reshape(mat, length(mat)))
    end
    s.travel = let
        mat = read(file, "travel")
        reinterpret(SVector{D, Int16}, reshape(mat, length(mat)))
    end
    s.wrapped = read(file, "wrapped")

    if haskey(file, "velocity")
        s.velocity = let
            mat = read(file, "velocity")
            reinterpret(SVector{D, F}, reshape(mat, length(mat)))
        end
    end
    if haskey(file, "force")
        s.velocity = let
            mat = read(file, "force")
            reinterpret(SVector{D, F}, reshape(mat, length(mat)))
        end
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
