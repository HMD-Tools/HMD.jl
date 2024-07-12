#####
##### Trajectory HDF5 types
#####

function serialize(vec::Vector{SVector{D, T}}, Tconvert=T) where {D, T<:Real}
    #return [vec[atom_id][dim] for dim in 1:D, atom_id in eachindex(vec)]
    spos = Matrix{Tconvert}(undef, D, length(vec))
    for atom_id in eachindex(vec)
        spos[:, atom_id] = vec[atom_id]
    end
    return spos
end


function _get_snap(
    props::HDF5.Dataset,
    i::Integer,     # snapshot index (1-origin)
    Natom::Integer, # number of atoms in each snapshots
    D::Integer,     # system dimension
    F::Type{<:AbstractFloat}, # system precision
    chunk::UnitRange{T} = 1:Natom
) where {T<:Integer}
    if D <= 0
        error("dimension D must be >= 1, found $D. ")
    elseif chunk.stop ≤ chunk.start
        error("chunk.stop must be >= chunk.start, found $chunk.")
    elseif chunk.start < 1 || chunk.stop > Natom
        error("chunk must be within 1:$Natom, found $chunk. ")
    end

    start = 1 + (D * (i-1) * Natom) + D * (chunk.start - 1)
    stop  = start + D * length(chunk) - 1
    #@assert (start - 1) % 3 == 0 "$start % 3 != 0"
    #@assert stop % 3 == 0 "$stop % 3 != 0"

    # Never remove ::Vector{F} because
    # reinterpret(::Type{SVector{D, F}}, arr::Vector{F1}) will
    # return invalid result if F1 != F
    # e.g. try `reinterpret(SVector{3, Int32}, [1,2,3])`
    arr::Vector{F} = props[start:stop, 1]

    # padding zeros if chunk is not full
    if chunk.start != 1
        prepend!(
            arr,
            zeros(F, D * (chunk.start - 1))
        )
    end
    if chunk.stop != Natom
        append!(
            arr,
            zeros(F, D * (Natom - chunk.stop))
        )
    end
    @assert length(arr) == D * Natom

    return reinterpret(SVector{D, F}, arr)
end



# TODO: IOの中間メモリ確保量の低減(Serialized typeがなくなってロジックだけになるので可読性確保の工夫が必要)
function _get_substructure(props::HDF5.Dataset, label::HLabel, Natom::Integer, D::Integer, F::Type{<:AbstractFloat})
    if D <= 0
        error("dimension D must be >= 1. ")
    end

    start = D * (i-1) * Natom + 1
    stop = D * i * Natom

    arr::Vector{F} = props[start:stop, 1]
    return reinterpret(SVector{D, F}, arr)
end

struct SerializedTopology
    num_node::Int64
    edges_org::Vector{Int64}
    edges_dst::Vector{Int64}
    denominator::Vector{BO_Precision}
    numerator::Vector{BO_Precision}
end

function serialize(topo::AbstractSimpleWeightedGraph)
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
    #topo = SimpleWeightedGraph{Int64, Rational{BO_Precision}}()
    topo = TopologyGraph{Int64, Rational{BO_Precision}}()
    add_vertices!(topo, ser_topo.num_node)
    for i in 1:length(ser_topo.edges_org)
        add_edge!(
            topo, ser_topo.edges_org[i],
            ser_topo.edges_dst[i],
            Rational{BO_Precision}(ser_topo.numerator[i], ser_topo.denominator[i])
        )
    end

    return topo
end

struct SerializedBoxes{F<:AbstractFloat}
    origins::Vector{F}
    axes::Vector{F}
end

function SerializedBoxes(origins::Matrix{F}, axes::Matrix{F}) where {F<:AbstractFloat}
    if size(origins, 2) == size(axes, 2) == 1
        return SerializedBoxes(
            reshape(origins, length(origins)),
            reshape(axes, length(axes))
        )
    else
        error(
            "origins and axes must be vector-like matrix. \n" *
            "origins: $(size(origins)) \n" *
            "axes: $(size(axes))"
        )
    end
end

function _get_snap(boxes::SerializedBoxes{F}, i::Integer, D::Integer) where {F<:AbstractFloat}
    if D <= 0
        error("dimension D must be >= 1. ")
    elseif length(boxes.origins) % D != 0 || length(boxes.axes) % (D*D) != 0
        error("length of boxes.origins and boxes.axes must be multiples of D and D*D, respectively. ")
    end

    #println(SMatrix{D, D, F, D*D}(boxes.axes[i] for i in D*D*(i-1)+1:D*D*i))
    return BoundingBox{D, F}(
        SVector{D, F}(boxes.origins[i] for i in D*(i-1)+1:D*i),
        SMatrix{D, D, F, D*D}(boxes.axes[i] for i in D*D*(i-1)+1:D*D*i)
    )
end

"""

structure:
    dimension
    precision
    system_type
    1d props
        times
        name1
        name2
    position
        step_i <- chunk?
    velocity
        step_i <- chunk?
    forces
        step_i <- chunk?
    box
        origin <- chunk?
        axis  <- chunk?
    topology
    element
    hierarchy
"""

mutable struct H5traj <: AbstractFileFormat
    file::HDF5.File
end

function h5traj(name::AbstractString, mode::AbstractString)
    if mode == "w" || (mode == "cw" && !isfile(name))
        error("if you want to create new HMD file, use h5traj(name, mode, template, traj_length)")
    end

    return H5traj(h5open(name, mode))
end

function h5traj(
    name::AbstractString,
    mode::AbstractString,
    template::System{D, F, S, L},
    traj_length::Integer,
    precision::Type{<:AbstractFloat} = F
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    file_handler = H5traj(h5open(name, mode))
    if mode in ("r", "r+") || (mode == "cw" && isfile(name))
        error("if you want to read or append existing HMD file, use h5traj(name, mode)")
    end
    file = get_file(file_handler)

    # metadata
    file["dimension"] = D
    file["precision"] = string(precision)
    file["system_type"] = string(S)
    file["natom"] = natom(template)

    # elements
    file["elements"] = all_elements(template)

    # scalar properties s.t. time, timestep, energy...
    # values of these properties are expected to be written after simulation done
    create_group(file, "1d_props")
    create_dataset(
        file,
        "times",
        datatype(precision),
        dataspace(traj_length, 1);
        chunk = (traj_length, 1)
    )

    # make chunks for each step for atomic & dynamic properties
    # chunk size is min(datasize, 16MB)
    chunk_size = min(D * natom(template) * traj_length, 16384000 ÷ sizeof(precision))
    for pname in ("position", "velocity", "force")
        isempty(getfield(template, Symbol(pname))) && continue
        create_dataset(
            file,
            "$pname",
            datatype(precision),
            dataspace(D*natom(template)*traj_length, 1);
            chunk = (chunk_size, 1)
        )
    end

    # simulation box
    create_dataset(
        file,
        "box_origin",
        datatype(precision),
        dataspace(D*traj_length, 1);
        chunk = (D*traj_length, 1)
    )
    create_dataset(
        file,
        "box_axis",
        datatype(precision),
        dataspace(D*D*traj_length, 1);
        chunk = (D*D*traj_length, 1)
    )

    # properties that changes only at reaction
    # dataset is created every reaction because the final size of the dataset is unknown
    create_group(file, "topology")
    create_group(file, "hierarchy")

    return file_handler
end

function close(file_handler::H5traj)
    close(get_file(file_handler))
end

function get_file(file_handler::H5traj)
    return file_handler.file
end

function add_snapshot!(
    file_handler::H5traj,
    s::System{D, F, SysType, L},
    index::Integer;
    reaction::Bool = false,
    unsafe::Bool = false
) where{D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if index < 1
        close(file_handler)
        error("index must be >= 1. ")
    end
    # metadata check
    if !unsafe
        _error_chk(file_handler, s)
    end

    file = get_file(file_handler)
    if isempty(file["topology"])
        reaction = true
    end

    change_wrap = false
    if wrapped(s)
        @warn "wrapping is not supported yet. unwrapping..."
        change_wrap = true
    end

    if reaction
        if nv(topology(s))==0 && isempty(all_elements(s))
            error("system's topology and elements are empty. ")
        end
        _add_topology!(file, topology(s), index)
        for hname in hierarchy_names(s)
            lh = hierarchy(s, hname)
            _add_hierarchy!(file, hname, lh, index)
        end
    end
    for pname in ("position", "velocity", "force")
        prop = getproperty(s, Symbol(pname))
        file_has_prop = haskey(file, pname)
        if isempty(prop) && !file_has_prop
            continue
        elseif !isempty(prop) && file_has_prop
            _add_atomprop!(file, pname, D, natom(s), prop, index)
        elseif isempty(prop) && file_has_prop
            error("file has $pname entry, but the system does not have $pname property. ")
        elseif !isempty(prop) && !file_has_prop
            error("file does not have $pname entry, but the system has $pname property. ")
        else
            error("unknown error occured. ")
        end
    end
    _add_box!(file, box(s), index, D)

    change_wrap && wrap!(s)

    return nothing
end

function _add_topology!(file, topo, index)
    if haskey(file, "topology/$index")
        close(file)
        error("topology at $index already exists. ")
    end

    stopo = serialize(topo)
    for fname in fieldnames(SerializedTopology)
        file["topology/$index/$(fname)"] = getfield(stopo, fname)
    end

    return nothing
end

function _add_hierarchy!(file, hname, lh, index)
    if !haskey(file, "hierarchy/$index")
        create_group(file, "hierarchy/$index")
    end

    ser_hierarchy = serialize(lh)
    for fname in fieldnames(PackedHierarchy)
        file["hierarchy/$index/$hname/$(fname)"] = getfield(ser_hierarchy, fname)
    end

    return nothing
end

function _add_atomprop!(file, pname, D, natom, prop, index)
    @assert length(prop) == natom
    @assert length(prop[1]) == D
    start = D * (index-1) * natom + 1
    stop = D * index * natom
    @assert length(start:stop) == D * natom
    file["$(pname)"][start:stop, 1] = reshape(stack(prop), D*natom)

    return nothing
end

function _add_box!(file, bbox, index, D)
    origin_range = D*(index-1)+1:D*index
    axis_range = D*D*(index-1)+1:D*D*index
    @assert length(origin_range) == D
    @assert length(axis_range) == D*D
    file["box_origin"][origin_range, 1] = bbox.origin
    file["box_axis"][axis_range, 1] = reshape(bbox.axis, D*D)

    return nothing
end

struct H5trajReader{D, F, S, L}
    file::HDF5.File
    times::Vector{F}
    elems::Vector{Atomic_Number_Precision}
    reactions::Vector{Int64}
    boxes::SerializedBoxes{F}
    Nsnap::Int64
    Natom::Int64
    has_velocity::Bool
    has_force::Bool
    hnames::Vector{String}
    atomprops::Dict{String, HDF5.Dataset}
end

function h5traj_reader(name::AbstractString, D::Integer, F::Type{<:AbstractFloat}, S::Type{<:AbstractSystemType})
    file_handler = h5traj(name, "r")
    file = get_file(file_handler)

    times = let
        x::Matrix{F} = read(file, "times")
        reshape(x, length(x))
    end
    #timesteps = let
    #    x::Matrix{Int64} = read(file, "timesteps")
    #    reshape(x, length(x))
    #end
    elems::Vector{Atomic_Number_Precision} = read(file, "elements")
    reactions::Vector{Int64} = let
        # if reaction is occured at step i, the dataset "topology/i" exists
        entries::Base.KeySet{String, Dict{String, Any}} = keys(read(file, "topology"))
        sort!([parse(Int64, str) for str in entries])
    end
    boxes = let
        origin::Matrix{F} = read(file, "box_origin")
        axis::Matrix{F} = read(file, "box_axis")
        SerializedBoxes(origin, axis)
    end
    Nsnap = length(times)
    Natom::Int64 = read(file, "natom")
    has_v = haskey(file, "velocity")
    has_f = haskey(file, "force")
    hnames = let
        entries::Base.KeySet{String, Dict{String, Any}} = keys(read(file, "hierarchy/1"))
        [str for str in entries]
    end
    atomprops = Dict{String, HDF5.Dataset}()
    for pname in ("position", "velocity", "force")
        if haskey(file, pname)
            atomprops[pname] = open_dataset(
                file,
                pname;
                chunk_cache = (521, UInt32(16384000), 0.75) # 16MB cache
            )
        end
    end

    #@assert length(timesteps) == length(times)
    @assert length(reactions) <= length(times)
    @assert allunique(reactions)

    return H5trajReader{D, F, S, D*D}(file, times, elems, reactions, boxes, Nsnap, Natom, has_v, has_f, hnames, atomprops)
end

eachindex(hmdfile::H5trajReader) = 1:hmdfile.Nsnap
length(hmdfile::H5trajReader) = hmdfile.Nsnap
natom(hmdfile::H5trajReader) = hmdfile.Natom
get_times(hmdfile::H5trajReader) = deepcopy(hmdfile.times) # must be deepcopy for compatibility with SubH5trajReader
firstindex(hmdfile::H5trajReader) = 1
lastindex(hmdfile::H5trajReader) = hmdfile.Nsnap


function read_traj(
    name::AbstractString,
    D::Integer, F::Type{<:AbstractFloat}, S::Type{<:AbstractSystemType}
)
    hmdfile = h5traj_reader(name, D, F, S)
    traj = Trajectory{D, F, S, D*D}()
    try
        read_traj!(traj, hmdfile)
    finally
        close(hmdfile.file)
    end

    return traj
end

function read_traj!(
    traj::Trajectory{D, F, S, L},
    hmdfile::H5trajReader{D, F, S, L}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    #traj.timesteps = hmdfile.timesteps
    traj.reactions = hmdfile.reactions
    traj.systems = [similar_system(traj) for _ in 1:hmdfile.Nsnap]
    # dynamic properties
    intval = hmdfile.Nsnap ÷ 100
    for i in 1:hmdfile.Nsnap
        if i % intval == 0
            print("progress: $(100*i÷hmdfile.Nsnap)%    \r")
        end
        set_time!(traj.systems[i], hmdfile.times[i])
        set_box!(traj.systems[i], _get_snap(hmdfile.boxes, i, D))
        traj.systems[i].travel = zeros(SVector{D, Int16}, hmdfile.Natom)
        traj.systems[i].wrapped = false
        traj.systems[i].position = _get_snap(hmdfile.atomprops["position"], i, hmdfile.Natom, D, F)
        if hmdfile.has_velocity
            traj.systems[i].velocity = _get_snap(hmdfile.atomprops["velocity"], i, hmdfile.Natom, D, F)
        end
        if hmdfile.has_force
            traj.systems[i].force = _get_snap(hmdfile.atomprops["force"], i, hmdfile.Natom, D, F)
        end
    end

    # static properties
    for rp in hmdfile.reactions
        traj.systems[rp].topology = _read_topology(hmdfile.file, rp)
        for hname in hmdfile.hnames
            traj.systems[rp].hierarchy[hname] = _read_hierarchy(hmdfile.file, hname, rp)
        end
        traj.systems[rp].element = hmdfile.elems
    end
    println()

    @assert length(traj.systems) == hmdfile.Nsnap
    @assert allunique(traj.reactions)
    @assert issorted(traj.reactions)

    return nothing
end


function Base.iterate(hmdfile::H5trajReader{D, F, S, L}) where {D, F, S, L}
    reader = System{D, F, S}()
    read_snapshot!(reader, hmdfile, 1)

    return reader, (2, reader)
end


function Base.iterate(
    hmdfile::H5trajReader{D, F, S, L},
    state::Tuple{Int64, System{D, F, S}}
) where {D, F, S, L}
    index, reader = state
    if index > hmdfile.Nsnap
        return nothing
    else
        read_snapshot!(reader, hmdfile, index)
        return reader, (index+1, reader)
    end
end


function read_snapshot!(
    reader::System{D, F, S, L},
    hmdfile::H5trajReader{D, F, S, L},
    index::Integer,
    chunk::UnitRange{T} = 1:natom(hmdfile)
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, T<:Integer}

    is_read_firsttime = if isempty(all_elements(reader)) && nv(topology(reader))==0 && isempty(hierarchy_names(reader))
        true
    elseif !isempty(all_elements(reader)) && nv(topology(reader))!=0 && !isempty(hierarchy_names(reader))
        false
    else
        error("some of static properties are empty. trajectory may be broken. ")
    end

    # static property
    atreaction = searchsorted(hmdfile.reactions, index)
    @assert length(atreaction) ≤ 1
    if atreaction == index:index || is_read_firsttime
        ii = if is_read_firsttime
            j = findfirst(ireac -> index < ireac, hmdfile.reactions)
            isnothing(j) ? hmdfile.reactions[end] : j
        else
            atreaction[1]
        end

        reader.topology = _read_topology(hmdfile.file, ii)
        for hname in hmdfile.hnames
            reader.hierarchy[hname] = _read_hierarchy(hmdfile.file, hname, ii)
        end
        reader.element = deepcopy(hmdfile.elems)
    end

    #dynamic property
    set_time!(reader, hmdfile.times[index])
    set_box!(reader, _get_snap(hmdfile.boxes, index, D))
    reader.travel = zeros(SVector{D, Int16}, hmdfile.Natom)
    reader.wrapped = false
    reader.position = _get_snap(hmdfile.atomprops["position"], index, hmdfile.Natom, D, F, chunk)
    if hmdfile.has_velocity
        reader.velocity = _get_snap(hmdfile.atomprops["velocity"], index, hmdfile.Natom, D, F, chunk)
    end
    if hmdfile.has_force
        reader.force = _get_snap(hmdfile.atomprops["force"], index, hmdfile.Natom, D, F, chunk)
    end

    return reader
end


struct SubH5trajReader{D, F, S, L, R<:OrdinalRange}
    h5traj::H5trajReader{D, F, S, L}
    time_range::R
end


function read_snapshot!(
    reader::System{D, F, S, L},
    hmdfile::SubH5trajReader{D, F, S, L, R},
    pseudo_idx::Integer,
    chunk::UnitRange{T} = 1:natom(hmdfile)
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, T<:Integer, R}
    real_idx = hmdfile.time_range[pseudo_idx]
    return read_snapshot!(reader, hmdfile.h5traj, real_idx, chunk)
end


function Base.getindex(hmdfile::H5trajReader{D, F, S, L}, time_range::R) where {D, F, S, L, R<:OrdinalRange}
    return SubH5trajReader{D, F, S, L, R}(hmdfile, time_range)
end

get_times(hmdfile::SubH5trajReader) = hmdfile.h5traj.times[hmdfile.time_range]
length(hmdfile::SubH5trajReader) = length(hmdfile.time_range)
natom(hmdfile::SubH5trajReader) = hmdfile.h5traj.Natom
time_range(hmdfile::SubH5trajReader) = hmdfile.time_range


function Base.iterate(hmdfile::SubH5trajReader{D,F,S,L,R}) where {D,F,S,L,R}
    reader = System{D, F, S}()
    read_snapshot!(reader, hmdfile.h5traj, first(hmdfile.time_range))

    return reader, (firstindex(hmdfile.time_range)+1, reader)
end


function Base.iterate(
    hmdfile::SubH5trajReader{D, F, S, L, R},
    state::Tuple{Int64, System{D, F, S}}
) where {D, F, S, L, R}
    # pusedo index is the index of hmdfile.time_range
    pseudo_idx, reader = state

    if pseudo_idx > lastindex(hmdfile.time_range)
        return nothing
    else
        # real index is the index of the whole trajectory
        real_idx = time_range(hmdfile)[pseudo_idx]

        read_snapshot!(reader, hmdfile.h5traj, real_idx)
        return reader, (pseudo_idx+1, reader)
    end
end


struct ChunkedH5trajReader{D, F, S, L, R}
    h5traj::SubH5trajReader{D, F, S, L, R}
    atoms::UnitRange{Int64}
end

function partial_read(hmdfile::SubH5trajReader{D, F, S, L, R}, atoms) where {D, F, S, L, R}
    return ChunkedH5trajReader(hmdfile, atoms)
end

function partial_read(hmdfile::H5trajReader{D, F, S, L}, atoms) where {D, F, S, L}
    sub_hmdfile = hmdfile[eachindex(hmdfile)]
    return ChunkedH5trajReader(sub_hmdfile, atoms)
end


function Base.iterate(chunked_traj::ChunkedH5trajReader{D, F, S, L, R}) where {D, F, S, L, R}
    reader = System{D, F, S}()
    hmdfile = chunked_traj.h5traj
    read_snapshot!(reader, hmdfile, first(time_range(hmdfile)), chunked_traj.atoms)

    return reader, (firstindex(hmdfile.time_range)+1, reader)
end


function Base.iterate(
    chunked_traj::ChunkedH5trajReader{D, F, S, L, R},
    state::Tuple{Int64, System{D, F, S}}
) where {D, F, S, L, R}
    hmdfile = chunked_traj.h5traj

    # pusedo index is the index of hmdfile.time_range
    pseudo_idx, reader = state

    if pseudo_idx > lastindex(hmdfile.time_range)
        return nothing
    else
        # real index is the index of the whole trajectory
        real_idx = time_range(hmdfile)[pseudo_idx]

        read_snapshot!(reader, hmdfile, real_idx, chunked_traj.atoms)
        return reader, (pseudo_idx+1, reader)
    end
end


function _import_dynamic!(
    reader::System{D, F, S, L},
    traj_file::H5traj,
    index::Integer,
    Natom::Integer;
    has_velocity::Bool = false,
    has_force::Bool = false,
    unsafe::Bool = false
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    if wrapped(reader)
        error("currently supports only unwrapped coordinates. ")
    end
    if !unsafe
        D_file, F_file, SysType_file = get_metadata(traj_file)
        if (D, F) != (D_file, F_file)
            println("(D, F) != (D_file, F_file): ($D, $F) != ($D_file, $F_file)")
            error("Dimension and precision of the reader and the trajectory file are different. ")
        end
    end

    file = get_file(traj_file)
    reader.position = _read_atomprop(file, "position", D, F, Natom, index)
    if has_velocity
        reader.velocity = _read_atomprop(file, "velocity", D, F, Natom, index)
    end
    if has_force
        reader.force = _read_atomprop(file, "force", D, F, Natom, index)
    end
    reader.travel = zeros(SVector{D, Int16}, Natom)
    reader.wrap = false
    # bounding box, time, and other 1d properties are not read here
    # because they should be read at once for all steps

    return nothing
end

function _import_static!(
    reader::System{D, F, S, L},
    traj_file::H5traj,
    index::Integer;
    unsafe::Bool = false
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    if !unsafe
        D_file, F_file, SysType_file = get_metadata(traj_file)
        if (D, F) != (D_file, F_file)
            println("(D, F) != (D_file, F_file): ($D, $F) != ($D_file, $F_file)")
            error("Dimension and precision of the reader and the trajectory file are different. ")
        end
    end

    file = get_file(traj_file)
    s.element = read(file, "element")
    reader.topology = _read_topology(file, index)
    for hname in hierarchy_names(reader)
        reader.hierarchy[hname] = _read_hierarchy(file, hname, index)
    end

    return nothing
end

function _read_topology(file, index)
    topo = SerializedTopology(
        read(file, "topology/$index/num_node"),
        read(file, "topology/$index/edges_org"),
        read(file, "topology/$index/edges_dst"),
        read(file, "topology/$index/denominator"),
        read(file, "topology/$index/numerator")
    )

    return deserialize(topo)
end

function _read_hierarchy(file, hname, index)
    ph = PackedHierarchy(
        read(file, "hierarchy/$index/$hname/num_node"),
        read(file, "hierarchy/$index/$hname/edges_org"),
        read(file, "hierarchy/$index/$hname/edges_dst"),
        read(file, "hierarchy/$index/$hname/label_ids"),
        read(file, "hierarchy/$index/$hname/chars"),
        read(file, "hierarchy/$index/$hname/bounds")
    )

    return deserialize(ph)
end

function get_metadata(traj_file::H5traj)
    file = get_file(traj_file)
    D = read(file, "dimension")
    F = read(file, "precision") |> Symbol |> eval
    SysType = read(file, "system_type")

    return D, F, SysType
end

function _error_chk(traj_file::H5traj, reader::System{D, F, S, L}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    D_file, F_file, S_file = get_metadata(traj_file)
    if read(file, "infotrype") != "Trajectory"
        error("The file is not a trajectory file. ")
    elseif D_file != D
        error("The dimension of the system is different from the dimension of the trajectory. ")
    elseif F_file != F
        error("The precision of the system is different from the precision of the trajectory. ")
    elseif S_file == S
        error("The system type of the system is different from the system type of the trajectory. ")
    end
end
