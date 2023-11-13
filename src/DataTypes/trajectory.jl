export Trajectory, SubTrajectory
export all_timesteps, get_timestep, is_reaction, get_system
export latest_reaction, latest_reaction_step, add!, update_reader!, add!
export setproperty!, length
export add_snapshot!, get_timesteps, get_reactions, get_metadata

Base.@kwdef mutable struct Trajectory{D, F<:AbstractFloat, SysType<:AbstractSystemType, L} <: AbstractTrajectory{D, F, SysType}
    # Vector of System. Properties which does not changes at evety timestep are empty except for reactions.
    systems::Vector{System{D, F, SysType, L}} = System{D, F, SysType, L}[]
    # indices corresponding to reactions (not timestep!)
    is_reaction::Vector{Int64} = Vector{Int64}(undef, 0)
    # index -> timestep
    timesteps::Vector{Int64} = Vector{Int64}(undef, 0)
end

function Trajectory(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return Trajectory{D, F, SysType, L}([s], [1], [1])
end

function Base.show(
    io::IO, ::MIME"text/plain", traj::Trajectory{D, F, S, L}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    "Trajectory{$D, $F, $S, $L}
        length: $(length(traj))
        time: $(time(traj.systems[1])) to $(time(traj.systems[end]))
    " |> println
end

function empty_trajectory(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return Trajectory{D, F, SysType, L}()
end

function is_reaction(traj::Trajectory{D, F, SysType, L}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return index ∈ traj.is_reaction
end

function get_system(traj::Trajectory{D, F, SysType, L}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return traj.systems[index]
end

function all_timesteps(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return traj.timesteps
end

function get_timestep(traj::Trajectory{D, F, SysType, L}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return all_timesteps(traj)[index]
end

function getindex(traj::AbstractTrajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if !(0 < index <= length(traj))
        throw(BoundsError(traj, index))
    end

    replica = similar_system(traj)
    # set properties that changes only at reaction
    rp = get_system(traj, latest_reaction(traj, index))
    replica.element = all_elements(rp) |> deepcopy
    replica.topology = topology(rp) |> deepcopy
    replica.hierarchy = deepcopy(rp.hierarchy)

    # set properties that changes at every step
    current = get_system(traj, index)
    set_time!(replica, time(current))
    set_box!(replica, deepcopy(box(current)))
    replica.position = all_positions(current) |> deepcopy
    replica.travel = deepcopy(current.travel)

    # others
    replica.wrapped = current.wrapped

    return replica
end

#function prop(traj::Trajectory, index::Integer, pname::AbstractString)
#    return prop(get_system(traj, index), pname)
#end

function Base.length(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return length(traj.systems)
end

function Base.iterate(traj::Trajectory{D, F, S, L}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    index = 1
    reader = similar_system(traj)
    import_dynamic!(reader, traj, index)
    import_static!(reader, traj, index)

    return (step=get_timestep(traj, index), snap=reader), (index+1, reader)
end

function Base.iterate(
    traj::Trajectory{D, F, SysType, L},
    state::Tuple{Int64, S}
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, S<:AbstractSystem{D, F, SysType}, L}
    index = state[1]
    if index <= length(traj)
        reader = state[2]
        rp = latest_reaction(traj, index)
        import_static!(reader, traj, rp)
        import_dynamic!(reader, traj, index)
        return (step=get_timestep(traj, index), snap=reader), (index+1, reader)
    else
        return nothing
    end
end

function add!(traj::Trajectory{D, F, SysType, L}, s::System{D, F, SysType, L}, timestep::Integer; reaction=false) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    @assert length(traj.systems) == length(traj.timesteps)
    if length(traj) > 0
        @assert traj.is_reaction[end] <= length(traj)
    else
        push!(traj.systems, s)
        push!(traj.is_reaction, 1)
        push!(traj.timesteps, timestep)
        return nothing
    end

    push!(traj.timesteps, timestep)
    if reaction
        if nv(topology(s))==0 && isempty(all_elements(s)) # && isempty(hierarchy(s))
            error("system's topology, hierarchy, and elements are empty. ")
        end
        push!(traj.systems, s)
        push!(traj.is_reaction, length(traj.systems))
    else
        replica = System{D, F, SysType}()
        set_box!(replica, box(s))
        set_time!(replica, time(s))
        replica.position = all_positions(s)
        replica.travel = s.travel
        replica.wrapped = s.wrapped
        push!(traj.systems, replica)
    end

    return nothing
end

function add!(
    traj::Trajectory{D, F, SysType, L}, addend::Trajectory{D, F, SysType, L}
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if topology(traj.systems[1]) != topology(addend.systems[1])
        error("topology of traj and addend are not compatible. ")
    end

    endstep = all_timesteps(traj)[end]
    nsnap = length(addend.systems)
    append!(traj.systems, view(addend.systems, 2:nsnap))
    append!(traj.is_reaction, addend.is_reaction[2:end])
    append!(traj.timesteps, addend.timesteps[2:end] .+ endstep)

    endtime = time(traj.systems[end])
    for i in nsnap+1:length(traj)
        traj.systems[i].time += endtime
    end

    return nothing
end

function import_dynamic!(reader::System{D, F, S1}, traj::Trajectory{D, F, S2}, index::Integer) where {D, F<:AbstractFloat, S1<:AbstractSystemType, S2<:AbstractSystemType}
    s = get_system(traj, index)
    set_time!(reader, time(s))
    set_box!(reader, box(s))
    reader.position = all_positions(s)
    reader.travel = s.travel
    reader.wrapped = s.wrapped

    return nothing
end

function import_static!(reader::System{D, F, S1}, traj::Trajectory{D, F, S2}, index::Integer) where {D, F<:AbstractFloat, S1<:AbstractSystemType, S2<:AbstractSystemType}
    s = get_system(traj, index)
    reader.topology = s.topology
    reader.hierarchy = s.hierarchy
    reader.element = s.element

    return nothing
end

function latest_reaction(traj::Trajectory{D, F, SysType, L}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if !(1 <= index <= length(traj))
        error("index $index ∉ [1, $(length(traj))]")
    end

    # find `i` s.t. traj.is_reaction[i] <= index < traj.is_reaction[i+1]
    ii = searchsortedlast(traj.is_reaction, index)
    return traj.is_reaction[ii]
end

#function is_reaction(s::System)
#    return all_elements(s) |> isempty
#end

function similar(
    traj::Trajectory{D, F, SysType, L};
    precision::Union{Type{<:AbstractFloat}, Nothing} = nothing
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    T = isnothing(precision) ? F : precision
    return Trajectory{D, T, SysType, L}()
end

function similar_system(
    traj::Trajectory{D, F, SysType, L};
    reserve_dynamic = false,
    reserve_static = false,
    precision::Union{Type{<:AbstractFloat}, Nothing} = nothing
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    T = isnothing(precision) ? F : precision
    s = System{D, T, SysType}()

    if reserve_dynamic
        import_dynamic!(s, traj, 1)
    elseif reserve_static
        import_static!(s, traj, 1)
    end

    return deepcopy(s)
end

function dimension(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return D
end

function precision(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return F
end

function system_type(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return SysType
end

function wrapped(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    is_wrapped = wrapped(traj.systems[1])
    @assert all(i -> wrapped(traj.systems[i]) == is_wrapped, 1:length(traj))
    return is_wrapped
end

function wrap!(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if !wrapped(traj)
        for i in 1:length(traj)
            wrap!(traj.systems[i])
        end
    end

    return nothing
end

function unwrap!(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if wrapped(traj)
        for i in 1:length(traj)
            unwrap!(traj.systems[i])
        end
    end

    return nothing
end



#####
##### SubTrajectory types
#####

struct SubTrajectory{D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange} <: AbstractTrajectory{D, F, S}
    traj::Trajectory{D, F, S, L}
    traj_range::R
end

function Base.show(
    io::IO, ::MIME"text/plain", st::SubTrajectory{D, F, S, L, R}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    if isempty(st)
        println("empty SubTrajectory{$D, $F, $S, $L}")
        return nothing
    end

    start = st.traj_range[1]
    final = st.traj_range[end]
    "SubTrajectory{$D, $F, $S, $L} with range $(st.traj_range)
        length: $(length(st))
        time: $(time(st.traj.systems[start])) to $(time(st.traj.systems[final]))
    " |> println

    return nothing
end

function Base.getindex(
    traj::Trajectory{D, F, S, L},
    traj_range::OrdinalRange{I, I}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, I<:Integer}
    if !(traj_range ⊆ 1:length(traj))
        error("traj_range must be a subset of $(1:length(traj)). found: $(traj_range)")
    end

    return SubTrajectory(traj, traj_range)
end

function Base.getindex(
    st::SubTrajectory{D, F, S, L, R},
    idx::Integer
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    real_index = st.traj_range[idx]
    return st.traj[real_index]
end

function Base.iterate(st::SubTrajectory{D, F, S, L, R}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    isempty(st.traj_range) && return nothing

    pseude_idx = 1
    real_idx = st.traj_range[pseude_idx]
    reader = similar_system(st.traj)
    import_dynamic!(reader, st.traj, real_idx)
    import_static!(reader, st.traj, real_idx)

    return (step=get_timestep(st.traj, real_idx), snap=reader), (pseude_idx+1, reader)
end

function Base.iterate(
    st::SubTrajectory{D, F, SysType, L, R},
    state::Tuple{Int64, S}
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, S<:AbstractSystem{D, F, SysType}, L, R<:OrdinalRange}
    pseude_idx = state[1]
    if pseude_idx <= length(st.traj_range)
        real_idx = st.traj_range[pseude_idx]
        reader = state[2]
        rp = latest_reaction(st.traj, real_idx)
        import_static!(reader, st.traj, rp)
        import_dynamic!(reader, st.traj, real_idx)
        return (step=get_timestep(st.traj, real_idx), snap=reader), (pseude_idx+1, reader)
    else
        return nothing
    end
end

function is_reaction(st::SubTrajectory{D, F, S, L, R}, index::Integer) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    return is_reaction(st.traj, st.traj_range[index])
end

function get_system(st::SubTrajectory{D, F, S, L, R}, index::Integer) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    return get_system(st.traj, st.traj_range[index])
end

function all_timesteps(st::SubTrajectory{D, F, S, L, R}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    return all_timesteps(st.traj)[st.traj_range]
end

function get_timestep(st::SubTrajectory{D, F, S, L, R}, index::Integer) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    return get_timestep(st.traj, st.traj_range[index])
end

function Base.length(st::SubTrajectory{D, F, S, L, R}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    return length(st.traj_range)
end

function wrapped(st::SubTrajectory{D, F, S, L, R}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    return wrapped(st.traj)
end



#####
##### Trajectory HDF5 types
#####

mutable struct H5traj <: AbstractFileFormat
    file::Union{HDF5.File, HDF5.Group}
end

function h5traj(name::AbstractString, mode::AbstractString)
    file_handler = H5traj(h5open(name, mode))
    file = get_file(file_handler)
    if mode != "w" && read(file, "infotype") != "Trajectory"
        close(file_handler)
        error("file $(name) is not a Trajectory file. ")
    end

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
    step::Int64,
    precision::Type{<:AbstractFloat} = F;
    reaction::Bool = false,
    unsafe::Bool = false,
) where{D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    file = get_file(file_handler)
    if haskey(file, "snapshots/$step")
        close(file)
        error("snapshot at $step already exists. ")
    end

    # construction
    if keys(file) |> isempty
        file["infotype"] = "Trajectory"
        file["dimension"] = dimension(s)
        file["precision"] = string(precision) #precision(s) |> string
        file["system_type"] = system_type(s) |> string
        file["wrapped"] = wrapped(s)
        create_group(file, "snapshots")
        create_group(file, "reactions")
    elseif !unsafe
        _error_chk(file, s; mode=mode)
    end

    # trajectory specific properties
    create_group(file, "snapshots/$step")
    snap = H5system(file["snapshots/$step"])
    if reaction
        if nv(topology(s))==0 && isempty(all_elements(s)) # && isempty(hierarchy(s))
            error("system's topology, hierarychy and elements are empty. ")
        end
        create_group(file, "reactions/$step")
        hmdsave(snap, s, precision)
    else
        hmdsave(
            snap,
            similar(s, reserve_dynamic=true, reserve_static=false),
            precision
        )
    end

    return nothing
end

function import_dynamic!(
    reader::System{D, F, S, L},
    traj_file::H5traj;
    index::Int64=typemin(Int64),
    step::Int64=typemin(Int64),
    unsafe = false
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    if index != typemin(Int64) && step == typemin(Int64)
        step = get_timesteps(traj_file)[index]
    elseif index == typemin(Int64) && step != typemin(Int64)
        # do nothing
    else
        error("""only one keyword argument "index" or "step" must be specified. """)
    end

    file = get_file(traj_file)
    snap = H5system(file["snapshots/$step"])
    if !unsafe
        D_file, F_file, SysType_file = get_metadata(snap)
        if (D, F) != (D_file, F_file)
            println("(D, F) != (D_file, F_file): ($D, $F) != ($D_file, $F_file)")
            error("Dimension and precision of the reader and the trajectory file are different. ")
        end
    end
    import_dynamic!(reader, snap)

    return nothing
end

function import_static!(
    reader::System{D, F, S, L},
    traj_file::H5traj;
    index::Int64=typemin(Int64),
    step::Int64=typemin(Int64),
    unsafe = false
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    if index != typemin(Int64) && step == typemin(Int64)
        step = get_timesteps(traj_file)[index]
    elseif index == typemin(Int64) && step != typemin(Int64)
        # do nothing
    else
        error("""only one keyword argument "index" or "step" must be specified. """)
    end

    file = get_file(traj_file)
    snap = H5system(file["snapshots/$step"])
    if !unsafe
        D_file, F_file, SysType_file = get_metadata(snap)
        if (D, F) != (D_file, F_file)
            println("(D, F) != (D_file, F_file): ($D, $F) != ($D_file, $F_file)")
            error("Dimension and precision of the reader and the trajectory file are different. ")
        end
    end
    import_static!(reader, snap)

    return nothing
end

function latest_reaction_step(traj_file::H5traj, current_step::Integer)
    # timestep at reaction (not index!)
    reaction_steps = get_reactions(traj_file)
    i = searchsortedlast(reaction_steps, current_step)

    return reaction_steps[i]
end

function get_timesteps(traj_file::H5traj)
    file = get_file(traj_file)
    return parse.(Int64, keys(file["snapshots"])) |> sort!
end

function get_reactions(traj_file::H5traj)
    file = get_file(traj_file)

    return parse.(Int64, keys(file["reactions"])) |> sort!
end

function get_metadata(traj_file::H5traj)
    file = get_file(traj_file)
    D = read(file, "dimension")
    F = read(file, "precision") |> Symbol |> eval
    SysType = read(file, "system_type")

    return D, F, SysType
end

function _error_chk(file, reader::System{D, F, SysType, L}; mode) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if read(file, "infotrype") != "Trajectory"
        error("The file is not a trajectory file. ")
    elseif read(file, "dimension") != dimension(reader)
        error("The dimension of the system is different from the dimension of the trajectory. ")
    elseif read(file, "precision") != string(precision(reader))
        error("The precision of the system is different from the precision of the trajectory. ")
    elseif read(file, "system_type") != string(system_type(reader))
        error("The system type of the system is different from the system type of the trajectory. ")
    elseif read(file, "mode") != mode
        error("mode mistmatch")
    end

    if wrapped(reader)
        error("currently supports only unwrapped coordinates. ")
    end
end

function is_reaction(traj_file::H5traj, index::Integer)
    file = get_file(traj_file)
    return read(file, "reactions/$index")
end

function Base.length(traj_file::H5traj)
    file = get_file(traj_file)
    return length(keys(file["snapshots"]))
end

function wrapped(traj_file::H5traj)
    file = get_file(traj_file)
    return read(file, "wrapped")
end

function similar_system(traj_file::H5traj)
    D, F, SysType = get_metadata(traj_file)
    SysType = SysType |> Symbol |> eval
    return System{D, F, SysType}()
end






#####
##### fast prototype Trajectory
#####
"""

structure:
    dimension
    precision
    system_type
    1d props
        times
        timesteps
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

export H5traj_new, h5traj_new
mutable struct H5traj_new <: AbstractFileFormat
    file::HDF5.File
end

function h5traj_new(name::AbstractString, mode::AbstractString)
    if mode == "w" || (mode == "cw" && !isfile(name))
        error("if you want to create new HMD file, use h5traj_new(name, mode, template, traj_length)")
    end

    return H5traj_new(h5open(name, mode))
end

function h5traj_new(
    name::AbstractString,
    mode::AbstractString,
    template::System{D, F, S, L},
    traj_length::Integer,
    precision::Type{<:AbstractFloat} = F
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    file_handler = H5traj_new(h5open(name, mode))
    if mode in ("r", "w+") || (mode == "cw" && isfile(name))
        error("if you want to read or append existing HMD file, use h5traj_new(name, mode)")
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
    create_dataset(
        file,
        "timesteps",
        datatype(Int64),
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

function close(file_handler::H5traj_new)
    close(get_file(file_handler))
end

function get_file(file_handler::H5traj_new)
    return file_handler.file
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

function _get_snap(props, i, Natom, D, F)
    if D <= 0
        error("dimension D must be >= 1. ")
    end

    start = D * (i-1) * Natom + 1
    stop = D * i * Natom
    return reinterpret(SVector{D, F}, props[start:stop, 1])
end

# hmdsave_new(traj) and read_traj_new()
function add_snapshot_new!(
    file_handler::H5traj_new,
    s::System{D, F, SysType, L},
    index::Integer;
    unsafe::Bool = false,
    reaction::Bool = false
) where{D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if index < 1
        close(file_handler)
        error("index must be >= 1. ")
    end

    file = get_file(file_handler)

    change_wrap = false
    if wrapped(s)
        @warn "wrapping is not supported yet. unwrapping..."
        change_wrap = true
    end

    # metadata
    if !unsafe
        #_error_chk(file, s; mode=mode)
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
    if haskey(file, "hierarchy/$index")
        close(file)
        error("hierarchy at $index already exists. ")
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

struct H5trajReader{F<:AbstractFloat}
    file::HDF5.File
    times::Vector{F}
    timesteps::Vector{Int64}
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

function h5traj_reader(name::AbstractString, F::Type{<:AbstractFloat})
    file_handler = h5traj_new(name, "r")
    file = get_file(file_handler)

    times = let
        x::Matrix{F} = read(file, "times")
        reshape(x, length(x))
    end
    timesteps = let
        x::Matrix{Int64} = read(file, "timesteps")
        reshape(x, length(x))
    end
    elems::Vector{Atomic_Number_Precision} = read(file, "elements")
    reactions::Vector{Int64} = let
        entries::Base.KeySet{String, Dict{String, Any}} = keys(read(file, "topology"))
        sort!([parse(Int64, str) for str in entries])
    end
    boxes = SerializedBoxes(
        read(file, "box_origin")::Matrix{F},
        read(file, "box_axis")::Matrix{F}
    )
    Nsnap = length(timesteps)
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

    @assert length(timesteps) == length(times)
    @assert length(reactions) < length(times)

    return H5trajReader{F}(file, times, timesteps, elems, reactions, boxes, Nsnap, Natom, has_v, has_f, hnames, atomprops)
end

function read_traj_new(
    name::AbstractString,
    D::Integer, F::Type{<:AbstractFloat}, S::Type{<:AbstractSystemType}
)
    hmdfile = h5traj_reader(name, F)
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
    hmdfile::H5trajReader{F}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    traj.timesteps = hmdfile.timesteps
    traj.is_reaction = hmdfile.reactions
    traj.systems = [similar_system(traj) for _ in 1:hmdfile.Nsnap]
    # dynamic properties
    for i in 1:hmdfile.Nsnap
        print("progress: $(100*i÷hmdfile.Nsnap)%    \r")
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
    for (i, rp) in enumerate(hmdfile.reactions)
        traj.systems[rp].topology = _read_topology(hmdfile.file, i)
        for hname in hmdfile.hnames
            traj.systems[rp].hierarchy[hname] = _read_hierarchy(hmdfile.file, hname, i)
        end
        traj.systems[rp].element = hmdfile.elems
        @assert traj.is_reaction[i] == hmdfile.reactions[i]
    end
    println()

    return nothing
end

function import_dynamic!(
    reader::System{D, F, S, L},
    traj_file::H5traj_new,
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

function import_static!(
    reader::System{D, F, S, L},
    traj_file::H5traj_new,
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

function get_metadata(traj_file::H5traj_new)
    file = get_file(traj_file)
    D = read(file, "dimension")
    F = read(file, "precision") |> Symbol |> eval
    SysType = read(file, "system_type")

    return D, F, SysType
end

#function _error_chk(file, reader::System{D, F, SysType, L}; mode) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
#    if read(file, "infotrype") != "Trajectory"
#        error("The file is not a trajectory file. ")
#    elseif read(file, "dimension") != dimension(reader)
#        error("The dimension of the system is different from the dimension of the trajectory. ")
#    elseif read(file, "precision") != string(precision(reader))
#        error("The precision of the system is different from the precision of the trajectory. ")
#    elseif read(file, "system_type") != string(system_type(reader))
#        error("The system type of the system is different from the system type of the trajectory. ")
#    elseif read(file, "mode") != mode
#        error("mode mistmatch")
#    end
#
#    if wrapped(reader)
#        error("currently supports only unwrapped coordinates. ")
#    end
#end
