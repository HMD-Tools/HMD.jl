export Trajectory
export all_timesteps, get_timestep, is_reaction, get_system
export latest_reaction, latest_reaction_step, add!, update_reader!, add!
export setproperty!, iterate, getindex, length
export add_snapshot!, get_timesteps, get_reactions, get_metadata

#use case
# 原子位置情報の追跡
# 時系列propの追跡

#設計
# 1. buffer::systemを与えて時系列データのポインタだけ差し替え
# 2. Systemのbase関数を多重ディスパッチでtrajctoryに拡張 read_onlyならいける?
# -> 2を採用(2の内部で結局1と同じことを実行する)


Base.@kwdef mutable struct Trajectory{D, F<:AbstractFloat, SysType<:AbstractSystemType} <: AbstractTrajectory{D, F, SysType}
    # Vector of System. Properties which does not changes at evety timestep are empty except for reactions.
    systems::Vector{System{D, F, SysType}} = System{D, F, SysType}[]
    # indices corresponding to reactions (not timestep!)
    is_reaction::Vector{Int64} = Vector{Int64}(undef, 0)
    # index -> timestep
    timesteps::Vector{Int64} = Vector{Int64}(undef, 0)
end

function Trajectory(s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return Trajectory{D, F, SysType}([s], [1], [1])
end

function empty_trajectory(s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return Trajectory{D, F, SysType}()
end

function is_reaction(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return index ∈ traj.is_reaction
end

function get_system(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return traj.systems[index]
end

function all_timesteps(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return traj.timesteps
end

function get_timestep(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return all_timesteps(traj)[index]
end

function prop(traj::Trajectory, index::Integer, pname::AbstractString)
    return prop(get_system(traj, index), pname)
end

function Base.length(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return length(traj.systems)
end

function add!(traj::Trajectory{D, F, SysType}, s::System{D, F, SysType}, timestep::Integer; reaction=false) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
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
        replica.props = s.props
        push!(traj.systems, replica)
    end

    return nothing
end

function add!(
    traj::Trajectory{D, F, SysType}, addend::Trajectory{D, F, SysType}
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
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
    reader.props = s.props

    return nothing
end

function import_static!(reader::System{D, F, S1}, traj::Trajectory{D, F, S2}, index::Integer) where {D, F<:AbstractFloat, S1<:AbstractSystemType, S2<:AbstractSystemType}
    s = get_system(traj, index)
    reader.topology = s.topology
    reader.hierarchy = s.hierarchy
    reader.element = s.element

    return nothing
end

function latest_reaction(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
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

function Base.similar(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return Trajectory{D, F, SysType}()
end

function similar_system(traj::Trajectory{D, F, SysType}; reserve_dynamic=false, reserve_static=false) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    s = System{D, F, SysType}()

    if reserve_dynamic
        import_dynamic!(s, traj, 1)
    elseif reserve_static
        import_static!(s, traj, 1)
    end

    return deepcopy(s)
end

function dimension(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return D
end

function precision(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return F
end

function system_type(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return SysType
end

function wrapped(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    is_wrapped = wrapped(traj.systems[1])
    @assert all(i -> wrapped(traj.systems[i]) == is_wrapped, 1:length(traj))
    return is_wrapped
end

function wrap!(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if !wrapped(traj)
        for i in 1:length(traj)
            wrap!(traj.systems[i])
        end
    end

    return nothing
end

function unwrap!(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if wrapped(traj)
        for i in 1:length(traj)
            unwrap!(traj.systems[i])
        end
    end

    return nothing
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

function add_snapshot!(file_handler::H5traj, s::System{D, F, SysType}, step::Int64; reaction::Bool=false, unsafe::Bool=false) where{D, F<:AbstractFloat, SysType<:AbstractSystemType}
    file = get_file(file_handler)
    # construction
    if keys(file) |> isempty
        file["infotype"] = "Trajectory"
        file["dimension"] = dimension(s)
        file["precision"] = precision(s) |> string
        file["system_type"] = system_type(s) |> string
        file["wrapped"] = wrapped(s)
        create_group(file, "snapshots")
        create_group(file, "reactions")
    elseif !unsafe
        _error_chk(file, s; mode=mode)
    end

    # trajectory specific properties
    if reaction
        if nv(topology(s))==0 && isempty(all_elements(s)) # && isempty(hierarchy(s))
            error("system's topology, hierarychy and elements are empty. ")
        end
        create_group(file, "reactions/$step")
    end
    create_group(file, "snapshots/$step")
    snap = H5system(file["snapshots/$step"])
    hmdsave(snap, s)

    return nothing
end

function import_dynamic!(reader::System{D, F, S}, traj_file::H5traj; index::Int64=typemin(Int64), step::Int64=typemin(Int64)) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    if index != typemin(Int64) && step == typemin(Int64)
        step = get_timesteps(traj_file)[index]
    elseif index == typemin(Int64) && step != typemin(Int64)
        # do nothing
    else
        error("""only one keyword argument "index" or "step" must be specified. """)
    end

    file = get_file(traj_file)
    snap = H5system(file["snapshots/$step"])
    D_file, F_file, SysType_file = get_metadata(snap)
    if (D, F) != (D_file, F_file)
        println("(D, F) != (D_file, F_file): ($D, $F) != ($D_file, $F_file)")
        error("Dimension and precision of the reader and the trajectory file are different. ")
    end
    import_dynamic!(reader, snap)

    return nothing
end

function import_static!(reader::System{D, F, S}, traj_file::H5traj; index::Int64=typemin(Int64), step::Int64=typemin(Int64)) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    if index != typemin(Int64) && step == typemin(Int64)
        step = get_timesteps(traj_file)[index]
    elseif index == typemin(Int64) && step != typemin(Int64)
        # do nothing
    else
        error("""only one keyword argument "index" or "step" must be specified. """)
    end

    file = get_file(traj_file)
    snap = H5system(file["snapshots/$step"])
    D_file, F_file, SysType_file = get_metadata(snap)
    if (D, F) != (D_file, F_file)
        println("(D, F) != (D_file, F_file): ($D, $F) != ($D_file, $F_file)")
        error("Dimension and precision of the reader and the trajectory file are different. ")
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

function _error_chk(file, reader::System{D, F, SysType}; mode) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
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
