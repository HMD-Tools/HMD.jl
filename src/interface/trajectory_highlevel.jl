"""
highlevelの引数をすべて具象型にすべきか？
    ・具象型で書き直すと新しい具象型を定義したときに再実装量が増える
    ・抽象型のままにしておく場合はinterfaceでNIを定義する意味が無い
iteratorをlowlevelに落とすことを検討
"""

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
    set_time!(replica, time(replica))
    set_box!(replica, deepcopy(box(current)))
    replica.position = all_positions(current) |> deepcopy
    replica.travel = deepcopy(current.travel)
    replica.props = deepcopy(current.props)

    # others
    replica.wrapped = current.wrapped

    return replica
end

function Base.iterate(traj::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    index = 1
    reader = similar_system(traj)
    DataTypes.import_dynamic!(reader, traj, index)
    DataTypes.import_static!(reader, traj, index)

    return (step=get_timestep(traj, index), reader=reader), index+1
end

function Base.iterate(traj::AbstractTrajectory{D, F, SysType}, state::Int64) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    index = state
    if index <= length(traj)
        reader = similar_system(traj)
        rp = latest_reaction(traj, index)
        DataTypes.import_static!(reader, traj, rp)
        DataTypes.import_dynamic!(reader, traj, index)
        return (step=get_timestep(traj, index), reader=reader), index+1
    else
        return nothing
    end
end

function firstindex(traj::AbstractTrajectory)
    return 1
end

function lastindex(traj::AbstractTrajectory)
    return length(traj)
end

#function Base.setproperty!(s::System{D, F, Immutable}, fieldname::Symbol) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#    error("""This type $(typeof(s)) is intended to be read-only. If you want to mutate some data in trajectory, "s = traj[i]" makes a deepcopy. """)
#end

#####
##### Trajectory HDF5 interface
#####

function hmdsave(name::AbstractString, traj::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    file_handler = h5traj(name, "w")
    index = 1
    for reader in traj
        add_snapshot!(file_handler, reader.reader, reader.step; reaction=is_reaction(traj, index), unsafe=true)
        index += 1
    end
    close(file_handler)
end

function read_traj(name::AbstractString, template::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    traj_file = h5traj(name, "r")

    D_file, F_file, SysType_file = get_metadata(traj_file)
    if (D_file, F_file) != (D, F)
        error("Trajectory file type ($D_file, $F_file) is not compatible with the template type ($D, $F).")
    end

    traj = similar(template)
    timesteps = get_timesteps(traj_file)
    reaction_points = get_reactions(traj_file)
    for (step, index) in zip(timesteps, 1:length(traj_file))
        s = snapshot(traj_file, index, similar_system(template))
        add!(traj, s, step; reaction=(step ∈ reaction_points))
    end

    return traj
end

function snapshot(traj_file::AbstractFileFormat, index::Integer, template::AbstractSystem{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    D_file, F_file, SysType_file = get_metadata(traj_file)
    if (D_file, F_file) != (D, F)
        error("Trajectory file $(name) is not compatible with the template $(template).")
    end
    step = get_timesteps(traj_file)[index]

    s = similar(template)
    DataTypes.import_dynamic!(s, traj_file; step=step)
    latest_react = latest_reaction_step(traj_file, step)
    DataTypes.import_static!(s, traj_file, step=latest_react)

    return s
end

function Base.iterate(traj_file::AbstractFileFormat)
    index = 1

    # timestep at reaction (not index!)
    reaction_steps = get_reactions(traj_file)
    timesteps = get_timesteps(traj_file)

    reader = similar_system(traj_file)
    DataTypes.import_static!(reader, traj_file, step=timesteps[index])
    static_cache = deepcopy(reader)
    DataTypes.import_dynamic!(reader, traj_file, step=timesteps[index])

    return (step=timesteps[index], reader=reader), (index+1, reaction_steps, timesteps, static_cache)
end

function Base.iterate(traj_file::AbstractFileFormat, state::Tuple{Int64, Vector{Int64}, Vector{Int64}, S}) where {S<:AbstractSystem}
    index, reaction_steps, timesteps, static_cache = state
    if index > length(traj_file)
        return nothing
    end

    reader = similar(static_cache)
    if timesteps[index] ∈ reaction_steps
        DataTypes.import_static!(reader, traj_file; step=timesteps[index])
        static_cache = deepcopy(reader)
    else
        reader.wrapped = static_cache.wrapped
        reader.element = static_cache.element
        reader.topology = static_cache.topology
        reader.hierarchy = static_cache.hierarchy
    end
    DataTypes.import_dynamic!(reader, traj_file; step=timesteps[index])

    return (step=timesteps[index], reader=reader), (index+1, reaction_steps, timesteps, static_cache)
end

# getindex?
#function slice(traj::AbstractTrajectory, index::Integer)
#
#end
#
#function slice(traj::AbstractTrajectory, time::AbstractFloat)
#
#end
#
#function to_system(traj::AbstractTrajectory)
#
#end
#
#function nearest_slice(traj::AbstractTrajectory, time::AbstractFloat)
#
#end
#
#function Base.push!(traj::AbstractTrajectory{D, F, SysType}, s::AbstractSystem{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#
#end
#
#function Base.append!(addend::AbstractTrajectory{D, F, SysType}, augend::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#
#end
#
#function Base.append!(addend::AbstractTrajectory{D, F, SysType}, augend::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#
#end
