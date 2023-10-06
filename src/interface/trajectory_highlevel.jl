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

function getindex(
    traj::AbstractTrajectory{D, F, SysType},
    rng::AbstractUnitRange{I}
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, I<:Integer}
    slice = empty_trajectory(traj[1])
    for (i, itr) in enumerate(traj)
        i ∉ rng && continue
        snap = deepcopy(itr.reader)
        add!(slice, snap, itr.step; reaction=is_reaction(traj, i))
    end

    return slice
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
    nsnap = length(traj)
    file_handler = h5traj(name, "w")
    index = 1
    for (i, reader) in enumerate(traj)
        print("progress: $(100*i÷nsnap)%    \r")
        add_snapshot!(file_handler, reader.reader, reader.step; reaction=is_reaction(traj, index), unsafe=true)
        index += 1
    end
    println()

    close(file_handler)
end

function read_traj(name::AbstractString, D::Integer, F::Type{<:AbstractFloat}, S::Type{<:AbstractSystemType})
    template = Trajectory{D, F, S, D*D}()
    return read_traj(name, template)
end

function read_traj(name::AbstractString, template::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    traj_file = h5traj(name, "r")
    nsnap = length(traj_file)

    D_file, F_file, SysType_file = get_metadata(traj_file)
    if (D_file, F_file, SysType_file) != (D, F, string(SysType))
        error("Trajectory file type ($D_file, $F_file, $SysType_file) is not compatible with the template type ($D, $F, $SysType).")
    end

    traj = similar(template)
    timesteps = get_timesteps(traj_file)
    reaction_points = get_reactions(traj_file)
    for (step, index) in zip(timesteps, 1:nsnap)
        print("progress: $(100*index÷nsnap)%    \r")
        s = snapshot(
            traj_file,
            similar_system(template);
            step = step,
            unsafe = true
        )
        add!(traj, s, step; reaction=(step ∈ reaction_points))
    end
    println()

    return traj
end

function snapshot(
    traj_file::AbstractFileFormat,
    template::AbstractSystem{D, F, SysType};
    index::Integer = typemin(Int64),
    step::Integer = typemin(Int64),
    unsafe::Bool = false
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    # check metadata
    if !unsafe
        D_file, F_file, SysType_file = get_metadata(traj_file)
        if (D_file, F_file) != (D, F)
            error("Trajectory file $(name) is not compatible with the template $(template).")
        end
    end

    # get timestep from index
    step = if step == typemin(Int64) && index != typemin(Int64)
        get_timesteps(traj_file)[index]
    elseif step != typemin(Int64) && index == typemin(Int64)
        step
    else
        error("Either index or step must be specified. ")
    end

    s = similar(template)
    DataTypes.import_dynamic!(s, traj_file; step=step, unsafe = unsafe)
    if step == latest_reaction_step(traj_file, step)
        DataTypes.import_static!(s, traj_file, step=step, unsafe = unsafe)
    end

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
