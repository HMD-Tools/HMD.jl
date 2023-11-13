"""
highlevelの引数をすべて具象型にすべきか？
    ・具象型で書き直すと新しい具象型を定義したときに再実装量が増える
    ・抽象型のままにしておく場合はinterfaceでNIを定義する意味が無い
iteratorをlowlevelに落とすことを検討
"""



#function Base.iterate(traj::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#    index = 1
#    reader = similar_system(traj)
#    DataTypes.import_dynamic!(reader, traj, index)
#    DataTypes.import_static!(reader, traj, index)
#
#    return (step=get_timestep(traj, index), snap=reader), (index+1, reader)
#end

#function Base.iterate(
#    traj::AbstractTrajectory{D, F, SysType},
#    state::Tuple{Int64, S}
#    #state::Int64
#) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, S<:AbstractSystem{D, F, SysType}}
#    index = state[1]
#    if index <= length(traj)
#        #reader = similar_system(traj)
#        reader = state[2]
#        rp = latest_reaction(traj, index)
#        DataTypes.import_static!(reader, traj, rp)
#        DataTypes.import_dynamic!(reader, traj, index)
#        return (step=get_timestep(traj, index), snap=reader), (index+1, reader)
#    else
#        return nothing
#    end
#end

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

function hmdsave(
    name::AbstractString,
    traj::AbstractTrajectory{D, F, SysType};
    precision = F
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if precision < F && precision != Float16
        @info "warning: saving precision is lower than the system precision. \n" *
            "This cause information loss."
    elseif precision == Float16
        @warn "For molecular dynamics, Float16 is not appropriate in most cases. \n"
    end

    change_wrap = false
    if wrapped(traj)
        @warn "Trajectory is wrapped. Saving unwrapped format..."
        change_wrap = true
        unwrap!(traj)
    end

    nsnap = length(traj)
    file_handler = h5traj(name, "w")
    index = 1
    for (i, reader) in enumerate(traj)
        print("progress: $(100*i÷nsnap)%    \r")
        add_snapshot!(
            file_handler,
            reader.snap,
            reader.step,
            precision;
            reaction = is_reaction(traj, index),
            unsafe = true,
        )
        index += 1
    end
    println()

    if change_wrap
        wrap!(traj)
    end

    close(file_handler)
end

function hmdsave_new(
    name::AbstractString,
    traj::AbstractTrajectory{D, F, SysType};
    precision = F
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if precision < F && precision != Float16
        @info "warning: saving precision is lower than the system precision. \n" *
            "This cause information loss."
    elseif precision == Float16
        @warn "For molecular dynamics, Float16 is not appropriate in most cases. \n"
    end

    change_wrap = false
    if wrapped(traj)
        @warn "Trajectory is wrapped. Saving unwrapped format..."
        change_wrap = true
        unwrap!(traj)
    end

    file_handler = h5traj_new(name, "w", traj[1], length(traj), precision)
    file = DataTypes.get_file(file_handler)
    nsnap = length(traj)
    try
        for (i, reader) in enumerate(traj)
            print("progress: $(100*i÷nsnap)%    \r")
            DataTypes.add_snapshot_new!(
                file_handler,
                reader.snap,
                i;
                reaction = is_reaction(traj, i),
                unsafe = true,
            )
        end
        file["/timesteps"][1:length(traj),1] = [get_timestep(traj, i) for i in 1:nsnap]
        file["/times"][1:length(traj),1] = [time(r.snap) for r in traj]
        println()
    finally
        close(file_handler)
        change_wrap && wrap!(traj)
    end

    return nothing
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

    return (step=timesteps[index], snap=reader), (index+1, reaction_steps, timesteps, static_cache)
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

    return (step=timesteps[index], snap=reader), (index+1, reaction_steps, timesteps, static_cache)
end
