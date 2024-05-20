export Trajectory, SubTrajectory
export is_reaction, get_system
export latest_reaction, add!
export length

Base.@kwdef mutable struct Trajectory{D, F<:AbstractFloat, SysType<:AbstractSystemType, L} <: AbstractTrajectory{D, F, SysType}
    # Vector of System. Properties which does not changes at evety timestep are empty except for reactions.
    systems::Vector{System{D, F, SysType, L}} = System{D, F, SysType, L}[]
    # indices corresponding to reactions (not timestep!)
    reactions::Vector{Int64} = Vector{Int64}(undef, 0)
    # index -> timestep
    #timesteps::Vector{Int64} = Vector{Int64}(undef, 0)
end

function Trajectory(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return Trajectory{D, F, SysType, L}([s], [1])
end

function Base.show(
    io::IO, ::MIME"text/plain", traj::Trajectory{D, F, S, L}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, L}
    if isempty(traj.systems)
        println("empty Trajectory{$D, $F, $S, $L}")
        return nothing
    end

    initial = @sprintf "%0.3e" time(traj.systems[1])
    final = @sprintf "%0.3e" time(traj.systems[end])
    "Trajectory{$D, $F, $S, $L}
        length: $(length(traj))
        time: $(initial) to $(final))
    " |> println
end

function empty_trajectory(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return Trajectory{D, F, SysType, L}()
end

function empty_trajectory(
    traj::Trajectory{D, F, SysType, L},
    precision::Type{<:AbstractFloat} = F
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return Trajectory{D, precision, SysType, L}()
end

function is_reaction(traj::Trajectory{D, F, SysType, L}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    #return index ∈ traj.reactions
    return insorted(index, traj.reactions)
end

function get_system(traj::Trajectory{D, F, SysType, L}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return traj.systems[index]
end

#function all_timesteps(traj::Trajectory{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
#    return traj.timesteps
#end

#function get_timestep(traj::Trajectory{D, F, SysType, L}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
#    return all_timesteps(traj)[index]
#end

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
    replica.velocity = deepcopy(current.velocity)
    replica.force = deepcopy(current.force)
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

    #return (step=get_timestep(traj, index), snap=reader), (index+1, reader)
    return reader, (index+1, reader)
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
        #return (step=get_timestep(traj, index), snap=reader), (index+1, reader)
        return reader, (index+1, reader)
    else
        return nothing
    end
end

function add_snapshot!(
    traj::Trajectory{D, F, SysType, L},
    s::System{D, F, SysType, L};
    reaction = false
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if length(traj) > 0
        @assert traj.reactions[end] <= length(traj)
    else
        push!(traj.systems, s)
        push!(traj.reactions, 1)
        #push!(traj.timesteps, timestep)
        return nothing
    end

    #push!(traj.timesteps, timestep)
    if reaction
        if nv(topology(s))==0 && isempty(all_elements(s)) # && isempty(hierarchy(s))
            error("system's topology, hierarchy, and elements are empty. ")
        end
        push!(traj.systems, s)
        push!(traj.reactions, length(traj.systems))
    else
        #replica = System{D, F, SysType}()
        #set_box!(replica, box(s))
        #set_time!(replica, time(s))
        #replica.position = all_positions(s)
        #replica.travel = s.travel
        #replica.velocity = s.velocity
        #replica.force = s.force
        #replica.wrapped = s.wrapped
        #push!(traj.systems, replica)
        empty!(s.topology)
        empty!(s.hierarchy)
        empty!(s.element)
        push!(traj.systems, s)
    end

    return nothing
end

# TODO: implement append!(::SubTrajectory)
function append!(
    traj::Trajectory{D, F, SysType, L},
    addend::Trajectory{D, F, SysType, L};
    include_head = false
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if topology(traj.systems[1]) != topology(addend.systems[1])
        error("topology of traj and addend are not compatible. ")
    end

    nsnap_orig = length(traj.systems)
    endtime = time(traj.systems[end])
    if include_head
        append!(traj.systems, addend.systems)
        append!(traj.reactions, addend.reactions)
        for i in (nsnap_orig+1) : length(traj.systems)
            traj.systems[i].time += endtime
        end
    else
        append!(traj.systems, @view addend.systems[2:end])
        append!(traj.reactions, addend.reactions[2:end])
        for i in (nsnap_orig+2) : length(traj.systems)
            traj.systems[i].time += endtime
        end
    end

    return nothing
end

function import_dynamic!(reader::System{D, F, S1}, traj::Trajectory{D, F, S2}, index::Integer) where {D, F<:AbstractFloat, S1<:AbstractSystemType, S2<:AbstractSystemType}
    s = get_system(traj, index)
    set_time!(reader, time(s))
    set_box!(reader, box(s))
    reader.position = all_positions(s)
    reader.velocity = s.velocity
    reader.force = s.force
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

    # find `i` s.t. traj.reactions[i] <= index < traj.reactions[i+1]
    ii = searchsortedlast(traj.reactions, index)
    return traj.reactions[ii]
end

function is_reaction(s::System)
    return all_elements(s) |> isempty
end

function similar_system(
    traj::Trajectory{D, F, SysType, L},
    precision::Type{<:AbstractFloat} = F;
    reserve_dynamic = false,
    reserve_static = false
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    s = System{D, precision, SysType}()
    return if !reserve_dynamic && !reserve_static
        s
    else
        reserve_dynamic && import_dynamic!(s, traj, 1)
        reserve_static && import_static!(s, traj, 1)
        deepcopy(s)
    end
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

    rp = latest_reaction(st.traj, real_idx)
    import_static!(reader, st.traj, rp)

    #return (step=get_timestep(st.traj, real_idx), snap=reader), (pseude_idx+1, reader)
    return reader, (pseude_idx+1, reader)
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
        #return (step=get_timestep(st.traj, real_idx), snap=reader), (pseude_idx+1, reader)
        return reader, (pseude_idx+1, reader)
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

#function all_timesteps(st::SubTrajectory{D, F, S, L, R}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
#    return all_timesteps(st.traj)[st.traj_range]
#end
#
#function get_timestep(st::SubTrajectory{D, F, S, L, R}, index::Integer) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
#    return get_timestep(st.traj, st.traj_range[index])
#end

function Base.length(st::SubTrajectory{D, F, S, L, R}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    return length(st.traj_range)
end

function wrapped(st::SubTrajectory{D, F, S, L, R}) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    return wrapped(st.traj)
end

function latest_reaction(st::SubTrajectory{D, F, S, L, R}, index::Integer) where {D, F<:AbstractFloat, S<:AbstractSystemType, L, R<:OrdinalRange}
    if index ∉ 1:length(st)
        error("index $index ∉ $(st.traj_range)")
    end

    # find `i` s.t. traj.reactions[i] <= index < traj.reactions[i+1]
    ii = searchsortedlast(st.traj.reactions, st.traj_range[index])
    return st.traj.reactions[ii]
end
