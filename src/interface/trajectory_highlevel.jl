function firstindex(traj::AbstractTrajectory)
    return 1
end

function lastindex(traj::AbstractTrajectory)
    return length(traj)
end

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
        error("saving trajectory with Float16 is not supported. ")
    end

    change_wrap = false
    if wrapped(traj)
        @warn "Trajectory is wrapped. Saving unwrapped format..."
        change_wrap = true
        unwrap!(traj)
    end

    file_handler = h5traj(name, "w", traj[1], length(traj), precision)
    file = get_file(file_handler)
    nsnap = length(traj)
    try
        for (i, reader) in enumerate(traj)
            print("progress: $(100*i÷nsnap)%    \r")
            add_snapshot!(
                file_handler,
                reader,
                i;
                reaction = is_reaction(traj, i),
                unsafe = true
            )
        end
        file["/times"][1:length(traj),1] = [time(r) for r in traj]
        println()
    finally
        close(file_handler)
        change_wrap && wrap!(traj)
    end

    return nothing
end
