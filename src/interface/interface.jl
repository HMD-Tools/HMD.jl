"""
    NotImplementedError{M}(m)
`Exception` thrown when a method from the `AbstractSystem` or `AbstractTrajectory` interface
is not implemented by a given AbstractSystem or Trajectory type.
"""
struct NotImplementedError{M} <: Exception
    m::M
    NotImplementedError(m::M) where {M} = new{M}(m)
end

function Base.showerror(io::IO, ie::NotImplementedError)
    return print(io, "method $(ie.m) not implemented.")
end

_NI(m) = throw(NotImplementedError(m))

# AbstractSystemのトポロジーは必ずAbstractGraphのsubtype
abstract type AbstractSystemType end
abstract type AbstractSystem{D, F<:AbstractFloat, SysType<:AbstractSystemType} end
abstract type AbstractTrajectory{D, F<:AbstractFloat, SysType<:AbstractSystemType} end
abstract type AbstractBbox{D, F<:AbstractFloat, L} end

# system core interface
dimension(s::AbstractSystem) = _NI("dimension")
precision(s::AbstractSystem) = _NI("precision")
system_type(s::AbstractSystem) = _NI("system_type")
similar(s::AbstractSystem; reserve_dynamic::Bool=false, reserve_static::Bool=false) = _NI("similar")
show(io::IO, ::MIME"text/plain", s::AbstractSystem{D, F, SysType}) where {D, F<:AbstractFloat, SysType <: AbstractSystemType} = _NI("show")
natom(s::AbstractSystem) = _NI("natom")
nbond(s::AbstractSystem) = _NI("nbond")
time(s::AbstractSystem) = _NI("time")
set_time!(s::AbstractSystem, time::AbstractFloat) = _NI("set_time!")
topology(s::AbstractSystem) = _NI("topology")
box(s::AbstractSystem) = _NI("box")
set_box!(s::AbstractSystem, box::AbstractBbox) = _NI("set_box!")
element_type(s::AbstractSystem) = _NI("element_type")
all_elements(s::AbstractSystem) = _NI("all_elements")
element(s::AbstractSystem, i::Integer) = _NI("element")
element(s::AbstractSystem, label::HLabel) = _NI("element")
add_element!(s::AbstractSystem, atomic_number::Integer) = NI("add_element!")
add_elements!(s::AbstractSystem, atomic_numbers::AbstractVector{<:Integer}) = NI("add_elements!")
set_element!(s::AbstractSystem, atom_id::Integer, atomic_number::Integer) = _NI("set_element!")
set_elements!(s::AbstractSystem, atom_ids::AbstractVector{<:Integer}, atomic_numbers::AbstractVector{<:Integer}) = _NI("set_elements!")
all_positions(s::AbstractSystem) = _NI("all_positions")
position(s::AbstractSystem, atom_id::Integer) = _NI("position")
add_position!(s::AbstractSystem, x::AbstractVector{<:AbstractFloat}) = NI("add_position!")
add_positions!(s::AbstractSystem, x::AbstractVector{AbstractVector{<:AbstractFloat}}) = NI("add_positions!")
set_position!(s::AbstractSystem, atom_id::Integer, x::AbstractVector{<:AbstractFloat}) = _NI("set_position!")
get_velocity(s::AbstractSystem, atom_id::Integer) = _NI("get_velocity")
all_velocities(s::AbstractSystem) = _NI("all_velocities")
add_velocity!(s::AbstractSystem, v::AbstractVector{<:Real}) = _NI("add_velocity!")
add_velocity!(s::AbstractSystem, v::AbstractVector{<:Unitful.Velocity}) = _NI("add_velocity!")
add_velocities!(s::AbstractSystem, v::AbstractVector{<:AbstractVector{<:Real}}) = _NI("add_velocities!")
add_velocities!(s::AbstractSystem, v::AbstractVector{<:AbstractVector{<:Unitful.Velocity}}) = _NI("add_velocities!")
set_velocity!(s::AbstractSystem, atom_id::Integer, v::AbstractVector{<:Real}) = _NI("set_velocity!")
set_velocity!(s::AbstractSystem, atom_id::Integer, v::AbstractVector{<:Unitful.Velocity}) = _NI("set_velocity!")
get_force(s::AbstractSystem, atom_id::Integer) = _NI("get_force")
all_forces(s::AbstractSystem) = _NI("all_forces")
add_force!(s::AbstractSystem, f::AbstractVector{<:Real}) = _NI("add_force!")
add_force!(s::AbstractSystem, f::AbstractVector{<:Unitful.Force}) = _NI("add_force!")
add_forces!(s::AbstractSystem, f::AbstractVector{<:AbstractVector{<:Real}}) = _NI("add_forces!")
add_forces!(s::AbstractSystem, f::AbstractVector{<:AbstractVector{<:Unitful.Force}}) = _NI("add_forces!")
set_force!(s::AbstractSystem, atom_id::Integer, f::AbstractVector{<:Real}) = _NI("set_force!")
set_force!(s::AbstractSystem, atom_id::Integer, f::AbstractVector{<:Unitful.Force}) = _NI("set_force!")
all_travels(s::AbstractSystem) = _NI("all_travels")
travel(s::AbstractSystem, atom_id::Integer) = _NI("travel")
set_travel!(s::AbstractSystem, atom_id::Integer, n::AbstractVector{<:Integer}) = _NI("set_travel!")
wrapped(s::AbstractSystem) = _NI("wrapped")
wrap!(s::AbstractSystem) = _NI("wrap!")
unwrap!(s::AbstractSystem) = _NI("unwrap!")
label2atom(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("label2atom")
merge!(augend::AbstractSystem, addend::AbstractSystem;
    augend_parents::Dict{String, HLabel}, addend_parents::Dict{String, HLabel},
    unsafe::Bool=false) = NI("merge!")

# system label manipulation
hierarchy_names(s::AbstractSystem) = _NI("hierarchy_names")
hierarchy(s::AbstractSystem, name::AbstractString) = _NI("hierarchy")
add_hierarchy!(s::AbstractSystem, name::AbstractString) = _NI("add_hierarchy!")
remove_hierarchy!(s::AbstractSystem, name::AbstractString) = _NI("remove_hierarchy!")
all_labels(s::AbstractSystem, hname::AbstractString) = _NI("all_labels")
all_labels(s::AbstractSystem, hname::AbstractString, label_type::AbstractString) = _NI("all_labels")
add_label!(s::AbstractSystem, hname::AbstractString, label::HLabel; unsafe::Bool=false) = _NI("add_label!")
#add_label!(s::AbstractSystem, hname::AbstractString, label_type::AbstractString) = _NI("add_label!")
add_labels!(s::AbstractSystem, hname::AbstractString, labels::AbstractVector{HLabel}) = _NI("add_labels!")
count_label(s::AbstractSystem, hname::AbstractString, label_type::String) = _NI("count_label")
add_relation!(s::AbstractSystem, hname::AbstractString; super::HLabel, sub::HLabel) = _NI("add_relation!")
add_relations!(s::AbstractSystem, hname::AbstractString; super::AbstractVector{HLabel}, sub::AbstractVector{HLabel}) = _NI("add_relations!")
insert_relation!(s::AbstractSystem, hname::AbstractString, label::HLabel; super::HLabel, sub::HLabel) = _NI("insert_relation!")
insert_relations!(s::AbstractSystem, hname::AbstractString, index::Integer; super::AbstractVector{HLabel}, sub::AbstractVector{HLabel}) = _NI("insert_relations!")
remove_label!(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("remove_label!")
remove_relation!(s::AbstractSystem, hname::AbstractString, super::HLabel, sub::HLabel) = _NI("remove_relation!")
label_unique(s::AbstractSystem) = _NI("label_unique")
label_nocycle(s::AbstractSystem) = _NI("label_nocycle")
label_connected(s::AbstractSystem) = _NI("label_connected")
contains(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("contains")
issuper(s::AbstractSystem, hname::AbstractString, label1::HLabel, label2::HLabel) = _NI("issuper")
issub(s::AbstractSystem, hname::AbstractString, label1::HLabel, label2::HLabel) = _NI("issub")
super(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("super")
sub(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("sub")

# system property interface
#prop_names(s::AbstractSystem) = NI("prop_names")
#prop(s::AbstractSystem, pname::AbstractString) = NI("prop")
#add_prop!(s::AbstractSystem, pname::AbstractString, p::AbstractArray) = NI("add_prop!")
#set_prop!(s::AbstractSystem, pname::AbstractString, p::AbstractArray) = NI("set_prop!")
#set_prop!(s::AbstractSystem, pname::AbstractString, p::Union{Float64, Float32}) = NI("set_prop!")

# system_io interface
abstract type AbstractFileFormat end
#close(file_handler::AbstractFileFormat) = _NI("close")

# trajectory interface
empty_trajectory(s::AbstractSystem) = _NI("empty_trajectory")
is_reaction(s::AbstractTrajectory, index::Integer) = _NI("is_reaction")
get_system(s::AbstractTrajectory, index::Integer) = _NI("get_system")
all_timesteps(traj::AbstractTrajectory) = _NI("all_timesteps")
get_timestep(traj::AbstractTrajectory, index::Integer) = _NI("get_timestep")
length(traj::AbstractTrajectory) = _NI("length")
add!(traj::AbstractTrajectory, s::AbstractSystem, timestep::Integer; reaction=false) = _NI("add!")
add!(traj::AbstractTrajectory, s::AbstractTrajectory) = _NI("add!")
import_dynamic!(reader::AbstractSystem, traj::AbstractTrajectory, index::Integer) = _NI("import_dynamic!")
import_static!(reader::AbstractSystem, traj::AbstractTrajectory, index::Integer) = _NI("import_static!")
latest_reaction(traj::AbstractTrajectory, index::Integer) = _NI("latest_reaction")
similar(traj::AbstractTrajectory) = _NI("similar")
similar_system(traj::AbstractTrajectory; reserve_dynamic::Bool=false, reserve_static::Bool=false) = _NI("similar_system")
dimension(traj::AbstractTrajectory) = _NI("dimension")
precision(traj::AbstractTrajectory) = _NI("precision")
system_type(traj::AbstractTrajectory) = _NI("system_type")
wrapped(traj::AbstractTrajectory) = _NI("wrapped")
wrap!(traj::AbstractTrajectory) = _NI("wrap!")
unwrap!(traj::AbstractTrajectory) = _NI("unwrap!")
prop(s::AbstractSystem, index::Integer, pname::AbstractString) = NI("prop")

# trajectory io interface
add_snapshot!(
    file_handler::AbstractFileFormat,
    s::AbstractSystem{D, F},
    index::Integer,
    reaction::Bool = false,
    unsafe::Bool = false
) where {D, F<:AbstractFloat} = _NI("add_snapshot!")
get_metadata(traj_file::AbstractFileFormat) = _NI("get_metadata")
