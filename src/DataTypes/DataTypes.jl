# presetとしてpolymer hierarchyを追加: polymer >: {monomer >: connection, endcap >: endatom}
# add_label, add, \oplus, add!を作成 テスト可能
# selectionなど読み込み機能を追加
#
#

# ここにはデータ型の内部を直接触る最小限の関数を書く
# HMD/interface/system.jlにAbstractSystem関数を作成することでDataTypesの内部構造を容易に変更できる

module DataTypes

using DataStructures
using Graphs
using HDF5
using LinearAlgebra
using MLStyle
using PeriodicTable
using Printf
using Reexport
using SimpleWeightedGraphs
using StaticArrays
using Unitful

using ..HierarchyLabels
using ..TopologyGraphs

@reexport import Base: *, +, -, /, <, <=, ==, >, >=, close, contains, convert,
    getindex, firstindex, lastindex, iterate, eachindex,
    length, position, precision, promote_rule, promote_type, setproperty!, show, similar,
    string, time, ∈, ∉, merge!, println, append!, replace!

@reexport import ..HMD: deserialize, serialize
@reexport import ..HMD:
    # system core interface
    AbstractBbox,
    AbstractSystem,
    AbstractSystemType,
    AbstractTrajectory,
    dimension,
    precision,
    system_type,
    similar,
    show,
    natom,
    eachatom,
    nbond,
    time,
    set_time!,
    topology,
    box,
    set_box!,
    element_type,
    all_elements,
    element,
    add_element!,
    add_elements!,
    set_element!,
    set_elements!,
    all_positions,
    position,
    add_position!,
    add_positions!,
    set_position!,
    get_velocity,
    all_velocities,
    add_velocity!,
    add_velocities!,
    set_velocity!,
    get_force,
    all_forces,
    add_force!,
    add_forces!,
    set_force!,
    all_travels,
    travel,
    set_travel!,
    wrapped,
    wrap!,
    unwrap!,
    label2atom,
    merge!,

    # system label manipulation
    hierarchy_names,
    hierarchy,
    add_hierarchy!,
    remove_hierarchy!,
    all_labels,
    add_label!,
    add_labels!,
    count_label,
    add_relation!,
    add_relations!,
    insert_relation!,
    remove_label!,
    remove_relation!,
    label_unique,
    label_nocycle,
    label_connected,
    contains,
    issuper,
    issub,
    super,
    sub,
    super_labels,
    sub_labels,

    # system io interface
    AbstractFileFormat,
    close,

    # trajectory interface
    empty_trajectory,
    get_system,
    all_timesteps,
    get_timestep,
    length,
    import_dynamic!,
    import_static!,
    latest_reaction,
    similar,
    similar_system,
    dimension,
    precision,
    system_type,
    wrapped,
    wrap!,
    unwrap!,
    is_reaction,
    length,
    getindex,
    lastindex,
    firstindex,
    wrapped,

    # trajectory io interface
    add_snapshot!,
    get_metadata,
    h5traj_reader


# core subtype signature
export Position, BoundingBox, HLabel, LabelHierarchy

# core immut signature
export GeneralSystem, System, print_to_string


# fileIO
export H5traj, SerializedTopology, PackedHierarchy
export h5traj, H5trajReader, SubH5trajReader, get_file, read_traj
export read_snapshot!, partial_read

#constants
export Entire_System, BO_Precision, atom_mass, Atom_Label, Atomic_Number_Precision

const Entire_System = HLabel("entire_system", 1)
const Atom_Label = ""
const BO_Precision = Int8

const Atomic_Number_Precision = Int8
const atom_mass = Dict(i => ustrip(elements[i].atomic_mass) for i in 1:length(elements))

include("position.jl")
include("boundingbox.jl")

#####
##### Type `System` definition
#####

struct GeneralSystem <: AbstractSystemType end

mutable struct System{D, F<:AbstractFloat, S<:AbstractSystemType, L} <: AbstractSystem{D, F, S}
    time::F
    #topology::SimpleWeightedGraph{Int64, Rational{BO_Precision}}
    topology::TopologyGraph{Int64, Rational{BO_Precision}}
    box::BoundingBox{D, F, L}

    # atom property
    position::Position{D, F}
    travel::Vector{SVector{D, Int16}}
    wrapped::Bool
    element::Vector{Atomic_Number_Precision}

    hierarchy::Dict{String, LabelHierarchy}

    # optional properties
    velocity::Vector{SVector{D, F}}
    force::Vector{SVector{D, F}}
end

function System{D, F, SysType}() where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    System{D, F, SysType, D*D}(
        zero(F),
        #SimpleWeightedGraph{Int64, Rational{BO_Precision}}(),
        TopologyGraph{Int64, Rational{BO_Precision}}(),
        BoundingBox{D, F}(),
        Position{D, F}(),
        Vector{SVector{D, Int16}}(undef, 0),
        false,
        Atomic_Number_Precision[],
        Dict{String, LabelHierarchy}(),
        Vector{SVector{D, F}}(undef, 0),
        Vector{SVector{D, F}}(undef, 0)
    )
end

function System{D, F}() where {D, F<:AbstractFloat}
    return System{D, F, GeneralSystem}()
end

function dimension(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return D
end

function precision(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return F
end

function system_type(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return SysType
end

function _mycopy(v::Vector{T}) where {T}
    x = Vector{T}(undef, length(v))
    @inbounds x .= v
    return x
end

function Base.similar(
    s::System{D, F, SysType, L},
    precision::Type{<:AbstractFloat} = F;
    reserve_dynamic::Bool = false,
    reserve_static::Bool = false
) where {D, F<:AbstractFloat, SysType, L}
    sim = System{D, precision, SysType}()
    if reserve_dynamic
        set_time!(sim, time(s))
        set_box!(sim, box(s))
        sim.position = deepcopy(s.position)
        sim.travel = deepcopy(s.travel)
        sim.wrapped = s.wrapped
        sim.velocity = deepcopy(s.velocity)
        sim.force = deepcopy(s.force)
    end
    if reserve_static
        sim.topology = deepcopy(s.topology)
        sim.hierarchy = deepcopy(s.hierarchy)
        sim.element = deepcopy(s.element)
    end
    return sim
end

function Base.show(io::IO, ::MIME"text/plain", s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    "System{$D, $F, $SysType}
        time: $(time(s))
        bbox: $(box(s))
        natoms: $(natom(s))
        hierarchy: $(hierarchy_names(s))
    " |> println
end

function natom(s::System)
    natm = length(s.position)
    return natm
end

eachatom(s::System) = 1:natom(s)

function nbond(s::System)
    return topology(s) |> ne
end

function time(s::System)
    s.time
end

function set_time!(s::System, time::AbstractFloat)
    s.time = time
end

function set_time!(s::System, time::Unitful.Time)
    s.time = uconvert(u"ns", time) |> ustrip
end

function topology(s::System)
    s.topology
end

function box(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    s.box
end

function set_box!(
    s::System{D, F, SysType, L},
    box::BoundingBox{D, T, L}
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L, T<:AbstractFloat}
    s.box = BoundingBox{D, F, L}(box.origin, box.axis)
end

function element_type(s::System)
    return Atomic_Number_Precision
end

function all_elements(s::System)
    s.element
end

function element(s::System, atom_id::Integer)
    s.element[atom_id]
end

function add_element!(s::System, atomic_number::Integer)
    push!(s.element, atomic_number)
end

function add_elements!(s::System, atomic_numbers::AbstractVector{<:Integer})
    append!(s.element, atomic_numbers)
end

function set_element!(s::System, atom_id::Integer, atomic_number::Integer)
    s.element[atom_id] = atomic_number
end

function set_elements!(s::System, atom_ids::AbstractVector{<:Integer}, atomic_numbers::AbstractVector{<:Integer})
    if length(atom_ids) != length(atomic_numbers)
        throw(DimensionMismatch("Dimension of atom_ids is $(length(atom_ids)) but atomic_numbers dimension is $(length(atomic_numbers))"))
    end
    s.element[atom_ids] .= atomic_numbers
end

function all_positions(s::System)
    s.position
end

function position(s::System, atom_id::Integer)
    s.position[atom_id]
end

function add_position!(s::System, x::AbstractVector{<:AbstractFloat})
    push!(s.position, x)
    if wrapped(s)
        error("atom coordinate is wrapped. Call unwrap!(s) before adding atom.")
    end
    push!(s.travel, SVector{3, Int16}(0, 0, 0))

    return nothing
end

function add_positions!(s::System, x::AbstractVector{<:AbstractVector{<:AbstractFloat}})
    append!(all_positions(s), x)
    if wrapped(s)
        error("atom coordinate is wrapped. Call unwrap!(s) before adding atom.")
    end
    append!(s.travel, [SVector{3, Int16}(0, 0, 0) for _ in 1:length(x)])

    return nothing
end

function set_position!(s::System, atom_id::Integer, x::AbstractVector{<:AbstractFloat})
    s.position[atom_id] = x
end

function set_position!(s::System, atom_id::Integer, x::AbstractVector{T}) where {T<:Unitful.Length}
    s.position[atom_id][1] = uconvert(u"Å", x[1]) |> ustrip
    s.position[atom_id][2] = uconvert(u"Å", x[2]) |> ustrip
    s.position[atom_id][3] = uconvert(u"Å", x[3]) |> ustrip
end

function get_velocity(s::System, atom_id::Integer)
    return s.velocity[atom_id]
end

function all_velocities(s::System)
    return s.velocity
end

function add_velocity!(s::System, x::AbstractVector{T}) where {T<:Real}
    push!(s.velocity, x)
    return nothing
end

function add_velocity!(
    s::System{D, F, S},
    x::AbstractVector{T}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, T<:Unitful.Velocity}
    push!(
        s.velocity,
        SVector{D, F}(uconvert(u"Å/ns", e).val for e in x)
    )
    return nothing
end

function add_velocities!(s::System, x::AbstractVector{AbstractVector{T}}) where {T<:Real}
    append!(s.velocity, x)
    return nothing
end

function add_velocities!(
    s::System{D, F, S},
    vecs::AbstractVector{AbstractVector{T}}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, T<:Unitful.Velocity}
    append!(
        s.velocity,
        [SVector{D, F}(uconvert(u"Å/ns", e).val for e in x) for x in vecs]
    )
    return nothing
end

function set_velocity!(s::System, atom_id::Integer, x::AbstractVector{T}) where {T<:Real}
    s.velocity[atom_id] = x
    return nothing
end

function set_velocity!(s::System, atom_id::Integer, x::AbstractVector{T}) where {T<:Unitful.Velocity}
    s.velocity[atom_id][1] = uconvert(u"Å/ns", x[1]) |> ustrip
    s.velocity[atom_id][2] = uconvert(u"Å/ns", x[2]) |> ustrip
    s.velocity[atom_id][3] = uconvert(u"Å/ns", x[3]) |> ustrip
    return nothing
end

function get_force(s::System, atom_id::Integer)
    return s.force[atom_id]
end

function all_forces(s::System)
    return s.force
end

function add_force!(s::System, x::AbstractVector{T}) where {T<:Real}
    push!(s.force, x)
    return nothing
end

function add_force!(
    s::System{D, F, S},
    x::AbstractVector{T}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, T<:Unitful.Force}
    push!(
        s.focrce[atom_id],
        SVector{D, F}(uconvert(u"eV/Å", e).val for e in x)
    )
    return nothing
end

function add_forces!(s::System, x::AbstractVector{AbstractVector{T}}) where {T<:Real}
    append!(s.force, x)
    return nothing
end

function add_forces!(
    s::System{D, F, S},
    vecs::AbstractVector{AbstractVector{T}}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, T<:Unitful.Force}
    append!(
        s.force,
        [SVector{D, F}(uconvert(u"eV/Å", e).val for e in x) for x in vecs]
    )
    return nothing
end

function set_force!(s::System, atom_id::Integer, x::AbstractVector{T}) where {T<:Real}
    s.force[atom_id] = x
    return nothing
end

function set_force!(s::System, atom_id::Integer, x::AbstractVector{T}) where {T<:Unitful.Force}
    s.force[atom_id][1] = uconvert(u"eV/Å", x[1]) |> ustrip
    s.force[atom_id][2] = uconvert(u"eV/Å", x[2]) |> ustrip
    s.force[atom_id][3] = uconvert(u"eV/Å", x[3]) |> ustrip
    return nothing
end

function all_travels(s::System)
    return s.travel
end

function travel(s::System, atom_id::Integer)
    return s.travel[atom_id]
end

function set_travel!(s::System, atom_id::Integer, n::AbstractVector{<:Integer})
    s.travel[atom_id] = n
end

function wrapped(s::System)
    s.wrapped
end

function _change_wrap!(s::System)
    s.wrapped = !(s.wrapped)
end

function _floor_rem(c::SVector{D, F}) where {D, F<:AbstractFloat}
    ipart = SVector{D, Int16}(floor(Int16, x) for x in c)
    fpart = c - ipart
    return ipart, fpart
end

function wrap!(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if wrapped(s)
        return nothing
    end

    axis = box(s).axis
    for id in 1:natom(s)
        # convert coordinates from Cartesian to box vectors i.e.
        # x = c[1] .* axis[:,1] .+ c[2] .* axis[:,2] .+ ...
        c = _box_coord(position(s, id), box(s))
        iparts, fparts = _floor_rem(c)
        pos = sum(fparts[dim] * axis[:,dim] for dim in 1:D)
        set_travel!(s, id, iparts)
        set_position!(s, id, pos)
    end
    _change_wrap!(s)

    return nothing
end

function unwrap!(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if !wrapped(s)
        return nothing
    end

    axis   = box(s).axis
    origin = box(s).origin
    for i in 1:natom(s)
        x = position(s, i) - origin
        n = travel(s, i)
        pos = x + sum(n[dim] * axis[:, dim] for dim in 1:D)
        set_position!(s, i, pos + origin)
        set_travel!(s, i, SVector{3, Int16}(0, 0, 0))
    end
    _change_wrap!(s)

    return nothing
end

"""
convert coordinates from Cartesian to box vectors.

Math equation is calculated by the below code with Symbolics.jl

# code
```julia-repl
    using Symbolics
    @variables p₁ p₂ p₃ # position
    @variables a₁ a₂ a₃ b₁ b₂ b₃ c₁ c₂ c₃ # box vectors
    @variables o₁ o₂ o₃  # box origin
    @variables α β γ # coefficients
    eq₁ = p₁ ~ α*a₁ + β*b₁ + γ*c₁ + o₁
    eq₂ = p₂ ~ α*a₂ + β*b₂ + γ*c₂ + o₂
    eq₃ = p₃ ~ α*a₃ + β*b₃ + γ*c₃ + o₃
    equations = Symbolics.solve_for([eq₁ eq₂ eq₃], [α, β, γ]; simplify=true)
```
"""
@inline function _box_coord(p::SVector{3, F}, box::BoundingBox{3, F, 9}) where {F<:AbstractFloat}
    a, b, c = box.axis[:,1], box.axis[:,2], box.axis[:,3]
    o = box.origin
    inv = 1 / (a[2]*b[3]*c[1] + a[3]*b[1]*c[2] + a[1]*b[2]*c[3] - a[1]*b[3]*c[2] - a[2]*b[1]*c[3] - a[3]*b[2]*c[1])

    α = inv * (b[2]*c[1]*o[3] + b[3]*c[2]*o[1] + b[3]*c[1]*p[2] + b[1]*c[2]*p[3] + b[1]*c[3]*o[2] + b[2]*c[3]*p[1] - b[2]*c[1]*p[3] - b[2]*c[3]*o[1] - b[3]*c[1]*o[2] - b[1]*c[2]*o[3] - b[1]*c[3]*p[2] - b[3]*c[2]*p[1])
    β = inv * (a[3]*c[1]*o[2] + a[1]*c[2]*o[3] + a[1]*c[3]*p[2] + a[3]*c[2]*p[1] + a[2]*c[1]*p[3] + a[2]*c[3]*o[1] - a[2]*c[1]*o[3] - a[3]*c[1]*p[2] - a[1]*c[2]*p[3] - a[1]*c[3]*o[2] - a[2]*c[3]*p[1] - a[3]*c[2]*o[1])
    γ = inv * (a[2]*b[1]*o[3] + a[3]*b[1]*p[2] + a[1]*b[3]*o[2] + a[2]*b[3]*p[1] + a[1]*b[2]*p[3] + a[3]*b[2]*o[1] - a[1]*b[2]*o[3] - a[3]*b[1]*o[2] - a[3]*b[2]*p[1] - a[2]*b[3]*o[1] - a[1]*b[3]*p[2] - a[2]*b[1]*p[3])

    return SVector{3, F}(α, β, γ)
end

function isatom(label::HLabel)
    return type(label) == ""
end

include("label_manipulation.jl")
#include("property.jl")
#include("system_io.jl")
include("trajectory.jl")
include("trajectory_HDF5.jl")
include("test.jl")

function merge!(
    augend::System{D, F, S1}, addend::System{D, F, S2};
    augend_parents::Dict{String, HLabel}, addend_parents::Dict{String, HLabel},
    unsafe::Bool = false
) where {D, F<:AbstractFloat, S1<:AbstractSystemType, S2<:AbstractSystemType}
    # unwrap if sysytem is wrapped
    change_augend = change_addend = false
    if wrapped(augend)
        unwrap(augend)
        change_augend = true
    end
    if wrapped(addend)
        unwrap(addend)
        change_addend = true
    end

    # check whther target hieararchies to merge exists in each systems
    if !(keys(augend_parents) ⊆ hierarchy_names(augend))
        error(
            "augend_parents must be a subset of augend hierarchy\n" *
            "   augend hierarchies: $(keys(augend_parents))\n" *
            "   targets: $(hierarchy_names(augend))"
        )
    elseif !(keys(addend_parents) ⊆ hierarchy_names(addend))
        error(
            "addend_parents must be a subset of addend hierarchy\n" *
            "   addend hierarchies: $(keys(addend_parents))\n" *
            "   targets: $(hierarchy_names(addend))"
        )
    end

    # check whether addend and aungend parents have the same keys
    if Set(keys(addend_parents)) != Set(keys(augend_parents))
        error(
            "augend_parents and addend_parents must have the same targets\n" *
            "   augend_parents: $(keys(augend_parents))\n" *
            "   addend_parents: $(keys(addend_parents))"
        )
    end

    #
    if !(hierarchy_names(addend) ⊆ hierarchy_names(augend))
        missings = filter(
            hname -> !(hname ∈ hierarchy_names(augend)),
            hierarchy_names(addend)
        )
        error(
            "addend hierarchy must be a subset of augend hierarchy\n" *
            "   augend hierarchies: $(hierarchy_names(augend))\n" *
            "   addend hierarchies: $(hierarchy_names(addend))" *
            "   missings: $(missings)"
        )
    end

    # add new atoms to augend
    atom_mapping = natom(augend)+1 : natom(augend)+natom(addend) # atom_mapping[addend_atom] == augend_atom
    append!(augend.position, all_positions(addend))
    append!(augend.travel, all_travels(addend))
    append!(augend.element, all_elements(addend))

    _merge_topology!(topology(augend), topology(addend), atom_mapping)

    for hname in keys(augend_parents)
        lh_augend = hierarchy(augend, hname)
        lh_addend = hierarchy(addend, hname)
        _merge_hierarchy!(
            lh_augend, lh_addend;
            augend_parent = augend_parents[hname],
            addend_parent = addend_parents[hname],
            atom_mapping,
            unsafe = unsafe
        )
    end

    # restore wrapping state
    change_augend && wrap!(augend)
    change_addend && wrap!(addend)

    return nothing
end

function _merge_topology!(
    augend::TopologyGraph{I, Rational{BO_Precision}},
    addend::TopologyGraph{I, Rational{BO_Precision}},
    id_mapping::AbstractVector{T}
) where {I<:Integer, T<:Integer}
    @assert id_mapping == nv(augend)+1 : nv(augend)+nv(addend)
    add_vertices!(augend, nv(addend))
    for edge in edges(addend)
        weight = get_weight(addend, src(edge), dst(edge))
        @assert add_edge!(augend, id_mapping[src(edge)], id_mapping[dst(edge)], weight)
    end

    return nothing
end

end #module
