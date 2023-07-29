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
using Reexport
using SimpleWeightedGraphs
using StaticArrays
using Unitful

using ..HierarchyLabels

@reexport import Base: *, +, -, /, <, <=, ==, >, >=, close, contains, convert, getindex,firstindex, lastindex, iterate,
    length, position, precision, promote_rule, promote_type, setproperty!, show, similar,
    string, time, ∈, ∉, merge!

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
    set_position!,
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
    insert_relations!,
    remove_label!,
    remove_relation!,
    contains,
    issuper,
    issub,
    super,
    sub,

    # system property interface
    prop_names,
    prop,
    set_prop!,

    # system io interface
    AbstractFileFormat,
    close,

    # trajectory interface
    empty_trajectory,
    get_system,
    all_timesteps,
    get_timestep,
    length,
    add!,
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
    add_snapshot!,
    latest_reaction_step,
    get_timesteps,
    get_reactions,
    get_metadata,
    is_reaction,
    length,
    getindex,
    lastindex,
    firstindex,
    wrapped,

    # trajectory io interface
    add_snapshot!,
    import_dynamic!,
    import_static!,
    latest_reaction_step,
    get_timesteps,
    get_reactions,
    get_metadata,
    is_reaction,
    length,
    getindex,
    lastindex,
    firstindex,
    wrapped

# core subtype signature
export Position, BoundingBox, HLabel, LabelHierarchy

# core immut signature
export GeneralSystem, System, print_to_string

# fileIO
export H5system, H5traj, SerializedTopology, PackedHierarchy
export h5system, h5traj, get_file

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

mutable struct System{D, F<:AbstractFloat, SysType<:AbstractSystemType, L} <: AbstractSystem{D, F, SysType}
    time::F
    topology::SimpleWeightedGraph{Int64, Rational{BO_Precision}}
    box::BoundingBox{D, F, L}

    # atom property
    position::Position{D, F}
    travel::Vector{SVector{D, Int16}}
    wrapped::Bool
    element::Vector{Atomic_Number_Precision}

    hierarchy::Dict{String, LabelHierarchy}
    props::Dict{String, Array}
end

function System{D, F, SysType}() where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    System{D, F, SysType, D*D}(
        zero(F),
        SimpleWeightedGraph{Int64, Rational{BO_Precision}}(),
        BoundingBox{D, F}(),
        Position{D, F}(),
        Vector{SVector{D, Int16}}(undef, 0),
        false,
        Atomic_Number_Precision[],
        Dict{String, LabelHierarchy}(),
        Dict{String, Array}()
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

function Base.similar(s::System{D, F, SysType, L}; reserve_dynamic::Bool=false, reserve_static::Bool=false) where {D, F, SysType, L}
    sim = System{D, F, SysType}()
    if reserve_dynamic
        set_time!(sim, time(s))
        set_box!(sim, box(s))
        sim.position = s.position
        sim.travel = s.travel
        sim.wrapped = s.wrapped
        sim.props = s.props
    end
    if reserve_static
        sim.topology = s.topology
        sim.hierarchy = s.hierarchy
        sim.element = s.element
    end
    return deepcopy(sim)
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
    #@assert natm == length(s.element)
    return length(s.position)
end

function nbond(s::System)
    return topology(s) |> ne
end

function time(s::System)
    s.time
end

function set_time!(s::System, time::AbstractFloat)
    s.time = time
end

function topology(s::System)
    s.topology
end

function box(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    s.box
end

function set_box!(s::System, box::BoundingBox)
    s.box = box
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
        error("atom addition with wrapped coordinates is not supprted. ")
    end
    push!(s.travel, zeros(Int16, 3))

    return nothing
end

function add_positions!(s::System, x::AbstractVector{<:AbstractVector{<:AbstractFloat}})
    append!(all_positions(s), x)
    if wrapped(s)
        error("atom addition with wrapped coordinates is not supprted. ")
    end
    append!(s.travel, [zeros(Int16, 3) for _ in 1:length(x)])

    return nothing
end

function set_position!(s::System, atom_id::Integer, x::AbstractVector{<:AbstractFloat})
    s.position[atom_id] = x
end

function set_position!(s::System, label::HLabel, x::AbstractVector{<:AbstractFloat})
    if !is_atom(label)
        error("label $label is not for atom. ")
    end
    set_position!(s, label, x)

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

function _frac_int_vector(mods::SVector{D, Tuple{F, F}}) where {D, F<:AbstractFloat}
    return SVector{D, F}(mods[i][1] for i in 1:D), SVector{D, Int16}(Int16(mods[i][2]) for i in 1:D)
end

function wrap!(s::System{D, F, SysType, L}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    if wrapped(s)
        return nothing
    end

    axis = box(s).axis
    origin = box(s).origin
    ## x = c[1] .* axis[:,1] .+ c[2] .* axis[:,2] .+ ...
    e_i_e_j = SMatrix{D, D, F, D*D}(
        dot(axis[:,i], axis[:,j]) for i in 1:D, j in 1:D
    ) #|> Symmetric
    for id in 1:natom(s)
        x = position(s, id) - origin
        c = e_i_e_j \ SVector{D, F}(dot(x, axis[:,dim]) for dim in 1:D)
        fparts, iparts = _frac_int_vector(modf.(c))
        pos = map(1:D) do dim
            @inbounds fparts[dim] * axis[:,dim]
        end |> p->reduce(+, p)
        set_travel!(s, id, iparts)
        set_position!(s, id, pos + origin)
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
        pos = x + mapreduce(dim -> n[dim] * axis[:, dim], +, 1:D)
        set_position!(s, i, pos + origin)
        set_travel!(s, i, zeros(Int16, 3))
    end
    _change_wrap!(s)

    return nothing
end

function label2atom(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    labels = _labels(lh)
    atom_ids = Int64[]

    stack = [_get_nodeid(lh, label)]
    if is_atom(labels[stack[1]])
        return [stack[1]]
    end

    while !isempty(stack)
        current_node = popfirst!(stack)
        next_nodes = _sub_id(lh, current_node)
        for node in next_nodes
            label = labels[node]
            if is_atom(label)
                push!(atom_ids, id(label))
            else
                pushfirst!(stack, node)
            end
        end
    end

    return atom_ids
end

function is_atom(label::HLabel)
    return type(label) == ""
end

include("label_manipulation.jl")
include("property.jl")
include("system_io.jl")
include("trajectory.jl")
include("test.jl")

function merge!(
    augend::System{D, F, SysType1}, addend::System{D, F, SysType2};
    augend_parent::HLabel, addend_parent::HLabel, unsafe::Bool=false
) where {D, F<:AbstractFloat, SysType1<:AbstractSystemType, SysType2<:AbstractSystemType}
    # wrap check
    wrapped(augend) != wrapped(addend) && error("augend and addend must have same wrap status. ")

    # add new atoms to augend
    append!(augend.position, all_positions(addend))
    append!(augend.travel, all_travels(addend))
    append!(augend.element, all_elements(addend))

    merge_topology!(topology(augend), topology(addend))

    for hname in hierarchy_names(augend)
        lh_augend = hierarchy(augend, hname)
        lh_addend = hierarchy(addend, hname)
        _merge_hierarchy!(
            lh_augend, lh_addend;
            augend_parent = augend_parent,
            addend_parent = addend_parent,
            unsafe = unsafe
        )
    end

    return nothing
end

function merge_topology!(augend::SimpleWeightedGraph, addend::SimpleWeightedGraph)
    id_mapping = nv(augend)+1 : nv(augend)+nv(addend)
    add_vertices!(augend, nv(addend))
    for edge in edges(addend)
        weight = get_weight(addend, src(edge), dst(edge))
        @assert add_edge!(augend, id_mapping[src(edge)], id_mapping[dst(edge)], weight)
    end

    return nothing
end

end #module
