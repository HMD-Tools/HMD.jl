module HMD

using Reexport

@reexport using Graphs
using HDF5
using LinearAlgebra
@reexport using PeriodicTable
@reexport using SimpleWeightedGraphs
using StaticArrays
@reexport using Unitful

@reexport import Base: getindex, firstindex, lastindex, setproperty!, iterate, length, precision, close, string, show, showerror
@reexport import Base: >, <, >=, <=, +, -, *, /, ==, position
@reexport import Base: time, contains, promote_type, promote_rule, similar, merge!
@reexport import Graphs: neighbors
@reexport import Base: sort, in, ∈

# utility
export oblique_coord, atom_label

# LabelHierarchy
export HLabel, LabelHierarchy, id, type

export
    # System constants
    Entire_System,
    Atom_Label,

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
    bond_order,
    box,
    set_box!,
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
    add_prop!,
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
    wrapped,

    # trajectory io interface
    add_snapshot!,
    get_metadata,

    # System highlevel interfaces
    is_atom,
    atom_label,
    add_atom!,
    add_atoms!,
    add_bond!,
    add_bonds!,
    bond_order,
    set_bondorder!,
    valence,
    neighbors,
    l2a,
    super_labels,
    sub_labels,
    hierarchy_leaf_isatom,

    # Trajectory highlevel interfaces
    length,
    getindex,
    firstindex,
    lastindex,
    getindex,
    iterate,
    hmdsave,
    read_traj,
    snapshot

include("util.jl")

include("HierarchyLabels/HierarchyLabels.jl")
using .HierarchyLabels

include("interface/interface.jl")

include("DataTypes/DataTypes.jl")
@reexport using .DataTypes

include("interface/system_highlevel.jl")
include("interface/trajectory_highlevel.jl")

#include("GenericIO/GenericIO.jl")
#using .GenericIO

end
