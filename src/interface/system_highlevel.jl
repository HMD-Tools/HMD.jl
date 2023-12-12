@reexport import .DataTypes: super, sub

const Entire_System = HLabel("entire_system", 1)

#const atom_mass = Dict{String, Float64}(
#    elements[:H ].symbol => 1.008,
#    elements[:C ].symbol => 12.012,
#    elements[:N ].symbol => 14.007,
#    elements[:O ].symbol => 16.000,
#    elements[:F ].symbol => 19.000,
#    elements[:Si].symbol => 28.086,
#    elements[:S ].symbol => 32.067,
#    elements[:Cl].symbol => 35.453
#)

function add_atom!(s::AbstractSystem, x::AbstractVector{<:AbstractFloat}, elem::Integer; super::HLabel)
    atom_id = natom(s) + 1
    add_position!(s, x)
    add_element!(s, elem)
    @assert add_vertex!(topology(s))
    for hname in hierarchy_names(s)
        add_label!(s, hname, atom_label(atom_id))
        add_relation!(s, hname; super=super, sub=atom_label(atom_id))
    end

    return nothing
end

function add_atoms!(s::AbstractSystem, x::AbstractVector{<:AbstractVector{<:AbstractFloat}}, elem::AbstractVector{<:Integer}; super::HLabel)
    if length(elem) != length(x)
        error("length of elem and x must be same. ")
    end

    front = natom(s) + 1
    back = natom(s) + length(x)
    atom_labels = atom_label.(front:back)
    add_positions!(s, x)
    add_elements!(s, elem)
    @assert add_vertices!(topology(s), length(x))
    for hname in hierarchy_names(s)
        add_labels!(s, hname, atom_labels)
        add_relations!(s, hname; super=super, subs=atom_labels)
    end

    return nothing
end

function add_bond!(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer; bond_order::Rational=1//1)
    if !(0//1 <= bond_order <= 8//1)
        error("bond order must be in [0, 8] ")
    end
    topo = topology(s)
    @assert add_edge!(topo, atom_id1, atom_id2, bond_order)

    return nothing
end

#function add_bonds!(s::AbstractSystem, pair::AbstractVector{Tuple{<:Integer, <:Integer, Rational{<:Integer}}})
function add_bonds!(s::AbstractSystem, pair)
    if any(p -> !(0//1 <= p[3] <= 8//1), pair)
        error("bond order must be in [0, 8] ")
    end

    topo = topology(s)
    for (atom_id1, atom_id2, bond_order) in pair
        @assert add_edge!(topo, atom_id1, atom_id2, bond_order)
    end
end

function add_bond!(s::AbstractSystem, label1::HLabel, label2::HLabel; bond_order::Rational=1//1)
    if !isatom(label1) || !isatom(label2)
        error("label is not for atom. ")
    end

    atom_id1 = convert(Int64, id(label1))
    atom_id2 = convert(Int64, id(label2))
    add_bond!(s, atom_id1, atom_id2; bond_order=bond_order)

    return nothing
end

# core apiに移管
function bond_order(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer)
    topo = topology(s)
    if !has_edge(topo, atom_id1, atom_id2)
        error("There is no bond beteen atoms $(atom_id1) and $(atom_id2). ")
    end

    return get_weight(topo, atom_id1, atom_id2)
end

# core apiに移管
function bond_order(s::AbstractSystem, label1::HLabel, label2::HLabel)
    if !isatom(label1) || !isatom(label2)
        error("label is not for atom. ")
    end
    atom_id1 = convert(Int64, id(label1))
    atom_id2 = convert(Int64, id(label2))

    return bond_order(s, atom_id1, atom_id2)
end

# core apiに移管
# topologyがAbstractGraphに従う限り現在のcore apiが動けばhighlevel apiも動くので移管しなくてもよい？
function set_bondorder!(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer, bo::Rational{<:Integer})
    topo = topology(s)
    @assert rem_edge!(topo, atom_id1, atom_id2)
    @assert add_edge!(topo, atom_id1, atom_id2, bo)
    return nothing
end

function valence(s::AbstractSystem, atom_id::Integer)
    topo = topology(s)
    valence = 0//1
    for neigh_id in all_neighbors(topo, atom_id)
        valence += get_weight(topo, atom_id, neigh_id)
    end

    return valence
end

function valence(s::AbstractSystem, label::HLabel)
    if !isatom(label)
        error("label $label is not for atom. ")
    end

    return valence(s, convert(Int64, id(label1)))
end

function atom_label(atom_id::Integer)
    return HLabel("", atom_id)
end

function isatom(label::HLabel)
    return type(label) == ""
end

function neighbors(s::AbstractSystem, atom_id::Integer)
    topo = topology(s)
    return all_neighbors(topo, atom_id)
end

function neighbors(s::AbstractSystem, label::HLabel)
    if !isatom(label)
        error("label $label is not for atom. ")
    end
    topo = topology(s)

    return all_neighbors(topo, convert(Int64, id(label1)))
end

function l2a(s::AbstractSystem, hname::AbstractString, label::HLabel)
    sub_labels = [label]
    atom_ids = Vector{Int64}(undef, 0)

    while any(!isatom, sub_labels)
        sub_next = Vector{HLabel}(undef, 0)
        for l in sub_labels
            if isatom(l)
                push!(atom_ids, convert(Int64, id(l)))
            else
                append!(sub_next, sub(s, hname, l))
            end
        end
        sub_labels = sub_next
    end
    append!(atom_ids, convert.(Int64, id.(sub_labels)))

    return atom_ids
end

function contains(atom_id::Integer, s::AbstractSystem, hname::AbstractString, label::HLabel)
    return contains(atom_label(atom_id), s, hname, label)
end

function contains(atom::HLabel, s::AbstractSystem, hname::AbstractString, label::HLabel)
    if isatom(atom)
        error("expected atom label, found $(atom). ")
    end

    return label ∈ super_labels(s, hname, atom)
end

function super_labels(s::AbstractSystem, hname::AbstractString, label::HLabel)
    return _traverse_from(s, hname, label, super)
end

function sub_labels(s::AbstractSystem, hname::AbstractString, label::HLabel)
    return _traverse_from(s, hname, label, sub)
end

function _traverse_from(s::AbstractSystem, hname::AbstractString, label::HLabel, func::Function)
    labels = Vector{HLabel}(undef, 0)

    # func is either super or sub
    stack = [label]
    while !isempty(stack)
        current = popfirst!(stack)
        next = func(s, hname, current)
        prepend!(stack, next)
        push!(labels, current)
    end

    return labels
end

function hierarchy_leaf_isatom(s::AbstractSystem)
    for hname in hierarchy_names(s)
        for label in all_labels(s, hname)
            if !isatom(label) && isempty(sub(s, hname, label))
                return false
            elseif isatom(label) && !isempty(sub(s, hname, label))
                return false
            end
        end
    end

    return true
end
