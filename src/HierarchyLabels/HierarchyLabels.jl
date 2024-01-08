module HierarchyLabels

using Graphs
using MLStyle
using Reexport
using StaticArrays

@reexport import Base: getindex, ==, ∈, ∋, ∉, show, println
@reexport import ..HMD: serialize, deserialize

export HLabel, id, type
export LabelResult, Label_Missing, Label_Occupied, Label_Duplication, Relation_Missing, Relation_Occupied, Success
export  LabelHierarchy, _labels, _add_label!, _add_labels!, _add_relation!, _remove_label!, _remove_relation!
export _label_unique, _label_nocycle, _label_connected
export _label2node, _contains, ∈, ∋, ∉, _has_relation ,_get_nodeid, getindex
export _issuper, _issub, _super_id, _sub_id, _super, _sub
export _root, _merge_hierarchy!
export PackedHierarchy

#####
##### HLabel definition
#####

struct HLabel
    type::String
    id::Int64
end

function HLabel(type::AbstractString, id::Integer)
    if !isascii(type)
        error("label type must be ascii string. ")
    end
    return HLabel(type, convert(Int64, id))
end

function id(label::HLabel)
    return label.id
end

function type(label::HLabel)
    return label.type
end

@data LabelResult begin
    Label_Missing
    Label_Occupied
    Label_Duplication
    Relation_Missing
    Relation_Occupied
    Cycle_Found
    Success
end

function ==(lhs::HLabel, rhs::HLabel)
    id(lhs) == id(rhs) && type(lhs) == type(rhs)
end

function Base.show(io::IO, ::MIME"text/plain", label::HLabel)
    i = id(label)
    t = type(label)
    print(io, "HLabel(\"$t\", $i)")
end

function Base.show(label::HLabel)
    i = id(label)
    t = type(label)
    print("HLabel(\"$t\", $i)")
end

function Base.print_to_string(label::HLabel)
    i = id(label)
    t = type(label)
    return "HLabel(\"$t\", $i)"
end

#####
##### LabelHiraraichy definition
#####

Base.@kwdef struct _Nodes
    super::Vector{Int64} = Int64[]
    sub::Vector{Set{Int64}} = Set{Int64}[]
end

"""

    LabelHierarchy

This type signifies a hierarchical structure of labels. This hierarchy is depicted through a
directed graph, with labels being maintained in a vector. The node id of the graph corresponds
to the index of the vector. The LabelHierarchy type includes a dictionary feature for
transforming a label into its corresponding node id.
"""
Base.@kwdef mutable struct LabelHierarchy
    g::_Nodes = _Nodes()
    labels::Vector{HLabel} = Vector{HLabel}()
    label2node::Dict{HLabel, Int64} = Dict{HLabel, Int64}()
end

function Base.length(nodes::_Nodes)
    @assert length(nodes.super) == length(nodes.sub)
    return length(nodes.super)
end

function Base.show(io::IO, ::MIME"text/plain", lh::LabelHierarchy)
    subtree(lh, label, level, indent) = begin
        label == _root(lh) && println("$label")
        for s in _sub(lh, label)
            str = join(fill(" ", level * indent)) * "$s"
            println(io, str)
            subtree(lh, s, level + 1, indent)
        end
    end
    if isempty(lh.labels)
        print()
    else
        subtree(lh, _root(lh), 1, 4)
    end
end

function Base.println(lh::LabelHierarchy)
    subtree(lh, label, level, indent) = begin
        label == _root(lh) && println("$label")
        labels = _sub(lh, label)
        nlines = length(labels)
        for (i, s) in enumerate(labels)
            bc = if _has_relation(lh, label, s)
                i == length(labels) ? "└" : "├"
            else
                " "
            end
            str = join(fill(" ", level * indent)) * bc * "$s\n"
            if degree(_nodes(lh), _get_nodeid(lh, s)) == 0
                printstyled(str; color=:red)
            else
                print(str)
            end
            subtree(lh, s, level + 1, indent)
        end
        #nlines
    end
    subtree(lh, _root(lh), 1, 4)
end

function _nodes(lh::LabelHierarchy)
    return lh.g
end

function _label2node(lh::LabelHierarchy)
    return lh.label2node
end

function _labels(lh::LabelHierarchy)
    return lh.labels
end

function _nv(lh::LabelHierarchy)
    @assert length(_nodes(lh)) == length(_labels(lh)) == length(_label2node(lh))
    return length(_nodes(lh))
end

function _nv(nodes::_Nodes)
    return length(nodes)
end

function _add_edge!(nodes::_Nodes, super::Integer, sub::Integer)
    if super == sub
        return false
    end
    if sub ∈ nodes.sub[super] # edge already exists
        return false
    end
    if nodes.super[sub] != 0 # edge already exists
        return false
    end

    push!(nodes.sub[super], sub)
    nodes.super[sub] = super

    return true
end

function _add_nodes!(nodes::_Nodes, n::Integer)
    @assert typemax(Int64) - length(nodes) ≥ n
    append!(nodes.super, zeros(Int64, n))
    append!(nodes.sub, [Set{Int64}() for _ in 1:n])
    return true
end

function _rem_edge!(nodes::_Nodes, super::Integer, sub::Integer)
    if super == sub
        println("\n==\n")
        return false
    end
    if sub ∉ nodes.sub[super] # edge does not exist
        println("\n$(sub)\n")
        return false
    end
    if nodes.super[sub] == 0 # edge does not exist
        println("\nsuper\n")
        return false
    end

    nodes.super[sub] = 0
    delete!(nodes.sub[super], sub)

    return true
end

function _edges(nodes::_Nodes)
    return ((super, id) for (id, super) in enumerate(nodes.super) if super != 0)
end

function _src(edge::Tuple{Int64, Int64}) # super
    return edge[1]
end

function _dst(edge::Tuple{Int64, Int64}) # sub
    return edge[2]
end

function _ne(nodes::_Nodes)
    return _nv(nodes) - 1
end

function _degree(nodes::_Nodes, id::Integer)
    return length(nodes.sub[id]) + (nodes.super[id] != 0 ? 1 : 0)
end

function _vertices(nodes::_Nodes)
    return eachindex(nodes.super)
end

"""

    _get_nodeid(lh, label)

Return the node id of the label in the LabelHierarchy.
"""
function _get_nodeid(lh::LabelHierarchy, label::HLabel)
    if _nv(_nodes(lh)) == 0
        error("There is no label in LabelHierarchy. ")
    end
    return _label2node(lh)[label]
end

function _get_label(lh::LabelHierarchy, id::Integer)
    if _nv(_nodes(lh)) == 0
        error("There is no label in LabelHierarchy. ")
    end
    return _labels(lh)[id]
end

function _set_label!(lh::LabelHierarchy, label::HLabel, id::Integer)
    if !_contains(lh, label)
        return Label_Missing
    end

    _labels(lh)[id] = label
    _label2node(lh)[label] = id

    return Success
end

function _add_label!(lh::LabelHierarchy, label::HLabel, unsafe::Bool=false)
    if !unsafe
        if _contains(lh, label)
            return Label_Occupied
        end
    end

    id = _nv(lh) + 1
    @assert _add_nodes!(_nodes(lh), 1)
    push!(_labels(lh), label)
    push!(_label2node(lh), label => id)

    return Success
end

function _add_labels!(lh::LabelHierarchy, labels::AbstractVector{HLabel})
    nv = _nv(lh)

    new_ids = (nv + one(nv)) : (nv + length(labels))
    @assert _add_nodes!(lh.g, length(labels))
    append!(lh.labels, labels)
    for (label, node_id) in zip(labels, new_ids)
        _label2node(lh)[label] = node_id
    end

    # if labels is not unique, unrecoverable error catched here
    if length(_label2node(lh)) != _nv(lh.g)
        error("labels are not unique. Hierarchy corrupted. ")
    else
        return Success
    end
end

function _add_relation!(lh::LabelHierarchy; super::HLabel, sub::HLabel, unsafe::Bool=false)
    if !unsafe
        if !_contains(lh, super) || !_contains(lh, sub)
            return Label_Missing
        elseif super == sub
            return Label_Duplication
        elseif _has_relation(lh, super, sub)
            return Relation_Occupied
        end
    end

    g, super_id, sub_id = _nodes(lh), _get_nodeid(lh, super), _get_nodeid(lh, sub)
    @assert _add_edge!(g, super_id, sub_id)

    if !unsafe && count(==(0), lh.g.super) < 1
        _rem_edge!(g, super_id, sub_id)
        return Cycle_Found
    end

    return Success
end

function _remove_label!(lh::LabelHierarchy, label::HLabel)
    if !_contains(lh, label)
        return Label_Missing
    end
    g, id = _nodes(lh), _get_nodeid(lh, label)
    @assert _nv(g) == length(_labels(lh)) == length(_label2node(lh))

    # when removing the ith node, the last node (== nv(g)) moves to the ith node
    super_node = g.super[id]
    sub_nodes = g.sub[id]
    delete!(_label2node(lh), label)
    if id == length(g)
        pop!(g.super); pop!(g.sub); pop!(_labels(lh))
        super_node != 0 && delete!(g.sub[super_node], id)
        for sid in sub_nodes
            g.super[sid] = 0
        end
    else
        tail = length(g)
        # update super node of tail
        if g.super[tail] != 0
            replace!(g.sub[g.super[tail]], tail => id)
        end
        # update sub nodes of tail
        for sid in g.sub[tail]
            @assert g.super[sid] == tail
            g.super[sid] = id
        end

        # update super node of id
        super_node != 0 && delete!(g.sub[super_node], id)
        # update sub nodes of id
        for sid in sub_nodes
            g.super[sid] = 0
        end

        lastnode = (pop!(g.super), pop!(g.sub))
        lastlabel = pop!(_labels(lh))
        g.super[id] = lastnode[1]
        g.sub[id] = lastnode[2]
        _labels(lh)[id] = lastlabel
        _label2node(lh)[lastlabel] = id
        @assert all(x -> tail ∉ x, g.sub)
    end

    return Success
end

function _remove_relation!(lh::LabelHierarchy, label1::HLabel, label2::HLabel)
    if !_contains(lh, label1) || !_contains(lh, label2)
        return Relation_Missing
    end

    g = _nodes(lh)
    n1, n2 = _get_nodeid(lh, label1), _get_nodeid(lh, label2)

    @assert _rem_edge!(g, n1, n2) || _rem_edge!(g, n2, n1)

    return Success
end

function _label_unique(lh::LabelHierarchy)
    @assert _label2node(lh) |> values |> allunique
    return allunique(_labels(lh))
end

function _label_nocycle(lh::LabelHierarchy)
    g = _nodes(lh)
    return count(==(0), g.super) == 1
end

function _label_connected(lh::LabelHierarchy)
    g = _nodes(lh)
    return all(id -> _degree(g, id) != 0, _vertices(g))
end

function _contains(lh::LabelHierarchy, label::HLabel)
    return haskey(_label2node(lh), label)
end

function ∈(label::HLabel, lh::LabelHierarchy)
    return _contains(lh, label)
end

function ∋(lh::LabelHierarchy, label::HLabel)
    return _contains(lh, label)
end

function ∉(label::HLabel, lh::LabelHierarchy)
    return !_contains(lh, label)
end

function _has_relation(lh::LabelHierarchy, label1::HLabel, label2::HLabel)
    if !(_contains(lh, label1) && _contains(lh, label2))
        error("labels not found in LabelHierarchy. ")
    end
    return _issuper(lh, label1, label2) || _issub(lh, label1, label2)
end

function _super_id(lh::LabelHierarchy, id::Integer)
    return _nodes(lh).super[id]
end

function _sub_id(lh::LabelHierarchy, id::Integer)
    return _nodes(lh).sub[id]
end

function _super_id(lh::LabelHierarchy, label::HLabel)
    id = _get_nodeid(lh, label)
    return _nodes(lh).super[id]
end

function _sub_id(lh::LabelHierarchy, label::HLabel)
    id = _get_nodeid(lh, label)
    return _nodes(lh).sub[id]
end

function _super(lh::LabelHierarchy, label::HLabel)
    super_id = _super_id(lh, label)
    return if super_id == 0
        nothing
    else
        _labels(lh)[super_id]
    end
end

function _sub(lh::LabelHierarchy, label::HLabel)
    sub_ids = _sub_id(lh, label)
    return [_labels(lh)[id] for id in sub_ids]
end

function _issuper(lh::LabelHierarchy, lhs::HLabel, rhs::HLabel)
    return _get_nodeid(lh, lhs) == _super_id(lh, rhs)
end

function _issub(lh::LabelHierarchy, lhs::HLabel, rhs::HLabel)
    return _get_nodeid(lh, lhs) ∈ _sub_id(lh, rhs)
end

function ==(lhs::LabelHierarchy, rhs::LabelHierarchy)
    if Set(_labels(lhs)) != Set(_labels(rhs))
        return false
    end

    lg= _nodes(lhs)
    for lhs_id in _vertices(lg)
        if lg.super[lhs_id] == 0 # root node
            continue
        end
        lhs_label = _get_label(lhs, lhs_id)
        _super_labels = _super(lhs, lhs_label)
        _sub_labels = _sub(lhs, lhs_label) |> Set
        @assert lhs_label ∈ rhs
        if _super_labels != _super(rhs, lhs_label) || _sub_labels != Set(_sub(rhs, lhs_label))
            return false
        end
    end

    return true
end

function _root_id(lh::LabelHierarchy)
    g = _nodes(lh)
    id = findall(==(0), g.super)
    @assert length(id) == 1

    return id[1]
end

function _root(lh::LabelHierarchy)
    i = _root_id(lh)
    return _get_label(lh, i)
end


function _merge_hierarchy!(
    augend::LabelHierarchy,
    addend::LabelHierarchy;
    augend_parent::HLabel,
    addend_parent::HLabel,
    unsafe::Bool = false
)
    # addend_parent自身とその上位ノードはマージから除外する
    exception = _find_exception(addend, addend_parent)
    @assert allunique(exception)
    @assert !isnothing(exception)

    # addend node id => augend node id
    node_mapping = Dict{Int64, Int64}()
    forward_shift = _nv(augend)
    back_shift = 0
    for id in 1:_nv(addend)
        if id ∈ exception
            back_shift += 1
            continue
        end
        node_mapping[id] = id + forward_shift - back_shift
    end

    # augendに既に存在するラベル数を記録
    # augendになくaddendにあるラベルのcounterは0
    counter = Dict{String, Int64}() #
    for label in _labels(augend)
        if !haskey(counter, type(label))
            counter[type(label)] = 0
        end
        counter[type(label)] += 1
    end
    for label in _labels(addend)
        if !haskey(counter, type(label))
            counter[type(label)] = 0
        end
    end

    # merge vertices and labels
    g_augend, g_addend = _nodes(augend), _nodes(addend)
    _add_nodes!(g_augend, _nv(g_addend) - length(exception))
    resize!(_labels(augend), _nv(g_augend))

    # merge edges
    for edge in _edges(g_addend)
        if _src(edge) ∉ exception && _dst(edge) ∉ exception
            _add_edge!(g_augend, node_mapping[_src(edge)], node_mapping[_dst(edge)])
        end
    end
    # connect top of addend hierarchy and augend_parent
    # number of addend part top may be >= 2
    addend_part_top = [node_mapping[id] for id in _sub_id(addend, addend_parent)]
    augend_parent_id = _get_nodeid(augend, augend_parent)
    for id in addend_part_top
        _add_edge!(
            g_augend,
            augend_parent_id, # super
            id, # sub
        )
    end

    # label
    for (old_id, new_id) in sort(collect(pairs(node_mapping)), by=x->x[2])
        ltype = _get_label(addend, old_id) |> type
        counter[ltype] += 1
        _labels(augend)[new_id] = HLabel(ltype, counter[ltype])
        push!(_label2node(augend), HLabel(ltype, counter[ltype]) => new_id)
    end

    if !unsafe
        _label_unique(augend) || error("label is not unique. ")
        _label_nocycle(augend) || error("label hierarchy has cycle. ")
        _label_connected(augend) || error("isolated label found. ")
    end

    return nothing
end

"""
木構造addendのうち、setdiff(addend, addend_parentより先の枝)のnode idを返す
"""
function _find_exception(addend, addend_parent)
    addend_parent_id = _get_nodeid(addend, addend_parent)
    # depth first search from addend root
    # skips sub branch of addend_parent
    buffer = [_root_id(addend)]
    exception = Set{Int64}()
    while !isempty(buffer)
        id = popfirst!(buffer)
        push!(exception, id)
        if id == addend_parent_id
            # added parentより上位のノードで別の枝が分岐している場合はその枝すべてが除外対象となる
            # これを保証するためにここではbreakせず、別の枝の探索を続ける
            continue
        else
            prepend!(buffer, _sub_id(addend, id))
        end
    end

    return exception
end

#####
##### hierarchy serialization for file storage
#####

struct PackedHierarchy
    num_node::Int64
    edges_org::Vector{Int64} # super
    edges_dst::Vector{Int64} # sub
    label_ids::Vector{Int64}
    # splitting Vector{String}
    chars::Vector{UInt8}
    bounds::Vector{Int64}
end

function serialize(lh::LabelHierarchy)
    # split Vecotr{HLabel} to Vector{Int64} and Vector{String}
    label_ids = [id(label) for label in _labels(lh)]
    chars, bounds = serialize([type(label) for label in _labels(lh)])

    g = _nodes(lh)
    num_node = _nv(g)
    edges_org, edges_dst = let
        orig, dest = Vector{Int64}(undef, _ne(g)), Vector{Int64}(undef, _ne(g))
        for (i, edge) in enumerate(_edges(g))
            orig[i], dest[i] = _src(edge), _dst(edge)
        end
        orig, dest
    end

    return PackedHierarchy(num_node, edges_org, edges_dst, label_ids, chars, bounds)
end

function deserialize(ph::PackedHierarchy)
    g = begin
        g = _Nodes()
        _add_nodes!(g, ph.num_node)
        for (o, d) in zip(ph.edges_org, ph.edges_dst)
            _add_edge!(g, o, d)
        end
        g
    end

    label_ids = ph.label_ids
    label_types = deserialize(ph.chars, ph.bounds)
    labels = [HLabel(type, id) for (type, id) in zip(label_types, label_ids)]

    lh = LabelHierarchy()
    lh.g = g
    lh.labels = labels
    lh.label2node = Dict(labels .=> 1:length(labels))

    return lh
end



include("test.jl")

end #module
