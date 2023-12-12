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
export _root, _depth, _merge_hierarchy!
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

"""

    LabelHierarchy

This type signifies a hierarchical structure of labels. This hierarchy is depicted through a
directed graph, with labels being maintained in a vector. The node id of the graph corresponds
to the index of the vector. The LabelHierarchy type includes a dictionary feature for
transforming a label into its corresponding node id.
"""
Base.@kwdef mutable struct LabelHierarchy
    g::DiGraph{Int64} = DiGraph{Int64}()
    labels::Vector{HLabel} = Vector{HLabel}()
    label2node::Dict{HLabel, Int64} = Dict{HLabel, Int64}()
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
            if degree(_hierarchy(lh), _get_nodeid(lh, s)) == 0
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

function _hierarchy(lh::LabelHierarchy)
    return lh.g
end

function _label2node(lh::LabelHierarchy)
    return lh.label2node
end

function _labels(lh::LabelHierarchy)
    return lh.labels
end

"""

    _get_nodeid(lh, label)

Return the node id of the label in the LabelHierarchy.
"""
function _get_nodeid(lh::LabelHierarchy, label::HLabel)
    if length(_hierarchy(lh)) == 0
        error("There is no label in LabelHierarchy. ")
    end
    return _label2node(lh)[label]
end

function _get_label(lh::LabelHierarchy, id::Integer)
    if length(_hierarchy(lh)) == 0
        error("There is no label in LabelHierarchy. ")
    end
    return _labels(lh)[id]
end

function _set_label!(lh::LabelHierarchy, label::HLabel, id::Integer)
    _labels(lh)[id] = label
    _label2node(lh)[label] = id
end

function _add_label!(lh::LabelHierarchy, label::HLabel, unsafe::Bool=false)
    g = _hierarchy(lh)

    if !unsafe
        if _contains(lh, label)
            return Label_Occupied
        end
    end

    @assert add_vertex!(g)

    current_id = nv(g)
    push!(_labels(lh), label)
    push!(_label2node(lh), label => current_id)

    return Success
end

function _add_labels!(lh::LabelHierarchy, labels::AbstractVector{HLabel})
    g = _hierarchy(lh)

    new_range = (nv(g) + one(nv(g))) : (nv(g) + length(labels))
    @assert add_vertices!(g, length(labels)) == length(labels)
    append!(lh.labels, labels)
    for (label, node_id) in zip(labels, new_range)
        _label2node(lh)[label] = node_id
    end

    if length(_label2node(lh)) != nv(g)
        return Label_Occupied
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

    g, super_id, sub_id = _hierarchy(lh), _get_nodeid(lh, super), _get_nodeid(lh, sub)
    @assert add_edge!(g, sub_id, super_id)

    if !unsafe && is_cyclic(g)
        rem_edge!(g, sub_id, super_id)
        return Cycle_Found
    end

    return Success
end

function _remove_label!(lh::LabelHierarchy, label::HLabel)
    if !_contains(lh, label)
        return Label_Missing
    end
    g, id = _hierarchy(lh), _get_nodeid(lh, label)

    # when removing the ith node, the last node (== nv(g)) moves to the ith node
    @assert rem_vertex!(g, id)
    end_label = pop!(_labels(lh))
    delete!(_label2node(lh), label)
    if id != nv(g)
        _set_label!(lh, end_label, id)
    end

    return Success
end

function _remove_relation!(lh::LabelHierarchy, label1::HLabel, label2::HLabel)
    if !_contains(lh, label1) || !_contains(lh, label2)
        return Relation_Missing
    end

    g = _hierarchy(lh)
    n1, n2 = _get_nodeid(lh, label1), _get_nodeid(lh, label2)

    @assert has_edge(g, n2, n1)
    @assert rem_edge!(g, n1, n2) || rem_edge!(g, n2, n1)

    return Success
end

function _label_unique(lh::LabelHierarchy)
    @assert _label2node(lh) |> values |> allunique
    return allunique(_labels(lh))
end

function _label_nocycle(lh::LabelHierarchy)
    return !is_cyclic(_hierarchy(lh))
end

function _label_connected(lh::LabelHierarchy)
    return is_connected(_hierarchy(lh))
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
    return outneighbors(_hierarchy(lh), id)
end

function _sub_id(lh::LabelHierarchy, id::Integer)
    return inneighbors(_hierarchy(lh), id)
end

function _super_id(lh::LabelHierarchy, label::HLabel)
    id = _get_nodeid(lh, label)
    return outneighbors(_hierarchy(lh), id)
end

function _sub_id(lh::LabelHierarchy, label::HLabel)
    id = _get_nodeid(lh, label)
    return inneighbors(_hierarchy(lh), id)
end

function _super(lh::LabelHierarchy, label::HLabel; recurse::Bool=false)
    super_ids = _super_id(lh, label)
    return _labels(lh)[super_ids]
end

function _sub(lh::LabelHierarchy, label::HLabel; recurse::Bool=false)
    sub_ids = _sub_id(lh, label)
    return _labels(lh)[sub_ids]
end

function _issuper(lh::LabelHierarchy, lhs::HLabel, rhs::HLabel)
    return _get_nodeid(lh, lhs) ∈ _super_id(lh, rhs)
end

function _issub(lh::LabelHierarchy, lhs::HLabel, rhs::HLabel)
    return _get_nodeid(lh, lhs) ∈ _sub_id(lh, rhs)
end

function ==(lhs::LabelHierarchy, rhs::LabelHierarchy)
    if Set(_labels(lhs)) != Set(_labels(rhs))
        println("labels")
        return false
    end

    lg, rg = _hierarchy(lhs), _hierarchy(rhs)
    for lhs_id in vertices(lg)
        lhs_label = _get_label(lhs, lhs_id)
        super_labels = _super(lhs, lhs_label) |> Set
        sub_labels = _sub(lhs, lhs_label) |> Set
        if lhs_label ∉ rhs
            return false
        end
        rhs_label = _get_label(rhs, _get_nodeid(rhs, lhs_label))
        if super_labels != Set(_super(rhs, rhs_label)) || sub_labels != Set(_sub(rhs, rhs_label))
            return false
        end
    end

    return true
end

function _root_id(lh::LabelHierarchy)
    g = _hierarchy(lh)
    root_id = filter(i -> isempty(_super_id(lh, i)), 1:nv(g))

    #@assert length(root_id) == 1
    return root_id[1]
end

function _root(lh::LabelHierarchy)
    i = _root_id(lh)
    return _get_label(lh, i)
end

# treeではないのでDFSは使えない
function _depth(lh::LabelHierarchy)
    sub_ids = _sub_id(lh, _root_id(lh))
    depth = 0
    while !isempty(sub_ids)
        depth += 1
        sub_ids = mapreduce(i -> _sub_id(lh, i), append!, sub_ids)
    end

    # excluding atom label
    return depth - 1
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

    # addend node id => augend node id
    node_mapping = Dict{Int64, Int64}()
    forward_shift = nv(_hierarchy(augend))
    back_shift = 0
    for id in 1:nv(_hierarchy(addend))
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
    g_augend, g_addend = _hierarchy(augend), _hierarchy(addend)
    add_vertices!(g_augend, nv(g_addend) - length(exception))
    resize!(_labels(augend), nv(g_augend))

    # merge edges
    for edge in edges(g_addend)
        if src(edge) ∉ exception && dst(edge) ∉ exception
            add_edge!(g_augend, node_mapping[src(edge)], node_mapping[dst(edge)])
        end
    end
    # connect top of addend hierarchy and augend_parent
    # number of addend part top may be >= 2
    addend_part_top = [node_mapping[id] for id in _sub_id(addend, addend_parent)]
    augend_parent_id = _get_nodeid(augend, augend_parent)
    for id in addend_part_top
        add_edge!(
            g_augend,
            id, # sub
            augend_parent_id # super
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
    exception = Int64[]
    while !isempty(buffer)
        id = popfirst!(buffer)
        pushfirst!(exception, id)
        if id == addend_parent_id
            # added parentより上位のノードで別の枝が分岐している場合はその枝すべてが除外対象となる
            # これを保証するためにここではbreakせず、別の枝の探索を続ける
            continue
        else
            prepend!(buffer, _sub_id(addend, id))
        end
    end

    if length(exception) > 50
        @warn "nodes excluded from merge is larger tha 50.\n" *
            "This may cause performance issue."
    end

    return SVector{length(exception), Int64}(exception)
end



#####
##### hierarchy serialization for file storage
#####

struct PackedHierarchy
    num_node::Int64
    edges_org::Vector{Int64}
    edges_dst::Vector{Int64}
    label_ids::Vector{Int64}
    # splitting Vector{String}
    chars::Vector{UInt8}
    bounds::Vector{Int64}
end

function serialize(lh::LabelHierarchy)
    # split Vecotr{HLabel} to Vector{Int64} and Vector{String}
    label_ids = [id(label) for label in _labels(lh)]
    chars, bounds = serialize([type(label) for label in _labels(lh)])

    g = _hierarchy(lh)
    num_node = nv(g)
    edges_org, edges_dst = let
        orig, dest = Vector{Int64}(undef, ne(g)), Vector{Int64}(undef, ne(g))
        for (i, edge) in enumerate(edges(g))
            orig[i], dest[i] = src(edge), dst(edge)
        end
        orig, dest
    end

    return PackedHierarchy(num_node, edges_org, edges_dst, label_ids, chars, bounds)
end

function deserialize(ph::PackedHierarchy)
    g = begin
        g = DiGraph()
        add_vertices!(g, ph.num_node)
        for (o, d) in zip(ph.edges_org, ph.edges_dst)
            add_edge!(g, o, d)
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
