function hierarchy_names(s::System)
    s.hierarchy |> keys |> collect
end

function hierarchy(s::System, hname::AbstractString)
    s.hierarchy[hname]
end

function add_hierarchy!(s::System, hname::AbstractString)
    if hname in hierarchy_names(s)
        error("hierarchy $(hname) already exists. ")
    end
    push!(s.hierarchy, hname => LabelHierarchy())
    #_add_label!(hierarchy(s, hname), Entire_System)
    add_label!(s, hname, Entire_System)

    return nothing
end

function add_hierarchy!(s::System, hnames::Union{Tuple, Base.Generator, V}) where {V<:AbstractVector{<:AbstractString}}
    for hname in hnames
        add_hierarchy!(s, hname)
    end

    return nothing
end

function remove_hierarchy!(s::System, hname::AbstractString)
    delete!(s.hierarchy, hname)
end

function all_labels(s::System, hname::AbstractString)
    lh = hierarchy(s, hname)
    return deepcopy(_labels(lh))
end

function all_labels(s::System, hname::AbstractString, label_type::AbstractString)
    #labels = hierarchy(s, hname) |> HierarchyLabels._label2node |> keys |> collect
    labels = hierarchy(s, hname) |> HierarchyLabels._labels

    return filter(label -> type(label)==label_type, labels)
end

function add_label!(s::System, hname::AbstractString, label::HLabel; unsafe::Bool=false)
    lh = hierarchy(s, hname)

    result = _add_label!(lh, label, unsafe)
    @match result begin
        Label_Occupied => error("label $(label) already exists in $(hname). ")
        Success        => return nothing
        _              => error("fatal error")
    end
end

function replace!(s::System, hname::AbstractString, old_new::Pair{HLabel, HLabel})
    lh = hierarchy(s, hname)

    result = _replace!(lh, old_new)
    @match result begin
        Label_Missing     => error("label $(old_new[1]) not found in $(hname). ")
        Label_Duplication => return s
        Label_Occupied    => error("label $(old_new[2]) already exists in $(hname). ")
        Success           => return s
        _                 => error("fatal error")
    end
end


#function add_label!(s::System, hname::AbstractString, label_type::AbstractString)
#    lh = hierarchy(s, hname)
#
#    n = count_label(s, hname, label_type)
#    addend = HLabel(label_type, n+1)
#    result = _add_label!(lh, addend)
#    @match result begin
#        Label_Occupied => error("label $(label) already exists. ")
#        Success        => return addend
#        _              => error("fatal error")
#    end
#end

function add_labels!(s::System, hname::AbstractString, labels::Union{Tuple, Base.Generator, V}) where {V<:AbstractVector{HLabel}}
    lh = hierarchy(s, hname)

    _add_labels!(lh, labels)

    return nothing
end

function count_label(s::System, hname::AbstractString, label_type::String)
    lh = hierarchy(s, hname)
    labels = _label2node(lh) |> keys # uniqueness is guaranteed
    return count(l -> type(l)==label_type, labels)
end

function add_relation!(
    s::System,
    hname::AbstractString;
    super::HLabel,
    sub::HLabel,
    unsafe::Bool = false
)
    lh = hierarchy(s, hname)

    result = _add_relation!(lh; super=super, sub=sub, unsafe=unsafe)
    @match result begin
        Label_Missing     => error("Super or sub not found. ")
        Label_Duplication => error("super and sub are equal. ")
        Relation_Occupied => error("""There in already relation between super and sub labels. Please use "insert_relation!()" instead. """)
        success           => return nothing
        _                 => error("fatal error")
    end
end

function add_relations!(s::System, hname::AbstractString; super::HLabel, subs::Union{Tuple, Base.Generator, V}) where {V<:AbstractVector{HLabel}}
    lh = hierarchy(s, hname)
    for sub in subs
        _add_relation!(lh; super=super, sub=sub, unsafe=true)
    end

    return nothing
end

function insert_relation!(
    s::System,
    hname::AbstractString,
    label::HLabel;
    super::HLabel,
    sub::HLabel,
    unsafe::Bool = false
)
    lh = hierarchy(s, hname)

    result = _remove_relation!(lh, super, sub)
    @match result begin
        Relation_Missing => error("relation between super and sub not found. ")
        Success          => nothing
    end
    add_label!(s, hname, label)
    add_relation!(s, hname; super=super, sub=label, unsafe=unsafe)
    add_relation!(s, hname; super=label, sub=sub, unsafe=unsafe)

    return nothing
end

function remove_label!(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)

    if !_contains(lh, label)
        error("label $(label) not found. ")
    end
    _remove_label!(lh, label)

    return nothing
end

function remove_relation!(s::System, hname::AbstractString; super::HLabel, sub::HLabel)
    lh = hierarchy(s, hname)

    if !_has_relation(lh, super, sub)
        error("relation betewwn $(super) and $(sub) not found. ")
    end
    _remove_relation!(lh, super, sub)

    return nothing
end

function label_unique(s::System)
    for hname in hierarchy_names(s)
        lh = hierarchy(s, hname)
        if !_label_unique(lh)
            error("label is not unique in hierarchy $(hname). ")
        end
    end

    return true
end

function label_nocycle(s::System)
    for hname in hierarchy_names(s)
        lh = hierarchy(s, hname)
        if !_label_nocycle(lh)
            error("label has cycle in hierarchy $(hname). ")
        end
    end

    return true
end

function label_connected(s::System)
    for hname in hierarchy_names(s)
        lh = hierarchy(s, hname)
        if !_label_connected(lh)
            error("isolated label found in hierarchy $(hname). ")
        end
    end

    return true
end

function contains(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    return _contains(lh, label)
end

function issuper(s::System, hname::AbstractString, label1::HLabel, label2::HLabel)
    lh = hierarchy(s, hname)
    return _issuper(lh, label1, label2)
end

function issub(s::System, hname::AbstractString, label1::HLabel, label2::HLabel)
    lh = hierarchy(s, hname)
    return _issub(lh, label1, label2)
end

function super(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    return _super(lh, label)
end

function sub(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    return _sub(lh, label)
end

function super_labels(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    labels = HLabel[]
    current = _get_nodeid(lh, label)
    while current != 0 # root node id
        push!(labels, _labels(lh)[current])
        current = _super_id(lh, current)
    end

    return labels
end

super_labels(s::System, hname::AbstractString, atom::Integer) = super_labels(s, hname, HLabel("", atom))

function sub_labels(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    current = _get_nodeid(lh, label)
    stack = [current]

    labels = HLabel[]
    while !isempty(stack)
        current = popfirst!(stack)
        append!(stack, _sub_id(lh, current))
        push!(labels, _labels(lh)[current])
    end

    return labels
end

sub_labels(s::System, hname::AbstractString, atom::Integer) = sub_labels(s, hname, HLabel("", atom))

function label2atom(s::System, hname::AbstractString, label::HLabel)
    if isatom(label)
        return [id(label)]
    end

    lh = hierarchy(s, hname)
    labels = _labels(lh)
    atoms = Int64[]

    stack = [_get_nodeid(lh, label)]
    while !isempty(stack)
        current = popfirst!(stack)
        for node_id in _sub_id(lh, current)
            l = labels[node_id]
            if isatom(l)
                push!(atoms, id(l))
            else
                push!(stack, node_id)
            end
        end
    end

    return atoms
end
