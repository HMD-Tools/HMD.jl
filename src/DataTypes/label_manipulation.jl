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

function remove_hierarchy!(s::System, hname::AbstractString)
    delete!(s.hierarchy, hname)
end

function all_labels(s::System, hname::AbstractString)
    lh = hierarchy(s, hname)
    return _labels(lh)
end

function all_labels(s::System, hname::AbstractString, label_type::AbstractString)
    labels = hierarchy(s, hname) |> HierarchyLabels._label2node |> keys |> collect

    return filter!(label -> type(label)==label_type, labels)
end

function add_label!(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)

    result = _add_label!(lh, label)
    @match result begin
        Label_Occupied => error("label $(label) already exists. ")
        Success        => return nothing
        _              => error("fatal error")
    end
end

function add_label!(s::System, hname::AbstractString, label_type::AbstractString)
    lh = hierarchy(s, hname)

    n = count_label(s, hname, label_type)
    addend = HLabel(label_type, n+1)
    result = _add_label!(lh, addend)
    @match result begin
        Label_Occupied => error("label $(label) already exists. ")
        Success        => return addend
        _              => error("fatal error")
    end
end

function add_labels!(s::System, hname::AbstractString, labels::AbstractVector{HLabel})
    lh = hierarchy(s, hname)

    _add_labels!(lh, labels)

    return nothing
end

function count_label(s::System, hname::AbstractString, label_type::String)
    lh = hierarchy(s, hname)
    labels = _label2node(lh) |> keys
    return count(l -> type(l)==label_type, labels)
end

function add_relation!(s::System, hname::AbstractString; super::HLabel, sub::HLabel)
    lh = hierarchy(s, hname)

    result = _add_relation!(lh; super=super, sub=sub)
    @match result begin
        Label_Missing     => error("Super or sub not found. ")
        Label_Duplication => error("super and sub are equal. ")
        Relation_Occupied => error("""There in already relation between super and sub labels. Please use "insert_relation!()" instead. """)
        success           => return nothing
        _                 => error("fatal error")
    end
end

function add_relations!(s::System, hname::AbstractString; super::HLabel, subs::AbstractVector{HLabel})
    lh = hierarchy(s, hname)
    for sub in subs
        _add_relation!(lh; super=super, sub=sub, unsafe=true)
    end

    return nothing
end

function insert_relation!(s::System, hname::AbstractString, label::HLabel; super::HLabel, sub::HLabel)
    lh = hierarchy(s, hname)

    add_label!(s, hname, label)
    add_relation!(s, hname; super=super, sub=label)
    add_relation!(s, hname; super=label, sub=sub)
    _remove_relation!(lh, super, sub)

    return nothing
end

function remove_label!(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)

    if !_contains!(lh, label)
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

function super(s::System, hname::AbstractString, label::HLabel; recurse::Bool=false)
    lh = hierarchy(s, hname)
    return _super(lh, label)
end

function sub(s::System, hname::AbstractString, label::HLabel; recurse::Bool=false)
    lh = hierarchy(s, hname)
    return _sub(lh, label)
end
