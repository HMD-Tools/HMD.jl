Base.@kwdef struct TopologyGraph{T<:Integer, U<:Real} <: AbstractSimpleWeightedGraph{T, U}
    adjlist::Vector{Vector{T}} = Vector{Vector{T}}()
    weights::Vector{Vector{U}} = Vector{Vector{U}}()
end

function TopologyGraph{T, U}(g::AbstractGraph{T}) where {T, U}
    tg = TopologyGraph{T, U}()
    add_vertices!(tg, nv(g))
    for e in edges(g)
        add_edge!(tg, src(e), dst(e), one(U))
    end

    return tg
end

function empty!(g::TopologyGraph{T, U}) where {T<:Integer, U<:Real}
    empty!(g.adjlist)
    empty!(g.weights)
end

function ==(g1::TopologyGraph{T, U}, g2::TopologyGraph{T, U}) where {T<:Integer, U<:Real}
    return all(eachindex(g1.adjlist)) do i
        g1.adjlist[i] == g2.adjlist[i] && g1.weights[i] == g2.weights[i]
    end
end

Base.zero(g::TopologyGraph{T, U}) where {T<:Integer, U<:Real} = TopologyGraph{T, U}()
Base.zero(::Type{TopologyGraph{T, U}}) where {T<:Integer, U<:Real} = TopologyGraph{T, U}()

function SimpleWeightedGraphs.edges(g::TopologyGraph{T, U}) where {T<:Integer, U<:Real}
    return (
        TopologyEdge{T,U}(v, g.adjlist[v][i], g.weights[v][i])
        for v in eachindex(g.adjlist) for i in eachindex(g.adjlist[v])
        if v ≤ g.adjlist[v][i]
    )
end

function Base.getindex(
    g::TopologyGraph{T, U},
    e::Tuple{T, T, U}
) where {T<:Integer, U<:Real}
    return get_weight(g, e[1], e[2])
end

function Base.getindex(
    g::TopologyGraph{T, U},
    e::Tuple{TopologyEdge{T, U}, U}
) where {T<:Integer, U<:Real}
    return get_weight(g, e[1][1], e[1][2])
end

SimpleWeightedGraphs.edgetype(g::TopologyGraph{T, U}) where {T<:Integer, U<:Real} = TopologyEdge{T, U}

function SimpleWeightedGraphs.has_edge(
    g::TopologyGraph{T, U},
    v::Integer,
    w::Integer
) where {T<:Integer, U<:Real}
    #if isempty(g.adjlist)
    #    return false
    #elseif has_vertex(g, v) && has_vertex(g, w)
    #    isempty(g.adjlist[v]) && return false
    #    hasedge = w ∈ g.adjlist[v]
    #    @assert hasedge == (v ∈ g.adjlist[w])
    #    return hasedge
    #else
    #    return false
    #end
    if has_vertex(g, v) && has_vertex(g, w)
        hasedge = w ∈ g.adjlist[v]
        @assert hasedge == (v ∈ g.adjlist[w])
        return hasedge
    else
        return false
    end
end

function SimpleWeightedGraphs.has_vertex(
    g::TopologyGraph{T, U},
    v::Integer
) where {T<:Integer, U<:Real}
    return v ≤ length(g.adjlist)
end

function SimpleWeightedGraphs.inneighbors(
    g::TopologyGraph{T, U},
    v::Integer
) where {T<:Integer, U<:Real}
    return _neighbors(g, v)
end

function SimpleWeightedGraphs.is_directed(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return false
end

function SimpleWeightedGraphs.is_directed(
    ::Type{TopologyGraph{T, U}}
) where {T<:Integer, U<:Real}
    return false
end

function SimpleWeightedGraphs.ne(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    n = zero(T)
    for v in vertices(g)
        nbr = _neighbors(g, v)
        n += length(nbr) + ifelse(v ∈ nbr, 1, 0)
    end
    @assert iseven(n)

    return n ÷ 2
end

function SimpleWeightedGraphs.nv(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return length(g.adjlist)
end

function SimpleWeightedGraphs.outneighbors(
    g::TopologyGraph{T, U},
    v::Integer
) where {T<:Integer, U<:Real}
    return _neighbors(g, v)
end

function SimpleWeightedGraphs.vertices(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return eachindex(g.adjlist)
end

function SimpleWeightedGraphs.add_vertex!(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return add_vertices!(g, 1)
end

function SimpleWeightedGraphs.add_vertices!(
    g::TopologyGraph{T, U},
    n::Integer
) where {T<:Integer, U<:Real}
    if length(g.adjlist) > typemax(T) - n
        return false
    end

    append!(g.adjlist, (T[] for _ in 1:n))
    append!(g.weights, (U[] for _ in 1:n))
    @assert length(g.adjlist) == length(g.weights)

    return true
end

function SimpleWeightedGraphs.add_edge!(
    g::TopologyGraph{T, U},
    u::Integer,
    v::Integer,
    weight::U
) where {T<:Integer, U<:Real}
    if !has_vertex(g, u) || !has_vertex(g, v) || has_edge(g, u, v)
        return false
    end
    ulist = g.adjlist[u]
    vlist = g.adjlist[v]
    uweight = g.weights[u]
    vweight = g.weights[v]

    @assert (v ∈ ulist) == (u ∈ vlist)
    @assert length(ulist) == length(uweight)
    @assert length(vlist) == length(vweight)

    iv = searchsortedfirst(vlist, u)
    insert!(vlist, iv, u)
    insert!(vweight, iv, weight)
    @assert issorted(vlist)
    @assert allunique(vlist)

    if u != v
        iu = searchsortedfirst(ulist, v)
        insert!(ulist, iu, v)
        insert!(uweight, iu, weight)
        @assert issorted(ulist)
        @assert allunique(vlist)
    end

    return true
end

function SimpleWeightedGraphs.add_edge!(
    g::TopologyGraph{T, U},
    e::TopologyEdge{T, U}
) where {T<:Integer, U<:Real}
    v, w, _weight = src(e), dst(e), weight(e)
    return add_edge!(g, v, w, _weight)
end

function SimpleWeightedGraphs.add_edge!(
    g::TopologyGraph{T, U},
    u::Integer,
    v::Integer
) where {T<:Integer, U<:Real}
    return add_edge!(g, u, v, one(U))
end

function SimpleWeightedGraphs.rem_edge!(
    g::TopologyGraph{T, U},
    u::Integer,
    v::Integer,
) where {T<:Integer, U<:Real}
    if !has_vertex(g, u) || !has_vertex(g, v) || !has_edge(g, u, v)
        return false
    end

    ulist = g.adjlist[u]
    vlist = g.adjlist[v]
    uweight = g.weights[u]
    vweight = g.weights[v]
    @assert (v ∈ ulist) == (u ∈ vlist)
    @assert length(ulist) == length(uweight)
    @assert length(vlist) == length(vweight)

    iu = findfirst(==(v), ulist)
    deleteat!(ulist, iu)
    deleteat!(uweight, iu)
    @assert v ∉ ulist
    @assert allunique(ulist)
    @assert issorted(ulist)

    if u != v
        iv = findfirst(==(u), vlist)
        deleteat!(vlist, iv)
        deleteat!(vweight, iv)
    end
    @assert u ∉ vlist
    @assert allunique(vlist)
    @assert issorted(vlist)

    return true
end

function SimpleWeightedGraphs.rem_edge!(
    g::TopologyGraph{T, U},
    e::TopologyEdge{T, U}
) where {T<:Integer, U<:Real}
    return rem_edge!(g, src(e), dst(e))
end

function SimpleWeightedGraphs.rem_vertex!(
    g::TopologyGraph{T, U},
    v::Integer,
) where {T<:Integer, U<:Real}
    if !has_vertex(g, v)
        return false
    end

    # remove edges containing v
    nbrs = deepcopy(_neighbors(g, v))
    for nbr in nbrs
        @assert rem_edge!(g, v, nbr)
    end
    @assert isempty(_neighbors(g, v))

    # the last vertex is moved to the position of v
    lastvertex = nv(g)
    if v == lastvertex
        pop!(g.adjlist)
        pop!(g.weights)
        return true
    end

    # if lastvertex has self loop, _neighbors(g, lastvertex) is changed during for loop
    nbrs = deepcopy(_neighbors(g, lastvertex))
    for nbr in nbrs
        replace!(g.adjlist[nbr], lastvertex => v)
        sort!(g.adjlist[nbr]; alg=InsertionSort)
        @assert allunique(g.adjlist[nbr])
    end
    g.adjlist[v] = pop!(g.adjlist)
    g.weights[v] = pop!(g.weights)

    return true
end

#function SimpleWeightedGraphs.all_neighbors(
#    g::TopologyGraph{T, U},
#    v::Integer
#) where {T<:Integer, U<:Real}
#    return _neighbors(g, v)
#end

#function SimpleWeightedGraphs.common_neighbors(
#    g::TopologyGraph{T, U},
#    u::Integer,
#    v::Integer
#) where {T<:Integer, U<:Real}
#    nbr = _neighbors(g, v)
#    return [w for w in _neighbors(g, u) if w ∈ nbr]
#end

function SimpleWeightedGraphs.degree(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return [degree(g, v) for v in vertices(g)]
end

function SimpleWeightedGraphs.degree(
    g::TopologyGraph{T, U},
    v::Integer
) where {T<:Integer, U<:Real}
    #return count(!=(v), _neighbors(g, v))
    return length(_neighbors(g, v))
end

#function SimpleWeightedGraphs.density(
#    g::TopologyGraph{T, U}
#) where {T<:Integer, U<:Real}
#    n = nv(g) * (nv(g) - 1) ÷ 2
#    return ne(g) / n
#end

function SimpleWeightedGraphs.has_self_loops(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return any(vertices(g)) do v
        v ∈ _neighbors(g, v)
    end
end

function SimpleWeightedGraphs.indegree(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return degree(g)
end

function SimpleWeightedGraphs.indegree(
    g::TopologyGraph{T, U},
    v::Integer
) where {T<:Integer, U<:Real}
    return degree(g, v)
end

function _neighbors(
    g::TopologyGraph{T, U},
    v::Integer
) where {T<:Integer, U<:Real}
    return g.adjlist[v]
end

#noallocextreme(f, comparison, initial, g)

function SimpleWeightedGraphs.num_self_loops(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return count(vertices(g)) do v
        v ∈ _neighbors(g, v)
    end
end

function SimpleWeightedGraphs.outdegree(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    return degree(g)
end

function SimpleWeightedGraphs.outdegree(
    g::TopologyGraph{T, U},
    v::Integer
) where {T<:Integer, U<:Real}
    return degree(g, v)
end

function SimpleWeightedGraphs.weights(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    mat = spzeros(U, nv(g), nv(g))
    for u in vertices(g)
        for v in _neighbors(g, u)
            mat[u, v] = get_weight(g, u, v)
        end
    end
    return mat
end

function SimpleWeightedGraphs.get_weight(
    g::TopologyGraph{T, U},
    u::Integer,
    v::Integer
) where {T<:Integer, U<:Real}
    # if u or v is out of range, throw error
    if v > length(g.adjlist)
        throw(BoundsError(g.adjlist, v))
    end

    nbrs = _neighbors(g, u)
    for i in eachindex(nbrs)
        if nbrs[i] == v
            return g.weights[u][i]
        end
    end

    return zero(U)
end

function SimpleWeightedGraphs.induced_subgraph(
    g::TopologyGraph{T, U},
    vlist::AbstractVector{V}
) where {T<:Integer, U<:Real, V<:Integer}
    if !allunique(vlist)
        throw(ArgumentError("vlist must be unique. "))
    end

    node_mapping = sort(vlist)
    if !has_vertex(g, node_mapping[end])
        error("at least one vertex $v out of range. ")
    end

    elist = Vector{TopologyEdge{T, U}}()
    for v in vlist
        for nbr in _neighbors(g, v)
            !insorted(nbr, node_mapping) && continue
            v ≤ nbr && push!(elist, TopologyEdge{T, U}(v, nbr, get_weight(g, v, nbr)))
        end
    end
    sort!(
        elist;
        lt = (e1, e2) -> if src(e1) != src(e2)
            src(e1) < src(e2)
        else
            dst(e1) < dst(e2)
        end
    )
    @assert issorted(elist; by=src)
    @assert allunique(elist)

    return _induced_subgraph(elist, node_mapping)
end

function SimpleWeightedGraphs.induced_subgraph(
    g::TopologyGraph{T1, U1},
    elist::AbstractVector{TopologyEdge{T2, U2}}
) where {T1<:Integer, U1<:Real, T2<:Integer, U2<:Real}
    if !allunique(elist)
        throw(ArgumentError("elist must be unique. "))
    end

    _elist = sort(elist; by=src)
    sort!(_elist; by=dst)

    node_mapping = T1[]
    for e in _elist
        if !has_vertex(g, src(e)) || !has_vertex(g, dst(e))
            error("at least one vertex $(src(e)) or $(dst(e)) out of range. ")
        end
        i = searchsortedfirst(node_mapping, src(e))
        node_mapping[i] != src(e) && insert!(node_mapping, i, src(e))
        i = searchsortedfirst(node_mapping, dst(e))
        node_mapping[i] != dst(e) && insert!(node_mapping, i, dst(e))
    end
    @assert issorted(node_mapping)
    @assert allunique(node_mapping)

    return _induced_subgraph(_elist, node_mapping)
end

function _induced_subgraph(
    elist::AbstractVector{TopologyEdge{T, U}},
    node_mapping::AbstractVector{V}
) where {T<:Integer, U<:Real, V<:Integer}
    sg = TopologyGraph{T, U}()
    add_vertices!(sg, length(node_mapping))
    for e in elist
        s = searchsortedfirst(node_mapping, src(e))
        d = searchsortedfirst(node_mapping, dst(e))
        @assert add_edge!(sg, s, d, weight(e))
        @assert node_mapping[s] == src(e)
        @assert node_mapping[d] == dst(e)
    end

    return sg, node_mapping
end

function SimpleWeightedGraphs.degree_matrix(
    g::TopologyGraph{T, U},
    UU::Type{<:Real} = U;
    dir::Symbol = :out
) where {T<:Integer, U<:Real}
    if dir ∉ (:in, :out, :both)
        error("dir must be either :out or :in or :both. ")
    end

    mat = spzeros(U, nv(g), nv(g))
    for v in vertices(g)
        mat[v, v] = UU(sum(get_weight(g, v, nbr) for nbr in _neighbors(g, v)))
    end

    return mat
end

function SimpleWeightedGraphs.adjacency_matrix(
    g::TopologyGraph{T, U},
    TT::Type{<:Integer} = T;
    dir::Symbol = :out
) where {T<:Integer, U<:Real}
    if dir ∉ (:in, :out, :both)
        error("dir must be either :out or :in or :both. ")
    end

    return TT.(weights(g))
end

function SimpleWeightedGraphs.connected_components(
    g::TopologyGraph{T, U}
) where {T<:Integer, U<:Real}
    ccs = Vector{Vector{T}}()
    unvisited = Set{T}(vertices(g))
    while !isempty(unvisited)
        # BFS since graph is sparse
        buffer = [first(unvisited)]
        cc = Set{T}(buffer[1])
        while !isempty(buffer)
            current = popfirst!(buffer)
            for nbr in _neighbors(g, current)
                nbr ∈ cc && continue
                push!(cc, nbr)
                push!(buffer, nbr)
            end
            delete!(unvisited, current)
        end
        push!(ccs, (sort!∘collect)(cc))
    end

    return ccs
end

function _bfs_shortestpath(
    g::TopologyGraph{T, U},
    start::Integer,
    finish::AbstractVector{T1},
    area::AbstractVector{T2} = vertices(g)
) where {T<:Integer, U<:Real, T1 <:Integer, T2 <:Integer} # -> path::Vector{Vector{T}}
    if isempty(finish)
        error("empty target atoms ")
    elseif isempty(area)
        error("empty area. ")
    elseif !has_vertex(g, start)
        error("start node $start is out of range. ")
    elseif !_unique_and_sorted(area)
        error("area must be unique and sorted. ")
    elseif !_unique_and_sorted(finish)
        error("target atoms for shortest path must be unique and sorted. ")
    elseif !_issubset_sorted(area, vertices(g))
        error("area must be a subset of vertices.")
    elseif !_issubset_sorted(finish, area)
        error("target atoms for shortest path must be a subset of area.")
    elseif (length(finish)==1 && start == finish[1]) || nv(g) == 1
        return [[T(start)]]
    end

    path = Vector{Vector{T}}()

    buffer = T[start]
    visited = Set{T}(start)
    elist = Vector{Tuple{T, T}}()
    while !isempty(buffer) && length(path) ≤ length(finish)
        v = popfirst!(buffer)
        for nbr in _neighbors(g, v)
            if nbr ∉ visited
                push!(buffer, nbr)
                push!(elist, (v, nbr))
                # visited is extended here, so edges in elist never merges
                push!(visited, nbr)
                insorted(nbr, finish) && push!(path, _path(elist, start, nbr))
            end
        end
    end

    # if start ∈ finish, the path [start] is not added in the BFS search above
    i = searchsortedfirst(finish, start)
    if finish[i] == start # equal to (start ∈ finish)
        #insert!(path, i, [start])
        pushfirst!(path, [start])
    end

    return path
end

# elist backtracking
function _path(elist::AbstractVector{Tuple{T,T}}, start::Integer, finish::Integer) where {T<:Integer}
    @assert elist[end][2] == finish

    # define support function to avoid runtime dispatch caused by captured variable
    findprev_equal(equalto, elist, i) = findprev(e -> e[2] == equalto, elist, i)

    path = T[elist[end][2]]
    current = lastindex(elist)
    head, tail = elist[current]
    while head != start
        @assert head != tail
        pushfirst!(path, head)
        # 呼び出し元のBFSに伴ってedgeが分岐することはあっても合流することはないのでcurrentは一意に定まる
        current = findprev_equal(head, elist, current)
        head, tail = elist[current]
    end

    return pushfirst!(path, head)
end


struct BFShotestState{T<:Integer}
    start::T
    distance::Dict{T, T} # distance[v] == distance from `start` to v
    predecessor::Dict{T, T} # predecessor[v] == predecessor of v in BFS
end

function getdist(state::BFShotestState{T}, v::Integer) where {T<:Integer}
    if !haskey(state.distance, v)
        error("vertex $v is not in the shortest path. ")
    end

    return state.distance[v]
end

function getpath(state::BFShotestState{T}, v::Integer) where {T<:Integer}
    if !haskey(state.predecessor, v)
        error("vertex $v is not in the shortest path. ")
    end

    path = Vector{T}(undef, state.distance[v]+1)
    path[end] = v
    for i in reverse(1:length(path)-1)
        path[i] = state.predecessor[path[i+1]]
    end
    @assert path[1] == state.start

    return path
end

function farthest(state::BFShotestState{T}) where {T<:Integer}
    return argmax(state.distance)
end

function bfs_shortestpath(
    g::TopologyGraph{T, U},
    start::Integer,
    area::AbstractVector{T1} = vertices(g)
) where {T<:Integer, U<:Real, T1 <:Integer}
    if isempty(area)
        error("empty area. ")
    elseif !has_vertex(g, start)
        error("start node $start is out of range. ")
    elseif !_unique_and_sorted(area)
        error("area must be unique and sorted. ")
    elseif !_issubset_sorted(area, vertices(g))
        error("area must be a subset of vertices.")
    elseif !insorted(start, area)
        error("start vertex must be in area. ")
    elseif length(area) == 1 || nv(g) == 1
        return BFShotestState{T}(
            start,
            Dict{T, T}(start => 0),
            Dict{T, T}(start => start)
        )
    end

    # dists[v] == distance from `start` to each v in `area`
    dists = Dict{T, T}(start => 0)

    # predecessor[v] == predecessor of each v in `area` in BFS
    predecessor = Dict{T, T}(start => start)

    buffer = T[start]
    while !isempty(buffer) && length(predecessor) < length(area)
        v = popfirst!(buffer)
        next_depth = dists[v] + 1

        for nbr in _neighbors(g, v)
            if nbr ∉ keys(predecessor) && insorted(nbr, area)
                push!(buffer, nbr)
                dists[nbr] = next_depth
                predecessor[nbr] = v
            end
        end
    end

    return BFShotestState{T}(start, dists, predecessor)
end


function neighborhood(
    g::TopologyGraph{T, U},
    start::Integer,
    d::Integer
) where {T<:Integer, U<:Real}
    nbr = neighborhood_dists(g, start, d)

    return [node for (node, dist) in nbr]
end

function neighborhood_dists(
    g::TopologyGraph{T, U},
    start::Integer,
    d::Integer
) where {T<:Integer, U<:Real}
    if !has_vertex(g, start)
        error("start node $start is out of range. ")
    elseif d < 0
        error("d must be non-negative. ")
    elseif d == 0
        return [(start, 0)]
    end

    visited = pushfirst!(
        [(nbr, 1) for nbr in _neighbors(g, start)],
        (start, 0)
    )
    d == 1 && return visited

    for depth in 2:d
        @assert length(visited) == depth-1
        ith_nbr = outer_neighbors(g, visited)
        @assert allunique(ith_nbr)
        append!(visited, ith_nbr)
    end
    @assert visited[end][2] == d

    return visited
end

function outer_neighbors(
    g::TopologyGraph{T, U},
    visited::AbstractVector{Tuple{T1, T1}}
) where {T<:Integer, U<:Real, T1<:Integer}
    if length(visited) < 2
        error("given visited area must have ≥2 length, found $(length(visited)). ")
    end

    outer = Vector{Tuple{T, T}}()
    maxdepth = visited[end][2]
    # visited node at maxdepth and maxdepth-1
    surface = view(
        visited,
        findlast(p -> p[2] == maxdepth-2, visited) + 1 : lastindex(visited) # equalto findall(p[2] > maxdepth-2)
    )
    for (node, depth) in surface
        for nbr in _neighbors(g, node)
            if all(nbr != p[1] for p in surface) && all(nbr != p[1] for p in outer)
                push!(outer, (nbr, maxdepth+1))
            end
        end
    end

    return outer
end



function _insertsorted!(vec::AbstractVector{T}, x::U) where {T<:Real, U<:Real}
    i = searchsortedfirst(vec, x)
    insert!(vec, i, x)

    return vec
end

function _deletesorted!(vec::AbstractVector{T}, x::U) where {T<:Real, U<:Real}
    i = searchsortedfirst(vec, x)
    if i > length(vec)
        return nothing
    elseif vec[i] == x
        deleteat!(vec, i)
    end

    return nothing
end

function _setdiffsorted!(vec::AbstractVector{T}, rem::AbstractVector{U}) where {T<:Real, U<:Real}
    for x in rem
        i = searchsortedfirst(vec, x)
        x != vec[i] && continue
        deleteat!(vec, i)
    end

    return nothing
end

function _unique_and_sorted(x::AbstractVector{<:Real})
    result = true
    for i in eachindex(x)[2:end]
        @inbounds result &= (x[i] > x[i-1])
    end
    return result
end

function _issubset_sorted(x::AbstractVector{T}, y::AbstractVector{U}) where {T<:Real, U<:Real}
    result = true
    for e in x
        result &= insorted(e, y)
    end
    return result
end
