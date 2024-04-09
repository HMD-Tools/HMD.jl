function test()

randweight(rng) = rand(rng, (1//1, 2//1, 3//1, 3//2))

function _path_graph(n, rng)
    g = TopologyGraph{Int64, Rational{Int8}}()
    add_vertices!(g, n)
    for i in 1:nv(g)-1
        add_edge!(g, i, i+1, randweight(rng))
    end
    return g
end

function _testgraph_shortest_path(rng)
    n = 30
    m = 10
    @assert m < 2n
    @assert m > 5
    # return graph with 2 path from node 1 to n
    # shorter path crosses to longer mainpath like (n+m÷2) -- (n÷2) -- (n+m÷2+1)
    # long main path: 1 -> 2 -> ... -> n
    # shortcut: 1 -> n+1 -> ... -> (n+m÷2) -> (n÷2) -> (n+m÷2+1) -> ... -> n+m -> n
    cross_shortcut = n+m÷2
    cross_mainpath = n ÷ 2
    shortest_path = vcat(
        1,
        n+1:cross_shortcut-1,
        [cross_shortcut, cross_mainpath, cross_shortcut+1],
        cross_shortcut+2:n+m,
        n
    )

    g = _path_graph(n, rng) # main path graph
    # add shortcut path
    add_vertices!(g, m)
    @assert sort(filter(i -> degree(g, i)==0, vertices(g))) == n+1:(n+m)
    # shortcut begin
    add_edge!(g, 1, n+1, randweight(rng))
    for i in n+1:cross_shortcut-1
        @assert add_edge!(g, i, i+1, randweight(rng))
    end
    @assert add_edge!(g, cross_shortcut, cross_mainpath  , randweight(rng))
    @assert add_edge!(g, cross_mainpath, cross_shortcut+1, randweight(rng))
    for i in cross_shortcut+1:(n+m)-1
        @assert add_edge!(g, i, i+1, randweight(rng))
    end
    # end shortcut
    add_edge!(g, n+m, n, randweight(rng))

    return shortest_path, g
end

function random_graph(nv, max_degree, rng)
    topo = TopologyGraph{Int64, Rational{Int8}}()
    wg = SimpleWeightedGraph{Int64, Rational{Int8}}()
    add_vertices!(topo, nv)
    add_vertices!(wg, nv)
    for i in 1:nv
        for _ in 1:rand(rng, 1:max_degree)
            nbr = rand(rng, 1:nv)
            weight = rand(rng, (1//1, 2//1, 3//1, 3//2))
            @assert has_edge(topo, i, nbr) == has_edge(wg, i, nbr)
            @assert get_weight(topo, i, nbr) == get_weight(wg, i, nbr)
            if !has_edge(topo, i, nbr)
                add_edge!(topo, i, nbr, weight)
                add_edge!(wg, i, nbr, weight)
                @assert get_weight(topo, i, nbr) == get_weight(wg, i, nbr) == weight
            end
            @assert has_edge(topo, i, nbr) == has_edge(wg, i, nbr) == true
        end
    end
    return topo, wg
end

function strict_eq(topo, wg)
    nv(topo) != nv(wg) && return false
    for i in 1:nv(topo)-1
        !has_vertex(wg, i) && return false
        for j in i+1:nv(topo)
            has_edge(topo, i, j) != has_edge(wg, i, j) && return false
            get_weight(topo, i, j) != get_weight(wg, i, j) && return false
        end
    end
    return true
end

for _ in 1:5
    topo, wg = random_graph(500, 6, Xoshiro(1))
    @test all(v -> allunique(inneighbors(topo, v)), vertices(topo))
    @test all(v -> issorted(inneighbors(topo, v)), vertices(topo))
    @test strict_eq(topo, wg)

    @test eltype(topo) == eltype(wg)
    @test is_directed(topo) == is_directed(wg)
    @test is_directed(typeof(topo)) == is_directed(typeof(wg))
    @test has_self_loops(topo) == has_self_loops(wg)
    @test num_self_loops(topo) == num_self_loops(wg)
    @test density(topo) == density(wg)
    @test degree(topo) == degree(wg)
    @test indegree(topo) == indegree(wg)
    @test outdegree(topo) == outdegree(wg)
    @test weights(topo) == weights(wg)
    @test Δ(topo) == Δ(wg)
    @test Δin(topo) == Δin(wg)
    @test Δout(topo) == Δout(wg)
    @test δ(topo) == δ(wg)
    @test δin(topo) == δin(wg)
    @test δout(topo) == δout(wg)

    # edge interface
    @test Set(Tuple(e) for e in edges(topo)) == Set(Tuple(e) for e in edges(wg))
    @test edgetype(topo) == TopologyEdge{Int64, Rational{Int8}}
    @test ne(topo) == ne(wg)
    for e in edges(topo)
        _src, _dst, _weight = Tuple(e)
        @test src(e) == _src
        @test dst(e) == _dst
        @test weight(e) == _weight
        @test get_weight(topo, _src, _dst) == _weight
        @test e == e
    end

    # vertex interface
    @test vertices(topo) == vertices(wg)
    @test nv(topo) == nv(wg)
    for v in shuffle(vertices(topo))
        @test issorted(all_neighbors(topo, v))
        @test has_vertex(topo, v) == has_vertex(wg, v)
        @test inneighbors(topo, v) == inneighbors(wg, v)
        @test outneighbors(topo, v) == outneighbors(wg, v)
        @test all_neighbors(topo, v) == all_neighbors(wg, v)
        @test neighbors(topo, v) == neighbors(wg, v)
        @test degree(topo, v) == degree(wg, v)
        @test indegree(topo, v) == indegree(wg, v)
        @test outdegree(topo, v) == outdegree(wg, v)
    end
    @test has_vertex(topo, nv(topo)+1) == has_vertex(wg, nv(wg)+1)

    # double vertex interface
    for v in shuffle(vertices(topo))
        for w in shuffle(vertices(topo))
            @test has_edge(topo, v, w) == has_edge(wg, v, w)
            @test get_weight(topo, v, w) == get_weight(wg, v, w)
            @test common_neighbors(topo, v, w) == common_neighbors(wg, v, w)
        end
    end

    # equality on edge removal
    _topo = deepcopy(topo)
    _wg = deepcopy(wg)
    _edges = (shuffle∘collect∘edges)(topo)
    @test allunique(_edges)
    for e in _edges
        s, d, _weight = Tuple(e)
        @test rem_edge!(_topo, s, d) == rem_edge!(_wg, s, d)
        @test strict_eq(_topo, _wg)
    end

    # equality on vertex removal
    _topo = deepcopy(topo)
    _wg = deepcopy(wg)
    rng = Xoshiro(1)
    while nv(_topo) > 1
        v = rand(rng, vertices(_topo))
        @assert rem_vertex!(_topo, v) == rem_vertex!(_wg, v)
    end

    # hihg level interfaces in SimpleWeightedGraphs
    @testset "connected_components" begin
        #@test connected_components(topo) == connected_components(wg)

    end
    @testset "degree matrix" begin
        @test degree_matrix(topo) == degree_matrix(wg)
    end
    @testset "a_star" begin
        u = rand(rng, 1:nv(topo))
        v = rand(rng, 1:nv(topo))
        @test a_star(topo, u, v) == a_star(wg, u, v)
    end
    @testset "dijkstra" begin
        v = rand(rng, 1:nv(topo))
        td = dijkstra_shortest_paths(topo, v)
        wd = dijkstra_shortest_paths(wg, v)
        for fname in fieldnames(typeof(td))
            @test getfield(td, fname) == getfield(wd, fname)
        end
    end
    @testset "induced_subgraph" begin
        a, b = let
            a = rand(rng, 1:nv(topo))
            b = rand(rng, 1:nv(topo))
            extrema((a, b))
        end
        st, tmap = induced_subgraph(topo, a:b)
        sw, wmap = induced_subgraph(wg, a:b)
        @test tmap == wmap
        @test strict_eq(st, sw)
    end

end # for loop


@testset "TopologyGraph" begin


rng = Xoshiro(1)
@testset "bfs_shortestpath_path" begin
    shortest_path, g = _testgraph_shortest_path(rng)
    @test_throws ErrorException bfs_shortestpath(g, 1, Int64[])
    @test_throws ErrorException bfs_shortestpath(g, 1, [1,4,6,8,105])
    @test_throws ErrorException bfs_shortestpath(g, 1, [20,1,4,6,8])
    @test_throws ErrorException bfs_shortestpath(g, 102, [1,4,6,8,20])

    @testset let x = bfs_shortestpath(g, 1, [1])
        @test getdist(x, 1) == 0
        @test getpath(x, 1) == [1]
    end
    @testset let x = bfs_shortestpath(g, 1, [1, 2])
        @test getdist(x, 1) == 0
        @test getpath(x, 1) == [1]
        @test getdist(x, 2) == 1
        @test getpath(x, 2) == [1,2]
    end
    @testset let x = bfs_shortestpath(g, 2, [1, 2])
        @test getdist(x, 1) == 1
        @test getpath(x, 1) == [2, 1]
        @test getdist(x, 2) == 0
        @test getpath(x, 2) == [2]
    end
    @testset let x = bfs_shortestpath(g, 2, [1, 2, 3])
        @test getdist(x, 1) == 1
        @test getpath(x, 1) == [2, 1]
        @test getdist(x, 2) == 0
        @test getpath(x, 2) == [2]
        @test getdist(x, 3) == 1
        @test getpath(x, 3) == [2, 3]
    end

    result = bfs_shortestpath(g, 1)
    @testset let
        for (i, v) in enumerate(shortest_path)
            @test getpath(result, v) == shortest_path[1:i]
            @test getdist(result, v) == i-1
        end
    end
    result_noncennect = let
        _g = deepcopy(g)
        add_vertices!(_g, 10)
        bfs_shortestpath(_g, 1, vertices(g))
    end
    @test result.distance == result_noncennect.distance
    @test result.predecessor == result_noncennect.predecessor

    sg = let
        sg = SimpleWeightedGraph{Int64, Rational{Int8}}()
        add_vertices!(sg, nv(g))
        for e in edges(g)
            # edge weight is 1//1 b.c. weight (= bond order) is not related to
            # topological distance on molecular graph
            add_edge!(sg, src(e), dst(e), 1//1)
        end
        sg
    end
    for (i, finish) in enumerate(shortest_path)
        i == 1 && continue # a_star returns empty edge list
        elist = a_star(sg, 1, finish)
        path = push!([src(e) for e in elist], dst(elist[end]))
        @test getpath(result, finish) == path
    end
    @inferred bfs_shortestpath(g, Int16(1), Int32[1, 4, 6, 8, 20])
end



@testset "cycle detection" begin
    wg = SimpleWeightedGraph{Int64, Rational{Int8}}(complete_graph(30))
    topo = let
        g = TopologyGraph{Int64, Rational{Int8}}()
        add_vertices!(g, nv(wg))
        for e in edges(wg)
            add_edge!(g, src(e), dst(e), weight(e))
        end
        g
    end
    @assert strict_eq(topo, wg)
    @test is_cyclic(topo) == is_cyclic(wg)
    tc = simplecycles_limited_length(topo, 10)
    wc = simplecycles_limited_length(wg, 10)
    @test length(tc) == length(wc)
    @test tc == wc
end


end # testset

end # test
