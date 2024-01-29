function test()

@testset "TopologyGraph" begin

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

for _ in 1:1
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
        # dijkstra_shortest_paths(topo, v)
        # dijkstra_shortest_paths(wg, v)
        #@test td == wd
        for fname in fieldnames(typeof(td))
            @test getfield(td, fname) == getfield(wd, fname)
        end
    end
    @testset "bfs_shortestpath_path" begin
        _path = TopologyGraph{Int64, Rational{Int8}}()
        add_vertices!(_path, 10)
        for i in 1:nv(_path)-1
            add_edge!(_path, i, i+1, rand(rng, (1//1, 2//1, 3//1, 3//2)))
        end
        result = bfs_shortestpath(_path, 1, 10)
        oracle = let
            elist = a_star(path_graph(10), 1, 10)
            push!([src(e) for e in elist], dst(elist[end]))
        end
        println(result)
        println(oracle)
        @test result == oracle
    end
    return nothing
    @testset "induced_subgraph" begin
        a, b = let
            a = rand(rng, 1:nv(topo))
            b = rand(rng, 1:nv(topo))
            extrema((a, b))
        end
        st = induced_subgraph(topo, a:b)
        sw = induced_subgraph(wg, a:b)
        return nothing
    end

end

end # testset

end # test
