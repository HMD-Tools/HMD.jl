using HMD
using Test
using Distributions
using Random

const sample = let
    sample = System{3, Float64}()
    # polymeric hierarchy
    add_hierarchy!(sample, "polymeric")
    add_labels!(
        sample,
        "polymeric",
        [HLabel("sublabel", 1), HLabel("molecule", 1)]
    )
    add_relation!(sample, "polymeric"; super=Entire_System, sub=HLabel("molecule", 1))
    add_relation!(sample, "polymeric"; super=HLabel("molecule", 1), sub=HLabel("sublabel", 1))

    # ruinfo hierarchy
    add_hierarchy!(sample, "ruinfo")
    add_labels!(
        sample,
        "ruinfo",
        [HLabel("ru1", 1), HLabel("molecule", 1)]
    )
    add_relation!(sample, "ruinfo"; super=Entire_System, sub=HLabel("molecule", 1))
    add_relation!(sample, "ruinfo"; super=HLabel("molecule", 1), sub=HLabel("ru1", 1))
    #atoms
    positions = [
        [-4.952, -0.895, 2.069],
        [-4.415, -1.252, 0.690],
        [-4.240, 0.045, -0.079],
        [-3.699, -0.289, -1.475],
        [-3.525, 0.979, -2.235],
        [-2.624, 1.971, -1.630],
        [-1.186, 1.756, -1.463],
        [-0.735, 0.615, -0.650],
        [0.824, 0.637, -0.507],
        [1.137, -0.571, 0.337],
        [2.622, -0.724, 0.591],
        [3.216, 0.440, 1.301],
        [4.683, 0.206, 1.529],
        [5.442, 0.027, 0.244],
        [6.918, -0.186, 0.566],
        [-5.412, -1.782, 2.558],
        [-5.740, -0.104, 1.919],
        [-4.103, -0.478, 2.637],
        [-5.142, -1.863, 0.140],
        [-3.435, -1.752, 0.805],
        [-5.208, 0.568, -0.157],
        [-3.555, 0.708, 0.473],
        [-2.857, -0.955, -1.377],
        [-4.553, -0.855, -2.001],
        [-4.590, 1.403, -2.326],
        [-3.216, 0.733, -3.296],
        [-2.788, 2.932, -2.265],
        [-3.096, 2.343, -0.640],
        [-0.734, 1.691, -2.522],
        [-0.699, 2.719, -1.094],
        [-0.846, -0.317, -1.294],
        [-1.187, 0.454, 0.309],
        [1.267, 0.622, -1.502],
        [1.029, 1.589, 0.024],
        [0.818, -1.466, -0.239],
        [0.551, -0.535, 1.253],
        [3.152, -0.887, -0.369],
        [2.762, -1.630, 1.213],
        [2.766, 0.604, 2.309],
        [3.145, 1.384, 0.698],
        [4.805, -0.743, 2.116],
        [5.147, 0.989, 2.163],
        [5.307, 0.907, -0.386],
        [5.132, -0.921, -0.273],
        [7.441, 0.788, 0.689],
        [6.991, -0.839, 1.450],
        [7.381, -0.702, -0.303]
    ]
    elems = vcat(fill(6, 15), fill(1, 32))
    @assert length(positions) == length(elems)
    add_atoms!(
        sample, positions, elems;
        super = Dict("polymeric" => HLabel("sublabel" , 1), "ruinfo" => HLabel("ru1" , 1))
    )

    #bonds
    add_bond!(sample,  1,  2; bond_order=1//1)
    add_bond!(sample,  1, 16; bond_order=1//1)
    add_bond!(sample,  1, 17; bond_order=1//1)
    add_bond!(sample,  1, 18; bond_order=1//1)
    add_bond!(sample,  2,  3; bond_order=1//1)
    add_bond!(sample,  2, 19; bond_order=1//1)
    add_bond!(sample,  2, 20; bond_order=1//1)
    add_bond!(sample,  3,  4; bond_order=1//1)
    add_bond!(sample,  3, 21; bond_order=1//1)
    add_bond!(sample,  3, 22; bond_order=1//1)
    add_bond!(sample,  4,  5; bond_order=1//1)
    add_bond!(sample,  4, 23; bond_order=1//1)
    add_bond!(sample,  4, 24; bond_order=1//1)
    add_bond!(sample,  5,  6; bond_order=1//1)
    add_bond!(sample,  5, 25; bond_order=1//1)
    add_bond!(sample,  5, 26; bond_order=1//1)
    add_bond!(sample,  6,  7; bond_order=1//1)
    add_bond!(sample,  6, 27; bond_order=1//1)
    add_bond!(sample,  6, 28; bond_order=1//1)
    add_bond!(sample,  7,  8; bond_order=1//1)
    add_bond!(sample,  7, 29; bond_order=1//1)
    add_bond!(sample,  7, 30; bond_order=1//1)
    add_bond!(sample,  8,  9; bond_order=1//1)
    add_bond!(sample,  8, 31; bond_order=1//1)
    add_bond!(sample,  8, 32; bond_order=1//1)
    add_bond!(sample,  9, 10; bond_order=1//1)
    add_bond!(sample,  9, 33; bond_order=1//1)
    add_bond!(sample,  9, 34; bond_order=1//1)
    add_bond!(sample, 10, 11; bond_order=1//1)
    add_bond!(sample, 10, 35; bond_order=1//1)
    add_bond!(sample, 10, 36; bond_order=1//1)
    add_bond!(sample, 11, 12; bond_order=1//1)
    add_bond!(sample, 11, 37; bond_order=1//1)
    add_bond!(sample, 11, 38; bond_order=1//1)
    add_bond!(sample, 12, 13; bond_order=1//1)
    add_bond!(sample, 12, 39; bond_order=1//1)
    add_bond!(sample, 12, 40; bond_order=1//1)
    add_bond!(sample, 13, 14; bond_order=1//1)
    add_bond!(sample, 13, 41; bond_order=1//1)
    add_bond!(sample, 13, 42; bond_order=1//1)
    add_bond!(sample, 14, 15; bond_order=1//1)
    add_bond!(sample, 14, 43; bond_order=1//1)
    add_bond!(sample, 14, 44; bond_order=1//1)
    add_bond!(sample, 15, 45; bond_order=1//1)
    add_bond!(sample, 15, 46; bond_order=1//1)
    add_bond!(sample, 15, 47; bond_order=1//1)
    sample
end


function strict_eq(lh1::LabelHierarchy, lh2::LabelHierarchy)
    if length(lh1.g) != length(lh2.g)
        error("label hierarchy node num mismatch. \n" *
            "    l1: $(length(lh1.g))" *
            "    l2: $(length(lh2.g))"
        )
    end

    ne(g) = HMD.HierarchyLabels._ne(g)
    if ne(lh1.g) != ne(lh2.g)
        error("label hierarchy edge num mismatch. \n" *
            "    l1: $(ne(l1))" *
            "    l2: $(ne(l2))"
        )
    end

    if lh1.labels != lh2.labels
        error("labels in hierarchy mistmatch (hierarchy broken)")
    end

    for p in lh1.label2node
        if p[2] != lh2.label2node[p[1]]
            error("label2node mismatch (hierarchy broken)")
        end
    end

    edges(g) = HMD.HierarchyLabels._edges(g)
    src(e) = HMD.HierarchyLabels._src(e)
    dst(e) = HMD.HierarchyLabels._dst(e)
    for e in edges(lh1.g)
        lh1.labels[src(e)] != lh2.labels[src(e)] && return false
        lh1.labels[dst(e)] != lh2.labels[dst(e)] && return false
    end

    return true
end

function strict_eq(s1::System, s2::System; dynamic_only=false)
    if any(x->abs(x)>1e-5, box(s1).origin - box(s2).origin)
        error("box origin mismatch. \n" *
            "    s1: $(box(s1).origin)" *
            "    s2: $(box(s2).origin)"
        )
    end

    if any(x->abs(x)>1e-5, box(s1).axis - box(s2).axis)
        error("box axis mismatch. \n" *
            "    s1: $(box(s1).axis)" *
            "    s2: $(box(s2).axis)"
        )
    end

    if !isapprox(time(s1), time(s2); atol=1e-5)
        error("time mismatch. \n" *
            "    s1: $(time(s1))" *
            "    s2: $(time(s2))"
        )
    end

    if dynamic_only
        if length(s1.position) != length(s2.position)
            error("natom mismatch. \n" *
                "    s1: $(length(s1.position))" *
                "    s2: $(length(s2.position))"
            )
        end
    elseif natom(s1) != natom(s2)
        error("natom mismatch. \n" *
            "    s1: $(natom(s1))" *
            "    s2: $(natom(s2))"
        )
    end

    for i in 1:length(s1.position)
        if any(x->abs(x)>1e-5, position(s1, i) - position(s2, i))
            error("atom $i position mismatch. \n" *
                "    s1: $(position(s1, i))" *
                "    s2: $(position(s2, i))"
            )
        end
        if travel(s1, i) != travel(s2, i)
            error("atom $i travel mismatch. \n" *
                "    s1: $(travel(s1, i))" *
                "    s2: $(travel(s2, i))"
            )
        end
    end

    for (e1, e2) in zip(all_elements(s1), all_elements(s2))
        if e1 != e2
            error("element mismatch. \n" *
                "    s1: $e1" *
                "    s2: $e2"
            )
        end
    end

    if wrapped(s1) != wrapped(s2)
        error("wrap mismatch. \n" *
            "    s1: $(wrapped(s1))" *
            "    s2: $(wrapped(s2))"
        )
    end

    if hierarchy_names(s1) != hierarchy_names(s2)
        error("hierarchy names mismatch. \n" *
            "    s1: $(hierarchy_names(s1))" *
            "    s2: $(hierarchy_names(s2))"
        )
    end
    for hname in hierarchy_names(s1)
        lh1 = hierarchy(s1, hname)
        lh2 = hierarchy(s2, hname)
        if !strict_eq(lh1, lh2)
            println("\nlh1: ")
            println(lh1)
            println("\nlh2: ")
            println(lh2)
            println()
            error("hierarchy $hname mismatch")
        end
    end

    t1, t2 = topology(s1), topology(s2)
    if nv(t1) != nv(t1)
        error("topology node num mismatch. \n" *
            "    s1: $(nv(t1))" *
            "    s2: $(nv(t2))"
        )
    end
    if ne(t1) != ne(t2)
        error("topology edge num mismatch. \n" *
            "    s1: $(ne(t1))" *
            "    s2: $(ne(t2))"
        )
    end
    for bond in edges(t1)
        head = src(bond)
        tail = dst(bond)
        if get_weight(t1, head, tail) != get_weight(t2, head, tail)
            error("bond order mismatch between $head and $tail. \n" *
                "    s1: $(get_weight(t1, head, tail))" *
                "    s2: $(get_weight(t2, head, tail))"
            )
        end
        for id in (head, tail)
            if any(x->abs(x)>1e-5, position(s1, id) - position(s2, id)) ||
                element(s1, id) != element(s2, id) ||
                travel(s1, id) != travel(s2, id)
                error("topology node id and atom index mismatch. \n" *
                    "    node id: $id" *
                    "    s1: position: $(position(s1, id))" *
                    "        element : $(element(s1, id))" *
                    "        travel  : $(travel(s1, id))" *
                    "    s2: position: $(position(s2, id))" *
                    "        element : $(element(s2, id))" *
                    "        travel  : $(travel(s2, id))"
                )
            end
        end
    end

    return true
end

function strict_eq(t1::Trajectory, t2::Trajectory; dynamic_only=false)
    if t1.reactions != t2.reactions
        error("is_reaction mismatch. \n" *
            "    t1: $(t1.reactions)" *
            "    t2: $(t2.reactions)"
        )
    end

    for i in 1:length(t1)
        strict_eq(t1.systems[i], t2.systems[i]; dynamic_only=dynamic_only)
    end

    return true
end


@testset "HMD.jl" begin
    @testset "TopologyGraphs" begin
        #HMD.TopologyGraphs.test()
    end

    @testset "HierarchyHLabels" begin
        HMD.DataTypes.HierarchyLabels.test()
    end

    @testset "label manipulation" begin
        include("label_manipulation.jl")
    end

    @testset "trajectory" begin
        #include("trajectory.jl")
    end

    @testset "merge system" begin
        s1 = deepcopy(sample)
        s2 = deepcopy(sample)
        merge!(
            s1, s2;
            augend_parents = Dict("polymeric" => Entire_System, "ruinfo" => Entire_System),
            addend_parents = Dict("polymeric" => Entire_System, "ruinfo" => Entire_System),
            unsafe = false
        )

        @test Set(all_labels(s1, "polymeric", "molecule")) == Set(HLabel("molecule", i) for i in 1:2)
        @test Set(all_labels(s1, "ruinfo"   , "molecule")) == Set(HLabel("molecule", i) for i in 1:2)

        @test Set(all_labels(s1, "polymeric", "sublabel")) == Set(HLabel("sublabel", i) for i in 1:2)
        @test Set(all_labels(s1, "ruinfo"   , "ru1"     )) == Set(HLabel("ru1"     , i) for i in 1:2)

        @test Set(sub(s1, "polymeric", Entire_System)) == Set([HLabel("molecule", i) for i in 1:2])
        @test Set(sub(s1, "ruinfo"   , Entire_System)) == Set([HLabel("molecule", i) for i in 1:2])

        @test sub(s1, "polymeric", HLabel("molecule", 1)) == [HLabel("sublabel", 1)]
        @test sub(s1, "polymeric", HLabel("molecule", 2)) == [HLabel("sublabel", 2)]

        @test sub(s1, "ruinfo", HLabel("molecule", 1)) == [HLabel("ru1", 1)]
        @test sub(s1, "ruinfo", HLabel("molecule", 2)) == [HLabel("ru1", 2)]

        @test sort(label2atom(s1, "polymeric", HLabel("molecule", 1))) == [i for i in 1:47]
        @test sort(label2atom(s1, "ruinfo"   , HLabel("molecule", 1))) == [i for i in 1:47]
        @test sort(label2atom(s1, "polymeric", HLabel("molecule", 2))) == [i for i in 48:94]
        @test sort(label2atom(s1, "ruinfo"   , HLabel("molecule", 2))) == [i for i in 48:94]
    end


    @testset "atom addition" begin
        @testset "add_label()-add_labels() equality" begin
            s1, s2 = deepcopy(sample), deepcopy(sample)
            pos = [rand(3) for _ in 1:100]
            elm = rand(DiscreteUniform(1, 90), 100)
            for i in 1:100
                add_atom!(
                    s1, pos[i], elm[i];
                    super = Dict("polymeric" => HLabel("sublabel" , 1), "ruinfo" => HLabel("ru1" , 1))
                )
            end
            add_atoms!(
                s2, pos, elm;
                super = Dict("polymeric" => HLabel("sublabel" , 1), "ruinfo" => HLabel("ru1" , 1))
            )
            strict_eq(s1, s2)
        end
    end


end
