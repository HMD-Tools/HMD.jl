using Test

function _random_tree(ndepth, maxwidth, rng)
    Entire_System = HLabel("Entire_System", 0)

    lh = LabelHierarchy()
    _add_label!(lh, Entire_System)

    last_labels = [Entire_System]
    for i in 2:ndepth
        next_labels = shuffle!(
            [HLabel(randstring(rng, rand(1:20)), i) for i in 1:rand(rng, 1:maxwidth)]
        )
        _add_labels!(lh, next_labels)
        for label in next_labels
            _add_relation!(
                lh;
                super = rand(rng, last_labels),
                sub = label
            )
        end
        last_labels = next_labels
    end

    return lh
end

function test()

@testset "HierarchyHLabels" begin
    l1, l2, l3 = HLabel("l1", 1), HLabel("l2", 2), HLabel("l3", 3)
    noexist    = HLabel("not_exist_in_h", -1)
    handmade = let
        lh = LabelHierarchy()
        @assert _add_nodes!(lh.g, 1); push!(lh.labels, l1); push!(lh.label2node, l1 => 1)
        @assert _add_nodes!(lh.g, 1); push!(lh.labels, l2); push!(lh.label2node, l2 => 2)
        @assert _add_nodes!(lh.g, 1); push!(lh.labels, l3); push!(lh.label2node, l3 => 3)
        _add_edge!(lh.g, 1, 2)
        _add_edge!(lh.g, 2, 3)
        lh
    end

    # check the order of add and insert matters
    addrel1 = let
        lh = LabelHierarchy()
        _add_label!(lh, l1, false)
        _add_label!(lh, l2, false)
        _add_label!(lh, l3, false)
        _add_relation!(lh; super = l1, sub = l2, unsafe=false)
        _add_relation!(lh; super = l2, sub = l3, unsafe=false)
        lh
    end
    addrel2 = let
        lh = LabelHierarchy()
        @test _add_label!(lh, l2, false) == Success
        @test _add_label!(lh, l1, false) == Success
        @test _add_label!(lh, l3, false) == Success
        @test _add_relation!(lh; super = l2, sub = l3, unsafe=false) == Success
        @test _add_relation!(lh; super = l1, sub = l2, unsafe=false) == Success
        lh
    end
    multiadd = let
        lh = LabelHierarchy()
        @test _add_labels!(lh, [l1, l2, l3]) == Success
        @test _add_relation!(lh; super = l2, sub = l3, unsafe=false) == Success
        @test _add_relation!(lh; super = l1, sub = l2, unsafe=false) == Success
        @test _label_unique(lh)
        @test _label_connected(lh)
        lh
    end
    @test addrel1 == addrel2 == handmade == multiadd

    @test _contains(addrel1, l1)
    @test !_contains(addrel1, noexist)
    @test _get_nodeid(addrel1, l1) == 1 && _get_nodeid(addrel1, l2) == 2 && _get_nodeid(addrel1, l3) == 3
    @test _super(addrel1, l3) == l2 && _super(addrel1, l2) == l1
    @test _sub(addrel1, l1)   == [l2] && _sub(addrel1, l2)   == [l3]
    @test _issuper(addrel1, l1, l2) && _issuper(addrel1, l2, l3)
    @test _issub(addrel1, l2, l1)   && _issub(addrel1, l3, l2)

    @inferred _contains(addrel1, l1)
    @inferred _get_nodeid(addrel1, l1)
    @inferred _sub(addrel1, l1)
    @inferred _issuper(addrel1, l1, l2)
    @inferred _issub(addrel1, l2, l1)

    @test_throws KeyError _get_nodeid(addrel1, noexist)
    @test_throws ErrorException _has_relation(addrel1, noexist, l2)
    @test_throws ErrorException _get_label(LabelHierarchy(), 0)
    @test_throws BoundsError _get_label(addrel1, -1)

    @test _add_label!(addrel1, l1, false) == Label_Occupied
    @test _add_relation!(addrel1; super=noexist, sub=l3, unsafe=false) == Label_Missing
    @test _add_relation!(addrel1; super=l3, sub=l3, unsafe=false)      == Label_Duplication
    @test _add_relation!(addrel1; super=l2, sub=l3, unsafe=false)      == Relation_Occupied
    @test _add_relation!(addrel1; super=l2, sub=l3, unsafe=false)      == Relation_Occupied
    @test addrel1 == multiadd

    # add_labels!() duplication detection
    let
        lh = LabelHierarchy()
        @test_throws ErrorException _add_labels!(lh, [l1, l2, l2])
    end

    # cycle detection in LabelHierarchy graph
    let
        lh = deepcopy(addrel1)
        @test _add_relation!(lh; super = l3, sub = l1, unsafe=false) == Cycle_Found
        @test lh == addrel1
    end

    # unsafe label construction
    let
        lh = LabelHierarchy()
        @test _add_labels!(lh, [l1, l2, l3]) == Success
        @test _add_relation!(lh; super = l2, sub = l3, unsafe=true) == Success
        @test _add_relation!(lh; super = l1, sub = l2, unsafe=true) == Success
        # make cycle
        cycled = deepcopy(lh)
        @test _add_relation!(cycled; super = l3, sub = l1, unsafe=true) == Success
        @test _label_nocycle(cycled) == false
        # add isolated label
        iso = deepcopy(lh)
        @test _label_connected(iso) == true
        @test _add_label!(iso, HLabel("isolate", 1), true) == Success
        @test _label_connected(iso) == false
    end

    addrel1 = let
        lh = LabelHierarchy()
        _add_label!(lh, l1, false)
        _add_label!(lh, l2, false)
        _add_label!(lh, l3, false)
        _add_relation!(lh; super = l1, sub = l2, unsafe=false)
        _add_relation!(lh; super = l2, sub = l3, unsafe=false)
        lh
    end
    lh = deepcopy(addrel1)
    @test lh == addrel1
    @test _remove_relation!(lh, noexist, l3) == Relation_Missing
    @test _remove_relation!(lh, l2, noexist) == Relation_Missing
    @test lh == addrel1
    @test _remove_label!(lh, l1) == Success
    @test lh != addrel1
    @test _remove_label!(lh, l1) == Label_Missing
    @test _remove_relation!(lh, l2, l3) == Success
    @test _remove_relation!(lh, l2, l1) == Relation_Missing

    @testset "Hierarchy Merge" begin
        # TODO: add atom label to test samples
        # hierarchy merge test
        # small piece
        lh1 = let
            lh = LabelHierarchy()
            _add_labels!(
                lh,
                [
                    HLabel("root", 1),
                    HLabel("parent", 1), HLabel("parent", 2), # parent layer
                    HLabel("child", 1), HLabel("child", 2),
                    HLabel("child", 4), HLabel("child", 5), # skip 3
                    HLabel("child", 6), HLabel("child", 10) #jump to 10
                ]
            )
            _add_relation!(lh; super = HLabel("root", 1), sub = HLabel("parent", 1), unsafe=false)
            _add_relation!(lh; super = HLabel("root", 1), sub = HLabel("parent", 2), unsafe=false)
            # under parent 1
            for i in (1,2,4)
                _add_relation!(lh; super = HLabel("parent", 1), sub = HLabel("child", i), unsafe=false)
            end
            # under parent 2
            for i in (5,6)
                _add_relation!(lh; super = HLabel("parent", 2), sub = HLabel("child", i), unsafe=false)
            end
            # under child 6
            _add_relation!(lh; super = HLabel("child", 6), sub = HLabel("child", 10), unsafe=false)
            lh
        end
        lh2 = let
            # replacing id 1 => -1, 10 => 9 from lh1
            lh = deepcopy(lh1)
            _replace!(lh, HLabel("child", 1 ) => HLabel("child", -1))
            _replace!(lh, HLabel("child", 10) => HLabel("child", 9 ))
            _add_label!(lh, HLabel("child", 11))
            _add_relation!(lh; super = HLabel("child", 6), sub = HLabel("child", 11), unsafe=false)
            lh
        end

        # test oracle
        oracle_p1p1 = let
            # oracle for
            #   _merge_hierarchy!(
            #       lh1,
            #       lh2;
            #       augend_parent = HLabel("parent", 1),
            #       addend_parent = HLabel("parent", 1)
            #   )
            lh = deepcopy(lh1)
            for i in (-1, 11, 12)
                @assert _add_label!(lh, HLabel("child", i)) == Success
                @assert _add_relation!(lh; super = HLabel("parent", 1), sub = HLabel("child", i), unsafe=false) == Success
            end
            @test _label_connected(lh)
            @test _label_nocycle(lh)
            @test _label_unique(lh)
            lh
        end
        oracle_p2p1 = let
            # oracle for
            #   _merge_hierarchy!(
            #       lh1,
            #       lh2;
            #       augend_parent = HLabel("parent", 2),
            #       addend_parent = HLabel("parent", 1)
            #   )
            lh = deepcopy(lh1)
            for i in (-1, 11, 12)
                @assert _add_label!(lh, HLabel("child", i)) == Success
                @assert _add_relation!(lh; super = HLabel("parent", 2), sub = HLabel("child", i), unsafe=false) == Success
            end
            @test _label_connected(lh)
            @test _label_nocycle(lh)
            @test _label_unique(lh)
            lh
        end
        oracle_p2p2 = let
            # oracle for
            #   _merge_hierarchy!(
            #       lh1,
            #       lh2;
            #       augend_parent = HLabel("parent", 2),
            #       addend_parent = HLabel("parent", 2)
            #   )
            lh = deepcopy(lh1)
            for i in (11, 12, 9, 13)
                @assert _add_label!(lh, HLabel("child", i)) == Success
                if i ∈ (11, 12)
                    @assert _add_relation!(
                        lh; super = HLabel("parent", 2),
                        sub = HLabel("child", i), unsafe=false
                    ) == Success
                else
                    _add_relation!(
                        lh; super = HLabel("child", 12),
                        sub = HLabel("child", i), unsafe=false
                    ) == Success
                end
            end
            @test _label_connected(lh)
            @test _label_nocycle(lh)
            @test _label_unique(lh)
            lh
        end

        tl = deepcopy(lh1)
        _merge_hierarchy!(
            tl,
            lh2;
            augend_parent = HLabel("parent", 1),
            addend_parent = HLabel("parent", 1)
        )
        @test tl == oracle_p1p1

        tl = deepcopy(lh1)
        _merge_hierarchy!(
            tl,
            lh2;
            augend_parent = HLabel("parent", 2),
            addend_parent = HLabel("parent", 1)
        )
        @test tl == oracle_p2p1

        tl = deepcopy(lh1)
        _merge_hierarchy!(
            tl,
            lh2;
            augend_parent = HLabel("parent", 2),
            addend_parent = HLabel("parent", 2)
        )
        @test tl == oracle_p2p2

    end # testset

end #testset

end # function
