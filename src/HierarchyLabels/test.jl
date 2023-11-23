using Test

function test()

@testset "HierarchyHLabels" begin
    l1, l2, l3 = HLabel("l1", 1), HLabel("l2", 2), HLabel("l3", 3)
    noexist    = HLabel("not_exist_in_h", -1)
    handmade = let
        lh = LabelHierarchy()
        @assert add_vertex!(_hierarchy(lh)); push!(lh.labels, l1); push!(lh.label2node, l1 => 1)
        @assert add_vertex!(_hierarchy(lh)); push!(lh.labels, l2); push!(lh.label2node, l2 => 2)
        @assert add_vertex!(_hierarchy(lh)); push!(lh.labels, l3); push!(lh.label2node, l3 => 3)
        add_edge!(_hierarchy(lh), 2, 1)
        add_edge!(_hierarchy(lh), 3, 2)
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
    @test _super(addrel1, l3) == [l2] && _super(addrel1, l2) == [l1]
    @test _sub(addrel1, l1)   == [l2] && _sub(addrel1, l2)   == [l3]
    @test _issuper(addrel1, l1, l2) && _issuper(addrel1, l2, l3)
    @test _issub(addrel1, l2, l1)   && _issub(addrel1, l3, l2)

    @inferred _contains(addrel1, l1)
    @inferred _get_nodeid(addrel1, l1)
    @inferred _sub(addrel1, l1)
    @inferred _super(addrel1, l3)
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
        @test _add_labels!(lh, [l1, l2, l2]) == Label_Occupied
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
        @test _add_relation!(lh; super = l3, sub = l1, unsafe=true) == Success
        @test _label_nocycle(lh) == false
        # add isolated label
        @test _label_connected(lh) == true
        @test _add_label!(lh, HLabel("isolate", 1), true) == Success
        @test _label_connected(lh) == false
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
        # hierarchy merge test
        # small piece
        lh1 = let
            lh = LabelHierarchy()
            _add_labels!(lh, [HLabel("lh1", i) for i in 1:6])
            # under 1
            for i in [2,3]
                _add_relation!(lh; super = HLabel("lh1", 1), sub = HLabel("lh1", i), unsafe=false)
            end
            # under 2
            _add_relation!(lh; super = HLabel("lh1", 2), sub = HLabel("lh1", 4), unsafe=false)
            # under 3
            for i in [5,6]
                _add_relation!(lh; super = HLabel("lh1", 3), sub = HLabel("lh1", i), unsafe=false)
            end
            lh
        end
        lh2 = let
            lh = LabelHierarchy()
            _add_labels!(lh, [HLabel("lh2", i) for i in 1:5])
            # under 1
            for i in [2,3]
                _add_relation!(lh; super = HLabel("lh2", 1), sub = HLabel("lh2", i), unsafe=false)
            end
            # under 2
            _add_relation!(lh; super = HLabel("lh2", 2), sub = HLabel("lh2", 4), unsafe=false)
            # under 3
            _add_relation!(lh; super = HLabel("lh2", 3), sub = HLabel("lh2", 5), unsafe=false)
            lh
        end
        lh3 = let
            lh = LabelHierarchy()
            _add_labels!(lh, [HLabel("lh3", i) for i in 1:5])
            # under 1
            for i in [2,3]
                _add_relation!(lh; super = HLabel("lh3", 1), sub = HLabel("lh3", i), unsafe=false)
            end
            # under 3
            for i in [4,5]
                _add_relation!(lh; super = HLabel("lh3", 3), sub = HLabel("lh3", i), unsafe=false)
            end
            lh
        end
        lh12 = let
            lh = LabelHierarchy()
            _add_labels!(lh, [HLabel("lh1", i) for i in 1:5])
            # under 1
            for i in [2,3]
                _add_relation!(lh; super = HLabel("lh1", 1), sub = HLabel("lh1", i), unsafe=false)
            end
            # under 2
            _add_relation!(lh; super = HLabel("lh1", 2), sub = HLabel("lh1", 4), unsafe=false)
            # under 3
            _add_relation!(lh; super = HLabel("lh1", 3), sub = HLabel("lh1", 5), unsafe=false)
            lh
        end
        lh13 = let
            lh = LabelHierarchy()
            _add_labels!(lh, [HLabel("lh1", i) for i in 1:5])
            # under 1
            for i in [2,3]
                _add_relation!(lh; super = HLabel("lh1", 1), sub = HLabel("lh1", i), unsafe=false)
            end
            # under 3
            for i in [4,5]
                _add_relation!(lh; super = HLabel("lh1", 3), sub = HLabel("lh1", i), unsafe=false)
            end
            lh
        end


        # test oracle
        oracle1 = let
            # oracle for
            #   _merge_hierarchy!(
            #       tl,
            #       lh2;
            #       augend_parent = HLabel("lh1", 1),
            #       addend_parent = HLabel("lh2", 1)
            #   )
            lh = deepcopy(lh1)
            # HLabel("lh2", 1) in lh2 is removed, then label ids in lh2 are shifted: 2:5 -> 1:4
            @test _add_labels!(lh, [HLabel("lh2", i) for i in (2-1):(5-1)]) == Success
            # under ("l1", 1)
            @test _add_relation!(lh; super=HLabel("lh1", 1), sub=HLabel("lh2", 2-1), unsafe=false) == Success
            @test _add_relation!(lh; super=HLabel("lh1", 1), sub=HLabel("lh2", 3-1), unsafe=false) == Success
            # under ("lh2", 2), ("lh2", 3)
            @test _add_relation!(lh; super=HLabel("lh2", 2-1), sub=HLabel("lh2", 4-1), unsafe=false) == Success
            @test _add_relation!(lh; super=HLabel("lh2", 3-1), sub=HLabel("lh2", 5-1), unsafe=false) == Success
            @test _label_connected(lh)
            @test _label_nocycle(lh)
            @test _label_unique(lh)
            lh
        end
        oracle2 = let
            # oracle for
            #   _merge_hierarchy!(
            #       tl,
            #       lh3;
            #       augend_parent = HLabel("lh1", 5),
            #       addend_parent = HLabel("lh3", 3)
            #   )
            lh = deepcopy(lh1)
            @test _add_labels!(lh, [HLabel("lh3", i) for i in (4-3):(5-3)]) == Success
            # under ("l1", 3)
            @test _add_relation!(lh; super=HLabel("lh1", 3), sub=HLabel("lh3", 4-3), unsafe=false) == Success
            @test _add_relation!(lh; super=HLabel("lh1", 3), sub=HLabel("lh3", 5-3), unsafe=false) == Success
            @test _label_connected(lh)
            @test _label_nocycle(lh)
            @test _label_unique(lh)
            lh
        end
        oracle11 = let
            # oracle for
            #   _merge_hierarchy!(
            #       tl,
            #       lh12;
            #       augend_parent = HLabel("lh1", 1),
            #       addend_parent = HLabel("lh1", 1)
            #   )
            lh = deepcopy(lh1)
            # HLabel("lh1", 1) in lh2 is removed, then label ids in lh2 are shifted: 2:5 -> 1:4
            @test _add_labels!(lh, [HLabel("lh1", i) for i in 7:10]) == Success
            # under ("l1", 1)
            @test _add_relation!(lh; super=HLabel("lh1", 1), sub=HLabel("lh1", 7), unsafe=false) == Success
            @test _add_relation!(lh; super=HLabel("lh1", 1), sub=HLabel("lh1", 8), unsafe=false) == Success
            # under ("lh1", 2), ("lh1", 3)
            @test _add_relation!(lh; super=HLabel("lh1", 7), sub=HLabel("lh1", 9),  unsafe=false) == Success
            @test _add_relation!(lh; super=HLabel("lh1", 8), sub=HLabel("lh1", 10), unsafe=false) == Success
            @test _label_connected(lh)
            @test _label_nocycle(lh)
            @test _label_unique(lh)
            lh
        end
        oracle12 = let
            # oracle for
            #   _merge_hierarchy!(
            #       tl,
            #       lh13;
            #       augend_parent = HLabel("lh1", 3),
            #       addend_parent = HLabel("lh1", 3)
            #   )
            lh = deepcopy(lh1)
            @test _add_labels!(lh, [HLabel("lh1", i) for i in 7:8]) == Success
            # under ("l1", 3)
            @test _add_relation!(lh; super=HLabel("lh1", 3), sub=HLabel("lh1", 7), unsafe=false) == Success
            @test _add_relation!(lh; super=HLabel("lh1", 3), sub=HLabel("lh1", 8), unsafe=false) == Success
            @test _label_connected(lh)
            @test _label_nocycle(lh)
            @test _label_unique(lh)
            lh
        end

        tl = deepcopy(lh1)
        _merge_hierarchy!(
            tl,
            lh2;
            augend_parent = HLabel("lh1", 1),
            addend_parent = HLabel("lh2", 1)
        )
        @test tl == oracle1

        tl = deepcopy(lh1)
        _merge_hierarchy!(
            tl,
            lh3;
            augend_parent = HLabel("lh1", 3),
            addend_parent = HLabel("lh3", 3)
        )
        @test tl == oracle2

        tl = deepcopy(lh1)
        _merge_hierarchy!(
            tl,
            lh12;
            augend_parent = HLabel("lh1", 1),
            addend_parent = HLabel("lh1", 1)
        )
        @test tl == oracle11

        tl = deepcopy(lh1)
        _merge_hierarchy!(
            tl,
            lh13;
            augend_parent = HLabel("lh1", 3),
            addend_parent = HLabel("lh1", 3)
        )
        @test tl == oracle12
    end # testset

end #testset

end # function
