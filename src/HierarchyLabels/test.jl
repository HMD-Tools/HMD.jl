using Test

function test()

@testset "HierarchyHLabels" begin
    l1, l2, l3 = HLabel("l1", 1), HLabel("l2", 2), HLabel("l3", 3)
    noexist    = HLabel("not_exist_in_h", -1)
    handmade = begin
        lh = LabelHierarchy()
        @assert add_vertex!(_hierarchy(lh)); push!(lh.labels, l1); push!(lh.label2node, l1 => 1)
        @assert add_vertex!(_hierarchy(lh)); push!(lh.labels, l2); push!(lh.label2node, l2 => 2)
        @assert add_vertex!(_hierarchy(lh)); push!(lh.labels, l3); push!(lh.label2node, l3 => 3)
        add_edge!(_hierarchy(lh), 2, 1)
        add_edge!(_hierarchy(lh), 3, 2)
        lh
    end

    # addの順番違いとinsert，エラー発生時の不変性をテスト
    addrel1 = begin
        lh = LabelHierarchy()
        _add_label!(lh, l1)
        _add_label!(lh, l2)
        _add_label!(lh, l3)
        _add_relation!(lh; super = l1, sub = l2)
        _add_relation!(lh; super = l2, sub = l3)
        lh
    end
    addrel2 = begin
        lh = LabelHierarchy()
        @test _add_label!(lh, l2) == Success
        @test _add_label!(lh, l1) == Success
        @test _add_label!(lh, l3) == Success
        @test _add_relation!(lh; super = l2, sub = l3) == Success
        @test _add_relation!(lh; super = l1, sub = l2) == Success
        lh
    end
    multiadd = begin
        lh = LabelHierarchy()
        @test _add_labels!(lh, [l1, l2, l3]) == Success
        @test _add_relation!(lh; super = l2, sub = l3) == Success
        @test _add_relation!(lh; super = l1, sub = l2) == Success
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

    @test_throws KeyError _get_nodeid(addrel1, noexist)
    @test_throws ErrorException _has_relation(addrel1, noexist, l2)
    @test_throws ErrorException _get_label(LabelHierarchy(), 0)
    @test_throws BoundsError _get_label(addrel1, -1)

    @test _add_label!(addrel1, l1) == Label_Occupied
    @test _add_relation!(addrel1; super=noexist, sub=l3) == Label_Missing
    @test _add_relation!(addrel1; super=l3, sub=l3)      == Label_Duplication
    @test _add_relation!(addrel1; super=l2, sub=l3)      == Relation_Occupied
    @test _add_relation!(addrel1; super=l2, sub=l3)      == Relation_Occupied

    addrel1 = begin
        lh = LabelHierarchy()
        _add_label!(lh, l1)
        _add_label!(lh, l2)
        _add_label!(lh, l3)
        _add_relation!(lh; super = l1, sub = l2)
        _add_relation!(lh; super = l2, sub = l3)
        lh
    end
    lh = deepcopy(addrel1)
    @test lh == addrel1
    @test _remove_relation!(lh, noexist, l3) == Relation_Missing
    @test _remove_relation!(lh, l2, noexist) == Relation_Missing
    @test lh == addrel1
    @test _remove_label!(lh, l1) == Label_Missing
    @test lh == addrel1
    @test _remove_label!(lh, l1) == Success
    @test _remove_relation!(lh, l2, l3) == Success
end

end #test
