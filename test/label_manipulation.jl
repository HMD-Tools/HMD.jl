let

function _all_relations(s, hname)
    rels = Set{Tuple{HLabel, HLabel}}()
    stack = [Entire_System]
    while !isempty(stack)
        v = popfirst!(stack)
        nbr = sub(s, hname, v)
        union!(rels, Set((v, u) for u in nbr))
        append!(stack, nbr)
    end

    return rels
end

@testset "label insertion" begin
    rels = _all_relations(sample, "polymeric")
    insert_point = rand(rels)

    inserted = deepcopy(sample)
    insert_relation!(
        inserted, "polymeric",
        HLabel("inserted", 1);
        super = insert_point[1],
        sub = insert_point[2]
    )
    rels_ins = _all_relations(inserted, "polymeric")

    @test setdiff(rels_ins, rels) == Set(
        [(insert_point[1], HLabel("inserted", 1)),
        (HLabel("inserted", 1), insert_point[2])]
    )
    @test setdiff(rels, rels_ins) == Set((insert_point,))
end





end # let
