let



reaction_flags = [Bool(rand(DiscreteUniform(0, 1))) for _ in 1:100]
@testset "invariance on trajectory save" begin
    traj = let
        t = Trajectory(sample)
        for i in 2:101
            s = similar(sample; reserve_dynamic=true, reserve_static=true)
            add_snapshot!(t, s; reaction = reaction_flags[i-1])
        end
        t
    end
    hmdsave("test.hmd", traj)
    t1 = read_traj("test.hmd", 3, Float64, GeneralSystem)
    for i in 2:length(traj)
        if nv(topology(traj.systems[i])) != nv(topology(t1.systems[i]))
            println(i, " ", reaction_flags[i-1])
        end
    end
    @test strict_eq(traj[1], t1[1])
    @test strict_eq(traj, t1; dynamic_only=false)
end




end # let
