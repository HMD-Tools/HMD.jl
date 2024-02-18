let



reaction_flags = [Bool(rand(DiscreteUniform(0, 1))) for _ in 1:100]
@testset "invariance on trajectory save" begin
    traj = let
        t = Trajectory(sample)
        for i in 2:101
            add_snapshot!(t, sample; reaction = reaction_flags[i-1])
        end
        t
    end
    hmdsave("test.hmd", traj)
    t1 = read_traj("test.hmd", 3, Float64, GeneralSystem)
    strict_eq(traj, t1; dynamic_only=true)
end




end # let
