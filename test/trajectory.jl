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
    f1 = h5traj_reader("test.hmd", 3, Float64, GeneralSystem)
    for i in 2:length(traj)
        if nv(topology(traj.systems[i])) != nv(topology(t1.systems[i]))
            println(i, " ", reaction_flags[i-1])
        end
    end
    @test strict_eq(traj[1], t1[1])
    @test strict_eq(traj[1], f1[1])
    @test strict_eq(traj, t1; dynamic_only=false)
    #@test strict_eq(traj, f1; dynamic_only=false)

    # slice test
    n = 0
    result = true
    s = System{3, Float64, GeneralSystem}()
    read_snapshot!(s, file, 1)
    for (i, r) in enumerate(traj[1:5:length(traj)])
        n += 1
        read_snapshot!(s, file, 5i-4)
        result &= strict_eq(r, s)
    end
    @test result
    @test n == length(1:5:length(t1))
end




end # let
