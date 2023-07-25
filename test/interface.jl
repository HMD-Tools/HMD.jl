s = System{3, Float64}()


@testset "system core interface" begin
    @test dimension(s) == 3
    @test precision(s) == Float64
    @test system_type(s) == GeneralSystem
    @test begin
        s2 = similar(s)
        dimension(s2) == dimension(s) && \
        precision(s2) == precision(s) && \
        system_type(s2) == system_type(s) 
    end

    

end