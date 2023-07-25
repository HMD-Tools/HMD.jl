using Test

function test()

@testset "DataTypes" begin
    @testset "Position" begin
        @test Position{3, Float64}() == Position(3, Float64, 0) == Position()
        @test Position{3, Float64}(10) == Position(3, Float64, 10) == Position(10)
        p = Position()
        for _ in 1:10
            push!(p, zeros(Float64, 3))
        end
        @test p == Position(10)
    end

    @testset "BoundingBox" begin
        @test BoundingBox{3, Float64}() == BoundingBox(3, Float64, zeros(Float64, 3), Matrix{Float64}(I, 3, 3))
        @test BoundingBox{3, Float64}() == BoundingBox(zeros(Float64, 3), Matrix{Float64}(I, 3, 3))
        @test_throws ErrorException BoundingBox(3, Float64, zeros(Float64, 3), zeros(3, 3))
        @test_throws ErrorException BoundingBox(-3, Float64, zeros(Float64, 3), zeros(3, 3))
        @test_throws ErrorException BoundingBox(3, Float64, zeros(Float64, 2), zeros(3, 3))
    end

    @testset "System" begin
        s = System{3, Float64}()

    end

    @testset "HLabel manipulation" begin

    end

    @testset "property" begin

    end
end

end #test
