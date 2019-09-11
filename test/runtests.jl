using Test, SpikeSynchrony

@testset "SPIKE distance" begin include("SPIKEdistance_tests.jl") end
@testset "vanRossum" begin include("vanRossum_tests.jl") end
