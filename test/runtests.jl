using linearstep
using Test


@time @testset "linearstep.jl" begin include("test_adam.jl") end
@time @testset "compare_adam.jl" begin include("test_compare_adam.jl") end
