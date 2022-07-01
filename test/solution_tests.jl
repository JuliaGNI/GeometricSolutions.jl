using GeometricSolutions
using Test


@testset "$(rpad("Interface Definition",80))" begin
    struct TestSolution{dType, tType, N} <: AbstractSolution{dType, tType, N} end
    test_sol = TestSolution{Float64, Float64, 2}()

    @test_throws ErrorException nsave(test_sol)
    @test_throws ErrorException ntime(test_sol)
    @test_throws ErrorException nsamples(test_sol)

    @test_throws ErrorException counter(test_sol)
    @test_throws ErrorException offset(test_sol)
    @test_throws ErrorException lastentry(test_sol)
    @test_throws ErrorException timesteps(test_sol)

    @test_throws ErrorException eachtimestep(test_sol)
    @test_throws ErrorException eachsample(test_sol)    
end
