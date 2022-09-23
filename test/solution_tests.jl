using GeometricBase
using GeometricEquations
using GeometricEquations.Tests
using GeometricSolutions
using Test


const step = 10


@testset "$(rpad("Interface Definition",80))" begin
    struct TestSolution{dType, tType} <: AbstractSolution{dType, tType} end
    test_sol = TestSolution{Float64, Float64}()

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


@testset "$(rpad("Geometric Solution",80))" begin

    prob = Tests.ExponentialGrowth.odeproblem()
    sol1 = GeometricSolution(prob)
    sol2 = GeometricSolution(prob; step = step)

    @test sol1.t == sol2.t
    @test sol1.s != sol2.s

    @test sol1.step == 1
    @test sol1.nstore == ntime(sol1.t)
    @test sol1.counter == 0

    @test sol2.step == 10
    @test sol2.nstore == div(ntime(sol2.t), step)
    @test sol2.counter == 0

    t = timesteps(sol1)
    q = initial_conditions(prob).q
    for i in eachtimestep(t)
        x = Tests.ExponentialGrowth.solution(t[i], q, zero(q), parameters(prob))
        sol1[i] = (q = x,)
        sol2[i] = (q = x,)
    end

    @test sol1.s.q[ntime(t)] == sol2.s.q[div(ntime(t), step)] == sol1.s.q[end] == sol2.s.q[end]

end


@testset "$(rpad("Ensemble Solution",80))" begin

end
