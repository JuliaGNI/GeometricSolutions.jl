using GeometricBase
using GeometricEquations
using GeometricEquations.Tests
using GeometricSolutions
using Test


const nstep = 10


@testset "$(rpad("Interface Definition",80))" begin
    struct TestSolution{dType, tType} <: AbstractSolution{dType, tType} end
    test_sol = TestSolution{Float64, Float64}()

    @test_throws ErrorException step(test_sol)
    @test_throws ErrorException nstore(test_sol)
    @test_throws ErrorException ntime(test_sol)

    # @test_throws ErrorException offset(test_sol)
    # @test_throws ErrorException lastentry(test_sol)
    @test_throws ErrorException timesteps(test_sol)
    @test_throws ErrorException eachtimestep(test_sol)
end


@testset "$(rpad("Geometric Solution",80))" begin

    prob = Tests.ExponentialGrowth.odeproblem()
    sol1 = GeometricSolution(prob)
    sol2 = GeometricSolution(prob; step = nstep)

    @test sol1.t == sol2.t
    @test sol1.s != sol2.s

    @test step(sol1) == sol1.step == 1
    @test ntime(sol1) == ntime(sol1.t)
    @test nstore(sol1) == sol1.nstore == ntime(sol1.t)
    @test timesteps(sol1) == collect(tbegin(prob):timestep(prob):tend(prob))

    @test step(sol2) == sol2.step == 10
    @test ntime(sol2) == ntime(sol2.t)
    @test nstore(sol2) == sol2.nstore == div(ntime(sol2), nstep)
    @test timesteps(sol2) == collect(tbegin(prob):timestep(prob):tend(prob))

    @test sol1[0].t == initial_conditions(prob).t
    @test sol1[0].q == initial_conditions(prob).q

    t = timesteps(sol1)
    q = initial_conditions(prob).q
    for i in eachtimestep(t)
        x = zero(q)
        Tests.ExponentialGrowth.solution(x, t[i], q, t[i-1], parameters(prob))
        sol1[i] = (q = x,)
        sol2[i] = (q = x,)
    end

    @test sol1.s.q[ntime(t)] == sol2.s.q[div(ntime(t), nstep)] == sol1.s.q[end] == sol2.s.q[end]

end


@testset "$(rpad("Ensemble Solution",80))" begin

end
