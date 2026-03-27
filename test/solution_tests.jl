using GeometricEquations
using GeometricEquations.Tests
using GeometricSolutions
using Test

using GeometricBase: eachsample
using GeometricSolutions: arrays

const nstep = 10

@testset "$(rpad("Geometric Solution",80))" begin
    prob = Tests.ExponentialGrowth.odeproblem()
    sol1 = GeometricSolution(prob)
    sol2 = GeometricSolution(prob, nstep)

    @test sol1.timeser == sol2.timeser
    @test sol1.dataser != sol2.dataser
    @test sol1.q == sol1[:q]

    @test step(sol1) == sol1.step == 1
    @test ntime(sol1) == ntime(sol1.timeser)
    @test nstore(sol1) == sol1.nstore == ntime(sol1.timeser)
    @test timesteps(sol1) == collect(initialtime(prob):timestep(prob):finaltime(prob))

    @test step(sol2) == sol2.step == 10
    @test ntime(sol2) == ntime(sol2.timeser)
    @test nstore(sol2) == sol2.nstore == div(ntime(sol2), nstep)
    @test timesteps(sol2) == collect(initialtime(prob):timestep(prob):finaltime(prob))

    @test sol1[0].t == initialstate(prob).t
    @test sol1[0].q == initialstate(prob).q

    @test sol1[1].t == initialstate(prob).t + timestep(prob)
    @test sol2[1].t == initialstate(prob).t + timestep(prob) * nstep

    t = timesteps(sol1)
    s = initialstate(prob)
    for i in eachtimestep(t)
        q̄ = copy(s.q)
        Tests.ExponentialGrowth.solution(s.q, t[i], q̄, t[i - 1], parameters(prob))
        copy!(sol1, s, i)
        copy!(sol2, s, i)
    end

    @test sol1.q[ntime(t)] ==
          sol2.q[div(ntime(t), nstep)] ==
          sol1.q[end] ==
          sol2.q[end]

    sol = GeometricSolution(prob)
    for n in eachtimestep(sol)
        Tests.ExponentialGrowth.solution(
            sol[n].q, sol[n].t, sol[0].q, sol[0].t, parameters(prob))
    end

    @test relative_maximum_error(sol, sol).q == 0

    # println(arrays(sol).q)

    # @test arrays(sol).q == hcat((vec(sol.q[n]) for n in eachtimestep(sol))...)
    @test arrays(sol).q == hcat((vec(q) for q in sol.q)...)
end

@testset "$(rpad("Ensemble Solution",80))" begin
    probs = Tests.ExponentialGrowth.odeensemble()

    esol1 = EnsembleSolution(probs)
    esol2 = EnsembleSolution(probs, nstep)

    @test esol1.t == esol2.t

    @test ntime(esol1) == ntime(esol1.t)
    @test timesteps(esol1) == collect(initialtime(probs):timestep(probs):finaltime(probs))

    @test ntime(esol2) == ntime(esol2.t)
    @test timesteps(esol2) == collect(initialtime(probs):timestep(probs):finaltime(probs))

    @test esol1[1] == solution(esol1, 1)
    @test esol1[begin] == solution(esol1, 1)
    @test esol1[:] == solution(esol1, :)

    sols = (solution(esol1, 1), solution(esol1, 2), solution(esol1, 3))

    for sol in esol1
        @test sol ∈ sols
    end

    sols = EnsembleSolution(probs)
    for (sol, prob) in zip(sols.s, probs)
        for n in eachtimestep(sol)
            Tests.ExponentialGrowth.solution(
                sol[n].q, sol[n].t, sol[0].q, sol[0].t, parameters(prob))
        end
    end

    @test relative_maximum_error(sols, sols).q == 0

    @test arrays(sols).q == hcat((arrays(sols[n]).q for n in eachsample(sols))...)
end
