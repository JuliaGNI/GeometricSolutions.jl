using GeometricEquations
using GeometricEquations.Tests
using GeometricSolutions
using HDF5
using Test

using GeometricBase: h5load, h5save

@test Base.get_extension(GeometricSolutions, :GeometricSolutionsHDF5Ext) !== nothing

const nstep = 5

# Every entry is distinct across members, steps and components, so a
# misplaced column or member cannot compare equal by accident.
function filled(problem, step)
    sol = EnsembleSolution(problem, step)
    for (j, s) in enumerate(sol), n in 1:nstore(s)

        s[n].q .= 1000j .+ 10n .+ reshape(1:length(s[n].q), size(s[n].q)) ./ 100
    end
    return sol
end

function roundtrip(sol, problem; path = "/")
    file = tempname() * ".h5"
    h5open(io -> h5save(io, sol; path), file, "w")
    return file, h5open(io -> h5load(EnsembleSolution, io, problem; path), file, "r")
end

@testset "$(rpad("EnsembleSolution with per-member parameters",80))" begin
    params = [(k = 0.5,), (k = 1.0,), (k = 2.0,)]
    problem = Tests.ExponentialGrowth.odeensemble(; parameters = params)
    sol = filled(problem, nstep)
    ns = nstore(sol[begin])

    file, sol2 = roundtrip(sol, problem)

    h5open(file, "r") do io
        q = read(io["q"])
        @test size(q) == (1, ns + 1, 3)
        for j in 1:3
            @test q[:, 1, j] == sol[j].q[0] == Tests.ExponentialGrowth.ics[j].q
            @test q[:, ns + 1, j] == sol[j].q[ns]
        end
        @test read(io["parameters/k"]) == [0.5, 1.0, 2.0]
        @test read(io["t"]) == collect(parent(sol[begin].t))
        @test read(attributes(io)["step"]) == nstep
    end

    @test sol2 isa EnsembleSolution
    @test parameters(sol2.problem) == params
    for j in 1:3, k in keys(sol)

        @test axes(sol2[j][k]) == (0:ns,)
        @test sol2[j][k] == sol[j][k]
    end
    @test sol2[2].q[0] == Tests.ExponentialGrowth.ics[2].q
    @test sol2[2][ns].q == sol[2][ns].q

    other = Tests.ExponentialGrowth.odeensemble(;
        parameters = [(k = 0.5,), (k = 1.5,), (k = 2.0,)])
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, other), file, "r")

    coarser = Tests.ExponentialGrowth.odeensemble(; parameters = params, Δt = 0.2)
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, coarser), file, "r")

    nullparams = ODEEnsemble(
        (v, t, x, params) -> (v .= x; nothing), (0.0, 10.0), 0.1,
        Tests.ExponentialGrowth.ics)
    @test_throws ArgumentError h5open(
        io -> h5load(EnsembleSolution, io, nullparams), file, "r")
end

@testset "$(rpad("EnsembleSolution with matrix states and array parameters",80))" begin
    vectorfield(v, t, x, params) = (v .= params.a .* x; nothing)
    ics = [(q = StateVariable(rand(2, 3)),), (q = StateVariable(rand(2, 3)),)]
    params = [(a = [1.0, 2.0],), (a = [3.0, 4.0],)]
    problem = ODEEnsemble(vectorfield, (0.0, 1.0), 0.1, ics; parameters = params)
    sol = filled(problem, 2)
    ns = nstore(sol[begin])

    file, sol2 = roundtrip(sol, problem; path = "ensemble")

    h5open(file, "r") do io
        @test size(read(io["ensemble/q"])) == (2, 3, ns + 1, 2)
        @test read(io["ensemble/q"])[:, :, 1, 2] == ics[2].q
        @test read(io["ensemble/parameters/a"]) == [1.0 3.0; 2.0 4.0]
    end

    for j in 1:2, k in keys(sol)

        @test axes(sol2[j][k]) == (0:ns,)
        @test sol2[j][k] == sol[j][k]
    end
end

@testset "$(rpad("EnsembleSolution without parameters",80))" begin
    problem = ODEEnsemble(
        (v, t, x, params) -> (v .= x; nothing), (0.0, 1.0), 0.1,
        Tests.ExponentialGrowth.ics)
    sol = filled(problem, 1)

    file, sol2 = roundtrip(sol, problem)

    @test !h5open(io -> haskey(io, "parameters"), file, "r")
    for j in 1:3
        @test sol2[j].q == sol[j].q
    end
end
