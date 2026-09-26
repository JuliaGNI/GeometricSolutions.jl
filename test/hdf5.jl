using GeometricEquations
using GeometricEquations.Tests
using GeometricSolutions
using HDF5
using Random
using Test

Random.seed!(1)

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
    for j in 1:3, k in keys(sol)

        @test axes(sol2[j][k]) == (0:ns,)
        @test sol2[j][k] == sol[j][k]
    end
    @test sol2[2].q[0] == Tests.ExponentialGrowth.ics[2].q
    @test sol2[2][ns].q == sol[2][ns].q

    other = Tests.ExponentialGrowth.odeensemble(;
        parameters = [(k = 0.5,), (k = 1.5,), (k = 2.0,)])
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, other), file, "r")

    extra = Tests.ExponentialGrowth.odeensemble(;
        parameters = [(k = 0.5, m = 1.0), (k = 1.0, m = 1.0), (k = 2.0, m = 1.0)])
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, extra), file, "r")
    extrafile, _ = roundtrip(filled(extra, nstep), extra)
    @test_throws ArgumentError h5open(
        io -> h5load(EnsembleSolution, io, problem), extrafile, "r")

    coarser = Tests.ExponentialGrowth.odeensemble(; parameters = params, Δt = 0.2)
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, coarser), file, "r")

    fewer = Tests.ExponentialGrowth.odeensemble(
        Tests.ExponentialGrowth.ics[1:2]; parameters = params[1:2])
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, fewer), file, "r")

    shorter = Tests.ExponentialGrowth.odeensemble(; parameters = params, tend = 5.0)
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, shorter), file, "r")

    # rand gives ics in [0, 1), so these initial conditions differ from every member's.
    moved = Tests.ExponentialGrowth.odeensemble(
        [(q = StateVariable([1.0 + j]),) for j in 1:3]; parameters = params)
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, moved), file, "r")

    pode = PODEEnsemble(
        (v, t, q, p, params) -> (v .= p; nothing),
        (f, t, q, p, params) -> (f .= q; nothing),
        (0.0, 10.0), 0.1,
        [(q = ic.q, p = StateVariable([0.0])) for ic in Tests.ExponentialGrowth.ics];
        parameters = params)
    @test_throws ArgumentError h5open(io -> h5load(EnsembleSolution, io, pode), file, "r")
    podefile, _ = roundtrip(filled(pode, nstep), pode)
    @test_throws ArgumentError h5open(
        io -> h5load(EnsembleSolution, io, problem), podefile, "r")

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

@testset "$(rpad("EnsembleSolution with a parameter HDF5 cannot store",80))" begin
    # HDF5 writes a Rational as a compound type and reads it back as a NamedTuple.
    for params in ([(k = "slow",), (k = "fast",)], [(k = 1 // 2,), (k = 1 // 3,)])
        problem = ODEEnsemble(
            (v, t, x, params) -> (v .= x; nothing), (0.0, 1.0), 0.1,
            Tests.ExponentialGrowth.ics[1:2]; parameters = params)
        sol = filled(problem, 1)

        file = tempname() * ".h5"
        h5open(file, "w") do io
            @test_throws ArgumentError h5save(io, sol; path = "ensemble")
            @test !haskey(io, "ensemble")
        end
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

    # The file-path forms come from GeometricBase and forward to the methods of this extension.
    fpath = tempname() * ".h5"
    h5save(fpath, sol; path = "ensemble")
    sol3 = h5load(EnsembleSolution, fpath, problem; path = "ensemble")
    for j in 1:3
        @test sol3[j].q == sol[j].q
    end
end
