module GeometricSolutionsHDF5Ext

using GeometricEquations: EnsembleProblem, NullParameters
using GeometricSolutions
using HDF5: H5DataStore, attributes, create_group

import GeometricBase: h5save, h5load, nsamples, parameters, timespan, timestep

# The file stores every state variable of every member in one dataset of size
# (size(x)..., nstore + 1, nsamples). Column n + 1 holds time index n, so the
# 0-based axis of a DataSeries maps to the 1-based axis of the file.

_group(h5::H5DataStore, path) = path == "/" ? h5 : create_group(h5, path)

_statekeys(sol::EnsembleSolution) = filter(!=(:t), keys(sol))

_stack(sol::EnsembleSolution, k) = stack(stack(parent(parent(s[k]))) for s in sol.s)

function _stack_parameters(ps::AbstractVector{<:NamedTuple})
    names = keys(ps[begin])
    all(p -> keys(p) == names, ps) ||
        throw(ArgumentError("the members' parameter sets have different names"))
    NamedTuple{names}(map(names) do k
        v = ps[begin][k]
        v isa Union{Number, AbstractArray{<:Number}} ||
            throw(ArgumentError("parameter $k is a $(typeof(v)), which HDF5 cannot store"))
        stack(p[k] for p in ps)
    end)
end

"""
    h5save(h5, sol::EnsembleSolution; path = "/")

Write the time series, every state variable and every member's parameters of `sol` into the
group `path` of `h5`.

A state variable `q` becomes the dataset `q` of size `(size(q)..., nstore + 1, nsamples)`, with
column `n + 1` holding time index `n`. A parameter `β` becomes `parameters/β` of size
`(size(β)..., nsamples)`. The time series is `t`, and `step`, `nstore`, `nsamples`, `timestep`
and `timespan` are attributes of the group. The equation is not stored:
[`h5load`](@ref GeometricBase.h5load) takes it from the problem it is given.
"""
function h5save(h5::H5DataStore, sol::EnsembleSolution; path::AbstractString = "/")
    # Stacked before the first write, so that a parameter HDF5 cannot store leaves no partial group.
    ps = parameters(sol.problem)
    sp = eltype(ps) <: NullParameters ? nothing : _stack_parameters(ps)

    g = _group(h5, path)
    s₀ = sol[begin]
    ns = nstore(s₀)

    attributes(g)["step"] = step(s₀)
    attributes(g)["nstore"] = ns
    attributes(g)["nsamples"] = nsamples(sol)
    attributes(g)["timestep"] = timestep(sol)
    attributes(g)["timespan"] = collect(timespan(sol))

    g["t"] = collect(parent(s₀[:t]))

    for k in _statekeys(sol)
        g[string(k)] = _stack(sol, k)
    end

    if !isnothing(sp)
        gp = create_group(g, "parameters")
        for (k, v) in pairs(sp)
            gp[string(k)] = v
        end
    end

    return h5
end

function _check(name, stored, expected)
    stored == expected ||
        throw(ArgumentError("the file has $name = $stored, the problem has $expected"))
end

"""
    h5load(EnsembleSolution, h5, problem::EnsembleProblem; path = "/")

Read the `EnsembleSolution` that [`h5save`](@ref GeometricBase.h5save) wrote into the group
`path` of `h5`.

HDF5 cannot hold the equation, so `problem` supplies it. The solution is built as
`EnsembleSolution(problem, step)` and filled from the file, which gives every `DataSeries` its
0-based time axis. The read throws an `ArgumentError` if the file does not belong to `problem`:
a different number of members, time step, time span or number of stored steps, or different
parameters for any member.
"""
function h5load(::Type{EnsembleSolution}, h5::H5DataStore, problem::EnsembleProblem;
        path::AbstractString = "/")
    g = path == "/" ? h5 : h5[path]
    attr(name) = read(attributes(g)[name])

    sol = EnsembleSolution(problem, attr("step"))
    ns = nstore(sol[begin])

    _check("nsamples", attr("nsamples"), nsamples(sol))
    _check("nstore", attr("nstore"), ns)
    _check("timestep", attr("timestep"), timestep(sol))
    _check("timespan", attr("timespan"), collect(timespan(sol)))

    ps = parameters(problem)
    if eltype(ps) <: NullParameters
        haskey(g, "parameters") &&
            throw(ArgumentError("the file has parameters, the problem has none"))
    else
        haskey(g, "parameters") ||
            throw(ArgumentError("the file has no parameters, the problem has some"))
        stored = g["parameters"]
        expected = _stack_parameters(ps)
        names = sort([string(k) for k in keys(expected)])
        _check("parameter names", sort(keys(stored)), names)
        for (k, v) in pairs(expected)
            _check("parameters $k", read(stored[string(k)]), v)
        end
    end

    for k in _statekeys(sol)
        A = read(g[string(k)])
        x₀ = sol[begin][k][0]
        _check("size of $k", size(A), (size(x₀)..., ns + 1, nsamples(sol)))
        for (j, s) in enumerate(sol.s), n in 0:ns

            s[k][n] = @view A[axes(x₀)..., n + 1, j]
        end
    end

    return sol
end

end
