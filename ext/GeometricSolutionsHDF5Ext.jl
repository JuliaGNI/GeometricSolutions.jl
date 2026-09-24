module GeometricSolutionsHDF5Ext

using GeometricBase: nsamples, parameters, timespan, timestep
using GeometricEquations: EnsembleProblem, NullParameters
using GeometricSolutions: EnsembleSolution, nstore
using HDF5: H5DataStore, attributes, create_group

import GeometricBase: h5save, h5load

# The file stores every state variable of every member in one dataset of size
# (size(x)..., nstore + 1, nsamples). Column n + 1 holds time index n, so the
# 0-based axis of a DataSeries maps to the 1-based axis of the file.

# The element types that HDF5 writes as a dataset and reads back as the same type.
const _Bits = Union{
    Bool, Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32,
    Float64}
const _Storable = Union{_Bits, Complex{<:_Bits}}

_group(h5::H5DataStore, path) = path == "/" ? h5 : create_group(h5, path)

_statekeys(sol::EnsembleSolution) = filter(!=(:t), keys(sol))

_stack(sol::EnsembleSolution, k) = stack(stack(parent(parent(s[k]))) for s in sol.s)

function _stack_parameters(ps::AbstractVector{<:NamedTuple})
    names = keys(ps[begin])
    all(p -> keys(p) == names, ps) ||
        throw(ArgumentError("the members' parameter sets have different names"))
    NamedTuple{names}(map(names) do k
        v = ps[begin][k]
        v isa Union{_Storable, AbstractArray{<:_Storable}} ||
            throw(ArgumentError("parameter $k is a $(typeof(v)), which HDF5 cannot store and read back"))
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

# Fills the series `x` of member `j` from `A`, behind a function barrier, since neither `read`
# nor the member of an `EnsembleSolution` infers a concrete type. Column 1 of `A` must equal the
# initial condition that `EnsembleSolution(problem, step)` took from the problem.
function _fill!(x, A::AbstractArray, k, j, ns)
    x₀ = x[0]
    _check("initial $k of member $j", A[axes(x₀)..., 1, j], x₀)
    for n in 1:ns
        x[n] = @view A[axes(x₀)..., n + 1, j]
    end
end

"""
    h5load(EnsembleSolution, h5, problem::EnsembleProblem; path = "/")

Read the `EnsembleSolution` that [`h5save`](@ref GeometricBase.h5save) wrote into the group
`path` of `h5`.

HDF5 cannot hold the equation, so `problem` supplies it. The solution is built as
`EnsembleSolution(problem, step)` and filled from the file, which gives every `DataSeries` its
0-based time axis. The read throws an `ArgumentError` if the file does not belong to `problem`:
a different number of members, time step, time span, number of stored steps or set of state
variables, or a different initial condition or different parameters for any member.
"""
function h5load(::Type{EnsembleSolution}, h5::H5DataStore, problem::EnsembleProblem;
        path::AbstractString = "/")
    g = path == "/" ? h5 : h5[path]
    attr(name) = read(attributes(g)[name])

    _check("nsamples", attr("nsamples"), nsamples(problem))
    _check("timestep", attr("timestep"), timestep(problem))
    _check("timespan", attr("timespan"), collect(timespan(problem)))

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

    sol = EnsembleSolution(problem, attr("step"))
    ns = nstore(sol[begin])
    _check("nstore", attr("nstore"), ns)

    datasets = filter(k -> k ∉ ("t", "parameters"), keys(g))
    _check("state variables", sort(datasets), sort([string(k) for k in _statekeys(sol)]))

    for k in _statekeys(sol)
        A = read(g[string(k)])
        _check("size of $k", size(A), (size(sol[begin][k][0])..., ns + 1, nsamples(sol)))
        for (j, s) in enumerate(sol.s)
            _fill!(s[k], A, k, j, ns)
        end
    end

    return sol
end

end
