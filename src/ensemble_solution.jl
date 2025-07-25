"""
`EnsembleSolution`: Collection of all solutions of an `EnsembleProblem`.
"""
struct EnsembleSolution{dType, tType, sType, probType} <: AbstractSolution{dType, tType}
    t::TimeSeries{tType}
    s::sType

    problem::probType

    function EnsembleSolution(problem::EnsembleProblem, step::Int = 1)
        t = TimeSeries(initialtime(problem), finaltime(problem), timestep(problem))
        s = [GeometricSolution(t, p, step) for p in problem]

        new{datatype(problem), timetype(problem), typeof(s), typeof(problem)}(t, s, problem)
    end
end

@inline Base.keys(sol::EnsembleSolution) = keys(sol[begin])

@inline GeometricBase.nsamples(sol::EnsembleSolution) = length(sol.s)
@inline GeometricBase.eachsample(sol::EnsembleSolution) = 1:nsamples(sol)

@inline GeometricBase.datatype(sol::EnsembleSolution{DT, TT}) where {DT, TT} = DT
@inline GeometricBase.timetype(sol::EnsembleSolution{DT, TT}) where {DT, TT} = TT

@inline GeometricBase.timespan(sol::EnsembleSolution) = timespan(sol.t)
@inline GeometricBase.timestep(sol::EnsembleSolution) = timestep(sol.t)

@inline GeometricBase.ntime(sol::EnsembleSolution) = ntime(sol.t)
@inline GeometricBase.timesteps(sol::EnsembleSolution) = sol.t
@inline GeometricBase.eachtimestep(sol::EnsembleSolution) = eachtimestep(sol.t)

@inline solution(sol::EnsembleSolution, i) = sol.s[i]

Base.getindex(sol::EnsembleSolution, i...) = getindex(sol.s, i...)
Base.firstindex(sol::EnsembleSolution) = firstindex(sol.s)
Base.lastindex(sol::EnsembleSolution) = lastindex(sol.s)

Base.length(sol::EnsembleSolution) = length(sol.s)
function Base.iterate(sol::EnsembleSolution, i = 1)
    i > length(sol) ? nothing : (solution(sol, i), i + 1)
end

function arrays(solutions::EnsembleSolution)
    z = (hcat((Array(sol[k]) for sol in solutions)...) for k in keys(solutions))
    NamedTuple{keys(solutions)}(z)
end

function relative_maximum_error(sols::EnsembleSolution, refs::EnsembleSolution)
    @assert nsamples(sols) == nsamples(refs)

    _errs = [relative_maximum_error(s...) for s in zip(sols, refs)]
    _keys = keys(_errs[begin])
    _vals = [getproperty.(_errs, i) for i in _keys]

    NamedTuple{_keys}(maximum.(_vals))
end
