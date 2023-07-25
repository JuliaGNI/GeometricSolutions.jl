"""
`EnsembleSolution`: Collection of all solutions of an `EnsembleProblem`.
"""
struct EnsembleSolution{dType, tType, sType, probType} <: AbstractSolution{dType,tType}
    t::TimeSeries{tType}
    s::sType

    problem::probType

    function EnsembleSolution(problem::GeometricEnsemble; step = 1)
        t = TimeSeries(tbegin(problem), tend(problem), tstep(problem))
        s = [GeometricSolution(t, p, step) for p in problem]
        
        new{datatype(problem), timetype(problem), typeof(s), typeof(problem)}(t, s, problem)
    end
end

@inline GeometricBase.nsamples(sol::EnsembleSolution) = length(sol.s)
@inline GeometricBase.eachsample(sol::EnsembleSolution) = 1:nsamples(sol)

@inline GeometricBase.tspan(sol::EnsembleSolution) = tspan(sol.t)
@inline GeometricBase.tstep(sol::EnsembleSolution) = tstep(sol.t)

@inline GeometricBase.ntime(sol::EnsembleSolution) = ntime(sol.t)
@inline GeometricBase.timesteps(sol::EnsembleSolution) = sol.t
@inline GeometricBase.eachtimestep(sol::EnsembleSolution) = eachtimestep(sol.t)

@inline solution(sol::EnsembleSolution, i) = sol.s[i]

Base.length(sol::EnsembleSolution) = length(sol.s)
Base.iterate(sol::EnsembleSolution, i=1) = i > length(sol) ? nothing : (solution(sol, i), i+1)

function relative_maximum_error(sols::EnsembleSolution, refs::EnsembleSolution)
    @assert nsamples(sols) == nsamples(refs)

    _errs = [relative_maximum_error(s...) for s in zip(sols, refs)]
    _keys = keys(_errs[begin])
    _vals = [getproperty.(_errs, i) for i in _keys]

    NamedTuple{_keys}(maximum.(_vals))
end
