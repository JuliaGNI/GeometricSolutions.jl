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

GeometricBase.nsamples(sol::EnsembleSolution) = length(sol.s)
GeometricBase.eachsample(sol::EnsembleSolution) = 1:nsamples(sol)
