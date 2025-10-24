"""
`EnsembleSolution`: Collection of solutions for an ensemble of geometric differential equations

Contains all fields necessary to store the solutions of an [`EnsembleProblem`](@extref GeometricEquations.EnsembleProblem),
which represents multiple instances of a geometric differential equation with different initial conditions
and/or parameters. Each solution in the ensemble is stored as a [`GeometricSolution`](@ref).

### Fields

* `t`: common time series shared across all solutions in the ensemble
* `s`: vector of [`GeometricSolution`](@ref) objects, one for each problem instance
* `problem`: the original [`EnsembleProblem`](@extref GeometricEquations.EnsembleProblem)

### Type Parameters

* `dType`: data type for solution values (e.g., `Float64`)
* `tType`: data type for time values (e.g., `Float64`)
* `sType`: type of the solution vector
* `probType`: type of the ensemble problem

### Constructor

```julia
EnsembleSolution(problem::EnsembleProblem, step::Int = 1)
```

Creates an `EnsembleSolution` from an [`EnsembleProblem`](@extref GeometricEquations.EnsembleProblem).
The optional parameter `step` determines the intervals for storing each solution,
i.e., if `step > 1` only every `step`'th solution point is stored for each ensemble member.

### Usage

An `EnsembleSolution` can be iterated over to access individual solutions:

```julia
for sol in ensemble_solution
    # sol is a GeometricSolution
    process_solution(sol)
end
```

Individual solutions can be accessed by index:

```julia
first_solution = ensemble_solution[1]
last_solution = ensemble_solution[end]
```

The number of solutions in the ensemble:

```julia
n_solutions = nsamples(ensemble_solution)
```

Convert all solutions to arrays for analysis:

```julia
solution_arrays = arrays(ensemble_solution)
```

### See Also

* [`GeometricSolution`](@ref): Individual solution type
* [`EnsembleProblem`](@extref GeometricEquations.EnsembleProblem): Ensemble problem type
* [`arrays`](@ref): Convert solutions to arrays for analysis
* [`relative_maximum_error`](@ref): Compare ensemble solutions
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
