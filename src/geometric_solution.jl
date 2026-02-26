"""
`GeometricSolution`: Solution of a geometric differential equation

Contains all fields necessary to store the solution of a GeometricProblem.

### Fields

* `t`:  time steps
* `s`:  NamedTuple of DataSeries for each solution component
* `step`: store every step'th time step (default: 1)
* `nstore`: number of time steps to store
* `offset`: offset of current time step

### Constructors

```julia
GeometricSolution(problem, step = 1)
```

The usual way to initialise a `Solution` is by passing a [`GeometricEquations.EquationProblem`](@extref), which
can for example be an [`GeometricEquations.ODEProblem`](@extref) or [`GeometricEquations.PODEProblem`](@extref).
The optional parameter `step` determines the intervals for storing the solution,
i.e., if `step > 1` only every `step`'th solution is actually stored.

"""
mutable struct GeometricSolution{
    dType, tType, tsType <: TimeSeries{tType}, dsType <: NamedTuple, probType, perType} <:
               AbstractSolution{dType, tType}
    t::tsType
    s::dsType

    problem::probType
    periodicity::perType

    step::Int
    nstore::Int
    offset::Int

    function GeometricSolution(t::TimeSeries, problem::GeometricProblem, step::Int = 1)
        @assert step ≥ 1
        nstore = div(ntime(t), step)
        s = NamedTuple{keys(problem.ics)}(Tuple(DataSeries(x, nstore) for x in problem.ics))
        period = _periodicity(s, periodicity(problem))
        sol = new{datatype(problem), timetype(problem),
            typeof(t), typeof(s), typeof(problem), typeof(period)}(
            t, s, problem, period, step, nstore, 0)
        sol[0] = initial_conditions(problem)
        return sol
    end
end

function GeometricSolution(problem::GeometricProblem, args...)
    t = TimeSeries(initialtime(problem), finaltime(problem), timestep(problem))
    GeometricSolution(t, problem, args...)
end

Base.keys(sol::GeometricSolution) = keys(sol.s)
Base.step(sol::GeometricSolution) = sol.step
nstore(sol::GeometricSolution) = sol.nstore

GeometricBase.datatype(sol::GeometricSolution{DT, TT}) where {DT, TT} = DT
GeometricBase.timetype(sol::GeometricSolution{DT, TT}) where {DT, TT} = TT

GeometricBase.solutions(sol::GeometricSolution) = sol.s
GeometricBase.timesteps(sol::GeometricSolution) = sol.t
GeometricBase.periodicity(sol::GeometricSolution) = sol.periodicity

GeometricBase.ntime(sol::GeometricSolution) = ntime(timesteps(sol))
GeometricBase.timespan(sol::GeometricSolution) = timespan(timesteps(sol))
GeometricBase.timestep(sol::GeometricSolution) = timestep(timesteps(sol))
GeometricBase.eachtimestep(sol::GeometricSolution) = eachtimestep(timesteps(sol))

Base.length(sol::GeometricSolution) = length(solutions(sol))
Base.iterate(sol::GeometricSolution, args...) = iterate(solutions(sol), args...)

function Base.hasproperty(
        ::GeometricSolution{DT, TT, tsType, dsType}, s::Symbol) where {
        DT, TT, tsType, dsType}
    hasfield(dsType, s) || hasfield(GeometricSolution, s)
end

function Base.getproperty(
        sol::GeometricSolution{DT, TT, tsType, dsType}, s::Symbol) where {
        DT, TT, tsType, dsType}
    if hasfield(dsType, s)
        return getfield(getfield(sol, :s), s)
    else
        return getfield(sol, s)
    end
end

Base.getindex(sol::GeometricSolution, ::Val{s}) where {s} = getindex(getfield(sol, :s), s)
Base.getindex(sol::GeometricSolution, s::Symbol) = getindex(sol, Val(s))

function Base.getindex(sol::GeometricSolution, n::Int)
    @assert n ≤ ntime(sol)
    NamedTuple{(:t, keys(sol)...)}((
        timesteps(sol)[n], (sol[Val(k)][n] for k in keys(sol))...))
end

function Base.setindex!(sol::GeometricSolution, s::NamedTuple, n::Int)
    @assert n ≤ ntime(sol)

    if mod(n, step(sol)) == 0
        for k in keys(sol) ∩ keys(s)
            sol[Val(k)][div(n, step(sol))] = s[k]
        end
    end

    return s
end

function Base.setindex!(sol::GeometricSolution, s::State, n::Int)
    @assert Val.(keys(sol)) ⊆ keys(s)
    @assert n ≤ ntime(sol)

    if mod(n, step(sol)) == 0
        for k in Val.(keys(sol))
            sol[k][div(n, step(sol))] = s[k]
        end
    end

    return s
end

function arrays(sol::GeometricSolution)
    NamedTuple{keys(sol)}((Array(sol[k]) for k in keys(sol)))
    # map(s -> Array(s), sol)
end

function relative_maximum_error(sol::GeometricSolution, ref::GeometricSolution)
    @assert keys(sol) == keys(ref)
    NamedTuple{keys(sol)}(relative_maximum_error(s...) for s in zip(sol.s, ref.s))
end
