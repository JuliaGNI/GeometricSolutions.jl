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
    dType, tType,
    tsType <: TimeSeries{tType},
    dsType <: NamedTuple,
    stType <: OffsetVector,
    probType, perType} <:
               AbstractSolution{dType, tType}
    timeser::tsType
    dataser::dsType
    statser::stType

    problem::probType
    periodicity::perType

    step::Int
    nstore::Int
    offset::Int

    function GeometricSolution(t::TimeSeries, problem::GeometricProblem, step::Int = 1)
        @assert step ≥ 1
        nstore = div(ntime(t), step)
        ics = initialstate(problem)

        states = [zero(ics) for _ in 0:nstore]
        statser = OffsetVector(states, 0:nstore)

        dats = Tuple(DataSeries([parent(states[i][k]) for i in eachindex(states)])
        for k in keys(ics))
        dataser = NamedTuple{keys(ics)}(dats)

        period = _periodicity(dataser, periodicity(problem))

        sol = new{datatype(problem), timetype(problem),
            typeof(t), typeof(dataser), typeof(statser), typeof(problem), typeof(period)}(
            t, dataser, statser, problem, period, step, nstore, 0)

        sol[0] = ics
        return sol
    end
end

function GeometricSolution(problem::GeometricProblem, args...)
    t = TimeSeries(initialtime(problem), finaltime(problem), timestep(problem))
    GeometricSolution(t, problem, args...)
end

Base.keys(sol::GeometricSolution) = keys(solutions(sol))
Base.step(sol::GeometricSolution) = sol.step
nstore(sol::GeometricSolution) = sol.nstore

GeometricBase.datatype(sol::GeometricSolution{DT, TT}) where {DT, TT} = DT
GeometricBase.timetype(sol::GeometricSolution{DT, TT}) where {DT, TT} = TT

states(sol::GeometricSolution) = sol.statser
GeometricBase.solutions(sol::GeometricSolution) = sol.dataser
GeometricBase.timesteps(sol::GeometricSolution) = sol.timeser
GeometricBase.periodicity(sol::GeometricSolution) = sol.periodicity

GeometricBase.ntime(sol::GeometricSolution) = ntime(timesteps(sol))
GeometricBase.timespan(sol::GeometricSolution) = timespan(timesteps(sol))
GeometricBase.timestep(sol::GeometricSolution) = timestep(timesteps(sol))
GeometricBase.eachtimestep(sol::GeometricSolution) = eachtimestep(timesteps(sol))

GeometricBase.variables(sol::GeometricSolution, n::Int) = variables(sol[n])

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
        return getfield(getfield(sol, :dataser), s)
    else
        return getfield(sol, s)
    end
end

state(sol::GeometricSolution, n::Int) = states(sol)[n]
solution(sol::GeometricSolution, s::Symbol) = solutions(sol)[s]
Base.time(sol::GeometricSolution, n::Int) = timesteps(sol)[n]

Base.getindex(sol::GeometricSolution, s::Symbol) = solution(sol, s)
Base.getindex(sol::GeometricSolution, n::Int) = state(sol, n)

function Base.setindex!(sol::GeometricSolution, s::NamedTuple, n::Int)
    @assert n ≤ nstore(sol)
    map(k -> sol[k][n] = s[k], keys(sol) ∩ keys(s))
    return s
end

function Base.setindex!(sol::GeometricSolution, st::State, n::Int)
    @assert keys(sol) == keys(st)
    @assert n ≤ nstore(sol)

    # dst = (sol[k][div(n, step(sol))] for k in keys(sol))
    # src = (s[k] for k in keys(sol))

    map((d, s) -> copy!(d, s), variables(sol, n), variables(st))

    # for k in Val.(keys(sol))
    #     sol[k][div(n, step(sol))] = s[k]
    # end

    return st
end

function Base.copy!(sol::GeometricSolution, src, n::Int)
    @assert n ≤ ntime(sol)
    if mod(n, step(sol)) == 0
        i = div(n, step(sol))
        sol[i] = src
    end
    return sol
end

function arrays(sol::GeometricSolution)
    NamedTuple{keys(sol)}((Array(sol[k]) for k in keys(sol)))
    # map(s -> Array(s), sol)
end

function relative_maximum_error(sol::GeometricSolution, ref::GeometricSolution)
    @assert keys(sol) == keys(ref)
    NamedTuple{keys(sol)}(relative_maximum_error(s...)
    for s in zip(solutions(sol), solutions(ref)))
end
