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

The usual way to initialise a `Solution` is by passing a [`GeometricProblem`](@ref), which
can for example be an [`ODEProblem`](@ref) or [`PODEProblem`](@ref).
The optional parameter `step` determines the intervals for storing the solution,
i.e., if `step > 1` only every `step`'th solution is actually stored.

"""
mutable struct GeometricSolution{dType, tType, dsType, probType, perType} <:
               AbstractSolution{dType, tType}
    t::TimeSeries{tType}
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
            typeof(s), typeof(problem), typeof(period)}(
            t, s, problem, period, step, nstore, 0)
        sol[0] = initial_conditions(problem)
        return sol
    end
end

function GeometricSolution(problem::GeometricProblem, args...)
    t = TimeSeries(tbegin(problem), tend(problem), tstep(problem))
    GeometricSolution(t, problem, args...)
end

@inline Base.keys(sol::GeometricSolution) = keys(sol.s)
@inline Base.step(sol::GeometricSolution) = sol.step
@inline nstore(sol::GeometricSolution) = sol.nstore

@inline GeometricBase.datatype(sol::GeometricSolution{DT, TT}) where {DT, TT} = DT
@inline GeometricBase.timetype(sol::GeometricSolution{DT, TT}) where {DT, TT} = TT

@inline GeometricBase.tspan(sol::GeometricSolution) = tspan(sol.t)
@inline GeometricBase.tstep(sol::GeometricSolution) = tstep(sol.t)

@inline GeometricBase.ntime(sol::GeometricSolution) = ntime(sol.t)
@inline GeometricBase.timesteps(sol::GeometricSolution) = sol.t
@inline GeometricBase.eachtimestep(sol::GeometricSolution) = eachtimestep(sol.t)
@inline GeometricBase.periodicity(sol::GeometricSolution) = sol.periodicity

@inline function Base.hasproperty(
        ::GeometricSolution{DT, TT, dsType}, s::Symbol) where {DT, TT, dsType}
    hasfield(dsType, s) || hasfield(GeometricSolution, s)
end

@inline function Base.getproperty(
        sol::GeometricSolution{DT, TT, dsType}, s::Symbol) where {DT, TT, dsType}
    if hasfield(dsType, s)
        return getfield(sol, :s)[s]
    else
        return getfield(sol, s)
    end
end

function Base.getindex(sol::GeometricSolution, s::Symbol)
    getfield(sol, :s)[s]
end

function Base.getindex(sol::GeometricSolution, n::Int)
    @assert n ≤ ntime(sol)
    NamedTuple{(:t, keys(sol.s)...)}((
        timesteps(sol)[n], (sol.s[k][n] for k in keys(sol.s))...))
end

function Base.setindex!(sol::GeometricSolution, s::NamedTuple, n::Int)
    # @assert keys(sol.s) ⊆ keys(s)
    @assert n ≤ ntime(sol)

    if mod(n, step(sol)) == 0
        for k in keys(sol.s) ∩ keys(s)
            sol.s[k][div(n, step(sol))] = s[k]
        end
    end

    return s
end

function arrays(sol::GeometricSolution)
    NamedTuple{keys(sol)}((Array(sol[k]) for k in keys(sol)))
end

function relative_maximum_error(sol::GeometricSolution, ref::GeometricSolution)
    @assert keys(sol.s) == keys(ref.s)
    NamedTuple{keys(sol.s)}(relative_maximum_error(s...) for s in zip(sol.s, ref.s))
end
