"""
`GeometricSolution`: Solution of an ordinary differential equation

Contains all fields necessary to store the solution of an ODE.

### Fields

* `nt`: number of time steps to store
* `t`:  time steps
* `s`:  NamedTuple of DataSeries for each solution component
* `ntime`: number of time steps to compute
* `nsave`: store every nsave'th time step (default: 1)

### Constructors

```julia
GeometricSolution(problem; nsave=DEFAULT_NSAVE)
GeometricSolution(t::TimeSeries, q::DataSeries, ntimesteps)
```

The usual way to initialise a `Solution` is by passing an equation, which for
`GeometricSolution` has to be an [`ODEProblem`](@ref) or [`SODEProblem`](@ref), a time step `Δt`
and the number of time steps `ntimesteps`. The optional parameters `nsave` and
`nwrite` determine the intervals for storing the solution and writing to file,
i.e., if `nsave > 1` only every `nsave`'th solution is actually stored, and
every `nwrite`'th time step the solution is stored to disk.

The other constructors, either passing a `TimeSeries` and a `DataSeries` or a
filename are used to read data from previous simulations.

"""
mutable struct GeometricSolution{dType, tType, dsType, probType, perType} <: AbstractSolution{dType,tType}
    t::TimeSeries{tType}
    s::dsType

    problem::probType
    periodicity::perType

    step::Int
    nstore::Int
    offset::Int

    function GeometricSolution(problem::GeometricProblem; step = 1)
        t = TimeSeries(tbegin(problem), tend(problem), tstep(problem))
        nstore = div(ntime(t), step)
        s = NamedTuple{keys(problem.ics)}(Tuple(DataSeries(x, nstore) for x in problem.ics))
        period = _periodicity(s, periodicity(problem))
        sol = new{datatype(problem), timetype(problem), typeof(s), typeof(problem), typeof(period)}(t, s, problem, period, step, nstore, 0)
        sol[0] = initial_conditions(problem)
        return sol
    end
end

@inline Base.step(sol::GeometricSolution) = sol.step
@inline nstore(sol::GeometricSolution) = sol.nstore

@inline GeometricBase.tspan(sol::GeometricSolution) = tspan(sol.t)
@inline GeometricBase.tstep(sol::GeometricSolution) = tstep(sol.t)

@inline GeometricBase.ntime(sol::GeometricSolution) = ntime(sol.t)
@inline GeometricBase.timesteps(sol::GeometricSolution) = sol.t
@inline GeometricBase.eachtimestep(sol::GeometricSolution) = eachtimestep(sol.t)
@inline GeometricBase.periodicity(sol::GeometricSolution) = sol.periodicity

@inline function Base.hasproperty(::GeometricSolution{DT,TT,dsType}, s::Symbol) where {DT,TT,dsType}
    hasfield(dsType, s) || hasfield(GeometricSolution, s)
end

@inline function Base.getproperty(sol::GeometricSolution{DT,TT,dsType}, s::Symbol) where {DT,TT,dsType}
    if hasfield(dsType, s)
        return getfield(sol, :s)[s]
    else
        return getfield(sol, s)
    end
end

function Base.getindex(sol::GeometricSolution, n::Int)
    @assert n ≤ ntime(sol)
    NamedTuple{(:t, keys(sol.s)...)}((timesteps(sol)[n], (sol.s[k][n] for k in keys(sol.s))...))
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
