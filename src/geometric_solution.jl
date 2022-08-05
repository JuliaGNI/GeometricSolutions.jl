
# _periodicity(::AT, periodicity::AT) where {AT <: AbstractArray} = periodicity
# _periodicity(q::AbstractArray, ::NullPeriodicity) = zero(q)

function _periodicity(s::NamedTuple, periodicity::NamedTuple)
    NamedTuple{keys(s)}(Tuple(haskey(periodicity, k) ? periodicity[k] : NullPeriodicity() for k in keys(s)))
end


"""
`GeometricSolution`: Solution of an ordinary differential equation

Contains all fields necessary to store the solution of an ODE.

### Fields

* `nt`: number of time steps to store
* `t`:  time steps
* `s`:  NamedTuple of DataSeries for each solution component
* `ntime`: number of time steps to compute
* `nsave`: store every nsave'th time step (default: 1)
* `counter`: counter for copied solution entries

### Constructors

```julia
GeometricSolution(problem; nsave=DEFAULT_NSAVE)
GeometricSolution(t::TimeSeries, q::DataSeries, ntimesteps)
```

The usual way to initialise a `Solution` is by passing an equation, which for
`GeometricSolution` has to be an [`ODE`](@ref) or [`SODE`](@ref), a time step `Δt`
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
    counter::Int

    function GeometricSolution(problem::GeometricProblem; step = 1)
        t = TimeSeries(tbegin(problem), tend(problem), tstep(problem))
        nstore = div(ntime(t), step)
        s = NamedTuple{keys(problem.ics)}(Tuple(DataSeries(x, nstore) for x in problem.ics))
        period = _periodicity(s, periodicity(problem))
        new{datatype(problem), timetype(problem), typeof(s), typeof(problem), typeof(period)}(t, s, problem, period, step, nstore, 0)
    end
end

@inline step(sol::GeometricSolution) = sol.step
@inline nstore(sol::GeometricSolution) = sol.nstore
@inline counter(sol::GeometricSolution) = sol.counter
@inline lastentry(sol::GeometricSolution) = sol.counter - 1

@inline GeometricBase.ntime(sol::GeometricSolution) = ntime(sol.t)
@inline GeometricBase.timesteps(sol::GeometricSolution) = sol.t
@inline GeometricBase.eachtimestep(sol::GeometricSolution) = eachtimestep(sol.t)
@inline GeometricBase.periodicity(sol::GeometricSolution) = sol.periodicity

function Base.setindex!(sol::GeometricSolution, s::NamedTuple, i::Int)
    @assert keys(sol.s) == keys(s)
    @assert i <= ntime(sol)

    if mod(i, step(sol)) == 0
        for k in keys(sol.s)
            sol.s[k][div(i, step(sol))] = s[k]
        end
        sol.counter += 1
    end
    
    return s
end
