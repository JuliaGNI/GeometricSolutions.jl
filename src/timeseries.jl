struct TimeSeries{T, N, OT <: OffsetVector{T}} <: AbstractVector{T}
    t::OT
    Δt::T

    function TimeSeries(ti::AbstractVector, Δt::T) where {T <: Real}
        n = length(ti) - 1
        t = OffsetVector(ti, 0:n)
        new{T, n, typeof(t)}(t, Δt)
    end
end

"""
    ntimesteps(tbegin::T, tend::T, Δt::T) where {T <: AbstractFloat}

The number of whole steps `Δt` from `tbegin` to `tend`, for floating-point times. The endpoints and
the step are first converted to `W = promote_type(T, Float64)`, so that `tend - tbegin` does not
round in a narrower type, and every expression below is evaluated in `W`. A count `n` whose end
`tbegin + n * Δt` lies within `8 * eps(T) * max(abs(tbegin), abs(tend))` of `tend` fills the span,
on either side of `tend`, so a quotient `(tend - tbegin) / Δt` that rounds just below an integer
loses no step. Any other count rounds down, and the last step ends before `tend`. A zero `Δt`
throws an `ArgumentError`.
"""
function ntimesteps(tbegin::T, tend::T, Δt::T) where {T <: AbstractFloat}
    iszero(Δt) && throw(ArgumentError("time step cannot be zero"))
    W = promote_type(T, Float64)
    t₀, t₁, h = W(tbegin), W(tend), W(Δt)
    q = (t₁ - t₀) / h
    n = round(q)
    fills = abs(t₀ + n * h - t₁) ≤ 8 * eps(T) * max(abs(t₀), abs(t₁))
    return Int(fills ? n : floor(q))
end

function TimeSeries(tbegin::T, tend::T, Δt::T) where {T <: AbstractFloat}
    n = ntimesteps(tbegin, tend, Δt)
    TimeSeries(range(tbegin; step = Δt, length = n + 1), Δt)
end

# Exact time types, such as integers and rationals, count their steps without rounding.
function TimeSeries(tbegin::T, tend::T, Δt::T) where {T <: Real}
    iszero(Δt) && throw(ArgumentError("time step cannot be zero"))
    TimeSeries(tbegin:Δt:tend, Δt)
end

function TimeSeries(n::Integer, Δt::T) where {T}
    TimeSeries(zero(Δt), n * Δt, Δt)
end

function TimeSeries(t::AbstractVector)
    @assert length(t) ≥ 2
    Δt = t[2] - t[1]
    return TimeSeries(t[begin], t[end], Δt)
end

@inline Base.parent(ts::TimeSeries) = ts.t
@inline Base.eltype(::TimeSeries{T}) where {T} = T
@inline Base.ndims(::TimeSeries) = 1

@inline Base.size(ts::TimeSeries, args...) = size(parent(ts), args...)

@inline Base.eachindex(ts::TimeSeries) = eachindex(parent(ts))
@inline Base.eachindex(ind::IndexCartesian, ts::TimeSeries) = eachindex(ind, parent(ts))
@inline Base.eachindex(ind::IndexLinear, ts::TimeSeries) = eachindex(ind, parent(ts))

@inline Base.firstindex(ts::TimeSeries, args...) = firstindex(parent(ts), args...)
@inline Base.lastindex(ts::TimeSeries, args...) = lastindex(parent(ts), args...)

@inline Base.axes(ts::TimeSeries, args...) = axes(parent(ts), args...)

@inline Base.strides(ts::TimeSeries) = strides(collect(parent(parent(ts))))

@inline Base.getindex(ts::TimeSeries, args...) = getindex(parent(ts), args...)

@inline GeometricBase.timespan(ts::TimeSeries) = (ts[begin], ts[end])
@inline GeometricBase.timestep(ts::TimeSeries) = ts.Δt

@inline GeometricBase.ntime(ts::TimeSeries{T, N}) where {T, N} = N
@inline GeometricBase.eachtimestep(ts::TimeSeries) = Base.OneTo(ntime(ts))

function Base.:(==)(ts1::TimeSeries, ts2::TimeSeries)
    (ts1.t == ts2.t && ts1.Δt == ts2.Δt)
end
function Base.:(==)(ts::TimeSeries{T1}, vec::AbstractVector{T2}) where {T1, T2}
    (T1 == T2 && collect(parent(parent(ts))) == vec)
end
Base.:(==)(vec::AbstractVector, ts::TimeSeries) = (ts == vec)
