struct TimeSeries{T, OT <: OffsetVector{T}} <: AbstractVector{T}
    n::Int
    t::OT
    Δt::T

    function TimeSeries(ti::AbstractVector, Δt::T) where {T <: Real}
        n = length(ti)-1
        t = OffsetVector(ti, 0:n)
        new{T, typeof(t)}(n, t, Δt)
    end
end

function TimeSeries(tbegin::T, tend::T, Δt::T) where {T <: Real}
    TimeSeries(tbegin:Δt:tend, Δt)
end

function TimeSeries(n::Integer, Δt::T) where {T}
    TimeSeries(zero(Δt), n*Δt, Δt)
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

@inline GeometricBase.ntime(ts::TimeSeries) = ts.n
@inline GeometricBase.eachtimestep(ts::TimeSeries) = Base.OneTo(ntime(ts))

Base.:(==)(ts1::TimeSeries, ts2::TimeSeries) = (ts1.n == ts2.n && ts1.t == ts2.t && ts1.Δt == ts2.Δt)
Base.:(==)(ts::TimeSeries{T1}, vec::AbstractVector{T2}) where {T1,T2} = (T1 == T2 && collect(parent(parent(ts))) == vec)
Base.:(==)(vec::AbstractVector, ts::TimeSeries) = (ts == vec)
