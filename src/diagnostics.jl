
function relative_norm_error(sol, ref, p = 2)
    norm(sol .- ref, p) / norm(ref, p)
end

function maximum_error(sol, ref)
    maximum(abs.(sol .- ref))
end

function relative_maximum_error(sol, ref)
    maximum_error(sol, ref) / maximum(abs.(ref))
end

"""
Computes the difference of two DataSeries.

Arguments: `(x::DataSeries{DT}, y::DataSeries{DT})`

Returns a DataSeries similar to `x` and `y` holding the time series of the difference between `x` and `y`.
"""
function compute_difference(x::DataSeries{DT}, y::DataSeries{DT}) where {DT}
    @assert axes(x) == axes(y)

    e = zero(x)
    parent(e) .= parent(x) .- parent(y)

    return e
end

"""
Takes a ScalarDataSeries holding an invariant and computes the relative error `(inv(t)-inv(0))/inv(0)`.

Returns a ScalarDataSeries similar to the argument holding the time series of the relativ errors.
"""
function compute_relative_error(invds::ScalarDataSeries{T}) where {T}
    errds = zero(invds)
    for i in eachindex(errds, invds)
        errds[i] = (invds[i] - invds[0]) / invds[0]
    end
    return errds
end

"""
Compute an invariant for the solution of an ODE or DAE system.

Arguments: `(t::TimeSeries, q::DataSeries{T}, invariant::Base.Callable)`

The `invariant` functions needs to take two arguments `(t,q)` and return the
corresponding value of the invariant.

Returns a ScalarDataSeries holding the time series of the invariant.
"""
function compute_invariant(
        t::TimeSeries, q::DataSeries{T}, invariant::Base.Callable) where {T}
    invds = DataSeries(T, ntime(q))
    try
        for i in eachindex(invds)
            invds[i] = invariant(t[i], q[i])
        end
    catch ex
        if isa(ex, DomainError)
            @warn("DOMAIN ERROR: Invariant diagnostics crashed.")
        else
            throw(ex)
        end
    end
    return invds
end

"""
Compute an invariant for the solution of a partitioned ODE or DAE system.

Arguments: `(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, invariant::Base.Callable)`

The `invariant` functions needs to take three arguments `(t,q,p)` and return the
corresponding value of the invariant.

Returns a ScalarDataSeries holding the time series of the invariant.
"""
function compute_invariant(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T},
        invariant::Base.Callable) where {T}
    invds = DataSeries(T, ntime(q))
    try
        for i in eachindex(invds)
            invds[i] = invariant(t[i], q[i], p[i])
        end
    catch ex
        if isa(ex, DomainError)
            @warn("DOMAIN ERROR: Invariant diagnostics crashed.")
        else
            throw(ex)
        end
    end
    return invds
end

"""
Compute the relative error of an invariant for the solution of an ODE or DAE system.

Arguments: `(t::TimeSeries, q::DataSeries{T}, invariant::Base.Callable)`

The `invariant` functions needs to take two arguments `(t,q)` and return the
corresponding value of the invariant.

Returns a tuple of two 1d DataSeries holding the time series of the invariant and the relativ error, respectively.
"""
function compute_invariant_error(t::TimeSeries, q::DataSeries, invariant::Base.Callable)
    invds = compute_invariant(t, q, invariant)
    errds = compute_relative_error(invds)
    (invds, errds)
end

"""
Compute the relative error of an invariant for the solution of a partitioned ODE or DAE system.

Arguments: `(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, invariant::Base.Callable)`

The `invariant` functions needs to take three arguments `(t,q,p)` and return the
corresponding value of the invariant.

Returns a tuple of two ScalarDataSeries holding the time series of the invariant and the relativ error, respectively.
"""
function compute_invariant_error(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T},
        invariant::Base.Callable) where {T}
    invds = compute_invariant(t, q, p, invariant)
    errds = compute_relative_error(invds)
    (invds, errds)
end

"""
Computes the drift in an invariant error.

Arguments: `(t::TimeSeries, invariant_error::DataSeries{T,1}, interval_length=100)`

The time series of the solution is split into intervals of `interval_length` time steps.
In each interval, the maximum of the absolute value of the invariant error is determined.
Returns a tuple of a TimeSeries that holds the centers of all intervals and a ScalarDataSeries
that holds the maxima.

This is useful to detect drifts in invariants that are not preserved exactly but whose error
is oscillating such as the energy error of Hamiltonian systems with symplectic integrators.
"""
function compute_error_drift(t::TimeSeries, invariant_error::ScalarDataSeries{T},
        interval_length = 100) where {T}
    @assert ntime(t) == ntime(invariant_error)
    @assert mod(t.n, interval_length) == 0

    nint = div(ntime(invariant_error), interval_length)
    Idrift = DataSeries(T, nint)
    Tdrift = [t[0]]

    for i in 1:nint
        i1 = interval_length * (i - 1) + 1
        i2 = interval_length * i
        Idrift[i] = maximum(abs.(invariant_error[i1:i2]))
        push!(Tdrift, (t[i1] + t[i2]) / 2)
    end

    (TimeSeries(Tdrift, (t[end] - t[begin]) / nint), Idrift)
end
