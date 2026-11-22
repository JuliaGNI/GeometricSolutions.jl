using GeometricSolutions
using OffsetArrays
using Test

@testset "$(rpad("Timeseries",80))" begin
    nt = 10
    ttbeg = 0.0
    ttend = 1.0
    Δt = 0.1
    ti = ttbeg:Δt:ttend

    ts = TimeSeries(ti, Δt)
    ts1 = TimeSeries(ttbeg, ttend, Δt)
    ts2 = TimeSeries(nt, Δt)
    ts3 = TimeSeries(collect(ti))

    @test ts == ts1 == ts2 == ts3

    @test typeof(ts) <: AbstractVector
    @test typeof(parent(ts)) <: OffsetVector
    @test parent(parent(ts)) == ti
    @test eltype(ts) == typeof(Δt)
    @test ndims(ts) == 1
    @test ntime(ts) == nt

    @test firstindex(ts) == 0
    @test firstindex(ts, 1) == firstindex(ts.t, 1)
    @test firstindex(ts, 2) == 1
    @test lastindex(ts) == nt
    @test lastindex(ts, 1) == lastindex(ts.t, 1)
    @test lastindex(ts, 2) == 1
    @test strides(ts) == (1,)
    @test stride(ts, 1) == 1
    @test axes(ts) == (0:nt,)
    @test axes(ts, 1) == 0:nt
    @test axes(ts, 2) == 1:1
    @test size(ts) == (nt+1,)
    @test size(ts.t) == (nt+1,)
    @test size(ts, 1) == nt+1
    @test size(ts, 2) == 1

    @test eachindex(ts) == 0:nt
    @test eachindex(IndexLinear(), ts) == 0:nt
    @test eachindex(IndexCartesian(), ts) == CartesianIndices((0:nt,))

    @test initialtime(ts) == ts[0] == ts[begin] == ttbeg
    @test finaltime(ts) == ts[nt] == ts[end] == ttend
    @test timespan(ts) == (ttbeg, ttend)
    @test timestep(ts) == Δt

    @test ts == vec(ti)
    @test vec(ti) == ts
end

# The wavelength of a 100 GHz wave in centimetres. Steps `λ/k` over spans `mλ` give quotients
# `(tend - tbegin)/Δt` that round just below or above the integer `m * k`; over the spans from zero
# the float range `tbegin:Δt:tend` drops the last step of 183 Float64 and 92 Float32 spans.
const λ = 2π / (2π * 100e9 / 2.99792458e10)

@testset "$(rpad("Timeseries of a span of whole steps",80))" begin
    # A span of the sweep that the float range shortens: 104 steps of 105.
    ts = TimeSeries(0.0, 5λ, λ / 21)
    @test ntime(ts) == 105
    @test isapprox(ts[end], 5λ; rtol = 8eps())

    for T in (Float64, Float32), t₀ in (zero(T), T(1.5))

        wrong = 0
        for k in 1:200, m in 1:10

            Δt = T(λ) / k
            t₁ = t₀ + m * T(λ)
            ts = TimeSeries(t₀, t₁, Δt)
            tol = 8 * eps(T) * max(abs(t₀), abs(t₁))
            wrong += !(ntime(ts) == m * k && ts[begin] == t₀ && abs(ts[end] - t₁) ≤ tol)
        end
        @test wrong == 0
    end

    # A time series from a step count holds that many steps of `Δt` from zero.
    for T in (Float64, Float32)
        @test all(
            ntime(TimeSeries(m * k, T(λ) / k)) == m * k &&
                TimeSeries(m * k, T(λ) / k)[1] == T(λ) / k for k in 1:200, m in 1:10)
    end
end

@testset "$(rpad("Timeseries of a span of partial steps",80))" begin
    # A span that is no whole number of steps ends at the last step before `tend`: the count
    # rounds down, also where it is nearer the next integer (1/0.35 ≈ 2.86).
    for (Δt, n) in ((0.3, 3), (0.35, 2))
        ts = TimeSeries(0.0, 1.0, Δt)
        @test ntime(ts) == n
        @test ts[end] ≈ n * Δt
        @test ts[end] < 1.0
    end

    # A negative time step runs backward from `tbegin`.
    ts = TimeSeries(1.0, 0.0, -0.1)
    @test ntime(ts) == 10
    @test ts[begin] == 1.0
    @test isapprox(ts[end], 0.0; atol = 8eps())

    # A zero time step has no step count.
    @test_throws r"^ArgumentError: time step cannot be zero" TimeSeries(0.0, 1.0, 0.0)
end
