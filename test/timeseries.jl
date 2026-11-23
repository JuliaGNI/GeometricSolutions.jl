using GeometricSolutions
import GeometricEquations
using OffsetArrays
using Random
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

@testset "$(rpad("Timeseries step count is public",80))" begin
    # Downstream packages count the steps of a run with `ntimesteps`, the rule of the time grid.
    @test Base.ispublic(GeometricSolutions, :ntimesteps)
    @test !Base.isexported(GeometricSolutions, :ntimesteps)
    @test GeometricSolutions.ntimesteps(0.0, 5λ, λ / 21) ==
          ntime(TimeSeries(0.0, 5λ, λ / 21))

    # The time grid counts its steps by the rule of the problem's step count, so a solution
    # holds as many steps as the problem takes.
    for T in (Float16, Float32, Float64), t₀ in (zero(T), T(1.5))

        wrong = 0
        for k in 1:200, m in 1:10

            Δt = T(λ) / k
            t₁ = t₀ + m * T(λ)
            wrong += GeometricSolutions.ntimesteps(t₀, t₁, Δt) !=
                     GeometricEquations.ntimesteps((t₀, t₁), Δt)
            wrong += GeometricSolutions.ntimesteps(t₁, t₀, -Δt) !=
                     GeometricEquations.ntimesteps((t₁, t₀), -Δt)
        end
        @test wrong == 0
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
    @test_throws r"^ArgumentError: time step cannot be zero" TimeSeries(10, 0.0)

    # A time step whose sign is opposite to that of `tend - tbegin` has no step count, also where
    # it is longer than the span.
    wrong_sign = r"^ArgumentError: time step .* does not point from"
    for (t₀, t₁, Δt) in ((0.0, 1.0, -0.3), (0.0, 1.0, -2.0),
        (1.0, 0.0, 0.1), (1.0f0, 0.0f0, 2.0f0))
        @test_throws wrong_sign TimeSeries(t₀, t₁, Δt)
        @test_throws wrong_sign GeometricSolutions.ntimesteps(t₀, t₁, Δt)
    end

    # An empty span has no steps, whatever the sign of the time step.
    @test ntime(TimeSeries(1.0, 1.0, 0.1)) == 0
    @test ntime(TimeSeries(1.0, 1.0, -0.1)) == 0
end

@testset "$(rpad("Timeseries step count in exact arithmetic",80))" begin
    # The counting rule in exact arithmetic on the binary values of the endpoints and the step.
    function exact_count(a::T, b::T, h::T) where {T}
        x₀, x₁, y = Rational{BigInt}(a), Rational{BigInt}(b), Rational{BigInt}(h)
        q = (x₁ - x₀) / y
        n = round(q)
        fills = abs(x₀ + n * y - x₁) ≤ 8 * Rational{BigInt}(eps(T)) * max(abs(x₀), abs(x₁))
        Int(fills ? n : floor(q))
    end

    # In Float16, `tend - tbegin = 1024.5` rounds to 1024, which is one step short.
    ts = TimeSeries(Float16(0.5), Float16(1025), Float16(0.5))
    @test ntime(ts) == 2049
    @test eltype(ts) == Float16
    @test ts[end] == Float16(1025)

    # Two Float32 spans whose rounded difference changes the count.
    for (a, b, h) in ((0.95651406f0, 8.94363f0, 0.001543704f0),
        (9.835699f0, 11.890267f0, 0.011607794f0))
        @test ntime(TimeSeries(a, b, h)) == exact_count(a, b, h)
    end

    Random.seed!(1)
    for T in (Float16, Float32, Float64)
        wrong = 0
        for _ in 1:10_000
            a = T(10 * rand())
            b = a + T(10 * rand())
            h = T(rand() / 10 + 1e-3)
            wrong += ntime(TimeSeries(a, b, h)) != exact_count(a, b, h)
        end
        @test wrong == 0
    end
end

@testset "$(rpad("Timeseries of exact time types",80))" begin
    # An integer time series keeps the time type in its grid.
    ts = TimeSeries(Int8(0), Int8(10), Int8(3))
    @test eltype(ts) == Int8
    @test ntime(ts) == 3
    @test ts[end] == Int8(9)

    # The span of integer endpoints near `typemax` does not wrap: one step of `3u` fits in `4u`.
    u = Int64(1) << 61
    ts = TimeSeries(-2u, 2u, 3u)
    @test ntime(ts) == 1
    @test ts[end] == u

    # A rational time series counts whole steps exactly, forward and backward.
    @test ntime(TimeSeries(0 // 1, 1 // 1, 3 // 10)) == 3
    @test ntime(TimeSeries(0 // 1, 1 // 1, 1 // 10)) == 10
    @test ntime(TimeSeries(1 // 1, 0 // 1, -1 // 10)) == 10
    zero_step = r"^ArgumentError: time step cannot be zero"
    @test_throws zero_step TimeSeries(0 // 1, 1 // 1, 0 // 1)
end
