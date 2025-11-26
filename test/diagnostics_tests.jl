using GeometricSolutions
using LinearAlgebra
using OffsetArrays
using Test

dt = Float64
nt = 10
nd = 2

invariant(t, q) = norm(q)
invariant(t, q, p) = norm(q) + norm(p)

q = OffsetArray([rand(dt, nd) for _ in 1:(nt + 1)], 0:nt)
p = OffsetArray([rand(dt, nd) for _ in 1:(nt + 1)], 0:nt)
p = OffsetArray([ones(dt, nd) for _ in 1:(nt + 1)], 0:nt)

# Error and Differences
qs = DataSeries(q)
ps = DataSeries(p)
df = compute_difference(qs, ps)

@test df == DataSeries(q .- p)

# Invariants for ODE solutions
ts = TimeSeries(0.0, 1.0, 1.0 / nt)
qs = DataSeries(q)
ti = compute_invariant(ts, qs, invariant)
tie, te = compute_invariant_error(ts, qs, invariant)
tdt, td = compute_error_drift(ts, te, 1)

@test typeof(ti) <: ScalarDataSeries{dt}
@test typeof(te) <: ScalarDataSeries{dt}
@test typeof(tie) <: ScalarDataSeries{dt}
@test typeof(td) <: ScalarDataSeries{dt}

@test parent(ti) == invariant.(parent(ts), parent(qs))
@test parent(ti) == parent(tie)
@test all(te .== (parent(ti) .- ti[begin]) ./ ti[begin])
@test parent(ts) == parent(tdt)
@test parent(td) == abs.(parent(te))

tdt, td = compute_error_drift(ts, te, 2)

@test all(parent(parent(tdt)) .≈ vcat(0.0, collect(0.15:0.2:0.95)))

# Invariants for PODE solutions
ts = TimeSeries(0.0, 1.0, 1.0 / nt)
qs = DataSeries(q)
ps = DataSeries(p)
ti = compute_invariant(ts, qs, ps, invariant)
tie, te = compute_invariant_error(ts, qs, ps, invariant)
tdt, td = compute_error_drift(ts, te, 1)

@test typeof(ti) <: ScalarDataSeries{dt}
@test typeof(te) <: ScalarDataSeries{dt}
@test typeof(tie) <: ScalarDataSeries{dt}
@test typeof(td) <: ScalarDataSeries{dt}

@test parent(ti) == invariant.(parent(ts), parent(qs), parent(ps))
@test parent(ti) == parent(tie)
@test all(te .== (parent(ti) .- ti[begin]) ./ ti[begin])
@test parent(ts) == parent(tdt)
@test parent(td) == abs.(parent(te))

tdt, td = compute_error_drift(ts, te, 2)

@test all(parent(parent(tdt)) .≈ vcat(0.0, collect(0.15:0.2:0.95)))
