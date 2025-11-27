using GeometricEquations
using GeometricSolutions
using LinearAlgebra
using OffsetArrays
using Test

dt = Float64
nt = 10
nd = 2

invariant(t, q, params) = norm(q)
invariant(t, q, p, params) = norm(q) + norm(p)

q = OffsetArray([rand(dt, nd) for _ in 1:(nt + 1)], 0:nt)
p = OffsetArray([rand(dt, nd) for _ in 1:(nt + 1)], 0:nt)
p = OffsetArray([ones(dt, nd) for _ in 1:(nt + 1)], 0:nt)
v = OffsetArray([ones(dt, nd) for _ in 1:(nt + 1)], 0:nt)
params = OffsetArray([NamedTuple() for _ in 1:(nt + 1)], 0:nt)

# Error and Differences
qs = DataSeries(q)
ps = DataSeries(p)
df = compute_difference(qs, ps)

@test df == DataSeries(q .- p)

# Invariants for ODE solutions
ts = TimeSeries(0.0, 1.0, 1.0 / nt)
qs = DataSeries(q)
ti = compute_invariant(ts, qs, NamedTuple(), invariant)
tie, te = compute_invariant_error(ts, qs, NamedTuple(), invariant)
tdt, td = compute_error_drift(ts, te, 1)

@test typeof(ti) <: ScalarDataSeries{dt}
@test typeof(te) <: ScalarDataSeries{dt}
@test typeof(tie) <: ScalarDataSeries{dt}
@test typeof(td) <: ScalarDataSeries{dt}

@test parent(ti) == invariant.(parent(ts), parent(qs), params)
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
ti = compute_invariant(ts, qs, ps, NamedTuple(), invariant)
tie, te = compute_invariant_error(ts, qs, ps, NamedTuple(), invariant)
tdt, td = compute_error_drift(ts, te, 1)

@test typeof(ti) <: ScalarDataSeries{dt}
@test typeof(te) <: ScalarDataSeries{dt}
@test typeof(tie) <: ScalarDataSeries{dt}
@test typeof(td) <: ScalarDataSeries{dt}

@test parent(ti) == invariant.(parent(ts), parent(qs), parent(ps), params)
@test parent(ti) == parent(tie)
@test all(te .== (parent(ti) .- ti[begin]) ./ ti[begin])
@test parent(ts) == parent(tdt)
@test parent(td) == abs.(parent(te))

tdt, td = compute_error_drift(ts, te, 2)

@test all(parent(parent(tdt)) .≈ vcat(0.0, collect(0.15:0.2:0.95)))

# Symplectic One-form
iode_ϑ! = (ϑ, t, q, v, params) -> ϑ .= [q[1]^2, q[2]^3]
iode_f! = (f, t, q, v, params) -> f .= [q[2], q[1]]
iode_g! = (g, t, q, v, λ, params) -> g .= [v[1], v[2]]
iode = IODEProblem(
    iode_ϑ!, iode_f!, iode_g!, (0.0, 1.0), 0.1, (q = [1.0, 0.0], p = [0.0, 1.0]))

sol = GeometricSolution(iode)

for i in eachindex(q)
    sol[i].q .= q[i]
    sol[i].p .= [q[i][1]^2, q[i][2]^3]
end

ϑs = compute_one_form(sol)
ϑc = OffsetArray([[_q[1]^2, _q[2]^3] for _q in q], 0:nt)

@test parent(ϑs) == ϑc

ϑe = compute_momentum_error(sol)

@test parent(ϑe) == OffsetArray([zeros(dt, nd) for _ in 1:(nt + 1)], 0:nt)
