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


# Norm and maximum errors

# relative_norm_error operates via `norm`, which recurses into arrays, so it
# accepts both plain numeric arrays and arrays of vector-valued points (one vector
# per time step). It does not support DataSeries.

# Identity: solution equals reference => error vanishes
nref = OffsetArray([dt[1.0, 2.0], dt[3.0, 4.0], dt[5.0, 6.0]], 0:2)

@test relative_norm_error(nref, nref) == 0

# Known nonzero value (differences chosen for an exact, checkable result)
nsol = OffsetArray([dt[1.0, 2.0], dt[3.0, 5.0], dt[7.0, 6.0]], 0:2)

@test relative_norm_error(nsol, nref) ≈ sqrt(5) / sqrt(91)

# p-norm argument (numeric vectors give the standard p-norms)
psol = dt[1.0, 2.0, 3.0]
pref = dt[1.0, 4.0, 3.0]

@test relative_norm_error(psol, pref) ≈ 2.0 / sqrt(26)     # default p = 2
@test relative_norm_error(psol, pref, 1) == 2.0 / 8.0      # p = 1
@test relative_norm_error(psol, pref, Inf) == 2.0 / 4.0    # p = ∞

# Solution points where all components are zero: some points of both sol and ref
# are the zero vector, but ref is not entirely zero, so the denominator stays nonzero.
znref = OffsetArray([dt[0.0, 0.0], dt[1.0, 2.0], dt[0.0, 0.0], dt[3.0, 4.0]], 0:3)
znsol = OffsetArray([dt[0.0, 0.0], dt[1.5, 2.0], dt[0.0, 0.0], dt[3.0, 4.0]], 0:3)

@test relative_norm_error(znsol, znref) ≈ 0.5 / sqrt(30)

# Zero points only in sol (ref has none): still no error, correct value
znref2 = OffsetArray([dt[1.0, 1.0], dt[2.0, 2.0], dt[4.0, 4.0]], 0:2)
znsol2 = OffsetArray([dt[0.0, 0.0], dt[2.0, 2.0], dt[4.0, 4.0]], 0:2)

@test relative_norm_error(znsol2, znref2) ≈ sqrt(2) / sqrt(42)


# maximum_error and relative_maximum_error operate via `abs.(sol .- ref)`, so the
# generic methods take plain numeric arrays (or a ScalarDataSeries); each element
# is treated as a scalar component.

# Identity: error vanishes
va = dt[1.0, 2.0, 3.0]

@test maximum_error(va, va) == 0
@test relative_maximum_error(va, va) == 0

# Known nonzero values
@test maximum_error(psol, pref) == 2.0
@test relative_maximum_error(psol, pref) == 2.0 / 4.0

# A matrix argument works the same way
@test maximum_error(dt[1 2; 3 4], dt[1 2; 3 6]) == 2.0
@test relative_maximum_error(dt[1 2; 3 4], dt[1 2; 3 6]) == 2.0 / 6.0

# Components equal to zero: some entries of sol and ref are zero, ref is not all zero
vzref = dt[0.0, 2.0, 0.0, 5.0]
vzsol = dt[0.0, 3.0, 0.0, 5.0]

@test maximum_error(vzsol, vzref) == 1.0                   # zero components contribute nothing
@test relative_maximum_error(vzsol, vzref) == 1.0 / 5.0    # divided by max abs component of ref

# ScalarDataSeries with points equal to zero. maximum_error uses the generic method,
# while relative_maximum_error dispatches to the DataSeries method that normalizes
# each point individually; a zero point matched exactly contributes 0 (not 0/0 = NaN).
sref = DataSeries(OffsetArray(dt[0.0, 2.0, 0.0, 3.0], 0:3))
ssol = DataSeries(OffsetArray(dt[0.0, 1.0, 0.0, 3.0], 0:3))

@test maximum_error(ssol, sref) == 1.0
@test relative_maximum_error(ssol, sref) == 1.0 / 2.0
@test !isnan(relative_maximum_error(ssol, sref))


# relative_maximum_error for a DataSeries of vector-valued points dispatches to the
# dedicated DataSeries method (src/dataseries.jl), which normalizes each point by
# its own maximum and takes the maximum of the per-point relative errors.
vref = OffsetArray([dt[1.0, 2.0], dt[3.0, 4.0], dt[5.0, 6.0]], 0:2)
vsol = OffsetArray([dt[1.0, 2.0], dt[3.0, 5.0], dt[7.0, 6.0]], 0:2)

@test relative_maximum_error(DataSeries(vsol), DataSeries(vref)) == 2.0 / 6.0

# All-zero solution points: a zero reference point matched exactly must contribute 0
# rather than 0/0 = NaN.
vzr = OffsetArray([dt[0.0, 0.0], dt[1.0, 2.0], dt[0.0, 0.0], dt[3.0, 4.0]], 0:3)
vzs = OffsetArray([dt[0.0, 0.0], dt[1.5, 2.0], dt[0.0, 0.0], dt[3.0, 4.0]], 0:3)

@test relative_maximum_error(DataSeries(vzs), DataSeries(vzr)) == 0.5 / 2.0
@test !isnan(relative_maximum_error(DataSeries(vzs), DataSeries(vzr)))


# relative_norm_error for two DataSeries returns a ScalarDataSeries holding the
# time trace of the per-step relative p-norm error.
rne = relative_norm_error(DataSeries(vsol), DataSeries(vref))

@test typeof(rne) <: ScalarDataSeries{dt}
@test axes(rne) == axes(DataSeries(vref))
@test rne[0] == 0.0                                        # sol == ref at step 0
@test rne[1] ≈ 1.0 / 5.0                                   # norm([0,1]) / norm([3,4])
@test rne[2] ≈ 2.0 / sqrt(61)                              # norm([2,0]) / norm([5,6])

# The p-norm argument is forwarded to each step
rne1 = relative_norm_error(DataSeries(vsol), DataSeries(vref), 1)

@test typeof(rne1) <: ScalarDataSeries{dt}
@test rne1[1] == 1.0 / 7.0                                 # norm([0,1],1) / norm([3,4],1)

# All-zero reference points matched exactly contribute 0, not 0/0 = NaN
rnez = relative_norm_error(DataSeries(vzs), DataSeries(vzr))

@test typeof(rnez) <: ScalarDataSeries{dt}
@test rnez[0] == 0.0
@test rnez[1] ≈ 0.5 / sqrt(5)
@test rnez[2] == 0.0
@test rnez[3] == 0.0
@test all(!isnan, parent(rnez))

# Also works for a ScalarDataSeries (scalar-valued points)
srne = relative_norm_error(ssol, sref)

@test typeof(srne) <: ScalarDataSeries{dt}
@test srne[1] == 1.0 / 2.0
@test all(!isnan, parent(srne))
