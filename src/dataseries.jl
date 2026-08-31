
const AbstractData{DT} = Union{DT, AbstractArray{DT}}

struct DataSeries{DT, AT <: AbstractData{DT}} <: AbstractVector{AT}
    d::OffsetVector{AT, Vector{AT}}

    function DataSeries(v::AbstractVector{AT}) where {DT, AT <: AbstractData{DT}}
        ntime = length(v) - 1
        new{DT, AT}(OffsetVector(v, 0:ntime))
    end
end

DataSeries(x::AbstractData, ntime::Integer) = DataSeries([zero(x) for _ in 0:ntime])
DataSeries(ds::DataSeries) = DataSeries(copy(parent(ds)))

const ScalarDataSeries{DT} = DataSeries{DT, DT}

@inline Base.parent(ds::DataSeries) = ds.d
@inline Base.eltype(::DataSeries{DT}) where {DT} = DT
@inline GeometricBase.arrtype(::DataSeries{DT, AT}) where {DT, AT} = AT
@inline Base.ndims(::DataSeries) = 1

@inline Base.size(ds::DataSeries, args...) = size(parent(ds), args...)

@inline Base.eachindex(ds::DataSeries) = eachindex(parent(ds))
@inline Base.eachindex(ind::IndexCartesian, ds::DataSeries) = eachindex(ind, parent(ds))
@inline Base.eachindex(ind::IndexLinear, ds::DataSeries) = eachindex(ind, parent(ds))

@inline Base.firstindex(ds::DataSeries, args...) = firstindex(parent(ds), args...)
@inline Base.lastindex(ds::DataSeries, args...) = lastindex(parent(ds), args...)

@inline Base.axes(ds::DataSeries, args...) = axes(parent(ds), args...)

@inline Base.strides(ds::DataSeries) = strides(parent(ds))

@inline Base.getindex(ds::DataSeries, args...) = getindex(parent(ds), args...)
@inline Base.setindex!(ds::DataSeries, args...) = setindex!(parent(ds), args...)

@inline Base.getindex(
    ds::DataSeries, ind::Union{
        Int, IndexLinear, AbstractRange}, i, args...) = getindex(
    parent(ds)[ind], i, args...)
@inline Base.setindex!(ds::DataSeries, x::AbstractArray, ind::Union{
    Int, IndexLinear, AbstractRange}) = copy!(
    parent(ds)[ind], x)

@inline function Base.getindex(ds::DataSeries, ::Colon, j::Union{Int, CartesianIndex})
    OffsetArray([ds[i][j] for i in eachindex(ds)], eachindex(ds))
end

@inline GeometricBase.ntime(ds::DataSeries) = lastindex(ds)

Base.:(==)(ds1::DataSeries, ds2::DataSeries) = parent(ds1) == parent(ds2)

GeometricBase.reset!(ds::DataSeries) = ds[begin] = ds[end]

function Base.show(
        io::IO, ds::DS) where {DT, AT <: AbstractArray{DT}, DS <: DataSeries{DT, AT}}
    print(io, "$(DS) with data type ", DT, " and array type ", AT, "\n")
    print(io, parent(parent(ds)))
end

function Base.show(io::IO, ds::DS) where {DT, DS <: DataSeries{DT, DT}}
    print(io, "$(DS) with data type ", DT, "\n")
    print(io, parent(parent(ds)))
end

function Base.zero(ds::DataSeries)
    DataSeries(zero(ds[begin]), ntime(ds))
end

_vec(x::AbstractArray) = vec(x)
_vec(x::Number) = [x]

function Base.Array(ds::DataSeries)
    nt = ntime(ds)
    nd = length(_vec(ds[begin]))
    z = zeros(eltype(ds), nd, nt + 1)
    for (i, elem) in enumerate(ds)
        z[:, i] .= _vec(elem)
    end
    return z
end

# Relative maximum error of a single point. If the point-wise absolute error
# vanishes (e.g. an all-zero reference point that is matched exactly) the relative
# error is zero, avoiding the 0/0 = NaN that dividing by a zero norm would produce.
function _relative_maximum_error(sol, ref)
    err = maximum_error(sol, ref)
    iszero(err) ? err : err / maximum(abs, ref)
end

"""
Computes the maximum over time of the relative maximum error between two DataSeries.

Arguments: `(ds::DataSeries, ref::DataSeries)`

For each time step the relative error `maximum_error(ds[i], ref[i]) / maximum(abs.(ref[i]))`
is evaluated, i.e. each point is normalized by its own maximum absolute value. Returns the
maximum of these per-step errors. An all-zero reference point matched exactly has a vanishing
point-wise error and therefore contributes `0` rather than `0/0 = NaN`.
"""
function relative_maximum_error(ds::DataSeries, ref::DataSeries)
    @assert axes(ds.d) == axes(ref.d)
    maximum(_relative_maximum_error.(ds.d, ref.d))
end

# Relative p-norm error of a single point. As above, a vanishing absolute norm
# error (e.g. an all-zero reference point matched exactly) yields zero instead of
# the 0/0 = NaN that dividing by a zero norm would produce.
function _relative_norm_error(sol, ref, p = 2)
    err = norm(sol .- ref, p)
    iszero(err) ? err : err / norm(ref, p)
end

"""
Computes the time trace of the relative `p`-norm error between two DataSeries.

Arguments: `(ds::DataSeries, ref::DataSeries, p = 2)`

For each time step the relative error `norm(ds[i] - ref[i], p) / norm(ref[i], p)`
is evaluated. Returns a `ScalarDataSeries` holding this time series.
"""
function relative_norm_error(ds::DataSeries, ref::DataSeries, p = 2)
    @assert axes(ds.d) == axes(ref.d)
    DataSeries([_relative_norm_error(d, r, p) for (d, r) in zip(ds.d, ref.d)])
end
