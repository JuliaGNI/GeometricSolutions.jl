
const AbstractData{DT} = Union{DT, AbstractArray{DT}}

struct DataSeries{DT, AT <: AbstractData{DT}} <: AbstractVector{AT}
    d::OffsetVector{AT, Vector{AT}}

    function DataSeries(v::AbstractVector{AT}) where {DT <: Real, AT <: AbstractData{DT}}
        ntime = length(v) - 1
        new{DT,AT}(OffsetVector(v, 0:ntime))
    end
end

DataSeries(x::AbstractData, ntime::Integer) = DataSeries([zero(x) for _ in 0:ntime])
DataSeries(ds::DataSeries) = DataSeries(copy(parent(ds)))

const ScalarDataSeries{DT} = DataSeries{DT,DT}


@inline Base.parent(ds::DataSeries) = ds.d
@inline Base.eltype(::DataSeries{DT}) where {DT} = DT
@inline arrtype(::DataSeries{DT,AT}) where {DT,AT} = AT
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

@inline Base.getindex(ds::DataSeries, ind::Union{Int,IndexLinear}, i, args...) = getindex(parent(ds)[ind], i, args...)

@inline GeometricBase.ntime(ds::DataSeries) = lastindex(ds)

Base.:(==)(ds1::DataSeries, ds2::DataSeries) = parent(ds1) == parent(ds2)

GeometricBase.reset!(ds::DataSeries) = ds[begin] = ds[end]

function Base.show(io::IO, ds::DS) where {DT, AT <: AbstractArray{DT}, DS <: DataSeries{DT,AT}}
    print(io, "$(DS) with data type ", DT, " and array type ", AT, "\n")
    print(io, parent(parent(ds)))
end

function Base.show(io::IO, ds::DS) where {DT, DS <: DataSeries{DT,DT}}
    print(io, "$(DS) with data type ", DT, "\n")
    print(io, parent(parent(ds)))
end

function Base.zero(ds::DataSeries)
    DataSeries(zero(ds[begin]), ntime(ds))
end
