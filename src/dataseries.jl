
using Base: Indices
using OffsetArrays: OffsetArray

abstract type AbstractDataSeries{T <: AbstractData, N, AT <: AbstractArray{T,N}} <: AbstractArray{T,N} end

function initialize!(ds::AbstractDataSeries{T,1}, q₀::T) where {T}
    for j in eachindex(ds)
        ds[j] = zero(q₀)
    end
    return ds
end

function initialize!(ds::AbstractDataSeries{T,1}, q₀::Vector{T}) where {T}
    @assert nsamples(ds) == length(q₀) == 1
    initialize!(ds, q₀[begin])
end

function initialize!(ds::AbstractDataSeries{T,2}, q₀::Vector{T}) where {T}
    @assert nsamples(ds) == length(q₀)
    for k in axes(ds,2)
        for j in axes(ds,1)
            ds[j,k] = zero(q₀[k])
        end
    end
    return ds
end


function DataSeriesConstructor(::Type{dsType}, ::Type{T}, nt::Int, ni::Int=1) where {dsType <: AbstractDataSeries, T <: Number}
    ds = dsType{T, ni == 1 ? 1 : 2}(nt, ni)
    initialize!(ds, zero(T))
end

function DataSeriesConstructor(::Type{dsType}, ::Type{T}, nd::Int, nt::Int, ni::Int) where {dsType <: AbstractDataSeries, T <: Number}
    ds = dsType{Vector{T}, ni == 1 ? 1 : 2}(nt, ni)
    initialize!(ds, [zeros(T, nd) for _ in 1:ni])
end

function DataSeriesConstructor(::Type{dsType}, q₀::T, nt::Int, ni::Int=1) where {dsType <: AbstractDataSeries, T <: AbstractData}
    @assert ni == 1
    ds = dsType{T,1}(nt, ni)
    initialize!(ds, q₀)
end

function DataSeriesConstructor(::Type{dsType}, q₀::Vector{T}, nt::Int, ni::Int=1) where {dsType <: AbstractDataSeries, T <: Number}
    @assert ni == 1 || ni == length(q₀)
    ds = ( ni == 1 ? dsType{Vector{T},1}(nt, ni) : dsType{T,2}(nt, ni) )
    initialize!(ds, q₀)
end

function DataSeriesConstructor(::Type{dsType}, q₀::Vector{T}, nt::Int, ni::Int=length(q₀)) where {dsType <: AbstractDataSeries, T <: AbstractArray{<:Number}}
    @assert ni == length(q₀)
    ds = dsType{T, ni == 1 ? 1 : 2}(nt, ni)
    initialize!(ds, q₀)
end

function DataSeriesConstructor(::Type{dsType}, d::AbstractArray{T,1}) where {dsType <: AbstractDataSeries, T}
    nt = size(d,1)-1
    ni = 1
    ds = dsType{T,1}(nt, ni)
    copy!(parent(ds), d)
    return ds
end

function DataSeriesConstructor(::Type{dsType}, d::AbstractArray{T,2}) where {dsType <: AbstractDataSeries, T}
    nt = size(d,1)-1
    ni = size(d,2)
    ds = dsType{T,2}(nt, ni)
    copy!(parent(ds), d)
    return ds
end


Base.parent(ds::DS) where {DS <: AbstractDataSeries} = error("parent() not implemented for $(DS)")

Base.eltype(ds::AbstractDataSeries{T,N}) where {T,N} = T
Base.ndims(ds::AbstractDataSeries{T,N}) where {T,N} = N

Base.:(==)(ds1::AbstractDataSeries, ds2::AbstractDataSeries) = (parent(ds1) == parent(ds2))

function Base.show(io::IO, ds::DS) where {T, N, DS <: AbstractDataSeries{T,N}}
    print(io, "$(DS) with data type ", T, " and ", N, " dimensions:\n")
    print(io, parent(ds))
end

GeometricBase.ntime(ds::AbstractDataSeries) = lastindex(ds.d, 1) - 1
GeometricBase.nsamples(ds::AbstractDataSeries) = size(ds.d, 2)

Base.size(ds::AbstractDataSeries) = size(ds.d)
Base.size(ds::AbstractDataSeries, d) = size(ds.d, d)


function Base.similar(ds::DS) where {T, DS <: AbstractDataSeries{T,1}}
    newds = DS(ntime(ds), nsamples(ds))
    for j in eachindex(ds)
        newds[j] = zero(ds[j])
    end
    return newds
end

function Base.similar(ds::DS) where {T, DS <: AbstractDataSeries{T,2}}
    newds = DS(ntime(ds), nsamples(ds))
    for k in axes(ds,2)
        for j in axes(ds,1)
            newds[j,k] = zero(ds[j,k])
        end
    end
    return newds
end


Base.eachindex(::IndexCartesian, ds::AbstractDataSeries) = CartesianIndices(axes(ds))
# Base.eachindex(::IndexLinear, ds::AbstractDataSeries{T,1}) where {T} = LinearIndices(axes(ds,1))

Base.firstindex(ds::AbstractDataSeries{T,1}) where {T} = 0
Base.firstindex(ds::AbstractDataSeries{T,N}) where {T,N} = firstindex(ds.d)

Base.firstindex(ds::AbstractDataSeries{T,1}, d) where {T} = d ≥ 1 && d ≤ 1 ? 0 : 1
Base.firstindex(ds::AbstractDataSeries{T,2}, d) where {T} = d ≥ 1 && d ≤ 2 ? (0, 1)[d] : 1

Base.lastindex(ds::AbstractDataSeries{T,1}) where {T} = ntime(ds)
Base.lastindex(ds::AbstractDataSeries{T,N}) where {T,N} = lastindex(ds.d)

Base.lastindex(ds::AbstractDataSeries{T,1}, d) where {T} = d ≥ 1 && d ≤ 1 ? ntime(ds) : 1
Base.lastindex(ds::AbstractDataSeries{T,2}, d) where {T} = d ≥ 1 && d ≤ 2 ? (ntime(ds), nsamples(ds))[d] : 1

@inline Base.axes(ds::AbstractDataSeries{T,1}) where {T} = (0:ntime(ds),)
@inline Base.axes(ds::AbstractDataSeries{T,2}) where {T} = (0:ntime(ds), 1:nsamples(ds))
@inline Base.axes(ds::AbstractDataSeries{T,N}, d) where {T,N} = d ≥ 1 && d ≤ N ? axes(ds)[d] : (1:1)

Base.strides(ds::AbstractDataSeries) = strides(ds.d)


function get_data!(ds::AbstractDataSeries{T,1}, x::T, n::Int, k::Int=1) where {T}
    @assert n ≥ 0 && n ≤ ntime(ds)
    @assert k == 1
    x .= ds.d[n+1]
end

function get_data!(ds::AbstractDataSeries{T,2}, x::T, n::Int, k::Int) where {T}
    @assert n ≥ 0 && n ≤ ntime(ds)
    @assert k ≥ 1 && k ≤ nsamples(ds)
    x .= ds.d[n+1,k]
end

function get_data!(ds::AbstractDataSeries{T,2}, x::Vector{T}, n::Int) where {T}
    @assert length(x) == nsamples(ds)
    @assert n ≥ 0 && n ≤ ntime(ds)
    @inbounds for k in eachindex(x)
        x[k] = ds.d[n+1,k]
    end
end

function set_data!(ds::AbstractDataSeries{T,1}, x::T, n::Int, k::Int=1) where {T}
    @assert n ≥ 0 && n ≤ ntime(ds)
    @assert k == 1
    if T <: Number
        ds.d[n+1] = x
    elseif T <: AbstractArray
        ds.d[n+1] .= x
    end
end

function set_data!(ds::AbstractDataSeries{T,1}, x::Vector{T}, n::Int, k::Int=1) where {T}
    @assert nsamples(ds) == length(x) == 1
    set_data!(ds, x[begin], n, k)
end

function set_data!(ds::AbstractDataSeries{T,2}, x::T, n::Int, k::Int) where {T}
    @assert n ≥ 0 && n ≤ ntime(ds)
    if T <: Number
        ds.d[n+1,k] = x
    elseif T <: AbstractArray
        ds.d[n+1,k] .= x
    end
end

function set_data!(ds::AbstractDataSeries{T,2}, x::Vector{T}, n::Int) where {T}
    @assert length(x) == nsamples(ds)
    @assert n ≥ 0 && n ≤ ntime(ds)
    @inbounds for k in axes(ds.d, 2)
        if T isa Number
            ds.d[n+1,k] = x[k]
        elseif T isa AbstractArray
            ds.d[n+1,k] .= x[k]
        end
    end
end

function GeometricBase.reset!(ds::AbstractDataSeries{T,1}) where {T}
    @inbounds ds[begin] = ds[end]
end

function GeometricBase.reset!(ds::AbstractDataSeries{T,2}) where {T}
    @inbounds for k in axes(ds,2)
        ds[begin,k] = ds[end,k]
    end
end


@inline function Base.getindex(ds::AbstractDataSeries{T,1}, I::Colon) where {T}
    OffsetArray([ds[i] for i in axes(ds,1)], axes(ds,1))
end

@inline function Base.getindex(ds::AbstractDataSeries{T,1}, j::Union{UnitRange,Int}) where {T}
    getindex(ds.d, j.+1)
end

@inline function Base.getindex(ds::AbstractDataSeries{T,1}, i::Union{Int,CartesianIndex}, j::Union{Int,CartesianIndex}) where {T}
    ds[j][i]
end

@inline function Base.getindex(ds::AbstractDataSeries{T,1}, i::Union{Int,CartesianIndex}, J::AbstractRange{Int}) where {T}
    [ds[j][i] for j in J]
end

@inline function Base.getindex(ds::AbstractDataSeries{T,1}, i::Union{Int,CartesianIndex}, J::Colon) where {T}
    OffsetArray([ds[j][i] for j in axes(ds,1)], axes(ds,1))
end

@inline Base.getindex(ds::AbstractDataSeries{T,1}, I, J::Colon) where {T} = getindex(ds, I, axes(ds,1))


@inline function Base.getindex(ds::AbstractDataSeries{T,2}, j::Union{UnitRange,Int}, k::Union{UnitRange,Int}) where {T}
    getindex(ds.d, j.+1, k)
end

@inline function Base.getindex(ds::AbstractDataSeries{T,2}, i::Union{Int,CartesianIndex}, j::Union{Int,CartesianIndex}, k::Union{Int,CartesianIndex}) where {T}
    ds[j,k][i]
end

@inline function Base.getindex(ds::AbstractDataSeries{T,2}, i::Union{Int,CartesianIndex}, J::AbstractRange{Int}, K::AbstractRange{Int}) where {T}
    [ds[j,k][i] for j in J, k in K]
end

@inline function Base.getindex(ds::AbstractDataSeries{T,2}, i::Union{Int,CartesianIndex}, J::Colon, K::AbstractRange{Int}) where {T}
    OffsetArray([ds[j,k][i] for j in axes(ds,1), k in K], axes(ds,1), K)
end

@inline function Base.getindex(ds::AbstractDataSeries{T,2}, i::Union{Int,CartesianIndex}, J::AbstractRange{Int}, k::Union{Int,CartesianIndex}) where {T}
    [ds[j,k][i] for j in J]
end

@inline function Base.getindex(ds::AbstractDataSeries{T,2}, i::Union{Int,CartesianIndex}, J::Colon, k::Union{Int,CartesianIndex}) where {T}
    OffsetArray([ds[j,k][i] for j in axes(ds,1)], axes(ds,1))
end

@inline function Base.getindex(ds::AbstractDataSeries{T,2}, i::Union{Int,CartesianIndex}, j::Union{Int,CartesianIndex}, K::AbstractRange{Int}) where {T}
    [ds[j,k][i] for k in K]
end

@inline Base.getindex(ds::AbstractDataSeries{T,2}, J::Colon, K::Colon) where {T} = getindex(ds, axes(ds,1), axes(ds,2))
@inline Base.getindex(ds::AbstractDataSeries{T,2}, J, K::Colon) where {T} = getindex(ds, J, axes(ds,2))
@inline Base.getindex(ds::AbstractDataSeries{T,2}, J::Colon, K) where {T} = getindex(ds, axes(ds,1), K)

@inline Base.getindex(ds::AbstractDataSeries{T,2}, I, J::Colon, K::Colon) where {T} = getindex(ds, I, axes(ds,1), axes(ds,2))
@inline Base.getindex(ds::AbstractDataSeries{T,2}, I, J::Colon, K) where {T} = getindex(ds, I, axes(ds,1), K)
@inline Base.getindex(ds::AbstractDataSeries{T,2}, I, J, K::Colon) where {T} = getindex(ds, I, J, axes(ds,2))


@inline function Base.setindex!(ds::AbstractDataSeries{T,1}, x, j::Union{UnitRange,Int}) where {T}
    setindex!(ds.d, x, j.+1)
end

@inline function Base.setindex!(ds::AbstractDataSeries{T,2}, x, j::Union{UnitRange,Int}, k::Union{UnitRange,Int}) where {T}
    setindex!(ds.d, x, j.+1, k)
end


# 
# @inline function Base.getindex(ds::AbstractDataSeries{T,2}, I::Union{Int,CartesianIndex}, J::Union{Int,CartesianIndex,AbstractRange{Int}}, K::Union{Int,CartesianIndex,AbstractRange{Int}}) where {T}
#     x = zeros(eltype(ds.d[begin]), (1 for i in 1:length(I))..., length(J))
#     for j in J
#         x[Tuple(I)..., j] = ds[j][I]
#     end
# end
# 
# similar(storagetype, shape)
# 


function fromarray(::Type{dsType}, d::AbstractArray{T,1}, ni::Int=1) where {dsType <: AbstractDataSeries, T <: Number}
    @assert ni == 1
    return dsType(d)
end

function fromarray(::Type{dsType}, d::AbstractArray{T,2}, ni::Int=1) where {dsType <: AbstractDataSeries, T <: Number}
    @assert ni == 1 || ni == size(d,2)
    
    if ni == 1
        ds = [Vector{T}(d[:,n]) for n in axes(d,2)]
    else
        ds = Array(d)
    end

    return dsType(ds)
end

function fromarray(::Type{dsType}, d::AbstractArray{T,N}, ni::Int=1) where {dsType <: AbstractDataSeries, T <: Number, N}
    @assert ni == 1 || ni == size(d)[end]
    @assert N ≥ 3

    if ni == 1
        AT = Array{T,ndims(d)-1}
        ds = [AT(d[axes(d)[1:end-1]..., n]) for n in axes(d)[end]]
    else
        AT = Array{T,ndims(d)-2}
        ds = [AT(d[axes(d)[1:end-2]..., n, k]) for n in axes(d)[end-1], k in axes(d)[end]]
    end

    return dsType(ds)
end



struct DataSeries{T <: AbstractData, N} <: AbstractDataSeries{T, N, Array{T,N}}
    d::Array{T,N}

    function DataSeries{T,N}(nt, ni) where {T <: AbstractData, N}
        @assert nt ≥ 0
        @assert ni > 0

        @assert N ∈ (1,2)

        if N == 1
            d = Array{T,N}(undef, nt+1)
        elseif N == 2
            d = Array{T,N}(undef, nt+1, ni)
        end

        new(d)
    end
end

Base.parent(ds::DataSeries) = ds.d

DataSeries(::Type{T}, nt::Int, ni::Int=1) where {T} = DataSeriesConstructor(DataSeries, T, nt, ni)
DataSeries(::Type{T}, nd::Int, nt::Int, ni::Int) where {T} = DataSeriesConstructor(DataSeries, T, nd, nt, ni)
DataSeries(q₀::T, nt::Int, ni::Int=1) where {T <: AbstractData} = DataSeriesConstructor(DataSeries, q₀, nt, ni)
DataSeries(q₀::Vector{T}, nt::Int, ni::Int=1) where {T <: Number} = DataSeriesConstructor(DataSeries, q₀, nt, ni)
DataSeries(q₀::Vector{T}, nt::Int, ni::Int=length(q₀)) where {T <: AbstractArray{<:Number}} = DataSeriesConstructor(DataSeries, q₀, nt, ni)
DataSeries(d::AbstractArray) = DataSeriesConstructor(DataSeries, d)
