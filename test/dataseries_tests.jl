using GeometricSolutions
using OffsetArrays
using Test

dt = Float64
nt = 10

@testset "$(rpad("Dataseries with scalar data",80))" begin
    ds = DataSeries(rand(dt), nt)
    @test typeof(ds) <: DataSeries{dt, dt}
    @test typeof(ds) <: ScalarDataSeries{dt}
    @test typeof(ds) <: AbstractVector{dt}
    @test ds == DataSeries(ds)
    @test ds == DataSeries(parent(ds))
    @test ds == DataSeries(parent(parent(ds)))
    @test firstindex(ds) == 0
    @test firstindex(ds, 1) == firstindex(ds.d, 1)
    @test firstindex(ds, 2) == 1
    @test lastindex(ds) == nt
    @test lastindex(ds, 1) == lastindex(ds.d, 1)
    @test lastindex(ds, 2) == 1
    @test strides(ds) == (1,)
    @test stride(ds, 1) == 1
    @test axes(ds) == (0:nt,)
    @test axes(ds, 1) == 0:nt
    @test axes(ds, 2) == 1:1
    @test size(ds.d) == (nt + 1,)
    @test size(ds.d) == size(ds)
    @test size(ds.d, 1) == nt + 1
    @test ndims(ds) == 1
    @test eltype(ds) == dt
    @test parent(ds) == ds.d
    @test typeof(parent(ds)) <: OffsetArray
    @test ntime(ds) == nt

    @test eachindex(ds) == 0:nt
    @test eachindex(IndexLinear(), ds) == 0:nt
    @test eachindex(IndexCartesian(), ds) == CartesianIndices((0:nt,))

    for i in 0:nt
        ds[i] = i
    end

    @test parent(ds)[0] == ds[0] == parent(ds)[begin] == ds[begin]
    @test parent(ds)[nt] == ds[nt] == parent(ds)[end] == ds[end]
    @test parent(ds.d) == collect(0:nt)
    @test ds[:] == parent(ds)[:]

    @test parent(ds) == OffsetVector([i for i in 0:nt], 0:nt)

    d = [x for x in ds]
    @test d == parent(ds)
    @test Array(ds) == parent(parent(ds))'

    d = rand(nt + 1)
    ds = DataSeries(d)
    @test parent(parent(ds)) == d

    @test parent(zero(ds)) == zero(parent(ds))

    @test relative_maximum_error(DataSeries(d), ds) == 0
end

@testset "$(rpad("Dataseries with vector-valued data",80))" begin
    nd = 2
    at = Vector{dt}

    ds = DataSeries(rand(dt, nd), nt)
    @test typeof(ds) <: DataSeries{dt, at}
    @test typeof(ds) <: AbstractVector{at}
    @test ds == DataSeries(ds)
    @test ds == DataSeries(parent(ds))
    @test ds == DataSeries(parent(parent(ds)))
    @test firstindex(ds) == 0
    @test firstindex(ds, 1) == firstindex(ds.d, 1)
    @test firstindex(ds, 2) == 1
    @test lastindex(ds) == nt
    @test lastindex(ds, 1) == lastindex(ds.d, 1)
    @test lastindex(ds, 2) == 1
    @test strides(ds) == (1,)
    @test stride(ds, 1) == 1
    @test axes(ds) == (0:nt,)
    @test axes(ds, 1) == 0:nt
    @test axes(ds, 2) == 1:1
    @test size(ds.d) == (nt + 1,)
    @test size(ds.d) == size(ds)
    @test size(ds.d, 1) == nt + 1
    @test ndims(ds) == 1
    @test eltype(ds) == dt
    @test arrtype(ds) == at
    @test parent(ds) == ds.d
    @test typeof(parent(ds)) <: OffsetArray
    @test ntime(ds) == nt

    @test eachindex(ds) == 0:nt
    @test eachindex(IndexLinear(), ds) == 0:nt
    @test eachindex(IndexCartesian(), ds) == CartesianIndices((0:nt,))

    ds2 = DataSeries(rand(dt, nd), nt)

    for i in 0:nt
        ds[i] = [i, i^2]
    end

    @test ds[1][1] == ds[1, 1]
    @test ds[1][2] == ds[1, 2]

    @test parent(ds)[0] == ds[0] == parent(ds)[begin] == ds[begin]
    @test parent(ds)[nt] == ds[nt] == parent(ds)[end] == ds[end]
    @test ds[:] == parent(ds)[:]

    @test parent(ds) == OffsetVector([Vector{dt}([i, i^2]) for i in 0:nt], 0:nt)

    d = [x for x in ds]
    @test d == parent(ds)

    d = rand(nt + 1)
    ds = DataSeries(d)
    @test parent(parent(ds)) == d

    @test parent(zero(ds)) == zero(parent(ds))

    @test relative_maximum_error(DataSeries(d), ds) == 0
end
