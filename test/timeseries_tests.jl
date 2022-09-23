using GeometricBase
using GeometricSolutions
using OffsetArrays
using Test


@testset "$(rpad("Timeseries",80))" begin
    nt = 10
    tbeg = 0.0
    tend = 1.0
    Δt = 0.1
    ti = tbeg:Δt:tend

    ts  = TimeSeries(ti, Δt)
    ts1 = TimeSeries(tbeg, tend, Δt)
    ts2 = TimeSeries(nt, Δt)
    ts3 = TimeSeries(collect(ti))

    @test ts == ts1 == ts2 == ts3

    @test typeof(ts) <: AbstractVector
    @test typeof(parent(ts)) <: OffsetVector
    @test parent(parent(ts)) == ti
    @test eltype(ts) == typeof(Δt)
    @test ndims(ts)  == 1
    @test ntime(ts)  == nt

    @test firstindex(ts)   == 0
    @test firstindex(ts,1) == firstindex(ts.t,1)
    @test firstindex(ts,2) == 1
    @test lastindex(ts)    == nt
    @test lastindex(ts,1)  == lastindex(ts.t,1)
    @test lastindex(ts,2)  == 1
    @test strides(ts)      == (1,)
    @test stride(ts,1)     == 1
    @test axes(ts)   == (0:nt,)
    @test axes(ts,1) == 0:nt
    @test axes(ts,2) == 1:1
    @test size(ts)   == (nt+1,)
    @test size(ts.t) == (nt+1,)
    @test size(ts,1) == nt+1
    @test size(ts,2) == 1

    @test eachindex(ts) == 0:nt
    @test eachindex(IndexLinear(), ts) == 0:nt
    @test eachindex(IndexCartesian(), ts) == CartesianIndices((0:nt,))

    @test ts[0] == ts[begin] == tbeg
    @test ts[nt] == ts[end] == tend

    @test ts == vec(ti)
    @test vec(ti) == ts
end
