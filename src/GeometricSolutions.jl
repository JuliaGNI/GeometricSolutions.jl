module GeometricSolutions

using Base: TwicePrecision
using OffsetArrays: OffsetVector
using Reexport
using Test

@reexport using GeometricBase


include("dataseries.jl")

export get_data!, set_data!
export AbstractDataSeries, DataSeries,
       DataSeriesConstructor


include("timeseries.jl")

export TimeSeries


include("abstract_solution.jl")

export AbstractSolution
export counter, offset, lastentry

end
