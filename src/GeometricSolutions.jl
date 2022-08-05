module GeometricSolutions

using Base: TwicePrecision
using OffsetArrays: OffsetVector
using Reexport
using Test

@reexport using GeometricBase
using GeometricEquations


include("dataseries.jl")

export DataSeries, ScalarDataSeries
export arrtype


include("timeseries.jl")

export TimeSeries


include("abstract_solution.jl")

export AbstractSolution
export counter, offset, lastentry


include("geometric_solution.jl")

export GeometricSolution

end
