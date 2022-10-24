module GeometricSolutions

using Base: TwicePrecision
using OffsetArrays: OffsetArray, OffsetVector
using Reexport
using Test

using GeometricBase
using GeometricEquations


include("dataseries.jl")

export DataSeries, ScalarDataSeries


include("timeseries.jl")

export TimeSeries


include("abstract_solution.jl")

export AbstractSolution
export nstore


include("geometric_solution.jl")

export GeometricSolution

end
