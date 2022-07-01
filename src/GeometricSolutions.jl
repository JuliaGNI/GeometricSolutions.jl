module GeometricSolutions

using Reexport
using Test

@reexport using GeometricBase


include("dataseries.jl")

export get_data!, set_data!
export AbstractDataSeries, DataSeries,
       DataSeriesConstructor


include("timeseries.jl")

export TimeSeries, compute_timeseries!


include("abstract_solution.jl")

export AbstractSolution
export counter, offset, lastentry

end
