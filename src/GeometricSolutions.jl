module GeometricSolutions

using GeometricBase
using GeometricEquations
using Test

using Base: TwicePrecision
using OffsetArrays: OffsetArray, OffsetVector

import GeometricBase: eachtimestep, timestep, timesteps

export arrtype, nstore, ntime, tbegin, tend
export eachtimestep, timestep, timesteps
export maximum_error, relative_maximum_error, relative_norm_error


include("utils.jl")


include("dataseries.jl")

export DataSeries, ScalarDataSeries


include("timeseries.jl")

export TimeSeries


include("abstract_solution.jl")

export AbstractSolution


include("geometric_solution.jl")

export GeometricSolution


include("ensemble_solution.jl")

export EnsembleSolution, solution

end
