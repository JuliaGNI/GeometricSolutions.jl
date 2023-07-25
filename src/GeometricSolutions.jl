module GeometricSolutions

using Base: TwicePrecision
using OffsetArrays: OffsetArray, OffsetVector
using Reexport
using Test

using GeometricBase
using GeometricEquations

import GeometricBase: eachtimestep, timestep, timesteps

export arrtype, nstore, ntime, tbegin, tend
export eachtimestep, timestep, timesteps
export relative_maximum_error


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
