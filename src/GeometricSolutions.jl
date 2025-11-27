module GeometricSolutions

using GeometricBase
using GeometricEquations
using Test

using Base: TwicePrecision
using OffsetArrays: OffsetArray, OffsetVector

import GeometricBase: AbstractSolution
import GeometricBase: initialtime, finaltime, timespan
import GeometricBase: eachtimestep, timestep, timesteps

export AbstractSolution
export arrtype, nstore, ntime
export initialtime, finaltime, timespan
export eachtimestep, timestep, timesteps
export maximum_error, relative_maximum_error, relative_norm_error

include("utils.jl")

include("dataseries.jl")

export DataSeries, ScalarDataSeries

include("timeseries.jl")

export TimeSeries

include("geometric_solution.jl")

export GeometricSolution

include("ensemble_solution.jl")

export EnsembleSolution, solution

include("diagnostics.jl")

export compute_invariant, compute_invariant_error, compute_momentum_error, compute_one_form
export compute_difference, compute_error_drift, compute_relative_error

const SolutionODE{dType, tType, dsType, probType, perType} = GeometricSolution{
    dType, tType, dsType, probType,
    perType} where {probType <: Union{ODEProblem, SODEProblem, SubstepProblem}}
const SolutionDAE{dType, tType, dsType, probType, perType} = GeometricSolution{
    dType, tType, dsType, probType, perType} where {probType <: DAEProblem}
const SolutionSDE{dType, tType, dsType, probType, perType} = GeometricSolution{
    dType, tType, dsType, probType, perType} where {probType <: SDEProblem}
const SolutionPODE{dType, tType, dsType, probType, perType} = GeometricSolution{
    dType, tType, dsType, probType,
    perType} where {probType <: Union{PODEProblem, HODEProblem, IODEProblem, LODEProblem}}
const SolutionPDAE{dType, tType, dsType, probType, perType} = GeometricSolution{
    dType, tType, dsType, probType,
    perType} where {probType <: Union{PDAEProblem, HDAEProblem, IDAEProblem, LDAEProblem}}
const SolutionPSDE{dType, tType, dsType, probType, perType} = GeometricSolution{
    dType, tType, dsType, probType, perType} where {probType <: PSDEProblem}

export SolutionODE, SolutionDAE, SolutionSDE
export SolutionPODE, SolutionPDAE, SolutionPSDE

end
