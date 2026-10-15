module GeometricSolutions

using GeometricBase
using GeometricEquations
using LinearAlgebra: norm
using Test

using OffsetArrays: OffsetArray, OffsetVector

import GeometricBase: AbstractSolution
import GeometricBase: initialstate, initialtime, finaltime, timespan
import GeometricBase: eachtimestep, timestep, timesteps, variables

export AbstractSolution
export arrtype, nstore, ntime
export initialtime, finaltime, timespan
export eachtimestep, timestep, timesteps
export maximum_error, relative_maximum_error, relative_norm_error
export states, solutions

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

const SolutionODE{dType, tType, tsType, dsType, stType,
    probType <: Union{ODEProblem, SODEProblem, SubstepProblem}, perType} = GeometricSolution{
    dType, tType, tsType, dsType, stType, probType, perType}
const SolutionDAE{dType, tType, tsType, dsType, stType,
    probType <: DAEProblem, perType} = GeometricSolution{
    dType, tType, tsType, dsType, stType, probType, perType}
const SolutionSDE{dType, tType, tsType, dsType, stType,
    probType <: SDEProblem, perType} = GeometricSolution{
    dType, tType, tsType, dsType, stType, probType, perType}
const SolutionPODE{dType, tType, tsType, dsType, stType,
    probType <: Union{PODEProblem, HODEProblem, IODEProblem, LODEProblem}, perType} = GeometricSolution{
    dType, tType, tsType, dsType, stType, probType, perType}
const SolutionPDAE{dType, tType, tsType, dsType, stType,
    probType <: Union{PDAEProblem, HDAEProblem, IDAEProblem, LDAEProblem}, perType} = GeometricSolution{
    dType, tType, tsType, dsType, stType, probType, perType}
const SolutionPSDE{dType, tType, tsType, dsType, stType,
    probType <: PSDEProblem, perType} = GeometricSolution{
    dType, tType, tsType, dsType, stType, probType, perType}

export SolutionODE, SolutionDAE, SolutionSDE
export SolutionPODE, SolutionPDAE, SolutionPSDE

end
