# Release Notes

All notable changes to GeometricSolutions.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31. 52 versions were released before it, the most recent tag
`v0.6.5`, and none of them are written up here: the record of that history is `git log` and the
tags. It is named as a gap rather than reconstructed, because a changelog assembled after the
fact loses exactly the reasoning that makes it worth keeping.

## [0.7.0]

### Breaking Changes

- **The time series of an `EnsembleSolution` has a concrete type.** The field `t` was declared
  `TimeSeries{tType}`, which leaves the other parameters of `TimeSeries` open, so every access
  to `sol.t` was inferred as an abstract type. `EnsembleSolution` takes a fifth type parameter,
  `tsType <: TimeSeries{tType}`, appended after `probType`, and `t::tsType`. A signature or an
  `isa` test that names the first four parameters or fewer, such as
  `EnsembleSolution{T, T1, Vector{ST}}`, matches as before. Only code that takes
  `EnsembleSolution{dType, tType, sType, probType}` for a concrete type, for example as a field
  type or in a `typeof` comparison, sees a change: that spelling is now a `UnionAll`.

### Changes

- **The `OffsetArrays` floor is 1.14** (`[compat] OffsetArrays = "1.14"`, from `"1"`).
  `OffsetArrays` 1.13 and earlier redefine `Base.fill(::Any, ::Tuple{Vararg{Union{Integer,
  AbstractUnitRange}}})`, and Julia 1.12 does not permit method overwriting during
  precompilation, so `OffsetArrays` and with it this package do not precompile at those versions.
- **The test environment declares `[compat] Aqua = "0.8.14"`.** Without an entry, a resolve at
  the lower bounds takes Aqua 0.1.0, whose `test_all` has no `ambiguities` keyword. Aqua 0.8.14 is
  the first release whose piracy check runs on Julia 1.12: earlier releases read the `mt` field of
  `Core.TypeName`, which Julia 1.12 removes.

## [0.6.8]

### Changes

- **Julia 1.12 is the floor** (`[compat] julia = "1.12"`, from 1.11). A patch release, because no
  API changes: a caller on Julia 1.11 can no longer install this version.
  No source file changes, and no workaround for the older floor existed to delete: the package
  has no version-conditional code.

## [0.6.7]

### Bug Fixes

- **`SolutionODE`, `SolutionDAE`, `SolutionSDE`, `SolutionPODE`, `SolutionPDAE` and
  `SolutionPSDE` match the solutions they name.** Each alias passed five type parameters to
  `GeometricSolution`, which has seven, so the problem-type bound fell on the `NamedTuple` slot of
  the data series and no solution matched any alias: `GeometricSolution(prob) isa SolutionPODE`
  was `false` for a `PODEProblem` or an `LODEProblem`. Each alias now has the seven parameters
  of `GeometricSolution`, with its bounds and the bound on the problem type, so each alias is a
  subtype of `GeometricSolution` and a method on an alias is more specific than one on
  `GeometricSolution`. A caller that branches on these
  aliases, such as `plot_energy_error` in GeometricProblems, now takes the partitioned branch
  for a partitioned solution
  ([JuliaGNI/NonlinearIntegrators.jl#117](https://github.com/JuliaGNI/NonlinearIntegrators.jl/issues/117)).

### Changes

- CI uploads coverage from the `Julia 1 - ubuntu-latest` job instead of `Julia min`, and a test
  job saves the Julia cache only when it succeeds.
- `test/solutions.jl` moved to `test/integration/solutions.jl`, and `test/hdf5.jl` to
  `test/integration/hdf5.jl`. The test convention keeps a test file at the top level of `test/`
  only where it mirrors `src/<name>.jl`, and there is no `src/solutions.jl` or `src/hdf5.jl`.
  `solutions.jl` tests `src/geometric_solution.jl` and `src/ensemble_solution.jl`, whose deepest
  common directory is `src/` itself. `hdf5.jl` tests the package extension
  `ext/GeometricSolutionsHDF5Ext.jl`, which has no mirror under `src/`.

## [0.6.6]

### New Features

- An `HDF5` package extension, `GeometricSolutionsHDF5Ext`, stores an `EnsembleSolution`.
  `GeometricBase.h5save(h5, sol; path)` writes it, and
  `GeometricBase.h5load(EnsembleSolution, h5, problem; path)` reads it back. A state variable `q`
  is one dataset of size `(size(q)..., nstore + 1, nsamples)`, where column `n + 1` holds time
  index `n`. Each parameter is a dataset under `parameters/`, of size `(size(p)..., nsamples)`.
  HDF5 cannot hold the equation, so the read takes the `EnsembleProblem` and rebuilds the
  solution from it, which gives every `DataSeries` its 0-based axis. The read throws an
  `ArgumentError` when the file does not belong to that problem: a different member count, time
  step, time span, stored-step count, set of state variables, initial condition, or any member's
  parameters.

  The two functions are GeometricBase's generics, and this package does not export them.
  `ReducedComplexityModeling` exports its own `h5save` and `h5load`, so an export here would make
  both names an `UndefVarError` for a caller that loads the two packages. GeometricBase 0.14.12
  is the first version with the generics.

### Changes

- The `[compat]` floors are raised to `GeometricBase = "0.15.0"`, `GeometricEquations = "0.21.5"`
  and `julia = "1.11"`, because GeometricBase 0.15 declares its stubs public and requires Julia
  1.11. Julia 1.10 users keep 0.6.5. The floor raise changes no source file.

### Tests

- A new `test/quality/explicit_imports.jl` runs ExplicitImports in the `core` group.
- The test suite follows the shared layout. The test dependencies are in `test/Project.toml`,
  and `Project.toml` has no `[extras]` or `[targets]`. `runtests.jl` runs the files through
  `@safetestset` in the `core` group. Each test file is named after the source file it tests:
  `dataseries.jl`, `timeseries.jl`, `diagnostics.jl`, `solutions.jl` for `GeometricSolution` and
  `EnsembleSolution`, and `hdf5.jl` for the HDF5 extension. The files that draw random numbers
  set a fixed seed.
- A new `test/quality/aqua.jl` runs Aqua. Its ambiguity check is marked broken with issue #33:
  two ambiguities between `==` on `TimeSeries` and GeometricBase's `==` on `AbstractVariable`.
  `Project.toml` gets the bound `Test = "1"`, which Aqua's compat check needs.
- The link to issue #33 is on the same line as the `broken = true` mark in
  `test/quality/aqua.jl`, as the shared `test-layout.jl --check` requires for every broken mark.
- `test/Project.toml` no longer has `[compat]` entries for `GeometricBase`, `GeometricEquations`,
  `HDF5`, `LinearAlgebra`, `OffsetArrays` and `Test`. These are dependencies of the root
  `Project.toml`, so the root's bounds apply to them in the test environment too, and a second
  bound can only duplicate or narrow the root's. The shared `test-layout.jl --check` reports such
  an entry.
- An unused binding goes out of `test/dataseries.jl`.
