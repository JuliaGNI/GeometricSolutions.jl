# Release Notes

All notable changes to GeometricSolutions.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31 and deliberately holds no entries. 52 versions were
released before it, the most recent tag `v0.6.5`, and none of them are written up here: the
record of that history is `git log` and the tags. It is named as a gap rather than
reconstructed, because a changelog assembled after the fact loses exactly the reasoning that
makes it worth keeping. The `[Unreleased]` target below is provisional — confirm it when the
first entry is written.

## [Unreleased] — targeting 0.7.0

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
  both names an `UndefVarError` for a caller that loads the two packages. The package now needs
  GeometricBase 0.14.12, the first version with the generics.

### Bug Fixes

### Breaking Changes

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
