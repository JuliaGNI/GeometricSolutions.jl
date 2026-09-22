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

## Open Issues

- `Project.toml` declares `version = "0.6.4"` while the newest tag is `v0.6.5`. Noticed
  2026-08-31 while seeding this file; the cause has not been established and nothing was
  changed. Resolve before the next release, since the close-out commit sets `version` by hand
  and would otherwise re-use a burnt number.
- `GeometricBase = "0.14.12"` names a version that is not registered yet. Until GeometricBase
  releases it, `Pkg.test()` and CI do not resolve; the suite runs only in an environment that
  develops a GeometricBase tree with `h5save` and `h5load`.
- A single `GeometricSolution` has no HDF5 methods. Only `EnsembleSolution` has them.

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
  step, time span, stored-step count, or any member's parameters.

  The two functions are GeometricBase's generics, and this package does not export them.
  `ReducedComplexityModeling` exports its own `h5save` and `h5load`, so an export here would make
  both names an `UndefVarError` for a caller that loads the two packages. The package now needs
  GeometricBase 0.14.12, the first version with the generics.

### Bug Fixes

### Breaking Changes
