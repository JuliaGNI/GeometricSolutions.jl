# Known issues

Defects found in review and not fixed in the change that found them.

### K2 · `Project.toml` declares `version = "0.6.4"` while the newest tag is `v0.6.5`.

- **location:** `Project.toml`
- **evidence:** Noticed
  2026-08-31 while seeding this file; the cause has not been established and nothing was
  changed. Resolve before the next release, since the close-out commit sets `version` by hand
  and would otherwise re-use a burnt number.
- **kind:** not verified
- **found:** 2026-08-31

### K3 · A single `GeometricSolution` has no HDF5 methods.

- **location:** —
- **evidence:** Only `EnsembleSolution` has them.
- **kind:** defect
- **found:** 2026-09-22

### K1 · No test in `test/dataseries.jl` catches two mutants of `src/dataseries.jl`.

- **location:** `src/dataseries.jl:55`, `src/dataseries.jl:90`
- **evidence:** two mutants survive `test/dataseries.jl`.
  - `GeometricBase.reset!(ds::DataSeries) = ds[begin] = ds[end]` changed to
    `ds[begin] = ds[begin]` survives every test file: no test calls `reset!`.
  - `err / maximum(abs, ref)` changed to `err / minimum(abs, ref)` in the pointwise relative
    error survives `test/dataseries.jl`. `test/diagnostics.jl:176` and `:183` catch it.

  Reproducer, with the mutation applied by hand to `src/dataseries.jl`:
  `julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/dataseries.jl")'`
  passes. Fix: add assertions for `reset!` and for the pointwise relative maximum error to
  `test/dataseries.jl`.
- **kind:** missing test
- **found:** 2026-09-26
