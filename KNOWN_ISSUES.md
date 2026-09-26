# Known issues

Defects found in review and not fixed in the change that found them.

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
