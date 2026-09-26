# Known issues

Defects found in review and not fixed in the change that found them.

## KI-1 · missing test · `src/dataseries.jl`

No test in `test/dataseries.jl` sees two mutants of `src/dataseries.jl`. The gap is also on the
earlier `test/dataseries_tests.jl`.

- `GeometricBase.reset!(ds::DataSeries) = ds[begin] = ds[end]` changed to `ds[begin] = ds[begin]`
  survives.
- `err / maximum(abs, ref)` changed to `err / minimum(abs, ref)` in the pointwise relative error
  survives.

Reproducer, with the mutation applied by hand to `src/dataseries.jl`:
`julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/dataseries.jl")'` passes.

Fix: add assertions for `reset!` and for the pointwise relative maximum error.
