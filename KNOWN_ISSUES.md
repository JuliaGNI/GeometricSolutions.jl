# Known issues

### K1 · `Project.toml` declares `version = "0.6.4"` while the newest tag is `v0.6.5`.

- **location:** `Project.toml`
- **evidence:** Noticed
  2026-08-31 while seeding this file; the cause has not been established and nothing was
  changed. Resolve before the next release, since the close-out commit sets `version` by hand
  and would otherwise re-use a burnt number.
- **kind:** not verified
- **found:** 2026-08-31

### K2 · A single `GeometricSolution` has no HDF5 methods.

- **location:** —
- **evidence:** Only `EnsembleSolution` has them.
- **kind:** defect
- **found:** 2026-09-22
