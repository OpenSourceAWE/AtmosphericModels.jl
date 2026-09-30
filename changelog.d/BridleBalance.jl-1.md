### Added
- `new_windfield` takes an `rng` keyword, so that a caller can generate independent
  turbulence realizations; the default, `StableRNG(1234)`, gives the field it always gave.
