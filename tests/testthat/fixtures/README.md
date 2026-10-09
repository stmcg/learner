# Fixed reference results for the extensions

`multisource-v030.rds` contains inputs, fitting arguments, and expected outputs
computed by Sean McGrath's base-R **learnerv2 0.3.0** implementation from the
`learner-paper` repository's simulations materials. The tests do not create
expected values by calling the `learner` implementation under test.

## Coverage

- Three simulated 9-by-7 source matrices with weights 2, 3, 5 and fixed rank 2.
- Common penalties and unequal U/V penalties, for complete and missing targets.
- Complete fitted matrices and all 15 objective values, rank, and stopping code.
- Separate-penalty CV with 2-by-2-by-2 candidates, three folds, and fixed RNG
  settings, for complete and missing targets. All scores and selected values
  are checked; each saved CV case has a unique minimum.
- The complete multiple-source D-LEARNER estimate, also checked at generation
  time against explicit dense weighted source operators.

The saved inputs are synthetic. They are generated with seed 405 and then stored
in full, rather than regenerated during tests. CV uses seed 514 with RNG kinds
`Mersenne-Twister`, `Inversion`, and `Rejection`; its test restores the incoming
RNG state. The fixture records the reference version, source-file MD5 hashes,
R version, and RNG kinds so the source of the expected values can be audited.

`expect_equal()` uses tolerance `1e-8` to allow ordinary floating-point and SVD
library differences. Matrix estimates and objective histories are checked, not
raw U/V factors, whose SVD signs are arbitrary. Rank is fixed to avoid depending
on future changes in automatic rank selection.

The 15-update examples use `threshold = 0` and a large growth limit to compare a
fixed optimization path. They do not demonstrate convergence or statistical
accuracy. They do not assert that the two implementations have identical stopping
behavior for all settings. The analytic tests in `test-stopping-behavior.R`
separately check the current package's documented stopping rules.

Normal tests only read this fixture and require neither `learnerv2` nor `foreach`.

## Deliberate regeneration

Run the generator manually from the package root, using an extracted, unmodified
copy of learnerv2 0.3.0. The generator requires `foreach`; it verifies the package
name and version and sources the reference R files without loading `learner`.

```sh
Rscript tests/testthat/fixtures/generate-reference.R \
  /path/to/learnerv2 \
  tests/testthat/fixtures/multisource-v030.rds
```

The generator refuses to overwrite an existing baseline unless `--overwrite` is
added. Do not regenerate the fixture merely to make a failing test pass. First
understand whether a change is an intended method/interface change, an optimization
difference, or a bug; review any new baseline independently and update this note.
Neither testthat nor R CMD check runs the generator automatically.
