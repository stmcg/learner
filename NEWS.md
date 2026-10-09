# learner 1.1.0

* Integrate multiple-source support and separate space penalties from the
  preliminary `learnerv2` 0.3.0 implementation into the package's R and C++
  implementation.
* `learner()`, `cv.learner()`, and `dlearner()` accept a list of aligned source
  matrices and optional `source_weights`. Nonnegative weights are normalized
  to sum to one and combine source projection operators. The common rank
  interface is retained: automatic rank selection uses the first source,
  and LEARNER initialization uses the first source even when rank is supplied.
* `learner()` accepts `lambda_1_row` and `lambda_1_col` to penalize the U and V
  source-space residuals separately. The original `lambda_1` interface,
  including positional calls, is retained. Equal separate coefficients reproduce
  the common-penalty calculation under the same settings.
* `cv.learner()` accepts separate row and column penalty grids and returns a
  three-dimensional validation-score array. Common-penalty calls retain their
  previous output structure and score scale (summed held-out squared errors).
* Cache source decompositions during cross-validation and use C++ data in
  parallel workers. Validate sources, weights, penalties, ranks, and fold counts.
* Expand argument documentation and examples for multiple sources, separate
  penalties, cross-validation, missing data, and the current stopping rules.
* Expand automated tests for compatibility, gradients, multiple-source fitting,
  cross-validation, and stopping behavior. Add fixed reference outputs from
  `learnerv2` 0.3.0 to detect unintended numerical changes.
* Add Jodie Li as a package author.

### learner version 1.0.0 (2025-03-02)

* Re-implemented most of `learner` and `cv.learner` in Rcpp, resulting in significant speed-ups
* Re-ordered matrix computations in `dlearner`, resulting in significant speed-ups for large matrices
* Fixed a bug in `cv.learner` when `Y_target` has missing entries
* Expanded unit testing

### learner version 0.1.0 (2025-01-08)

* First version released on GitHub (https://github.com/stmcg/learner) and CRAN (https://CRAN.R-project.org/package=learner)
