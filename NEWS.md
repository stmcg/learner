### learner version 1.1.0 (unreleased)

* Add support for multiple source populations and optional source weights in
  `learner()`, `dlearner()`, and `cv.learner()`.
* Add `dat_multisource`, an example dataset with one target and two source
  populations with independently generated noise.
* Add separate row- and column-space penalties in `learner()` and
  cross-validation over these penalties in `cv.learner()`.
* Add plotting functions for objective function values, cross-validation
  results, population latent spaces, and variable contributions.

### learner version 1.0.0 (2025-03-02)

* Re-implemented most of `learner` and `cv.learner` in Rcpp, resulting in significant speed-ups
* Re-ordered matrix computations in `dlearner`, resulting in significant speed-ups for large matrices
* Fixed a bug in `cv.learner` when `Y_target` has missing entries
* Expanded unit testing

### learner version 0.1.0 (2025-01-08)

* First version released on GitHub (https://github.com/stmcg/learner) and CRAN (https://CRAN.R-project.org/package=learner)
