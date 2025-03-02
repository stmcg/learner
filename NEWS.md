### learner version 1.0.0 (2025-03-02)

* Re-implemented most of `learner` and `cv.learner` in Rcpp, resulting in significant speed-ups
* Re-ordered matrix computations in `dlearner`, resulting in significant speed-ups for large matrices
* Fixed a bug in `cv.learner` when `Y_target` has missing entries
* Expanded unit testing

### learner version 0.1.0 (2025-01-08)

* First version released on GitHub (https://github.com/stmcg/learner) and CRAN (https://CRAN.R-project.org/package=learner)
