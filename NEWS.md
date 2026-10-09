# learner (development)

* Add `plot_objective()` for full and recent objective histories, with stopping
  codes and an explicitly labeled recorded minimum.
* Add `plot_cv()` for common and separate penalty grids. The ggplot2 heatmaps
  show relative error above the minimum using reversed magma colors, percentage
  labels, and a black selection box. Separate panels share one legend.
* Report absolute SSE or MSE ranges, exact ties, and boundary selections.
  Finite errors above an optional relative cutoff are labeled as high error,
  never optimizer divergence. Zero minima use absolute scores instead.
  Return ggplot objects for customization and `ggsave()` export.
* Add `plot_cv_profile()` for penalty sensitivity curves, including separate
  row and column penalties and optional base-10 horizontal coordinates.
* Add `plot_source_spaces()` for population projector submatrices and
  `plot_contributions()` for squared singular-vector contribution scores.
* Add `plot_simulation_results()` for replicate-level errors, with shared
  scales and optional explicitly labeled SD or SE bars.
* Include reproducible plotting examples and automated tests. These plots do
  not establish convergence, assign biological labels, or reproduce the
  published experiments.
* Add ggplot2 as a plotting dependency and require R >= 3.6.0 for the color
  palette functions used by the simulation plots.

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
