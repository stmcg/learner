# Developer-only script. It is never run by testthat or R CMD check.
# Usage: Rscript generate-reference.R /path/to/learnerv2-0.3.0 /path/to/output.rds
# Add --overwrite only after reviewing an intentional change to the baseline.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop('Supply the extracted reference package directory and output .rds path.')
reference_dir <- normalizePath(args[1], mustWork = TRUE)
output <- args[2]
if (file.exists(output) && !('--overwrite' %in% args)) {
  stop('Reference already exists; review the change before using --overwrite.')
}
description <- read.dcf(file.path(reference_dir, 'DESCRIPTION'))
stopifnot(description[1, 'Package'] == 'learnerv2', description[1, 'Version'] == '0.3.0')
# The reference CV implementation uses foreach even for sequential execution.
if (!requireNamespace('foreach', quietly = TRUE)) stop('Install foreach to generate this fixture.')
reference <- new.env(parent = asNamespace('foreach'))
reference_files <- file.path('R', c('learner.R', 'cv_learner.R', 'dlearner.R'))
for (f in reference_files) sys.source(file.path(reference_dir, f), envir = reference)

# Inputs are saved in the fixture, so tests do not regenerate them from the RNG.
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
set.seed(405)
sources <- list(matrix(rnorm(63), 9, 7), matrix(rnorm(63), 9, 7),
                matrix(rnorm(63), 9, 7))
target <- matrix(rnorm(63), 9, 7)
masked <- target
masked[c(3, 8, 19)] <- NA_real_
targets <- list(complete = target, missing = masked)
control <- list(max_iter = 15L, threshold = 0, max_value = 1e6)
rank <- 2L
weights <- c(2, 3, 5)
step <- 0.01
penalties <- list(common = list(lambda_1 = 3),
                  separate = list(lambda_1_row = 3, lambda_1_col = 0.25))
grids <- list(lambda_1_row_all = c(0.25, 3), lambda_1_col_all = c(0.5, 2),
              lambda_2_all = c(0, 1))
fit_cases <- list()
cv_cases <- list()
for (target_name in names(targets)) {
  for (penalty_name in names(penalties)) {
    fit_args <- c(list(Y_source = sources, Y_target = targets[[target_name]],
                       r = rank, source_weights = weights, lambda_2 = 1,
                       step_size = step, control = control), penalties[[penalty_name]])
    reference_args <- fit_args
    reference_args$r <- NULL
    reference_args$r_source <- rank
    reference_args$r_target <- rank
    expected <- do.call(reference$learner, reference_args)
    stopifnot(expected$convergence_criterion == 2L,
              length(expected$objective_values) == control$max_iter,
              all(is.finite(expected$objective_values)))
    fit_cases[[paste(target_name, penalty_name, sep = '_')]] <-
      list(args = fit_args, expected = expected)
  }
  cv_args <- c(list(Y_source = sources, Y_target = targets[[target_name]],
                    r = rank, source_weights = weights, step_size = step,
                    control = control, n_folds = 3L, n_cores = 1L), grids)
  reference_args <- cv_args
  reference_args$r <- NULL
  reference_args$r_source <- rank
  reference_args$r_target <- rank
  set.seed(514)
  expected <- do.call(reference$cv.learner, reference_args)
  stopifnot(all(is.finite(expected$mse_all)))
  sorted_scores <- sort(as.vector(expected$mse_all))
  stopifnot(sorted_scores[2] - sorted_scores[1] > 1e-5)
  cv_cases[[target_name]] <- list(args = cv_args, seed = 514L, expected = expected)
}
direct_args <- list(Y_source = sources, Y_target = target, r = rank,
                    source_weights = weights)
direct_expected <- do.call(reference$dlearner, direct_args)
# Independently verify D-LEARNER by forming the dense weighted operators.
bases <- lapply(sources, svd, nu = rank, nv = rank)
w <- weights / sum(weights)
pu <- Reduce(`+`, Map(function(b, wk) wk * tcrossprod(b$u), bases, w))
pv <- Reduce(`+`, Map(function(b, wk) wk * tcrossprod(b$v), bases, w))
stopifnot(isTRUE(all.equal(direct_expected$dlearner_estimate,
                          pu %*% target %*% pv, tolerance = 1e-12)))
fixture <- list(
  provenance = list(reference_package = 'learnerv2', reference_version = '0.3.0',
                    reference_md5 = setNames(as.character(tools::md5sum(
                      file.path(reference_dir, reference_files))), reference_files),
                    r_version = as.character(getRversion()),
                    rng_kind = RNGkind(), tolerance = 1e-8),
  fit = fit_cases, cv = cv_cases,
  direct = list(args = direct_args, expected = direct_expected))
saveRDS(fixture, output, version = 2)
cat('Wrote independently generated reference fixture:', output, '\n')
