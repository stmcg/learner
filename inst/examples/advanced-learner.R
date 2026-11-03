# Run from top to bottom with Source with Echo in RStudio, or:
# source("advanced-learner.R", echo = TRUE, max.deparse.length = Inf)
# Needs the development version of learner. Install it with:
# devtools::install_github("stmcg/learner")
# This compiles C++ code, so you need a compiler (Xcode Command Line Tools on
# macOS, Rtools on Windows). Restart R after installing.
# Resolve the first error before continuing.

# load-package
install_help <- paste(
  "Install the development version of learner with",
  "devtools::install_github(\"stmcg/learner\"), then restart R and rerun",
  "from the top.")
if (!requireNamespace("learner", quietly = TRUE)) {
  stop(install_help,
       call. = FALSE)
}
required_data <- c("dat_asymmetric", "dat_multisource")
available_data <- utils::data(package = "learner")$results[, "Item"]
missing_data <- setdiff(required_data, available_data)
if (length(missing_data)) {
  stop(paste0("Your installed learner package is missing: ",
              paste(missing_data, collapse = ", "),
              ". ", install_help), call. = FALSE)
}
library(learner)
data("dat_asymmetric", "dat_multisource", package = "learner")

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    utils::packageVersion("ggplot2") < "3.4.0") {
  stop('Install ggplot2 >= 3.4.0 with install.packages("ggplot2"), then rerun.',
       call. = FALSE)
}
# Check behavior rather than relying on a development version number.
compatible <- local({
  y <- diag(c(4, 2, 1))
  target <- diag(c(4, 0, 0))
  cv <- cv.learner(y, target, r = 1, lambda_1_all = 1, lambda_2_all = 1,
                   step_size = 0.001, n_folds = 2, control = list(max_iter = 2))
  all(c("lambda_1_all", "lambda_2_all", "n_observed") %in% names(cv))
})
if (!compatible) {
  stop(paste("This learner installation lacks the CV metadata used here.",
             install_help), call. = FALSE)
}


# asymmetric-data
asym <- dat_asymmetric

# asymmetric-cv
asym_space_grid <- c(0, 10, 100)
asym_balance_grid <- c(0.001, 0.01, 0.1, 1)
set.seed(815)
cv_separate <- cv.learner(
  asym$Y_source, asym$Y_target, r = 3,
  lambda_1_row_all = asym_space_grid,
  lambda_1_col_all = asym_space_grid,
  lambda_2_all = asym_balance_grid,
  n_folds = 3, step_size = 0.003, control = list(max_iter = 2000))
c(U = cv_separate$lambda_1_row_min,
  V = cv_separate$lambda_1_col_min,
  balance = cv_separate$lambda_2_min)

# asymmetric-refit
fit_separate <- learner(
  asym$Y_source, asym$Y_target, r = cv_separate$r,
  lambda_1_row = cv_separate$lambda_1_row_min,
  lambda_1_col = cv_separate$lambda_1_col_min,
  lambda_2 = cv_separate$lambda_2_min,
  step_size = 0.003, control = list(max_iter = 2000))
dim(fit_separate$learner_estimate)
fit_separate$convergence_criterion

# asymmetric-cv-plot
plot_cv(cv_separate)

# multiple-data
sources <- dat_multisource$Y_source
target <- dat_multisource$Y_target
weights <- c(0.7, 0.3)

# multiple-direct
fit_direct <- dlearner(sources, target, r = 3, source_weights = weights)
dim(fit_direct$dlearner_estimate)

# multiple-cv
set.seed(815)
cv_multi <- cv.learner(
  sources, target, r = 3, source_weights = weights,
  lambda_1_all = c(0, 1, 10), lambda_2_all = c(0.1, 1),
  n_folds = 3, step_size = 0.0003, control = list(max_iter = 2000))
c(space = cv_multi$lambda_1_min, balance = cv_multi$lambda_2_min)

# multiple-cv-plot
plot_cv(cv_multi)

# multiple-refit
fit_multi <- learner(
  sources, target, r = cv_multi$r, source_weights = weights,
  lambda_1 = cv_multi$lambda_1_min, lambda_2 = cv_multi$lambda_2_min,
  step_size = 0.0003, control = list(max_iter = 2000))
dim(fit_multi$learner_estimate)
fit_multi$convergence_criterion

# grid-boundaries
needed <- c("asym", "asym_space_grid", "asym_balance_grid",
            "cv_separate", "fit_separate")
if (!all(vapply(needed, exists, logical(1), inherits = TRUE))) {
  stop("Run Section 2 through the final refit before starting Section 4.",
       call. = FALSE)
}
selected <- c(U = cv_separate$lambda_1_row_min,
              V = cv_separate$lambda_1_col_min,
              balance = cv_separate$lambda_2_min)
grids <- list(cv_separate$lambda_1_row_all, cv_separate$lambda_1_col_all,
              cv_separate$lambda_2_all)
data.frame(parameter = names(selected), selected = unname(selected),
           lower = vapply(grids, min, numeric(1)),
           upper = vapply(grids, max, numeric(1)),
           on_boundary = vapply(seq_along(grids), function(i) {
             selected[i] %in% range(grids[[i]])
           }, logical(1)))

# asymmetric-objective
plot_objective(fit_separate, iterations = "near_optimum")

# sensitivity-checks
fit_summary <- function(cv, fit, step, budget) {
  data.frame(step_size = step, max_iter = budget,
             U = cv$lambda_1_row_min, V = cv$lambda_1_col_min,
             balance = cv$lambda_2_min,
             min_cv_mse = min(cv$mse_all) / cv$n_observed,
             updates = length(fit$objective_values),
             stop_code = fit$convergence_criterion)
}
check_setting <- function(step, budget) {
  set.seed(815)
  cv <- cv.learner(
    asym$Y_source, asym$Y_target, r = 3,
    lambda_1_row_all = asym_space_grid,
    lambda_1_col_all = asym_space_grid, lambda_2_all = asym_balance_grid,
    n_folds = 3, step_size = step, control = list(max_iter = budget))
  fit <- learner(
    asym$Y_source, asym$Y_target, r = cv$r,
    lambda_1_row = cv$lambda_1_row_min, lambda_1_col = cv$lambda_1_col_min,
    lambda_2 = cv$lambda_2_min, step_size = step,
    control = list(max_iter = budget))
  fit_summary(cv, fit, step, budget)
}
sensitivity <- rbind(
  fit_summary(cv_separate, fit_separate, 0.003, 2000),
  check_setting(0.003, 4000),
  check_setting(0.0015, 4000))
print(sensitivity, row.names = FALSE)

# inner-cv-stopping
set.seed(815)
indices <- sample(which(!is.na(asym$Y_target)))
n <- length(indices)
folds <- lapply(1:3, function(k) {
  indices[seq.int(floor((k - 1) * n / 3) + 1, floor(k * n / 3))]
})
inner_fits <- do.call(rbind, lapply(1:3, function(k) {
  training <- asym$Y_target
  training[folds[[k]]] <- NA
  fit <- learner(
    asym$Y_source, training, r = cv_separate$r,
    lambda_1_row = cv_separate$lambda_1_row_min,
    lambda_1_col = cv_separate$lambda_1_col_min,
    lambda_2 = cv_separate$lambda_2_min,
    step_size = 0.003, control = list(max_iter = 2000))
  data.frame(fold = k, updates = length(fit$objective_values),
             stop_code = fit$convergence_criterion)
}))
print(inner_fits, row.names = FALSE)
