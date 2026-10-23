test_that('CV heatmaps reject invalid titles and separate arguments in common mode', {
  args <- list(x = cv_plot_example(FALSE), lambda_1_all = c(10, 0),
               lambda_2_all = c(0, 1))
  for (title in list(NA_character_, character(), c('a', 'b'), 1)) {
    expect_error(do.call(plot_cv, c(args, list(main = title))),
                 'main must be one nonmissing character string', fixed = TRUE)
  }
  for (extra in list(list(lambda_1_row_all = c(1, 2)),
                    list(lambda_1_col_all = c(1, 2)),
                    list(lambda_2_index = 1))) {
    expect_error(do.call(plot_cv, c(args, extra)),
                 'Separate-penalty arguments require a 3D score array', fixed = TRUE)
  }
})

test_that('CV heatmaps explain how to install the optional plotting dependency', {
  # Older supported testthat versions do not provide local_mocked_bindings().
  skip_if_not_installed('testthat', minimum_version = '3.2.0')
  original_require <- base::requireNamespace
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, 'ggplot2')) return(FALSE)
      original_require(package, ...)
    }, .package = 'base')
  expect_error(plot_cv(cv_plot_example(FALSE), lambda_1_all = c(10, 0),
                       lambda_2_all = c(0, 1)),
               'Install it with install.packages("ggplot2")', fixed = TRUE)
})

test_that('CV heatmaps explain when the plotting dependency needs updating', {
  skip_if_not_installed('testthat', minimum_version = '3.2.0')
  # Simulate availability/version checks without changing installed packages.
  original_require <- base::requireNamespace
  original_version <- utils::packageVersion
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, 'ggplot2')) return(TRUE)
      original_require(package, ...)
    }, .package = 'base')
  testthat::local_mocked_bindings(
    packageVersion = function(pkg, ...) {
      if (identical(pkg, 'ggplot2')) return(package_version('3.3.6'))
      original_version(pkg, ...)
    }, .package = 'utils')
  expect_error(plot_cv(cv_plot_example(FALSE), lambda_1_all = c(10, 0),
                       lambda_2_all = c(0, 1)),
               'requires ggplot2 >= 3.4.0. Update it', fixed = TRUE)
})
