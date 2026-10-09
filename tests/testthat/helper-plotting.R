# Devices are temporary and closed even when an assertion fails.
plot_device <- function() {
  path <- tempfile(fileext = '.pdf')
  grDevices::pdf(path, width = 10, height = 5)
  path
}

cv_plot_example <- function(separate = TRUE) {
  if (separate) {
    scores <- array(c(8, 7, 6, 5, 4, 3, 2, 1), c(2, 2, 2))
    list(mse_all = scores, lambda_1_min = NA_real_, lambda_1_row_min = 0,
         lambda_1_col_min = 3, lambda_2_min = 1, r = 2)
  } else list(mse_all = matrix(c(8, 7, 6, 5), 2, 2),
              lambda_1_min = 0, lambda_2_min = 1, r = 2)
}
