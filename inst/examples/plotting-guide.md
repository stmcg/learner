# Plotting guide

These functions display existing fits or measured results. Plot appearance alone
cannot establish convergence or statistical superiority.

## Objective history

```r
plot_objective(fit)                                # All iterations
plot_objective(fit, iterations = "near_optimum")   # Around the smallest recorded value
plot_objective(fit, iterations = 1:500)            # A specific range
```

`"near_optimum"` shows up to 41 consecutive records centered on the first
smallest finite recorded objective value: 20 before and 20 after. At either
boundary the window shifts to stay within the history; histories shorter than
41 records are shown in full. Both axes adjust to this window. To choose a
different width, supply numeric iteration indices, such as `79:119`.
The shaded band spans actual adjacent objective values and the line shows their
pairwise average, positioned halfway between the two iterations. Missing values
interrupt the band and line; the band is not an uncertainty interval.
Numeric selections such as `1:1000` keep only recorded iterations. A selection
with no recorded iterations produces an error.

## CV heatmaps

Only the heatmap requires the optional package ggplot2 >= 3.4.0. Install it
with `install.packages("ggplot2")`. Estimation and all other plotting functions
work without ggplot2.

```r
plot_cv(cv, lambda_1_all = common_grid, lambda_2_all = balance_grid,
        n_observed = sum(!is.na(target)), log_labels = TRUE)
```

Use the exact grids supplied to `cv.learner()`. With independent penalties,
replace `lambda_1_all` with `lambda_1_row_all` and `lambda_1_col_all`.
Heatmaps then have one panel per balance penalty. Zero candidates require
`log_labels = FALSE`. Heatmap positions are equally spaced grid slots.

The CV object's `mse_all` contains pooled held-out squared-error sums. Supplying
the total number of observed target entries converts these to MSE. Without that
number the plots explicitly show SSE. Source entries are not in the denominator.
Heatmaps use a reversed magma palette: pale cells have lower error. Labels and
colors show percentage increases over the minimum, with a black selection box.
Finite errors more than 30% above the minimum are gray and labeled >30%, not
as optimizer divergence. Use `relative_limit = NULL` to color all finite errors
or `show_values = FALSE` to hide cell text. If the minimum is zero, the heatmap
uses absolute scores instead. Cell labels are rounded, so matching printed
percentages do not establish an exact tie.

The plots report the score range, minimum, exact ties, and boundary selections.
A boundary minimum deserves investigation; a small change across a parameter
is not an exact tie or evidence that the parameter is irrelevant.

For a customized export, capture the returned object:

```r
result <- plot_cv(cv, lambda_1_all = common_grid, lambda_2_all = balance_grid)
ggplot2::ggsave("cv-heatmap.png", plot = result$plot, width = 8, height = 5.5, dpi = 200)
```

For multiple pages, `result$plots` contains every ggplot page. Inspect exported
images at the intended size to check text and legend placement.

## Population latent spaces

```r
plot_source_spaces(target, list(A = source_A, B = source_B), r = 3)
```

For each population, compute an SVD and retain orthonormal bases B_U and B_V.
Plot B_U B_U^T for the row side and B_V B_V^T for the column side. These are
individual population projectors, not raw data matrices or the weighted average
used by the optimizer. Shared color limits within each side support comparisons.
Negative entries are red and positive entries blue.

SVD uses full matrices before display subsetting. `row_index`, `column_index`,
or `max_display` limit plotting cost without changing the fitted spaces. All
populations must have aligned variables in the same order. Labels default to
Row index and Column index; use biological labels only when appropriate.

## Variable contributions

```r
plot_contributions(fit, r = 3)
```

Recompute the SVD of the estimated signal. For component j, a row's score is
B_U[i,j]^2 and a column's score is B_V[i,j]^2. The scores sum to one per component
on each side over the full set of variables. They are not squared optimization
factors, which include scaling. The plot uses separate color scales on the two
sides and reports display clipping if explicit limits are supplied.
These component-specific scores depend on basis choice in repeated singular
subspaces; the function does not match components or assign biological meaning.
