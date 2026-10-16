#' Plot the objective function values of a LEARNER fit
#'
#' Plots the values of the objective function at each iteration of the
#' numerical optimization algorithm in \code{learner}.
#'
#' @param x A result returned by \code{learner}.
#' @param iterations Iterations to show: \code{"all"} (the default),
#' \code{"near_optimum"}, or a vector of iteration numbers such as
#' \code{1:1000}. \code{"near_optimum"} shows the iterations from the one
#' with the smallest objective value through the last iteration. If the
#' smallest value occurs within the last 10\% of the iterations (and at least
#' the last 10), the last 10\% of the iterations are shown instead. The
#' horizontal axis is the same as for \code{"all"}, except in that case.
#' @param main Title of the plot. By default, \code{"Objective function"},
#' with \code{"(near optimum)"} added when \code{iterations = "near_optimum"}.
#'
#' @return Invisibly, a data frame of the plotted iterations and objective
#' values.
#'
#' @details The red point marks the smallest objective value, which is also
#' given above the plot. The objective value at the initialization is not
#' included. With \code{iterations = "near_optimum"}, the plot shows a shaded
#' band spanning each pair of consecutive objective values, so that values
#' that alternate between iterations remain visible, and a line at their
#' average.
#'
#' @examples
#' fit <- learner(dat_highsim$Y_source, dat_highsim$Y_target,
#'                r = 3, lambda_1 = 1, lambda_2 = 1, step_size = 0.003,
#'                control = list(max_iter = 1000))
#' plot_objective(fit)
#' plot_objective(fit, iterations = "near_optimum")
#' plot_objective(fit, iterations = 1:100)
#'
#' @export
plot_objective <- function(x, iterations = "all", main = NULL) {
  objective <- x$objective_values
  if (!is.numeric(objective) || length(objective) == 0L ||
      !any(is.finite(objective))) {
    stop('x must be a result returned by learner().')
  }
  n <- length(objective)
  finite <- is.finite(objective)
  minimum <- which.min(replace(objective, !finite, Inf))
  near <- identical(iterations, 'near_optimum')
  tail_view <- FALSE
  if (identical(iterations, 'all')) {
    shown <- seq_len(n)
  } else if (near) {
    # If the minimum is the end of the history, show its final stretch instead.
    tail_length <- max(10L, ceiling(0.1 * n))
    start <- min(minimum, max(1L, n - tail_length + 1L))
    tail_view <- start != minimum
    shown <- seq.int(start, n)
  } else if (is.numeric(iterations) && length(iterations) > 0L &&
             all(iterations %in% seq_len(n))) {
    shown <- sort(unique(as.integer(iterations)))
  } else {
    stop("iterations must be 'all', 'near_optimum', or iteration numbers between 1 and ",
         n, '.')
  }
  if (!any(finite[shown])) {
    stop('The selected iterations contain no finite objective values; choose a different range.')
  }
  if (is.null(main)) {
    main <- if (near) 'Objective function (near optimum)' else 'Objective function'
  }
  values <- replace(objective[shown], !finite[shown], NA_real_)

  old <- graphics::par(mar = c(5, 4.5, 4, 1.5))
  on.exit(graphics::par(old), add = TRUE)
  # Near the optimum, keep the horizontal axis of the full history so that the
  # two plots line up; only the vertical axis is zoomed.
  xlim <- if (near && !tail_view) c(1, n) else range(shown)
  graphics::plot(shown, values, type = 'n', xlim = xlim, xlab = 'Iteration',
                 ylab = 'Objective value', bty = 'l')
  graphics::title(main, line = 2)
  graphics::mtext(paste0('Minimum: ', format(objective[minimum], digits = 6),
                         ' at iteration ', minimum),
                  side = 3, line = 0.6, cex = 0.8, col = 'gray30')
  if (near && length(shown) > 1L) {
    # The band spans each pair of consecutive values, so alternating values
    # appear as a band rather than as dense vertical lines.
    nxt <- c(values[-1], values[length(values)])  # the last value pairs with itself
    keep <- !is.na(values) & !is.na(nxt)
    graphics::polygon(c(shown[keep], rev(shown[keep])),
                      c(pmin(values, nxt)[keep], rev(pmax(values, nxt)[keep])),
                      col = grDevices::adjustcolor('#2166AC', 0.4), border = NA)
    graphics::lines(shown[keep], ((values + nxt) / 2)[keep], col = '#2166AC', lwd = 1.5)
  } else {
    graphics::lines(shown, values, col = '#2166AC', type = if (length(shown) == 1L) 'p' else 'l')
  }
  if (minimum %in% shown) {
    graphics::points(minimum, objective[minimum], pch = 16, cex = 0.9, col = '#B2182B')
  }
  # Symbol keys distinguish the band, mean line, and minimum without colored text.
  band_color <- grDevices::adjustcolor('#2166AC', 0.4)
  labels <- if (near && length(shown) > 1L)
    c('Consecutive-iteration range', 'Pairwise average') else 'Objective value'
  colors <- if (near && length(shown) > 1L) c(band_color, '#2166AC') else '#2166AC'
  symbols <- if (near && length(shown) > 1L) c(22, NA) else NA
  line_types <- if (near && length(shown) > 1L) c(0, 1) else 1
  if (minimum %in% shown) {
    labels <- c(labels, 'Minimum recorded objective')
    colors <- c(colors, '#B2182B')
    symbols <- c(symbols, 16)
    line_types <- c(line_types, 0)
  }
  graphics::legend('topright', inset = 0.02, legend = labels,
                   col = colors, text.col = 'gray30', pch = symbols,
                   pt.bg = band_color, pt.cex = 1.1, lty = line_types,
                   lwd = 1.5, cex = 0.8, y.intersp = 1.3, bty = 'n')
  invisible(data.frame(iteration = shown, objective = objective[shown]))
}

#' Plot the LEARNER cross-validation scores
#'
#' Draw a common-penalty heatmap or separate-penalty heatmaps, with a shared
#' color scale and a black box at the first minimum in R array order.
#'
#' @param x A result returned by \code{cv.learner()}.
#' @param lambda_1_all Common-penalty candidates in the exact order used for CV.
#' Required for a two-dimensional score array; omit for separate penalties.
#' @param lambda_2_all Balance-penalty candidates in the exact order used for CV.
#' @param lambda_1_row_all,lambda_1_col_all Separate-penalty candidates in the
#' exact order used for CV. Both are required for a three-dimensional array.
#' @param n_observed Optional positive integer giving the number of observed
#' target entries used in CV. If supplied, scores are divided by this count and
#' labeled Validation MSE. Otherwise the original sums are labeled Validation
#' SSE. Missing target entries must not be included in this count.
#' @param lambda_2_index Optional vector of distinct integer indices into
#' \code{lambda_2_all}, selecting panels to draw in separate-penalty mode.
#' \code{NULL} draws all panels. Not used in common-penalty mode.
#'
#' @param log_labels Logical scalar. If TRUE, display base-10 logarithms of
#' the candidates on the heatmap axes.
#' All candidates on those axes must be strictly positive. Log labels are
#' rounded to two decimal places. Defaults to FALSE.
#'
#' @param relative_limit Optional positive relative increase cutoff, default 0.3
#' (30 percent). Larger finite errors are gray and labeled above the cutoff,
#' never as diverged. NULL keeps all finite values in the color scale.
#' @param show_values Whether to print scores inside cells. Defaults to TRUE.
#' @param main Plot title.
#'
#' @return Invisibly, a list with \code{scores}, \code{scale}, \code{grids},
#' \code{best_index}, \code{shown_lambda_2_index}, \code{color_limits}, and
#' \code{unavailable_count}, \code{log_labels}, \code{score_range}, \code{minimum},
#' \code{exact_minima}, \code{boundary}, and \code{lambda_2_slice_range}.
#' Indices refer to the original candidate order. Boundary flags indicate a
#' selected numerical endpoint on a grid with at least two distinct values.
#' Also returns \code{display_data}, \code{display_limits}, \code{relative},
#' \code{above_limit_count}, \code{plot} (the first ggplot page), and \code{plots}
#' (all ggplot pages). The original absolute score diagnostics are unchanged.
#'
#' @details
#' Existing CV results do not store the full candidate grids, so callers must
#' provide them. Lengths and the selected penalty values are checked, but these
#' checks cannot detect every incorrectly supplied grid. Use the original grids,
#' without sorting or removing duplicates. Axes are equally spaced candidate
#' positions labeled with actual values, or their base-10 logarithms when
#' \code{log_labels=TRUE}. This changes labels only, not scores, candidate order,
#' or cell spacing. The reversed magma palette uses pale colors for lower error
#' and dark colors for higher error. Cell labels show increases over the minimum.
#'
#' Common mode shows \code{lambda_1_all} on the horizontal axis and
#' \code{lambda_2_all} on the vertical axis. Separate mode shows the U-space
#' penalty horizontally and the V-space penalty vertically, with one panel for
#' each balance penalty. All panels share limits computed from the full score
#' array, even when only some panels are selected. Up to six panels are drawn
#' per page. Gray cells have non-finite scores and are excluded from selection.
#' Only the first minimum is marked; a panel need not contain a black box.
#' Captions report the score range, exact ties, and boundary selections. These are
#' descriptive diagnostics, not significance or convergence tests.
#'
#' The stored \code{mse_all} values are sums of held-out squared errors. MSE
#' conversion requires the correct \code{n_observed}; the function cannot infer
#' this count from the scores. Neither this plot nor CV certifies convergence
#' of the individual fits. No fitting is performed and \code{x} is unchanged.
#' Requires optional ggplot2 >= 3.4.0; install it with
#' \code{install.packages("ggplot2")}. Estimation and other plots do not require it.
#' Draws with ggplot2 without modifying base graphics parameters. Each page shares
#' one legend. A sufficiently large device is recommended for multiple panels.
#' Relative increases are undefined when the minimum is zero; the display then
#' uses absolute scores and ignores relative_limit. High relative error alone
#' never establishes optimizer divergence. Non-finite scores are unavailable.
#' Use ggplot2::ggsave(plot = result$plot, ...) to export the first page, or
#' export each member of result$plots for multiple pages.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     utils::packageVersion("ggplot2") >= "3.4.0") {
#' row_grid <- c(0.1, 1)
#' col_grid <- c(0.1, 2)
#' balance_grid <- c(0.1, 1)
#' target <- dat_highsim$Y_target[1:12, 1:8]
#' set.seed(42)
#' cv <- cv.learner(dat_highsim$Y_source[1:12, 1:8], target,
#'                  r = 2, lambda_1_row_all = row_grid,
#'                  lambda_1_col_all = col_grid, lambda_2_all = balance_grid,
#'                  step_size = 0.0003, n_folds = 2,
#'                  control = list(max_iter = 20))
#' plot_cv(cv, lambda_1_row_all = row_grid, lambda_1_col_all = col_grid,
#'         lambda_2_all = balance_grid, n_observed = sum(!is.na(target)))
#'
#' }
#'
#' @export
plot_cv <- function(x, lambda_1_all = NULL, lambda_2_all = NULL,
                    lambda_1_row_all = NULL, lambda_1_col_all = NULL,
                    n_observed = NULL, lambda_2_index = NULL, log_labels = FALSE,
                    relative_limit = 0.3, show_values = TRUE,
                    main = 'Cross-validation error') {
  d <- prepare_cv_plot(x, lambda_1_all, lambda_2_all, lambda_1_row_all,
                       lambda_1_col_all, n_observed, lambda_2_index, log_labels)
  if (!is.null(relative_limit) && (!is.numeric(relative_limit) ||
      is.complex(relative_limit) || length(relative_limit) != 1L ||
      !is.finite(relative_limit) || relative_limit <= 0)) {
    stop('relative_limit must be NULL or a finite positive number.')
  }
  if (!is.logical(show_values) || length(show_values) != 1L || is.na(show_values)) {
    stop('show_values must be TRUE or FALSE.')
  }
  if (!is.character(main) || length(main) != 1L || is.na(main)) {
    stop('main must be one nonmissing character string.')
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop('plot_cv() requires the optional ggplot2 package. Install it with install.packages("ggplot2").',
         call. = FALSE)
  }
  if (utils::packageVersion("ggplot2") < "3.4.0") {
    stop('plot_cv() requires ggplot2 >= 3.4.0. Update it with install.packages("ggplot2").',
         call. = FALSE)
  }
  separate <- length(dim(d$scores)) == 3L
  indices <- as.data.frame(arrayInd(seq_along(d$scores), dim(d$scores)))
  names(indices) <- if (separate) c('x_index', 'y_index', 'slice') else c('x_index', 'y_index')
  if (!separate) indices$slice <- 1L
  indices$score <- as.vector(d$scores)
  indices$available <- is.finite(indices$score)
  # Relative increases are undefined when the best score is zero.
  relative <- d$minimum > 0
  indices$display <- if (relative) indices$score / d$minimum - 1 else indices$score
  indices$above_limit <- indices$available & relative & !is.null(relative_limit)
  if (!is.null(relative_limit)) indices$above_limit <- indices$above_limit & indices$display > relative_limit + 8 * .Machine$double.eps * max(1, relative_limit)
  indices$fill_value <- indices$display
  indices$fill_value[!indices$available | indices$above_limit] <- NA_real_
  finite_fill <- indices$fill_value[is.finite(indices$fill_value)]
  upper <- max(finite_fill)
  limits <- c(0, if (upper > 0) upper else if (relative) 0.01 else 1e-8)
  # One global scale across all slices, even if only a subset is displayed.
  label_digits <- if (limits[2] <= 0.05) 2L else 1L
  indices$label <- if (relative) sprintf(paste0('%.', label_digits, 'f%%'), 100 * indices$display) else
    formatC(indices$display, digits = 4, format = 'g')
  indices$label[!indices$available] <- 'NA'
  if (any(indices$above_limit)) indices$label[indices$above_limit] <-
    paste0('>', format(100 * relative_limit, trim = TRUE), '%')
  indices$text_color <- ifelse(is.finite(indices$fill_value) &
                                indices$fill_value > limits[2] * 0.55, 'white', '#202020')
  best_slice <- if (separate) d$best_index[3] else 1L
  indices$selected <- indices$x_index == d$best_index[1] &
    indices$y_index == d$best_index[2] & indices$slice == best_slice
  panels <- if (separate) d$shown_lambda_2_index else 1L
  labels_x <- if (separate) d$grids$row else d$grids$common
  labels_y <- if (separate) d$grids$col else d$grids$balance
  xlab <- if (separate) expression(lambda[1*','*plain(row)]) else expression(lambda[1])
  ylab <- if (separate) expression(lambda[1*','*plain(col)]) else expression(lambda[2])
  if (log_labels) {
    labels_x <- round(log10(labels_x), 2)
    labels_y <- round(log10(labels_y), 2)
    xlab <- bquote(log[10](.(xlab[[1]])))
    ylab <- bquote(log[10](.(ylab[[1]])))
  }
  metric <- if (is.null(n_observed)) 'SSE' else 'MSE'
  selected_values <- vapply(seq_along(d$grids), function(i) d$grids[[i]][d$best_index[i]], numeric(1))
  selected_names <- if (separate) c('row', 'col', 'balance') else c('lambda1', 'lambda2')
  subtitle <- paste0('Selected ', paste(paste0(selected_names, ' = ',
                    trimws(formatC(selected_values, digits = 4, format = 'g'))), collapse = ', '),
                    '  |  Minimum ', metric, ' = ', formatC(d$minimum, digits = 6, format = 'g'))
  caption <- paste0('Black box: selected parameters. Cell labels rounded.\n', metric, ' range: ',
                    formatC(d$score_range[1], digits = 6, format = 'g'), ' to ',
                    formatC(d$score_range[2], digits = 6, format = 'g'), '.')
  if (any(d$boundary)) caption <- paste0(caption, ' Selection on grid boundary.')
  if (d$exact_minima > 1) caption <- paste0(caption, '\n', d$exact_minima,
                               ' exact minima; first in CV order selected.')
  if (any(indices$above_limit)) caption <- paste0(caption, '\nGray >',
       format(100 * relative_limit, trim = TRUE), '%: high error, not a divergence diagnosis.')
  if (d$unavailable_count > 0) caption <- paste0(caption, '\nGray NA: unavailable score.')
  if (!relative) caption <- paste0(caption, '\nMinimum is zero; colors and labels show absolute ', metric, '.')
  breaks <- seq(0, limits[2], length.out = 5L)
  key_labels <- if (relative) sprintf(paste0('%.', label_digits, 'f%%'), 100 * breaks) else
    formatC(breaks, digits = 4, format = 'g')
  # Declare aesthetic bindings for package checks; ggplot evaluates them in data.
  x_index <- y_index <- fill_value <- label <- text_color <- NULL
  plots <- list()
  for (page in split(panels, ceiling(seq_along(panels) / 6))) {
    data <- indices[indices$slice %in% page, , drop = FALSE]
    facet_labels <- paste0('lambda[2] == ', format(d$grids$balance[page], digits = 4, trim = TRUE))
    if (anyDuplicated(facet_labels)) {
      facet_labels <- paste0('lambda[2] == ', format(d$grids$balance[page], digits = 17, trim = TRUE))
    }
    # Grid indices identify panels; rounded display labels must not merge them.
    data$panel <- factor(data$slice, levels = page)
    panel_labels <- stats::setNames(facet_labels, as.character(page))
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = x_index, y = y_index, fill = fill_value)) +
      ggplot2::geom_tile(color = 'white', linewidth = 0.35, width = 1, height = 1) +
      ggplot2::scale_fill_viridis_c(option = 'magma', direction = -1,
          na.value = 'grey85', limits = limits, breaks = breaks, labels = key_labels,
          name = if (relative) paste0('CV ', metric, '\nabove minimum') else paste('CV', metric)) +
      ggplot2::scale_x_continuous(breaks = seq_along(labels_x), labels = labels_x,
                                 expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_y_continuous(breaks = seq_along(labels_y), labels = labels_y,
                                 expand = ggplot2::expansion(mult = 0)) +
      ggplot2::labs(title = main, subtitle = subtitle, x = xlab, y = ylab, caption = caption) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(face = 'bold', size = 14),
          plot.subtitle = ggplot2::element_text(size = 10, margin = ggplot2::margin(b = 10)),
          plot.caption = ggplot2::element_text(size = 9, hjust = 0, margin = ggplot2::margin(t = 10)),
          legend.title = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = 9),
          axis.text = ggplot2::element_text(color = '#303030', size = 10),
          strip.text = ggplot2::element_text(size = 11),
          plot.margin = ggplot2::margin(12, 18, 10, 12))
    if (show_values) plot <- plot +
      ggplot2::geom_text(ggplot2::aes(label = label, color = text_color), size = 2.8) +
      ggplot2::scale_color_identity()
    plot <- plot + ggplot2::geom_tile(data = data[data$selected, , drop = FALSE],
                       fill = NA, color = 'black', linewidth = 0.9, width = 0.94, height = 0.94)
    if (separate) plot <- plot + ggplot2::facet_wrap(~panel, ncol = min(3L, length(page)),
                                                    labeller = ggplot2::as_labeller(panel_labels, default = ggplot2::label_parsed))
    print(plot)
    plots[[length(plots) + 1L]] <- plot
  }
  d$display_data <- indices
  d$display_limits <- limits
  d$relative <- relative
  d$above_limit_count <- sum(indices$above_limit)
  d$plots <- plots
  d$plot <- plots[[1L]]
  invisible(d)
}

prepare_cv_plot <- function(x, lambda_1_all = NULL, lambda_2_all = NULL,
                    lambda_1_row_all = NULL, lambda_1_col_all = NULL,
                    n_observed = NULL, lambda_2_index = NULL, log_labels = FALSE) {
  if (!is.list(x) || !is.numeric(x$mse_all) || is.complex(x$mse_all) ||
      !(length(dim(x$mse_all)) %in% c(2L, 3L)) || any(dim(x$mse_all) < 1L)) {
    stop('x must contain a nonempty two- or three-dimensional mse_all array.')
  }
  if (!is.logical(log_labels) || length(log_labels) != 1L || is.na(log_labels)) {
    stop('log_labels must be TRUE or FALSE.')
  }
  scores <- x$mse_all
  finite <- is.finite(scores)
  if (!any(finite)) stop('No finite validation scores are available to plot.')
  if (any(scores[finite] < 0)) stop('Validation squared-error scores cannot be negative.')
  separate <- length(dim(scores)) == 3L
  check_grid <- function(values, size, name) {
    validate_penalty_grid(values, name)
    if (!is.null(dim(values)) || length(values) != size) {
      stop(name, ' must match its score-array dimension in the original CV order.')
    }
  }
  if (separate) {
    if (!is.null(lambda_1_all)) stop('Do not supply lambda_1_all for separate penalties.')
    check_grid(lambda_1_row_all, dim(scores)[1], 'lambda_1_row_all')
    check_grid(lambda_1_col_all, dim(scores)[2], 'lambda_1_col_all')
    check_grid(lambda_2_all, dim(scores)[3], 'lambda_2_all')
    grids <- list(row = lambda_1_row_all, col = lambda_1_col_all,
                  balance = lambda_2_all)
    panels <- if (is.null(lambda_2_index)) seq_along(lambda_2_all) else lambda_2_index
    if (!is.numeric(panels) || length(panels) == 0L || any(!is.finite(panels)) ||
        any(panels != floor(panels)) || any(panels < 1 | panels > length(lambda_2_all)) ||
        anyDuplicated(panels)) {
      stop('lambda_2_index must contain distinct valid balance-grid indices.')
    }
  } else {
    if (!is.null(lambda_1_row_all) || !is.null(lambda_1_col_all) ||
        !is.null(lambda_2_index)) stop('Separate-penalty arguments require a 3D score array.')
    check_grid(lambda_1_all, dim(scores)[1], 'lambda_1_all')
    check_grid(lambda_2_all, dim(scores)[2], 'lambda_2_all')
    grids <- list(common = lambda_1_all, balance = lambda_2_all)
    panels <- 1L
  }
  axis_grids <- if (separate) c(grids$row, grids$col) else c(grids$common, grids$balance)
  if (log_labels && any(axis_grids <= 0)) {
    stop('log_labels requires strictly positive candidates on both heatmap axes.')
  }
  selectable <- scores
  selectable[!finite] <- Inf
  best <- as.integer(arrayInd(which.min(selectable), dim(scores)))
  expected <- if (separate) {
    c(lambda_1_row_min = grids$row[best[1]], lambda_1_col_min = grids$col[best[2]],
      lambda_2_min = grids$balance[best[3]])
  } else c(lambda_1_min = grids$common[best[1]], lambda_2_min = grids$balance[best[2]])
  for (name in names(expected)) {
    if (!is.numeric(x[[name]]) || length(x[[name]]) != 1L ||
        !isTRUE(all.equal(unname(x[[name]]), unname(expected[[name]]), tolerance = 1e-12))) {
      stop('Supplied grids do not match the selected penalties stored in x: ', name)
    }
  }
  scale <- 'Validation SSE'
  if (!is.null(n_observed)) {
    validate_count(n_observed, 'n_observed', upper = Inf)
    scores <- scores / n_observed
    scale <- 'Validation MSE'
  }
  limits <- range(scores[finite])
  if (diff(limits) == 0) {
    pad <- max(abs(limits[1]) * 0.01, 1e-8)
    limits <- c(max(0, limits[1] - pad), limits[2] + pad)
  }
  score_range <- range(scores[finite])
  minimum <- min(scores[finite])
  boundary <- vapply(seq_along(grids), function(i) {
    g <- grids[[i]]
    length(unique(g)) > 1L && g[best[i]] %in% range(g)
  }, logical(1))
  names(boundary) <- names(grids)
  slice <- if (separate) scores[best[1], best[2], ] else scores[best[1], ]
  list(scores = scores, scale = scale, grids = grids, best_index = best,
       shown_lambda_2_index = if (separate) panels else NULL,
       color_limits = limits, unavailable_count = sum(!finite), log_labels = log_labels,
       score_range = score_range, minimum = minimum, exact_minima = sum(scores[finite] == minimum),
       boundary = boundary, lambda_2_slice_range = range(slice[is.finite(slice)]))
}

# Caption reports numerical facts rather than inferring a flat direction from color.
plot_cv_caption <- function(d) {
  number <- function(z) formatC(z, digits = 6, format = 'fg', flag = '#')
  relative <- if (d$minimum > 0) paste0(' (span ', formatC(100 * diff(d$score_range) / d$minimum, digits = 2, format = 'f'), '% of minimum)') else ''
  line1 <- paste0('Range: ', number(d$score_range[1]), ' to ', number(d$score_range[2]), relative,
                   '; selected minimum = ', number(d$minimum), '.')
  boundary <- names(d$boundary)[d$boundary]
  line2 <- paste0(if(length(boundary)) paste0('Grid boundary: ', paste(boundary, collapse = ', '), '; ') else '',
                   'lambda_2 slice span at selected space penalty: ',
                   formatC(diff(d$lambda_2_slice_range), digits = 3, format = 'g'), '.')
  line3 <- if (d$exact_minima > 1) paste(d$exact_minima, 'exact minima; first in CV order selected.') else 'Unique minimum among evaluated candidates.'
  line3 <- paste('X: selected parameters.', line3,
                 if(d$unavailable_count > 0) 'Gray: unavailable.' else '')
  for (i in 1:3) graphics::mtext(c(line1,line2,line3)[i], outer = TRUE,
                                 side = 1, line = 0.4 + (i-1)*1.1, cex = 0.78)
}

# Five ticks include the numerical endpoints even for a narrow score range.
draw_color_key <- function(left, right, bottom, top, limits, colors) {
  edges <- seq(bottom, top, length.out=length(colors)+1L)
  graphics::rect(left,edges[-length(edges)],right,edges[-1],col=colors,border=NA)
  ticks <- seq(limits[1],limits[2],length.out=5)
  graphics::axis(4,pos=right,at=seq(bottom,top,length.out=5),
                 labels=formatC(ticks,format='g',digits=5),las=1,cex.axis=0.9,tck=-0.015)
}
