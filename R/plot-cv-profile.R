#' Plot cross-validation sensitivity curves
#'
#' Compare validation scores across the common space penalty, or across the
#' row-space penalty at each fixed column-space penalty. Each line corresponds
#' to one balance penalty, making nearly overlapping curves visible.
#'
#' @inheritParams plot_cv
#' @return Invisibly, a list of absolute-score diagnostics: \code{scores},
#' \code{scale}, \code{grids}, \code{best_index}, \code{shown_lambda_2_index},
#' \code{color_limits}, \code{unavailable_count}, \code{log_labels},
#' \code{score_range}, \code{minimum}, \code{exact_minima}, \code{boundary},
#' \code{lambda_2_slice_range}, and \code{x_order} (horizontal drawing order).
#' This base-graphics function does not return ggplot objects or the relative
#' display fields returned by \code{plot_cv()}.
#' @details
#' Supply the exact original grids. The horizontal candidates are sorted only
#' for drawing lines; arrays and selected indices retain CV order. Duplicate
#' horizontal candidates are rejected. In separate mode there is one panel per
#' column-space penalty; \code{lambda_2_index} selects the balance-penalty lines.
#' At most four panels are drawn per page. With \code{log_labels=TRUE}, the
#' horizontal coordinate is the actual base-10 logarithm, not a rounded label.
#' All vertical limits use the complete score array. Non-finite values break
#' lines. The graph does not quantify uncertainty or establish that a parameter
#' has no effect; exact ties and boundary selections are reported separately.
#' @examples
#' set.seed(42)
#' cv <- cv.learner(dat_highsim$Y_source[1:12, 1:8],
#'                  dat_highsim$Y_target[1:12, 1:8], r = 2,
#'                  lambda_1_all = c(0.1, 1, 10), lambda_2_all = c(0.1, 1),
#'                  step_size = 0.0003, n_folds = 2, control = list(max_iter = 20))
#' plot_cv_profile(cv, lambda_1_all = c(0.1, 1, 10),
#'                 lambda_2_all = c(0.1, 1), n_observed = 96, log_labels = TRUE)
#' @export
plot_cv_profile <- function(x, lambda_1_all = NULL, lambda_2_all = NULL,
                            lambda_1_row_all = NULL, lambda_1_col_all = NULL,
                            n_observed = NULL, lambda_2_index = NULL, log_labels = FALSE) {
  d <- prepare_cv_plot(x, lambda_1_all, lambda_2_all, lambda_1_row_all,
                       lambda_1_col_all, n_observed, lambda_2_index, log_labels)
  separate <- length(dim(d$scores)) == 3L
  horizontal <- if (separate) d$grids$row else d$grids$common
  if (anyDuplicated(horizontal)) stop('Profile horizontal candidates must be distinct.')
  order_x <- order(horizontal)
  horizontal <- horizontal[order_x]
  if(log_labels) horizontal <- log10(horizontal)
  balance <- if(separate) d$shown_lambda_2_index else seq_along(d$grids$balance)
  panels <- if(separate) seq_along(d$grids$col) else 1L
  colors <- grDevices::colorRampPalette(c('#2166AC','#1B9E77','#D95F02','#984EA3'))(length(d$grids$balance))
  limits <- d$color_limits
  old <- graphics::par(no.readonly=TRUE);on.exit(graphics::par(old),add=TRUE)
  for(page in split(panels,ceiling(seq_along(panels)/4))) {
    nc <- min(2,length(page))
    graphics::par(mfrow=c(ceiling(length(page)/nc),nc),mar=c(3.8,4.5,2.5,7),
                  oma=c(4.2,0,0,0),mgp=c(2.5,0.6,0),cex.axis=0.9,cex.lab=0.95)
    for(j in page) {
      label <- if(separate) expression(lambda[1*','*plain(row)]) else expression(lambda[1])
      if(log_labels)label<-bquote(log[10](.(label[[1]])))
      title <- if(separate) bquote('CV profile at '~lambda[1*','*plain(col)] == .(d$grids$col[j])) else 'Cross-validation sensitivity'
      graphics::plot(horizontal,rep(NA_real_,length(horizontal)),type='n',ylim=limits,
                     xlab=label,ylab=sub('Validation','Held-out',d$scale),main=title,cex.main=1.1)
      graphics::grid(col='gray90',lty=1)
      for(k in balance) {
        y <- if(separate) d$scores[order_x,j,k] else d$scores[order_x,k]
        y[!is.finite(y)]<-NA_real_
        graphics::lines(horizontal,y,type='b',pch=16,cex=0.55,col=colors[k],lty=1+(k-1)%%6,lwd=1.5)
      }
      if(!separate || (j==d$best_index[2] && d$best_index[3] %in% balance)) {
        selected_x<-if(separate)d$grids$row[d$best_index[1]] else d$grids$common[d$best_index[1]]
        if(log_labels)selected_x<-log10(selected_x)
        graphics::points(selected_x,d$minimum,pch=4,cex=1.3,lwd=2)
      }
      usr <- graphics::par('usr')
      graphics::legend(usr[2] + 0.06 * diff(usr[1:2]), usr[4], xjust=0, yjust=1,
                       xpd=NA,bty='n',title=expression(lambda[2]),
                       legend=format(d$grids$balance[balance],digits=3,trim=TRUE),
                       col=colors[balance],lty=1+(balance-1)%%6,pch=16,cex=0.8)
    }
    plot_cv_caption(d)
  }
  d$x_order<-order_x
  invisible(d)
}
