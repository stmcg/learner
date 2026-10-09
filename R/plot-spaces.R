# Internal helpers shared by the space and contribution plots.
check_plot_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || any(dim(x)==0) || any(!is.finite(x))) {
    stop(name, ' must be a complete finite real numeric matrix.')
  }
}

plot_indices <- function(index, n, max_display, name) {
  if(is.null(index)) return(unique(as.integer(round(seq(1,n,length.out=min(n,max_display))))))
  if(!is.numeric(index)||is.complex(index)||length(index)==0||any(!is.finite(index))||
     any(index!=floor(index))||any(index<1|index>n)||anyDuplicated(index)) {
    stop(name, ' must contain distinct valid integer indices.')
  }
  as.integer(index)
}

plot_svd <- function(x, r) {
  s<-svd(x,nu=r,nv=r)
  if(s$d[r] <= max(dim(x))*s$d[1]*.Machine$double.eps) {
    stop('r includes numerically zero singular directions; use a smaller rank.')
  }
  s
}

plot_limits_pair <- function(limits, automatic) {
  if(is.null(limits))return(automatic)
  if(!is.numeric(limits)||is.complex(limits)||length(limits)!=2L||
     any(!is.finite(limits))||any(limits<=0))stop('limits must contain two finite positive values: row then column.')
  stats::setNames(as.numeric(limits),c('row','column'))
}

# Use a raster to avoid drawing hundreds of thousands of individual rectangles.
draw_matrix_panel <- function(values, x_index, y_index, limits, colors, title, xlab, ylab) {
  nx<-nrow(values);ny<-ncol(values)
  clipped<-pmax(limits[1],pmin(limits[2],values))
  bins<-1L+floor((clipped-limits[1])/diff(limits)*(length(colors)-1L))
  pixels<-matrix(colors[pmax(1L,pmin(length(colors),bins))],nx,ny)
  graphics::plot.new()
  # Reserve a constant fraction of the panel width for the color key.
  graphics::plot.window(xlim=c(0.5,nx+0.5+nx*0.15),ylim=c(0.5,ny+0.5),xaxs='i',yaxs='i')
  graphics::rasterImage(grDevices::as.raster(t(pixels[,ny:1,drop=FALSE])),
                        0.5,0.5,nx+0.5,ny+0.5,interpolate=FALSE)
  graphics::rect(0.5,0.5,nx+0.5,ny+0.5,border='gray40',lwd=0.7)
  ix<-unique(as.integer(round(seq(1,nx,length.out=min(5,nx)))))
  iy<-unique(as.integer(round(seq(1,ny,length.out=min(8,ny)))))
  graphics::axis(1,at=ix,labels=x_index[ix],cex.axis=0.85)
  graphics::axis(2,at=iy,labels=y_index[iy],las=1,cex.axis=0.85)
  graphics::title(main=title,xlab=xlab,ylab=ylab,cex.main=1.05,cex.lab=0.9)
  draw_color_key(nx+0.5+nx*0.05,nx+0.5+nx*0.075,0.5,ny+0.5,limits,colors)
}

#' Compare target and source latent spaces
#'
#' Plot column-space and row-space projection submatrices for the target and
#' each source, following the layout of Figure 4 in the LEARNER paper.
#'
#' @param Y_target Complete finite target matrix, or a completed estimate.
#' Missing values are rejected; the function does not impute them.
#' @param Y_source Aligned complete source matrix or nonempty list of matrices.
#' @param r Positive integer giving the retained rank for every matrix.
#' @param row_index,column_index Optional distinct integer indices selecting
#' variables to display, in display order. The same indices are used in all
#' populations. SVD always uses the full matrices before subsetting.
#' @param max_display Maximum number of displayed variables on each side when
#' indices are omitted; default 200. Deterministic evenly spaced indices are
#' used, without changing the random-number state.
#' @param limits Optional two positive color limits, in row-side then column-side
#' order. Each side uses a symmetric scale about zero shared across populations.
#' Values beyond supplied limits are clipped for display and counts are returned.
#' @param row_label,column_label Axis labels. Use Variant index and Phenotype
#' index only when these accurately describe the input rows and columns.
#' @param target_name Panel label describing the supplied target matrix.
#' @return Invisibly, a list with \code{row_projectors}, \code{column_projectors},
#' \code{row_index}, \code{column_index}, \code{limits}, and \code{clipped_count}.
#' Each projector is the displayed submatrix, computed without allocating the
#' full row-by-row projector. Its population name identifies the input.
#' @details
#' For an SVD basis B, the displayed matrix is a submatrix of B B-transpose.
#' Column-side panels appear on top, row-side panels below. Red denotes negative
#' values and blue positive values. Color scales are shared within each side,
#' including across pages, but can differ between the two sides. At most three
#' populations appear per page. Source names are taken from the input list when
#' supplied. No name-based alignment or matching of factors is performed.
#' These are individual population projectors, not the weighted source operator.
#' A submatrix need not itself be an idempotent projector.
#' Numerically zero retained singular directions are rejected. With repeated
#' singular values at the truncation boundary the retained space is not unique.
#' The input order must already encode any desired phenotype or genomic ordering.
#' @examples
#' a <- dat_highsim$Y_source[1:20,1:10]
#' y <- dat_highsim$Y_target[1:20,1:10]
#' plot_source_spaces(y, list(A=a, B=a*1.1), r=2, max_display=20)
#' @export
plot_source_spaces <- function(Y_target, Y_source, r, row_index=NULL, column_index=NULL,
                               max_display=200L, limits=NULL, row_label='Row index',
                               column_label='Column index', target_name='Target population') {
  check_plot_matrix(Y_target,'Y_target')
  sources<-prepare_sources(Y_source,Y_target,allow_missing_target=FALSE)$matrices
  for(s in sources)check_plot_matrix(s,'Y_source')
  r<-resolve_rank(r,Y_target);validate_count(max_display,'max_display')
  ri<-plot_indices(row_index,nrow(Y_target),max_display,'row_index')
  ci<-plot_indices(column_index,ncol(Y_target),max_display,'column_index')
  labels<-names(sources)
  if(is.null(labels))labels<-rep('',length(sources))
  labels[is.na(labels)]<-''
  labels<-ifelse(nzchar(labels),paste('Source',labels),paste('Source',seq_along(sources)))
  populations<-c(list(Y_target),sources);names(populations)<-make.unique(c(target_name,labels))
  bases<-lapply(populations,plot_svd,r=r)
  up<-lapply(bases,function(s)tcrossprod(s$u[ri,,drop=FALSE]))
  vp<-lapply(bases,function(s)tcrossprod(s$v[ci,,drop=FALSE]))
  auto<-c(row=max(vapply(up,function(p)max(abs(p)),numeric(1))),
          column=max(vapply(vp,function(p)max(abs(p)),numeric(1))))
  auto[auto==0]<-1e-8
  cap<-plot_limits_pair(limits,auto)
  counts<-list(row=vapply(up,function(p)sum(abs(p)>cap[1]),integer(1)),
               column=vapply(vp,function(p)sum(abs(p)>cap[2]),integer(1)))
  colors<-grDevices::colorRampPalette(c('firebrick3','white','dodgerblue3'))(101)
  old<-graphics::par(no.readonly=TRUE);on.exit(graphics::par(old),add=TRUE)
  for(page in split(seq_along(populations),ceiling(seq_along(populations)/3))) {
    graphics::par(mfrow=c(2,length(page)),mar=c(3.5,3.8,2.5,4),oma=c(2,0,0,0),mgp=c(2.2,0.6,0))
    for(i in page)draw_matrix_panel(vp[[i]],ci,ci,c(-cap[2],cap[2]),colors,names(populations)[i],column_label,column_label)
    for(i in page)draw_matrix_panel(up[[i]],ri,ri,c(-cap[1],cap[1]),colors,names(populations)[i],row_label,row_label)
    caption<-paste('Red: negative; blue: positive. Shared scales within each row.',
                   if(sum(unlist(counts))>0)paste(sum(unlist(counts)),'displayed entries clipped across all panels.')else '')
    graphics::mtext(caption,outer=TRUE,side=1,line=0.5,cex=0.8)
  }
  invisible(list(row_projectors=up,column_projectors=vp,row_index=ri,column_index=ci,
                 limits=cap,clipped_count=counts))
}

#' Plot variable contributions to latent components
#'
#' Show squared singular-vector entries for a fitted target matrix, following
#' Figure 5 of the LEARNER paper. Higher contribution scores are darker red.
#'
#' @param x A complete estimated matrix, a \code{learner()} result, or a
#' \code{dlearner()} result. No rotation of the singular vectors is performed.
#' @param r Number of components. Required for a matrix; defaults to the stored
#' rank for a fitting result. Numerically zero directions are rejected.
#' @inheritParams plot_source_spaces
#' @param main Title describing the input estimate.
#' @return Invisibly, a list with full \code{row_scores} and \code{column_scores},
#' displayed \code{row_index} and \code{column_index}, \code{limits}, and
#' \code{clipped_count}. Score matrices have variables in rows and components
#' in columns; all values are returned without display clipping.
#' @details
#' An SVD of the estimated matrix supplies U and V. Row scores are U squared
#' elementwise; column scores are V squared elementwise. Scores sum to one per
#' component over all variables on the corresponding side. They are not squared
#' entries of the raw optimization factors, effect sizes, or correlations.
#' Subsetting is applied only for display and scores are not renormalized.
#' The left panel displays columns and the right panel rows. Color limits are
#' separate for the two sides and start at zero. Optional limits specify row
#' then column upper limits, with clipped-entry counts reported.
#' Signs do not affect squared scores, but rotations within repeated singular
#' subspaces can change component-specific scores. Components from different
#' populations are not automatically aligned or given biological labels.
#' @examples
#' fit <- dlearner(dat_highsim$Y_source,dat_highsim$Y_target,r=3)
#' plot_contributions(fit,max_display=50)
#' @export
plot_contributions <- function(x, r=NULL, row_index=NULL, column_index=NULL,
                               max_display=200L, limits=NULL, row_label='Row index',
                               column_label='Column index', main='Estimated target contributions') {
  if(is.list(x)) {
    estimate<-if(!is.null(x$learner_estimate))x$learner_estimate else x$dlearner_estimate
    if(is.null(r))r<-x$r
  } else estimate<-x
  check_plot_matrix(estimate,'x')
  if(is.null(r))stop('Supply r for a matrix without a stored rank.')
  r<-resolve_rank(r,estimate);validate_count(max_display,'max_display')
  ri<-plot_indices(row_index,nrow(estimate),max_display,'row_index')
  ci<-plot_indices(column_index,ncol(estimate),max_display,'column_index')
  s<-plot_svd(estimate,r)
  row_scores<-s$u^2;col_scores<-s$v^2
  cap<-plot_limits_pair(limits,c(row=max(row_scores[ri,,drop=FALSE]),column=max(col_scores[ci,,drop=FALSE])))
  cap[cap==0]<-1e-8
  counts<-c(row=sum(row_scores[ri,,drop=FALSE]>cap[1]),column=sum(col_scores[ci,,drop=FALSE]>cap[2]))
  colors<-grDevices::colorRampPalette(c('white','firebrick3'))(100)
  old<-graphics::par(no.readonly=TRUE);on.exit(graphics::par(old),add=TRUE)
  graphics::par(mfrow=c(1,2),mar=c(3.8,3.8,2.2,4),oma=c(2,0,2,0),mgp=c(2.2,0.6,0))
  draw_matrix_panel(col_scores[ci,,drop=FALSE],ci,seq_len(r),c(0,cap[2]),colors,'Column contributions',column_label,'Component')
  draw_matrix_panel(row_scores[ri,,drop=FALSE],ri,seq_len(r),c(0,cap[1]),colors,'Row contributions',row_label,'Component')
  graphics::mtext(main,outer=TRUE,side=3,line=0.4,font=2,cex=1.1)
  caption<-paste('Squared singular-vector entries; separate color scales.',
                 if(sum(counts)>0)paste(sum(counts),'displayed entries clipped.')else '')
  graphics::mtext(caption,outer=TRUE,side=1,line=0.4,cex=0.8)
  invisible(list(row_scores=row_scores,column_scores=col_scores,row_index=ri,column_index=ci,
                 limits=cap,clipped_count=counts))
}
