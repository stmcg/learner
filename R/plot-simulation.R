#' Compare methods across simulation settings
#'
#' Plot replicate means by noise-variance ratio, with panels for each scenario
#' and rank. The default method colors and line types follow the simulation
#' figures in the LEARNER paper. This function summarizes supplied results;
#' it does not run simulations or reproduce the paper's experiments.
#'
#' @param data Data frame with columns \code{scenario}, \code{rank},
#' \code{variance_ratio}, \code{method}, and \code{error}. Each row must be one
#' replicate result, not an already aggregated mean. Scenario and method are
#' nonempty character strings or factors; rank is a positive integer; the
#' variance ratio is positive; error is finite and nonnegative. Use the same
#' error definition across all methods and scenarios. For the paper's metric,
#' supply the Frobenius norm of estimated minus true signal.
#' @param uncertainty One of \code{"none"}, \code{"sd"}, or \code{"se"}.
#' Bars show one sample standard deviation or one standard error of the mean,
#' not a confidence interval. Bars are omitted for groups with one replicate.
#' @param main Overall title.
#' @param ylab Vertical label naming the supplied error metric.
#' @param xlab Horizontal label; defaults to the target/source noise variance ratio.
#' @return Invisibly, a data frame with grouping columns and \code{mean},
#' \code{sd}, \code{se}, and \code{n}. Missing groups are not filled with zero.
#' @details
#' Scenarios and methods follow their first appearance in the data; ranks and
#' ratios are sorted numerically. At most two ranks and three scenarios are
#' drawn per page. All pages share vertical limits and method styles. Missing
#' method-by-ratio groups break curves. Canonical method names are
#' \code{"Target-only SVD"}, \code{"D-LEARNER"}, and \code{"LEARNER"}.
#' Other names are permitted and receive distinct styles. Replicate identifiers
#' are not inferred: callers must avoid duplicate records or mixing experiments.
#' @examples
#' results <- expand.grid(scenario = c("High similarity", "Low similarity"),
#'                        rank = 2, variance_ratio = c(1, 3, 5),
#'                        method = c("Target-only SVD", "D-LEARNER", "LEARNER"),
#'                        replicate = 1:2, stringsAsFactors = FALSE)
#' # Illustrative values only; replace with measured errors from simulations.
#' results$error <- 3 + 1 / results$variance_ratio + results$replicate / 10
#' plot_simulation_results(results, main = "Illustrative input format")
#' @export
plot_simulation_results <- function(data, uncertainty=c('none','sd','se'),
                                    main='Simulation comparison', ylab='Estimation error',
                                    xlab=expression(sigma[0]^2/sigma[1]^2)) {
  uncertainty<-match.arg(uncertainty)
  required<-c('scenario','rank','variance_ratio','method','error')
  if(!is.data.frame(data)||!all(required %in% names(data))||nrow(data)==0L)
    stop('data must be a nonempty data frame with scenario, rank, variance_ratio, method, and error.')
  data<-data[,required,drop=FALSE]
  for(name in c('scenario','method')) {
    if(!(is.character(data[[name]])||is.factor(data[[name]]))||anyNA(data[[name]])||
       any(!nzchar(trimws(as.character(data[[name]])))))stop(name,' must contain nonempty labels.')
    data[[name]]<-as.character(data[[name]])
  }
  for(name in c('rank','variance_ratio','error')) {
    z<-data[[name]]
    if(!is.numeric(z)||is.complex(z)||any(!is.finite(z)))stop(name,' must be finite and numeric.')
  }
  if(any(data$rank<1|data$rank!=floor(data$rank)))stop('rank must contain positive integers.')
  if(any(data$variance_ratio<=0))stop('variance_ratio must be positive.')
  if(any(data$error<0))stop('error must be nonnegative.')
  # Integer codes avoid collisions when labels contain punctuation.
  keys<-lapply(data[,required[1:4]],function(z)match(z,unique(z)))
  groups<-split(seq_len(nrow(data)),interaction(keys,drop=TRUE))
  summary<-do.call(rbind,lapply(groups,function(i) {
    first<-data[i[1],required[1:4],drop=FALSE]
    n<-length(i);s<-if(n>1)stats::sd(data$error[i])else NA_real_
    cbind(first,mean=mean(data$error[i]),sd=s,se=s/sqrt(n),n=n)
  }))
  rownames(summary)<-NULL
  scenarios<-unique(data$scenario);ranks<-sort(unique(data$rank));methods<-unique(data$method)
  summary<-summary[order(match(summary$scenario,scenarios),summary$rank,
                         match(summary$method,methods),summary$variance_ratio),]
  rownames(summary)<-NULL
  canonical<-c('Target-only SVD','D-LEARNER','LEARNER')
  colors<-stats::setNames(grDevices::hcl.colors(length(methods),'Dark 3'),methods)
  types<-stats::setNames(1+(seq_along(methods)-1)%%6,methods)
  for(j in seq_along(canonical))if(canonical[j] %in% methods) {
    colors[canonical[j]]<-c('#E41A1C','#4DAF4A','#377EB8')[j]
    types[canonical[j]]<-c(1,2,6)[j]
  }
  bar<-if(uncertainty=='none')rep(0,nrow(summary))else summary[[uncertainty]]
  bar[is.na(bar)]<-0
  ylim<-c(min(0,summary$mean-bar),max(summary$mean+bar))
  if(diff(ylim)==0)ylim<-c(0,1)
  ylim[2]<-ylim[2]+0.05*diff(ylim)
  old<-graphics::par(no.readonly=TRUE);on.exit(graphics::par(old),add=TRUE)
  for(sp in split(scenarios,ceiling(seq_along(scenarios)/3))) {
    for(rp in split(ranks,ceiling(seq_along(ranks)/2))) {
      graphics::par(mfrow=c(length(sp),length(rp)),mar=c(3.5,4,2.8,1),
                    oma=c(4.5,0,2.5,0),mgp=c(2.3,0.6,0),cex.axis=0.9,cex.lab=0.95)
      for(sc in sp)for(rk in rp) {
        panel<-summary[summary$scenario==sc & summary$rank==rk,,drop=FALSE]
        ratios<-sort(unique(data$variance_ratio))
        graphics::plot(range(ratios),ylim,type='n',ylim=ylim,xlab=xlab,ylab=ylab,
                       main=paste0(sc,' (Rank = ',rk,')'),cex.main=1.05)
        graphics::grid(col='gray85',lty=1)
        for(m in methods) {
          z<-panel[panel$method==m,,drop=FALSE]
          pos<-match(ratios,z$variance_ratio)
          graphics::lines(ratios,z$mean[pos],type='b',pch=16,col=colors[m],lty=types[m],lwd=2)
          if(uncertainty!='none') {
            valid<-is.finite(z[[uncertainty]]) & z[[uncertainty]]>0
            if(any(valid))graphics::arrows(z$variance_ratio[valid],z$mean[valid]-z[[uncertainty]][valid],
                     z$variance_ratio[valid],z$mean[valid]+z[[uncertainty]][valid],
                     angle=90,code=3,length=0.035,col=colors[m])
          }
        }
      }
      # A full-device overlay places one shared legend beneath the panel grid.
      graphics::par(fig=c(0,1,0,1),new=TRUE,mar=rep(0,4),oma=rep(0,4),cex=1)
      graphics::plot.new();graphics::plot.window(c(0,1),c(0,1),xaxs='i',yaxs='i')
      graphics::legend('bottom',inset=0.022,legend=methods,col=colors,lty=types,pch=16,
                       ncol=min(3,length(methods)),bty='n',cex=0.9,xpd=NA)
      label<-if(uncertainty=='none') 'Replicate means; no uncertainty bars.' else
        paste0('Replicate means; bars: +/- 1 ',toupper(uncertainty),'. No bars for n = 1.')
      graphics::text(0.5,0.012,label,cex=0.75,xpd=NA)
      graphics::text(0.5,0.98,main,font=2,cex=1.05,xpd=NA)
    }
  }
  invisible(summary)
}
