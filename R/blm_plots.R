#### PLOTTING FUNCTIONS ####


#' Plot of a Bayesian Linear Model Fit
#' 
#' Plots the blm object along with the MAP fit line and quantiles for the fit.
#' 
#' @param x a blm object.
#' @param explanatory name of the explanatory variable to be plotted on the x 
#'   axis.
#' @param show_blm_interval display interval of the blm fitted values?
#' @param blm_interval_level level for the interval of the distribution of
#'   response to be shown.
#' @param expand_fit expand fit lines to the limits of the plot?
#' @param blm_col color for the blm fit lines.
#' @param blm_map_lty linetype for the blm MAP estimate line.
#' @param blm_interval_lty linetype for the blm quantiles line.
#' @param show_lm display the regular lm fit on the plot?
#' @param lm_col color for the lm fit line.
#' @param lm_lty linetype for the lm fit line.
#' @param caption title to be added to the plot.
#' @param add_formula append formula to the caption.
#' @param show_legend add a legend for the fit lines to the plot?
#' @param legend_parm list of additional arguments for the legend.
#' @param xlab label for the x axis.
#' @param ylab label for the y axis.
#' @param ... other arguments passed to plot.
#' 
#' @details Plots the data points used to create the \code{\link{blm}} object, 
#'   together with the blm fit, 95% quanitle of the blm fit and the 
#'   (frequentist) lm fit. 
#'   \itemize{ 
#'     \item \code{explanatory} can be used to specify what parameter is to be 
#'       plotted on the x axis. This can be useful ie when plotting models
#'       without a 'pure' explanatory ie y~sin(x); where untransformed 'x'
#'       itself is not part of the model. Providing \code{explanatory = 'x'}
#'       will produce the correct output for such models.
#'     \item \code{expand_fit} shows fit lines and extend to plot limits, 
#'       otherwise the fit is only shown for the range covered by the data. 
#'     \item \code{fit_legend_param} a list (of lists) of parameters to
#'      customize the legend. The first item in each slot targets the lm fit
#'      line (if present), the second paramter the blm fit line and the third
#'      paramter the interval lines for the blm fit.
#'   }
#'   
#' @examples
#' x <- rnorm(100)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3 ; w2 <- 10
#' 
#' y <- rnorm(100, mean = w0 + w1 * x + w2 *sin(x), sd = sqrt(1/b))
#' model <- blm(y ~ x + sin(x), prior = NULL, beta = b, 
#'                data = data.frame(x=x, y=y)) 
#' 
#' ##need to do more in case of complex model that does not contain a term 'x'
#' y <- rnorm(100, mean = w0 + w1 * cos(x), sd = sqrt(1/b))
#' model <- blm(y ~ cos(x), prior = NULL, beta = b, 
#'              data = data.frame(x=x, y=y))
#' 
#' \dontrun{
#' #plots y~cos(x) on the x-axis and this does probably not produce the plot we
#' want.
#' plot(model, xlim=c(-10,10))
#' }
#'
#' #one can specify that we want 'x' and not 'cos(x)' on the x axis by providing
#' #the name of the explanatory and we get the correct results 
#' plot(model, explanatory='x', xlim=c(-10,10))
#' 
#' 
#' #response and explanatory names do not matter
#' z_not_x <- rnorm(100)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3 ; w2 <- 10
#' 
#' means = w0 + w1 * z_not_x + w2 *sin(z_not_x)
#' v_not_y <- rnorm(100, means, sd = sqrt(1/b))
#' formula <- v_not_y ~ z_not_x + sin(z_not_x)
#' model <- blm(formula, prior = NULL, beta = b, 
#'                data = data.frame(z_not_x=z_not_x, v_not_y=v_not_y))
#' 
#' plot(model)
#' 
#' #customizing the plot, almost anything is possible
#' plot(model, col='purple', pch=19, cex=.3, type='p')
#' plot(model, col='purple', pch=19, cex=.3, xlim=c(-10,10), ylim=c(-50,50))
#' plot(model, col='purple', pch=19, cex=.3, xlim=c(-10,10), ylim=c(-50,50), 
#'      expand_fit=FALSE)
#' 
#' #customizing the plot lines
#'  plot(model, show_lm=FALSE)
#'  plot(model, show_lm=FALSE, blm_interval_lty=3, blm_col='red')
#'  plot(model, lm_col='blue', lm_lty=4, blm_interval_lty=3, blm_col='red')
#'  
#'  plot(model, show_legend=FALSE)
#'  plot(model, legend_parm=list(x='bottomright', cex=.5))
#'  
#' 
#' @export
plot.blm <- function(x, explanatory = NULL,
                     show_blm_interval = TRUE, 
                     blm_interval_level = .95,
                     expand_fit = TRUE,
                     blm_col = 'blue',
                     blm_map_lty = 1, 
                     blm_interval_lty = 2,
                     show_lm = TRUE, 
                     lm_col = 'red', 
                     lm_lty = 2,
                     caption = 'Bayesian lm fit',
                     add_formula = TRUE,
                     show_legend = TRUE, 
                     legend_parm = NULL,
                     xlab, ylab, ...){
  
  #want to keep argument name x to be consistent with generic fun
  object <- x
  
  if ( missing(explanatory)) {
    #take the first term of the formula per default
    frame <- object[['frame']]
    explanatory <- labels(terms(frame))[1]
  } else {
    frame <- expand.model.frame(object, explanatory, na.expand = FALSE)
  }
  
  x <- frame[,explanatory]
  response_lab <- colnames(frame)[attr(terms(frame), 'response')]
  response <- model.response(frame)
  
  if (missing('xlab')) xlab <- explanatory
  if (missing('ylab')) ylab <- response_lab
  caption <- ifelse(add_formula, 
                    paste(caption, deparse(object$formula), sep='\n'),
                    NULL)
  
  plot(x, response, 
       main=caption,
       xlab=xlab, ylab=ylab, ...)
  
  if ( expand_fit ) {
    plot_lims <- par('usr')[1:2]
    x_fit <- seq(plot_lims[1], plot_lims[2], length.out=1000)
    data <- data.frame(x=x_fit)
    colnames(data)[1] <- explanatory
  } else {
    x_order <- order(frame[,explanatory])
    x_fit <- frame[x_order,explanatory]
    data <- frame[x_order,]
  }
  
  if ( show_lm ) {
    #add lm fit
    y_lm <- predict(lm(object$formula, object$frame), data)
    lines(y_lm~x_fit, col = lm_col, lty= lm_lty)
  } 
  
  #add blm fit
  y_blm <- predict(object, data, se.fit=TRUE)
  
  #MAP fit lines
  y_blm_map <- y_blm$fit
  lines(y_blm_map~x_fit, col = blm_col, lty= blm_map_lty)
  
  #add quantiles of the fit
  if (show_blm_interval) {
    upper <- qnorm(blm_interval_level, mean = y_blm$fit, sd = y_blm$se.fit)
    lower <- qnorm(1-blm_interval_level, mean = y_blm$fit, sd = y_blm$se.fit)
    lines(lower~x_fit, col = blm_col, lty = blm_interval_lty)
    lines(upper~x_fit, col = blm_col, lty = blm_interval_lty)
  }
  
  ##add the legend
  if (show_legend) {
    if (missing(legend_parm)) {
      legend_parm <- list()
    }
    if (is.null(legend_parm[['x']])) {
      legend_parm[['x']] <- 'topleft'
    }
    if (is.null(legend_parm[['legend']])) {
      if (show_blm_interval) {
        legend_parm[['legend']] <- c('blm MAP', 
                                     paste(c('blm', 
                                             100*blm_interval_level, 
                                             '% quantile'), collapse=' '))
        legend_parm[['lty']] <- c(blm_map_lty, blm_interval_lty)
        legend_parm[['col']] <- rep(blm_col,2)
      } else {
        legend_parm[['legend']] <- 'blm MAP'
        legend_parm[['lty']] <- blm_map_lty
        legend_parm[['col']] <- blm_col
      }
    }
    if ( show_lm ) {
      legend_parm[['legend']] <- c(legend_parm[['legend']], 'lm fit' )
      legend_parm[['lty']] <- c(legend_parm[['lty']], lm_lty)
      legend_parm[['col']] <- c(legend_parm[['col']], lm_col)
    }
    do.call(legend, legend_parm)
  }
  
}



#' Sample Lines
#' 
#' Adds random sampled lines from a \code{blm} object to an existing plot.
#' 
#' @param object a \code{blm} object
#' @param n number of lines to show
#' @param sample_from_prior draw samples from the prior distribution?
#' @param ... additional arguments passed to lines
#' 
#' @details Draws random samples from the posterior (or prior if
#'   sample_from_prior = TRUE) distribution of \code{\link{blm}} object.
#'   Requires the mvrnorm from the MASS package for random sampling. NOTE: only
#'   works for straight lines, assuming first and second coefficient are
#'   intercept and slope.
#' 
#' @examples
#' x <- rnorm(100)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' 
#' y <- rnorm(100, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, 
#'                data = data.frame(x=x, y=y)) 
#' 
#' plot(model)
#' \dontrun{
#' add_sample_lines(model, col='green')                
#' }
#'
#' @export 
add_sample_lines <- function(object, n = 5, sample_from_prior = FALSE, ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop('add_sample_lines requires package MASS',
         call. = FALSE)
  }
  
  if (sample_from_prior) {
    d <- object$prior
  } else {
    d <- object$posterior
  }
  s <- MASS::mvrnorm(n, d$means, d$covar, empirical=TRUE)
  
  apply(s, 1, abline, ...)
}



#' Density Plot of Multivariate Normal Distribution
#' 
#' Plots a kernel density for a \code{mvnd} object.
#' 
#' @param object a \code{mvnd} object
#' @param n sampling depth used for kernel density estimation.
#' @param ... additional arguments passed to the plot function.
#' 
#' @details Draws random samples from a multivariate normal distribution of
#'   \code{\link{mvnd}} object and plots the kernel density of those.
#'   Requires the mvrnorm from the MASS package for random sampling. If present,
#'   the package LSD is used for kernel density coloring. ... argments are
#'   passed to generic plot or \code{\link[LSD]{heatscatter}} functions if LSD
#'   package is installed.
#'   
#' @examples
#' x <- rnorm(100)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' 
#' y <- rnorm(100, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, 
#'                data = data.frame(x=x, y=y)) 
#' 
#' \dontrun{
#' par(mfrow=c(1,2))
#' kernel_density(model$prior, xlim=c(-4,4), ylim=c(-4,4), main='prior')    
#' kernel_density(model$posterior, xlim=c(-4,4), ylim=c(-4,4), main='posterior')
#' par(mfrow=c(1,1)) 
#' } 
#'                
#' @export 
kernel_density <- function(object, n = 10000, ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop('Kernel density requires package MASS',
         call. = FALSE)
  }
  s <- MASS::mvrnorm(n, object$means, object$covar, empirical=TRUE)
  if (requireNamespace("LSD", quietly = TRUE)) {
    LSD::heatscatter(s[,2], s[,1], 
                     xlab=names(object$means)[2], 
                     ylab=names(object$means)[1], ...)
  } else {
    plot(s[,2], s[,1], 
         xlab=names(object$means)[2], 
         ylab=names(object$means)[1], ...)
  }
}



#' Plot of Residuals against Fitted Values
#' 
#' Plots residuals against fitted values for a \code{\link{blm}} object.
#' 
#' @param object a \code{\link{blm}} object.
#' @param show_var display errorbars showing standard deviation on the 
#'   residuals.
#' @param id.n number of points to be labelled in each plot, starting with the
#'   most extreme.
#' @param labels.id  vector of labels, from which the labels for extreme points
#'   will be chosen. NULL uses observation numbers.
#' @param cex.id  magnification of point labels.
#' @param ... other arguments passed to plot.
#'   
#' @details Standard deviation of the residuals is same as for the fitted
#'   values.
resid_vs_fitted_plot <- function(object, 
                                 show_var = FALSE, 
                                 id.n = 3, 
                                 labels.id = names(residuals(object)), 
                                 cex.id = 0.75,
                                 ...) {
  fits <- fitted(object, var = T)
  y <- resid(object)
  x <- fits$mean
  plot(y~x,  
       main = 'Residuals vs Fitted',
       xlab = paste('Fitted Values\n',
                    deparse(object$formula),sep=''),
       ylab = 'Residuals',
       ...)
  
  if (id.n > 0) {
    ids <- order(abs(y), decreasing = TRUE)[1:id.n]
    text(y[ids]~x[ids], label=labels.id[ids], cex=cex.id, adj=-.5)
  }
  
  if (show_var) {
    sd <- sqrt(fits$var)
    segments(x, y-sd, x, y+sd)
    x_range <- par('usr')[2] - par('usr')[1]
    epsilon = x_range/200
    segments(x-epsilon,y-sd,x+epsilon,y-sd)
    segments(x-epsilon,y+sd,x+epsilon,y+sd)
  }
}



#' Scale-Location Plot
#' 
#' Scale-Location Plot for a \code{\link{blm}} object.
#' 
#' @param object a \code{\link{blm}} object.
#' @param id.n number of points to be labelled in each plot, starting with the
#'   most extreme.
#' @param labels.id  vector of labels, from which the labels for extreme points
#'   will be chosen. NULL uses observation numbers.
#' @param cex.id  magnification of point labels.
#' @param ... other arguments passed to plot.
#'   
#' @details Scale-location plots the sqrt(| residuals |) against fitted values.
#'   This is supposed to diminish skewness. Note: In contrast to the
#'   \code{\link{plot.lm}} implementation, the residuals here are taken for
#'   sqrt(| residuals |) and not the standardized residuals.
scale_location_plot <- function(object,                                  
                                id.n = 3, 
                                labels.id = names(residuals(object)), 
                                cex.id = 0.75,
                                ...) {
  x <- fitted(object, var = FALSE)
  y <- sqrt(abs(resid(object)))
  plot(y~x, 
       main = 'Scale-Location',
       xlab = paste('Fitted Values\n',
                    deparse(object$formula),sep=''),
       ylab = expression(sqrt(abs(Residuals))),
       ...)
  
  if (id.n > 0) {
    ids <- order(y, decreasing = TRUE)[1:id.n]
    text(y[ids]~x[ids], label=labels.id[ids], cex=cex.id, adj=-.5)
  }
  
}


#' Q-Q Plot
#' 
#' Q-Q Plot for a \code{\link{blm}} object.
#' 
#' @param object a \code{\link{blm}} object.
#' @param qqline add qqline to the plot.
#' @param ... other arguments passed to qqplot.
#'   
#' @details Q-Q plots the sqrt(| residuals |) against fitted values. This
#'   transformation of the residuals is supposed to diminish skewness (see
#'   \code{\link{plot.lm}}. Note: In contrast to the \code{\link{plot.lm}}
#'   implementation, the residuals here are taken for sqrt(| residuals |) and
#'   not the standardized residuals.
qq_blm_plot <- function(object, qqline = TRUE, ...) {
  x <- fitted(object, var = FALSE)
  y <- resid(object)
  
  qqnorm(x, 
         main = 'Normal Q-Q',
         xlab = paste('Theoretical Quantiles\n',
                      deparse(object$formula),sep=''),
         ylab = 'Residuals',
         ...)
  if (qqline) qqline(x)
  
}



#' Diagnostic Plots
#' 
#' Diagnostic plots for a \code{\link{blm}} object.
#' 
#' @param object a \code{\link{blm}} object.
#' @param which subset of the plots to display.
#' @param show_var display errorbars showing standard deviation on the 
#'   residuals vs fitted plot.
#' @param id.n number of points to be labelled in each plot, starting with the
#'   most extreme. Applies to plots 1 and 2.
#' @param labels.id  vector of labels, from which the labels for extreme points
#'   will be chosen. NULL uses observation numbers. Applies to plots 1 and 2.
#' @param cex.id  magnification of point labels. Applies to plots 1 and 2.
#' @param qqline add qqline to the Q-Q plot.
#' @param ... other arguments passed to the plot functions.
#'   
#' @details The following plots are available: 
#'  \enumerate{ 
#'    \item \code{\link{resid_vs_fitted_plot}} 
#'    \item \code{\link{scale_location_plot}} 
#'    \item \code{\link{qq_blm_plot}}
#'  }
#'  
#' @export
diagnostic_plots <- function(object, 
                             which = c(1:3),
                             show_var = FALSE, 
                             id.n = 3, 
                             labels.id = names(residuals(object)), 
                             cex.id = 0.75,
                             qqline = TRUE, ...) {
  
  if ( !('blm' %in% class(object)) ) {
    stop('function needs an class blm object')
  }
  
  if ( 1 %in% which) {
    resid_vs_fitted_plot(object, show_var=show_var, 
                         id.n=id.n, labels.id=labels.id, 
                         cex.id=cex.id, ...)
  }
  
  if ( 2 %in% which) {
    scale_location_plot(object, id.n=id.n, labels.id=labels.id, 
                        cex.id=cex.id, ...)
  }
  
  if ( 3 %in% which) {
    qq_blm_plot(object, qqline=qqline, ...)
  }
  
}