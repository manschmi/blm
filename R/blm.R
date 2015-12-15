#' blm: A package for bayesian linear regression.
#'
#' The \code{blm} package is for construction and analysis of bayesian linear
#' regression models.
#' 
#' 
#' @section blm functions:
#' blm
#' update.blm
#' print.blm
#' summary.blm
#' fitted.blm
#' predict.blm
#' summary.blm
#' resid.blm
#' coef.blm
#' deviance.blm
#' plot.blm
#' 
#' 
#' @docType package
#' @name blm
#' 
NULL



#' Prior distribution
#' 
#' Creates a prior distribution for a blm object
#' 
#' @param mean mean or a vector of the mean(s) of the distribution.
#' @param alpha single value or vector for precision(s) of the distribution.
#' @param dim number of variables for the distribution.
#'   
#' @details Creates a multivariate normal distribution object of class 
#'   \code{\link{mv_dist}} that holds \code{means} and covariance \code{covar} 
#'   of a multivariate normal distribution. If number of means or covars is less
#'   than \code{dim}, the values provided are expanded to fill the multivariate 
#'   distribution.
#'   
#' @return An object of class \code{mv_dist} with provided means and covariances
#'   computed as 1/\code{alpha}.
#'   
#' @examples
#'  make_prior(0,1,2)
#'  
#'  make_prior(c(0,1),c(1,2))
#'  
#' @export
make_prior <- function(mean, alpha, dim=length(mean)){
  covar <- diag(1/alpha, dim, dim)
  if (length(mean) != dim) {
    mean <- rep(mean,dim)[1:dim]
  } 
  mv_dist(mean, covar)
}


#' Bayesian Linear Regression
#' 
#' Computes a bayesian linear regression fit. The current implementation 
#' requires knowledge about the precision \code{beta} of the observations.
#' 
#' @param formula an object of class "formula": a symbolic description of the 
#'   model to be fitted.
#' @param prior the prior distribution, either a \code{mv_dist} or a \code{blm} 
#'   object.
#' @param beta the precision of the data used for the model.
#' @param ... Additional arguments and values.
#'   
#' @details Models for \code{blm} are provided as for \code{\link{lm}}. If prior
#'   distribution \code{prior} of the wieghts is not provided the prior means 
#'   are set to 0 and the variances to 1.
#'   
#' @return A object of class \code{blm}. An object of class "lm" is a list 
#'   containing at least the following components: 
#'   \describe{ 
#'    \item{call}{the matched call} 
#'    \item{formula}{the formula used} 
#'    \item{frame}{the model frame used} 
#'    \item{matrix}{the model matrix used} 
#'    \item{beta}{the precision of the data} 
#'    \item{prior}{the prior distribution used} 
#'    \item{posterior}{the posterior distribution}
#'   }
#'   
#' @examples
#'  
#'    w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
#'    x <- rnorm(50)
#'    y <- rnorm(50, w1 * x + w0, 1/b)
#'    mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
#'    
#'    summary(mod)
#'    plot(mod)
#'    
#'    #use of a prior, typically from an existing model
#'    x2 <- rnorm(50)
#'    y2 <- rnorm(50, w1 * x2 + w0, 1/b)
#'    mod <- blm(y~x, prior=mod, beta=b, data=data.frame(x=x2, y=y2))
#'    mod
#'    
#'    #use with 2 explanatory variables
#'    w2 <- 3.3
#'    z <- rnorm(50)
#'    y <- rnorm(50, w2 * z + w1 * x + w0, 1/b)
#'    mod <- blm(y~x+z, beta=b, data=data.frame(x=x, y=y, z=z))
#'    mod
#'  
#' @export
blm <- function(formula, prior = NULL, beta = 1, ...) {
  
  frame <- model.frame(as.formula(formula), ...)
  matrix <- model.matrix(as.formula(formula), frame)
  mod_vars <- colnames(matrix)
  
  if (is.null(prior)) {
    n <- ncol(matrix)
    prior <- make_prior(mean = 0, alpha = 1, dim = n)
    prior <- mv_dist.set_var_names(prior, mod_vars)
  }
  
  prior <- distribution(prior)
  
  prior_vars <- mv_dist.var_names(prior)
  if ( !all(prior_vars %in% mod_vars) ) {
    #can occur ie when updating to a new formula with fewer terms
    if ( all(mod_vars %in% prior_vars) ){
      message('prior contains more variables than the model : variables not used in the model are ignored in the fit')
      prior <- mv_dist.subset(prior, colnames(matrix))
    } else {
      missing <- mod_vars[!(mod_vars %in% prior_vars)]
      stop(paste('model formula contains variable names not provided in the prior', missing ))
    }
  }
  
  prior_covar <- mv_dist.covar(prior)
  
  response <- model.response(frame)
  
  lhd_term <- beta * t(matrix) %*% matrix
  
  covar <- solve(prior_covar + lhd_term)
  
  means <- beta * covar %*% t(matrix) %*%  response
  
  posterior <- distribution(mv_dist(means, covar))
  
  mod <- structure(list(call = sys.call(),
                        formula = formula,
                        frame = frame,
                        matrix = matrix,
                        beta = beta,
                        prior = prior,
                        posterior = posterior),
                   class = 'blm')
  
}


#' Update a blm distribution
#' 
#' Updates a \code{blm} object, using new data, parameters 
#'  or a new model formula.
#' 
#' @param object a \code{blm} object.
#' @param formula a formula for the updated object.
#' @param beta precision estimate for the updated object.
#' @param prior prior distribution of the updated object.
#' @param data data of the updated object.
#' @param ... other arguments. These are passed to blm.
#' 
#' @details Updates a \code{\link{blm}} object, using features of the input object except, if provided otherwise. The prior for the updated object is typically the posterior distribution of the input object, but can be specified specifically (as for all other parameters. Importantly, new data will only be used when specified as named argument.
#' 
#' @return A object of class \code{blm}.
#' 
#' @examples
#'    w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
#'    x <- rnorm(50)
#'    y <- rnorm(50, w1 * x + w0, 1/b)
#'    mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
#'    mod
#'    
#'    #use of a prior, typically from an existing model
#'    x2 <- rnorm(50)
#'    y2 <- rnorm(50, w1 * x2 + w0, 1/b)
#'    
#'    #using posterior of mod as prior
#'    new_mod <- update(mod, data=data.frame(x=x2, y=y2)) 
#'    new_mod
#'    
#'    #using same prior for mod and new model
#'    new_mod2 <- update(mod, prior=mod$prior, data=data.frame(x=x2, y=y2)) 
#'    new_mod2
#' 
#' @export
update.blm <- function(object, 
                       formula = object$formula, 
                       beta = object$beta,
                       prior = object, 
                       data = object$frame, 
                       ...) {
  
  blm(formula=formula, prior=prior, beta=beta, data=data, ...)
  
}


#' Coefficients
#' 
#' Coefficients from the posterior distribution of a bayesian model.
#'
#' @param object a \code{blm} object.
#' @param ... other arguments (currently ignored).
#' 
#' @return A named vector for the maxium a posteriori likelihood of coefficients. 
#'  That is, means of the posterior distribution of the coefficients of the
#'   object. To get the covariance matrix for the coefficients use 
#'   \code{\link{covar.blm}}.
#'  
#' @examples
#' x <- rnorm(10)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' coef(model)
#' coefficients(model)
#' 
#' @export 
coef.blm <- function(object, ...){
  mv_dist.means(object$posterior)
}




#' Covariance of coefficients
#' 
#' Covariance of coefficients from the posterior distribution of a bayesian 
#'  model.
#'  
#' @param object a \code{blm} object.
#' @param ... other arguments (currently ignored).
#' 
#' @return Covariance matrix for the distribution of the coefficients of a 
#'  \code{\link{blm}} model. 
#' 
#' @examples
#' x <- rnorm(10)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' covar.blm(model)
#' 
#' @export 
covar.blm <- function(object, ...){
  mv_dist.covar(object$posterior)
}



#' Confidence Interval for Coefficients of a Bayesian Linear Model
#' 
#' Compute the confidence interval \code{\link{confint}} for a bayesian linear
#'  model.
#'  
#' @param object a \code{blm} object.
#' @param parm parameter for which confidence interval should be computed. Either a number or vector of numbers. If missig all parameters are considered.
#' @param level confidence level
#' @param ... other arguments (currently ignored).
#'
#' @details Deviance is the sum of squares of the residuals, extracted from the 
#'  maxium a posterior estimate of the \code{blm} object.
#'  
#' @return Pseudo confidence interval of parameter(s) of the \code{blm} object.
#' Returns the quantile specified by level for the distribution of \code{param}.
#'  
#' @examples
#' x <- rnorm(10)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' confint(model)
#' confint(model, 1)
#' confint(model, c(1,2))
#' 
#' @export 
confint.blm <- function (object, 
                         parm, 
                         level = .95, ...) {
  
  if (missing(parm)) {
    param <- seq_along(coef(object))
  }
  
  cf <- coef(object)[parm]
  cf_var <- diag(covar.blm(object))[parm]
  
  l <- (1-level)/2
  lower <- qnorm(l, cf, sqrt(cf_var))
  upper <- qnorm((1-l), cf, sqrt(cf_var))
  
  matrix(c(lower, upper),
         nrow=length(lower),
         dimnames=list(names(cf),
                       paste(100*c(l,1-l), '%')))
}




#' Predict 
#' 
#' Predict means and variance of the response variable for a \code{\link{blm}}
#'  object.
#'  
#' @param object a \code{blm} object.
#' @param newdata an optional data frame containing variables with which to 
#'  predict. If missing, values used for fitting will be extracted from the 
#'  object.
#' @param report.var report also variance for the predicted values.
#' @param ... other arguments (currently ignored).
#' 
#' @return A vector of predicted values. If report.var = TRUE, a list containing
#'  means and variances of the predicted values.
#'  
#' @examples
#' x <- rnorm(10)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' predict(model)
#' 
#' x <- rnorm(10) 
#' predict(model, data.frame(x=x))
#' 
#' 
#' @export
predict.blm <- function(object, newdata = NULL, report.var = FALSE, ...){
  
  if ( missing(newdata) ) {
    return( fitted(object, report.var=report.var) )   
  }
  
  responseless_formula <- delete.response(terms(object$formula))
  frame <- model.frame(responseless_formula, newdata)
  matrix <- model.matrix(responseless_formula, frame)
  
  means <- apply(matrix, 1, function(xi) sum(coef(object) * xi))
  
  if (report.var) {
    vars <- apply(matrix, 1, function(xi) 
      1/object$beta +  t(xi) %*% covar.blm(object) %*% xi)
    
    return(list(mean = means, var = vars))
  }
  
  means
}



#' Extract Model Fitted Values 
#' 
#' Extracts the fitted values for the response variable of a \code{\link{blm}}
#'  object.
#'  
#' @param object a \code{blm} object.
#' @param report.var Report also variance for the fitted values.
#' @param ... other arguments (currently ignored).
#'
#' @return A vector of fitted values. If report.var = TRUE, a list containing
#'  means and variances of the fitted values. Identical to \code{predict.blm}
#'  without newdata.
#'  
#' @examples
#' x <- rnorm(10)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' fitted(model)
#' 
#' @export
fitted.blm <- function(object, report.var = FALSE, ...){
  
  responseless_formula <- delete.response(terms(object$formula))
  matrix <- model.matrix(responseless_formula, object$frame)    
  
  means <- apply(matrix, 1, function(xi) sum(coef(object) * xi))
  
  if (report.var) {
    vars <- apply(matrix, 1, function(xi) 
      1/object$beta +  t(xi) %*% covar.blm(object) %*% xi)
    
    return(list(mean = means, var = vars))
  }
  
  means
}



#' Residuals of Bayesian Linear Model
#' 
#' \code{residuals.blm} returns \code{\link{residuals}} for a bayesian linear
#'  model.
#'  
#' @param object a \code{blm} object.
#' @param ... other arguments (currently ignored).
#' 
#' @return Residuals of the \code{blm} object.
#'  
#' @details Residuals are extracted from the maxium a posterior estimate of the 
#'  \code{blm} object.
#'  
#' @examples
#' x <- rnorm(10)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' residuals(model)
#' resid(model)
#' 
#' @export 
residuals.blm <- function(object, ...){
  frame <- object[['frame']]
  response_col <- attr(terms(frame), 'response')
  response <- frame[, response_col]
  
  prediction <- fitted.blm(object, report.var = F)
  
  response - prediction
}



#' Deviance of Bayesian Linear Model
#' 
#' \code{deviance.blm} returns \code{\link{deviance}} for a bayesian linear
#'  model.
#'  
#' @param object a \code{blm} object.
#' @param ... other arguments (currently ignored).
#'  
#' @details Deviance is the sum of squares of the residuals, extracted from the 
#'  maxium a posterior estimate of the \code{blm} object.
#'  
#' @return Deviance of the \code{blm} object.
#' 
#' @examples
#' x <- rnorm(10)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3
#' y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#' model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' deviance(model)
#' 
#' 
#' @export 
deviance.blm <- function(object, ...){
  sum( resid(object)^2 )
}



#' Print 
#' 
#' Print a blm object
#' 
#' @param x a \code{\link{blm}} object.
#' @param ... other arguments (currently ignored).
#' 
#' @export
print.blm <- function(x, ...){
  cat('Call:\n')
  print(x$call)
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
}



#' Summary 
#' 
#' Summary of a blm object
#' 
#' @param object a \code{\link{blm}} object.
#' @param ... other arguments (currently ignored).
#' 
#' @return An object of class \code{summary.blm} containing:
#'  \describe{
#'   \item{call}{the call of the object }
#'   \item{prior}{the prior distribution of the fit}
#'   \item{posterior}{the posterior distribution of the fit}
#'  }
#'   
#' @export
summary.blm <- function(object, ...){
  
  structure(list(
    call = object$call,
    prior = object$prior,
    posterior = object$posterior
  ),
  class = "summary.blm")
  
}



#' Plot of a Bayesian Linear Model Fit
#' 
#' Plots the blm object along with the MAP fit line and quantiles for the fit.
#' 
#' @param x a blm object.
#' @param explanatory name of the explanatory variable to be plotted on the x 
#'   axis.
#' @param se display variance of blm fit?
#' @param se_level level for the interval of the distribution of response to be
#'   shown.
#' @param show_lm display the regular lm fit on the plot?
#' @param expand_fit expand fit lines to the limits of the plot?
#' @param show_fit_legend add a legend for the fit lines to the plot?
#' @param fit_legend_param parameters for the fit legend?
#' @param xlab label for the x axis.
#' @param ylab label for the y axis.
#' @param ... other arguments passed to plot.
#' 
#' @details Plots the data points used to create the \code{\link{blm}} object,
#'   together with the blm fit, 95% quanitle of the blm fit and the
#'   (frequentist) lm fit. \code{explanatory} can be used to specify what
#'   parameter is to be plotted on the x axis. This can be useful ie when
#'   plotting models with a 'pure' explanatory ie y~sin(x); where 'x' is not
#'   part of the model. Providing \code{explanatory = 'x'} will produce the
#'   correct output for such models. \code{data} can be used to display
#'   additional data points on the plot. The data used to create the model will
#'   always be shown.
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
#' plot(model, xlim=c(-10,10))
#' 
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
#' #we can specify that we want 'x' and not 'cos(x)' by providing the name of
#' # the explanatory and we get the correct results
#' plot(model, explanatory='x', xlim=c(-10,10))
#' 
#' 
#' #response and explanatory names do not matter
#' z_not_x <- rnorm(100)
#' b <- 1.3
#' w0 <- 0.2 ; w1 <- 3 ; w2 <- 10
#' 
#' v_not_y <- rnorm(100, mean = w0 + w1 * z_not_x + w2 *sin(z_not_x), sd = sqrt(1/b))
#' model <- blm(v_not_y ~ z_not_x + sin(z_not_x), prior = NULL, beta = b, 
#'                data = data.frame(z_not_x=z_not_x, v_not_y=v_not_y))
#' 
#' plot(model)
#' 
#' 
#' @export
plot.blm <- function(x, explanatory = NULL, 
                     se=TRUE, se_level = .95, 
                     show_lm=TRUE, expand_fit=TRUE, 
                     show_fit_legend=TRUE, fit_legend_param=NULL,
                     xlab, ylab, ...){
  
  object <- x #want to keep argument name x to be consistent with generic fun
  
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
  
  plot(x, response, type = "p", pch=19, xlab=xlab, ylab=ylab, ...)
  
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
    lines(y_lm~x_fit, col='red', lty=2)
  } 
  
  #add blm fit
  y_blm <- predict(object, data, report.var=T)
  y_blm_map <- y_blm$mean
  lines(y_blm_map~x_fit, col='blue', lty=1)
  
  #add quantiles of the fit
  if (se) {
    upper <- qnorm(se_level, mean = y_blm$mean, sd = sqrt(y_blm$var))
    lower <- qnorm(1-se_level, mean = y_blm$mean, sd = sqrt(y_blm$var))
    
    lines(lower~x_fit, col='blue', lty=2)
    lines(upper~x_fit, col='blue', lty=2) 
  }

  ##add the legend
  if (show_fit_legend) {
    if (missing(fit_legend_param)) {
      fit_legend_param <- list()
    }
    if (is.null(fit_legend_param[['x']])) {
      fit_legend_param[['x']] <- 'topleft'
    }
    if (is.null(fit_legend_param[['legend']])) {
      if ( show_lm ) {
        fit_legend_param[['legend']] <- c('lm fit', 'blm MAP')
      } else {
        fit_legend_param[['legend']] <- 'blm MAP'
      }
      if (se) {
        fit_legend_param[['legend']] <- c(fit_legend_param[['legend']],
                                          paste(c('blm', 
                                                  100*se_level, 
                                                  '% quantile'), collapse=' '))
      }
    }
    if (is.null(fit_legend_param[['lty']])) {
      if ( show_lm ) {
        fit_legend_param[['lty']] <- c(2, 1)
      } else {
        fit_legend_param[['lty']] <- 1
      }
      if (se) {
        fit_legend_param[['lty']] <- c(fit_legend_param[['lty']], 2)
      }
    }
    if (is.null(fit_legend_param[['col']])) {
      fit_legend_param[['col']] <- c()
      if ( show_lm ) {
        fit_legend_param[['col']] <- c('red', 'blue')
      } else {
        fit_legend_param[['col']] <- 'blue'
      }
      if (se) {
        fit_legend_param[['col']] <- c(fit_legend_param[['col']], 'blue')
      }   
    }
    do.call(legend, fit_legend_param)
  }
  
}