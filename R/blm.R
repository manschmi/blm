#' blm: A package for bayesian linear regression.
#'
#' The \code{blm} package is for construction and analysis of bayesian linear
#'  regression models.
#' 
#' 
#' @section blm functions:
#' blm
#' update.blm
#' print.blm
#' summary.blm
#' resid.blm
#' coef.blm
#' deviance.blm
#' plot.blm
#' 
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
#' @return An object of class \code{mv_dist} containing a "means": vector of 
#' means and "covar": the covariance matrix.
make_prior <- function(mu, alpha, dim){
  covar <- diag(1/alpha, dim, dim)
  if (length(mu) != dim) {
    mu_vec <- rep(mu,dim)[1:dim]
  }
  mv_dist(mu_vec, covar)
}


#' Bayesian Linear Regression
#' 
#' Computes a bayesian linear regression fit. It only works for
#'  fixed variance beta.
#'  
#'  @param formula Object of class "formula": a symbolic description of the 
#'    model to be fitted.
#'  @param prior The prior distribution, either a mv_dist or a blm object.
#'  @param ... Additional arguments, or arguments with changed values
#'  
#'  @export
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
  if ( all(mod_vars %in% prior_vars) ) {
    #can occur ie when updating to a new formula with fewer terms
    if ( mod_vars %in% prior_vars ){
      warning('prior contains more variables than the model : variables not used in the model are ignored in the fit')
      prior <- subset(prior, colnames(matrix))
    } else {
      missing <- mod_vars[!(mod_vars %in% prior_vars)]
      stop(paste('model formula contains variable names not provided in the prior', missing ))
    }
  }
  
  print(prior)

  prior_covar <- mv_dist.covar(prior)
  
  response <- model.response(frame)
  
  lhd_term <- beta * t(matrix) %*% matrix
  
  print(lhd_term)
  
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
#' @export
update.blm <- function(model, formula = model$formula, prior=model, beta, ...) {
    
  if ( (formula != model$formula)  &  (prior == model) ) {
    #to do -> trim prior to fit formula
  }
  
  blm(formula, prior, beta, ...)

  }


#' Coefficients
#' 
#' Coefficients from the posterior distribution of a bayesian model.
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
#' Compute the confidence interval \code{\link{deviance}} for a bayesian linear
#'  model.
#'  
#' @param object a \code{blm} object.
#' @param param parameter for which confidence interval should be computed. Either a number or vector of numbers. If missig all parameters are considered.
#' @param level confidence level
#'
#' @return Pseudo confidence interval of parameter(s) of the \code{blm} object.
#' Returns the quantile specified by level for the distribution of \code{param}.
#'  
#' @details Deviance is the sum of squares of the residuals, extracted from the 
#'  maxium a posterior estimate of the \code{blm} object.
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
                         param, 
                         level = .05, ...) {

  if (missing(param)) {
    param <- seq_along(cf)
  }
  
  cf <- coef(object)[param]
  cf_var <- diag(cov.blm(object))[param]
  
  l <- level/2
  lower <- qnorm(l, cf, sqrt(cf_var))
  upper <- qnorm((1-l), cf, sqrt(cf_var))
  
  matrix(c(lower, upper),
         ncol=2,
         dimnames=list(names(cf),
                      paste(100 * c(l,1-l), '%')))
}




#' Predict 
#' 
#' Predict means and variance of the response variable for a \code{\link{blm}}
#'  object.
#'  
#' @param object a \code{blm} object.
#' @param newdata An optional data frame containing variables with which to 
#'  predict. If missing, values used for fitting will be extracted from the 
#'  object.
#' @param report.var Report also variance for the predicted values.
#'
#' @return A vector of predicted values. If report.var = TRUE, a list containing
#'  means and variances of the predicted values.
#'  
#' @export
predict.blm <- function(object, newdata = NULL, report.var = F, ...){
  responseless_formula <- delete.response(terms(object$formula))
  
  if ( missing(newdata) ) {
    matrix <- model.matrix(responseless_formula, object$frame)    
  } else {
    frame <- model.frame(responseless_formula, newdata)
    matrix <- model.matrix(responseless_formula, frame)
  }
  
  means <- apply(matrix, 1, function(xi) sum(coef(object) * xi))
  
  if (report.var) {
    vars <- apply(matrix, 1, function(xi) 
                    1/object$beta +  t(xi) %*% cov.blm(object) %*% xi)
    
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
  
  prediction <- predict.blm(object, report.var = F)
  
  response - prediction
}



#' Deviance of Bayesian Linear Model
#' 
#' \code{devaince.blm} returns \code{\link{devaince}} for a bayesian linear
#'  model.
#'  
#' @param object a \code{blm} object.
#'
#' @return Deviance of the \code{blm} object.
#'  
#' @details Deviance is the sum of squares of the residuals, extracted from the 
#'  maxium a posterior estimate of the \code{blm} object.
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
  sum( resid(model)^2 )
}



#' Print 
#' 
#' Print a blm object
#' @export
print.blm <- function(object, ...){
  cat('Call:\n')
  print(object$call)
  cat('\n')
  cat('Coefficients:\n')
  print(coef(object))
}



#' Summary 
#' 
#' Summary of a blm object
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
#' @param object a blm object.
#' @param explanatory a single number or vector for explanatory variables the 
#'  response will be plotted against.
#' 
#' @export
plot.blm <- function(object, explanatory = 1, ...){
  
  frame <- object[['frame']]
  for (exp in explanatory) {
    x_lab <- labels(terms(frame))[exp]
    x <- frame[,x_lab]
    
    response_lab <- attr(terms(frame), 'response')
    response <- model.response(frame)
    
    plot(x, response, type = "p", 
         xlab = x_lab, ylab = response_lab, ...)
    abline(lm(y~x), col = 'red') #regular lm fit
    abline(model, col='blue') #blm fit
    
    #ugly hack to get the 'se' of the fit plotted as line
    x_order <- order(frame[,x_lab])
    y_fit <- predict(object, report.var=T)
    upper <- qnorm(.95, mean = y_fit$mean, sd = sqrt(y_fit$var))
    lower <- qnorm(.05, mean = y_fit$mean, sd = sqrt(y_fit$var))
    lines(lower[x_order]~x[x_order], col='blue', lty=2)
    lines(upper[x_order]~x[x_order], col='blue', lty=2)
  }
  
}