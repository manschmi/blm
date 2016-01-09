#' blm: A package for bayesian linear regression.
#'
#' The \code{blm} package is for construction and analysis of bayesian linear
#' regression models. The core functions all operate on a S3 class \code{\link{blm}}.
#' 
#' 
#' @section Generic functions for blm objects:
#' \itemize{
#'   \item \link{update.blm}
#'   \item \link{print.blm}
#'   \item \link{summary.blm}
#'   \item \link{fitted.blm}
#'   \item \link{predict.blm}
#'   \item \link{summary.blm}
#'   \item \link{residuals.blm}
#'   \item \link{coef.blm}
#'   \item \link{deviance.blm}
#'   \item \link{plot.blm}
#' }
#' 
#' @section Other plots for analysis of blm objects:
#' \describe{
#'   \item{\link{add_sample_lines}}{Adds random sampled fit lines to a plot.blm.}
#'   \item{\link{parameter_kernel}}{Plot of mdistribution of coefficients.}
#'   \item{\link{diagnostic_plots}}{Diagnostic plots, similar to \code{\link{plot.lm}}.}
#' }
#' 
#' @section Other experimental functions for analysis of blm objects:
#' \describe{
#'   \item{\link{pblm}}{Total probability of the model.}
#'   \item{\link{bayes_factor}}{Bayes factor for comparison of 2 models.}
#'   \item{\link{bic}}{Bayesian information criterion (BIC) of the model.}
#' }
#' 
#' 
#' @docType package
#' @name blmr
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
#'   make_prior(0,1,2)
#'  
#'   make_prior(c(0,1),c(1,2))
#'  
#' @export
make_prior <- function(mean, alpha, dim = length(mean)){
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
#'   distribution \code{prior} of the weights is not provided the prior means 
#'   are set to 0 and the variances to 1.
#'   
#' @return A object of class \code{blm}. An object of class "lm" is a list 
#'   containing at least the following components: 
#'   \itemize{
#'     \item{call:        the matched call} 
#'     \item{formula:     the formula used} 
#'     \item{df.residual: the degrees of freedom of the model}
#'     \item{frame:       the model frame used} 
#'     \item{matrix:      the model matrix used} 
#'     \item{beta:        the precision of the data} 
#'     \item{prior:       the prior distribution used} 
#'     \item{posterior:   the posterior distribution}
#'   }
#'   
#' @examples
#'   w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3 
#'   x <- rnorm(50)
#'   y <- rnorm(50, w1 * x + w0, 1/b) 
#'   mod <- blm(y~x, data=data.frame(x=x, y=y))
#'   mod
#'   plot(mod)
#'   
#'   #use with known precision
#'   mod_det <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
#'   mod_det
#'   
#'   #use of a prior, typically from an existing model 
#'   x2 <- rnorm(50) 
#'   y2 <- rnorm(50, w1 * x2 + w0, 1/b) 
#'   mod2 <- blm(y~x, prior=mod, data=data.frame(x=x2, y=y2)) 
#'   mod2
#'   
#'   #use with 2 explanatory variables 
#'   w2 <- 3.3 
#'   z <- rnorm(50) 
#'   y <- rnorm(50, w2 * z + w1 * x + w0, 1/b) 
#'   mod <- blm(y~x+z, data=data.frame(x=x, y=y, z=z)) 
#'   mod
#'  
#' @export
blm <- function(formula, prior = NULL, beta, ...) {
  
  frame <- model.frame(as.formula(formula), ...)
  matrix <- model.matrix(as.formula(formula), frame)
  df <- nrow(matrix) - ncol(matrix)
  response <- model.response(frame)
  coef_names <- colnames(matrix)
  
  if (is.null(prior)) {
    n <- ncol(matrix)
    prior <- make_prior(mean = 0, alpha = 1, dim = n)
    prior <- mv_dist.set_var_names(prior, coef_names)
  }
  
  prior <- distribution(prior)
  
  prior_names <- mv_dist.var_names(prior)
  if ( length(prior_names) != length(coef_names) & !all(prior_names %in% coef_names) ) {
    #can occur ie when updating to a new formula with fewer terms
    if ( all(coef_names %in% prior_names) ){
      message('prior contains more variables than the model : variables not used in the model are ignored in the fit')
      prior <- mv_dist.subset(prior, colnames(matrix))
    } else {
      missing <- coef_names[!(coef_names %in% prior_names)]
      stop(paste('model formula contains variable names not provided in the prior', missing ))
    }
  }
  
  prior_means <- mv_dist.means(prior)
  prior_covar <- mv_dist.covar(prior)
  prior_precision <- solve(prior_covar)
  
  if (missing(beta)) {
    covar <- solve( prior_precision + t(matrix) %*% matrix )
    
    means <- covar %*% 
      drop(prior_precision %*% prior_means + 
             t(matrix) %*% matrix(response,ncol=1))
    
    beta <- 1/(sum((response - drop(matrix %*% means))^2)/df)
  } else {
    covar <- solve(prior_covar + beta * t(matrix) %*% matrix)
    means <- beta * covar %*% t(matrix) %*% response
  }
  
  posterior <- distribution(mv_dist(means, covar))
  
  structure(list(call = match.call(),
                        formula = formula,
                        frame = frame,
                        matrix = matrix,
                        df.residual = df,
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
#' @details Updates a \code{\link{blm}} object, using features of the input 
#'   object, except, if provided otherwise. The prior for the updated object is 
#'   typically the posterior distribution of the input object, but can be 
#'   specified specifically (as for all other parameters). Importantly, new data
#'   will only be used when specified as named argument. \cr The function can
#'   also be used to update the model formula, as described for
#'   \code{\link{update}}.
#'   
#' @return A object of class \code{blm}.
#'   
#' @examples
#'   w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3 
#'   x <- rnorm(50) 
#'   y <- rnorm(50, w1 * x + w0, 1/b) 
#'   mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y)) 
#'   mod
#'   
#'   #use of a prior, typically from an existing model 
#'   x2 <- rnorm(50) 
#'   y2 <- rnorm(50, w1 * x2 + w0, 1/b)
#'   
#'   #using posterior of mod as prior 
#'   new_mod <- update(mod, data=data.frame(x=x2, y=y2)) 
#'   new_mod
#'   
#'   #using same prior for mod and new model 
#'   new_mod2 <- update(mod, prior=mod$prior, data=data.frame(x=x2, y=y2)) 
#'   new_mod2
#'   
#'   #update model formula 
#'   new_mod2 <- update(mod, y~x+0, prior=mod$prior, 
#'                      data=data.frame(x=x2, y=y2)) 
#'   new_mod2
#'   
#'   #also works with standard R formula update semantics 
#'   new_mod2 <- update(mod, ~.+0, prior=mod$prior, data=data.frame(x=x2, y=y2)) 
#'   new_mod2
#'    
#' 
#' @export
update.blm <- function(object, 
                       formula = object$formula, 
                       beta = object$beta,
                       prior = object, 
                       data = object$frame, 
                       ...) {
  
  co <- getCall(object)
  cm <- match.call()
  
  updated_formula = update(object$formula, formula)
  co$formula <- updated_formula
  
  co$beta <- beta
  
  if ( missing(prior) ) {
    co$prior <- 'object$posterior'
  } else {
    co$prior <- cm$prior
  }
  
  if ( !missing(data) ) {
    co$data <- cm$data
  }

  mod <- blm(updated_formula, prior=prior, beta=beta, data=data,  ...) 
  mod$call <- co
  mod
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
#' @param parm parameter for which confidence interval should be computed.
#'   Either a number or vector of numbers. If missig all parameters are
#'   considered.
#' @param level confidence level
#' @param ... other arguments (currently ignored).
#'
#' @details Deviance is the sum of squares of the residuals, extracted from the 
#'  maxium a posterior estimate of the \code{blm} object.
#'  
#' @return Pseudo confidence interval of parameter(s) of the \code{blm} object. 
#'   Returns the quantile specified by level for the distribution of
#'   \code{param}.
#'  
#' @examples
#'   x <- rnorm(10)
#'   b <- 1.3
#'   w0 <- 0.2 ; w1 <- 3
#'   y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#'   model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#'   
#'   confint(model)
#'   confint(model, 1)
#'   confint(model, c(1,2))
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
#'   object.
#'  
#' @param object a \code{blm} object.
#' @param newdata an optional data frame containing variables with which to 
#'   predict. If missing, values used for fitting will be extracted from the 
#'   object.
#' @param report.var report also variance for the predicted values.
#' @param ... other arguments (currently ignored).
#' 
#' @return A vector of predicted values. If report.var = TRUE, a list containing
#'   means and variances of the predicted values.
#'  
#' @examples
#'   x <- rnorm(100)
#'   b <- 1.3
#'   w0 <- 0.2 ; w1 <- 3
#'   y <- rnorm(100, mean = w0 + w1 * x, sd = sqrt(1/b))
#'   model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#'   
#'   predict(model)
#'   
#'   x <- rnorm(10) 
#'   predict(model, data.frame(x=x))
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
    error_var <- 1/object$beta
    
    vars <- apply(matrix, 1, function(xi) 
      error_var +  t(xi) %*% covar.blm(object) %*% xi)
    
    return(list(mean = means, var = vars))
  }
  
  means
}



#' Extract Model Fitted Values 
#' 
#' Extracts the fitted values for the response variable of a \code{\link{blm}}
#'   object.
#'  
#' @param object a \code{blm} object.
#' @param report.var Report also variance for the fitted values.
#' @param ... other arguments (currently ignored).
#'
#' @return A vector of fitted values. If report.var = TRUE, a list containing 
#'   means and variances of the fitted values. Identical to \code{predict.blm} 
#'   without newdata.
#'  
#' @examples
#'   x <- rnorm(10)
#'   b <- 1.3
#'   w0 <- 0.2 ; w1 <- 3
#'   y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#'   model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#'   
#'   fitted(model)
#'   
#' @export
fitted.blm <- function(object, report.var = FALSE, ...){
  
  responseless_formula <- delete.response(terms(object$formula))
  matrix <- model.matrix(responseless_formula, object$frame)    
  
  means <- drop(matrix %*% coef(object))
  
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
#'   \code{blm} object.
#'  
#' @examples
#'   x <- rnorm(10)
#'   b <- 1.3
#'   w0 <- 0.2 ; w1 <- 3
#'   y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#'   model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#'   
#'   residuals(model)
#'   resid(model)
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
#'   model.
#'  
#' @param object a \code{blm} object.
#' @param ... other arguments (currently ignored).
#' 
#' @details Deviance is the sum of squares of the residuals, extracted from the 
#'   maxium a posterior estimate of the \code{blm} object.
#'  
#' @return Deviance of the \code{blm} object.
#' 
#' @examples
#'   x <- rnorm(10)
#'   b <- 1.3
#'   w0 <- 0.2 ; w1 <- 3
#'   y <- rnorm(10, mean = w0 + w1 * x, sd = sqrt(1/b))
#'   model <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#'   
#'   deviance(model)
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
#'   \describe{
#'     \item{call}{the call of the object }
#'     \item{prior}{the prior distribution of the fit}
#'     \item{posterior}{the posterior distribution of the fit}
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



#' Print 
#' 
#' Print a summary.blm object
#' 
#' @param x a \code{\link{summary.blm}} object.
#' @param ... other arguments (currently ignored).
#' 
#' @export
print.summary.blm <- function(x, ...){
  cat('Call:\n')
  print(x$call)
  cat('\n----\n')
  cat('Posterior Coefficients:\n')
  cat('\nEstimate:\n')
  print(coef(x))
  cat('\nCovariance:\n')
  print(covar.blm(x))
  cat('\n----\n')
  cat('Prior Coefficients:\n')
  cat('\nEstimate:\n')
  print(x$prior$means)
  cat('\nCovariance:\n')
  print(x$prior$covar)
}