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
#'   mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
#'   
#'   summary(mod) 
#'   
#'   plot(mod)
#'   
#'   #use of a prior, typically from an existing model 
#'   x2 <- rnorm(50) 
#'   y2 <- rnorm(50, w1 * x2 + w0, 1/b) 
#'   mod <- blm(y~x, prior=mod, beta=b, data=data.frame(x=x2, y=y2)) 
#'   mod
#'   
#'   #use with 2 explanatory variables 
#'   w2 <- 3.3 
#'   z <- rnorm(50) 
#'   y <- rnorm(50, w2 * z + w1 * x + w0, 1/b) 
#'   mod <- blm(y~x+z, beta=b, data=data.frame(x=x,y=y, z=z)) 
#'   mod
#'  
#' @export
blm <- function(formula, prior = NULL, beta = NULL, ...) {
  
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
  if ( length(prior_vars) != length(mod_vars) & !all(prior_vars %in% mod_vars) ) {
    #can occur ie when updating to a new formula with fewer terms
    if ( all(mod_vars %in% prior_vars) ){
      message('prior contains more variables than the model : variables not used in the model are ignored in the fit')
      prior <- mv_dist.subset(prior, colnames(matrix))
    } else {
      missing <- mod_vars[!(mod_vars %in% prior_vars)]
      stop(paste('model formula contains variable names not provided in the prior', missing ))
    }
  }
  
  prior_means <- mv_dist.means(prior)
  prior_covar <- mv_dist.covar(prior)
  
  response <- model.response(frame)
  
  covar <- solve(prior_covar + t(matrix) %*% matrix)
  
  means <- covar %*% 
    drop(prior_means %*% prior_covar + t(t(matrix) %*% matrix(response,ncol=1)))
  
  posterior <- distribution(mv_dist(means, covar))
  
  structure(list(call = match.call(),
                        formula = formula,
                        frame = frame,
                        matrix = matrix,
                        df.residual = nrow(matrix) - ncol(matrix),
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
    if (is.null(object$beta)) {
      error_var(object)
    } else {
      var <- 1/object$beta
    }
    
    vars <- apply(matrix, 1, function(xi) 
      var +  t(xi) %*% covar.blm(object) %*% xi)
    
    return(list(mean = means, var = vars))
  }
  
  means
}



#' Variance of the Error 
#' 
#' Variance of the error of a (bayesian) linear regression model.
#' 
#' @details Simply \code{deviance/df}. Required to compute the distribution of
#'   fitted values.
error_var <- function(object){
  n <- nrow(object$matrix)
  p <- ncol(object$matrix)
  deviance(object)/(n-p)
}



#' Variance of the Error 
#' 
#' Variance of the error of a (bayesian) linear regression model.
#' 
#' @param object a \code{\link{blm}} object.
#' 
#' @details Simply \code{deviance/df}. Required to compute the distribution of
#'   fitted values.
error_var <- function(object){
  n <- nrow(object$matrix)
  p <- ncol(object$matrix)
  deviance(object)/(n-p)
}



#' Probability of Response 
#' 
#' Probability of observing the response given a (bayesian) linear regression model.
#' 
#' @param object a \code{\link{blm}} object.
#' 
#' @details Computes as described in \url{https://en.wikipedia.org/wiki/Bayesian_linear_regression}.
pY <- function(object, prior_mean = 1, prior_var = 1){
  a_prior <- prior_mean/2
  b_prior <- (prior_mean * prior_var)/2
  
  resp <- model.response(object$frame)
  
  mu_prior <- object$prior$means
  cov_prior <- object$prior$covar
  
  mu_post <- object$posterior$means
  cov_post <- object$posterior$covar
  
  n <- nrow(object$matrix)
  
  a_post <- a_prior + n/2
  b_post <- b_prior + 0.5*( t(resp)%*%resp + 
                                  t(mu_prior) %*% cov_prior %*% mu_prior -
                                  t(mu_post) %*% cov_post %*% mu_post )
  
  1/((2*pi)^(n/2)) * 
    sqrt( det(cov_prior)/det(cov_post) ) * 
    (b_prior^a_prior) / (b_post^a_post) *
    gamma(a_post) / gamma(a_prior)
}



#' Bayes Factor
#' 
#'  Computes the Bayes factor for comparison of 2 \code{\link{blm}} objects.
#'  
#'  @param model1 first \code{blm} object
#'  @param model2 second \code{blm} object
#'  
#'  @details Returns the ratio of the probability of model1 divided by the probability of model2. A value of K > 1 means that M1 is more strongly supported by the data under consideration than M2 (\url{https://en.wikipedia.org/wiki/Bayes_factor}). Wikipedia mentions additional scales for interpretation, ie by Harold Jeffreys: 
#'    \itemize{
#'      \item K < 1: supports model2
#'      \item K > 10: strong support for model1
#'      \item K > 100: decisive support for model1
#'    }.
#'  An alternative table, widely cited, is provided by Kass and Raftery:
#'    \itemize{
#'      \item K 1-3: not worth more than a bare mention
#'      \item K 3-20: positive support for model1
#'      \item K 20-150: strong support for model1
#'      \item K >150: strong support for model1
#'    }.
#' 
bayes_factor <- function(model1, model2) {
  pY(model1) / pY(model2)
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




#' Bayesian information criterion (BIC)
#' 
#' Bayesian information criterion (BIC) for the fit of a \code{blm} object.
#' 
#' @param object a \code{blm} object.
#' @param ... other arguments (currently ignored).
#' 
#' @details Bayesian information criterion (BIC) of the \code{blm} object. 
#'   \code{BIC} = n * ln(RSS/n) + k * ln(n). Where n is the number of 
#'   observations, k is the number of free parameters to be estimated and RSS is
#'   the residuals sum of squares (ie the \code{\link{deviance}}). When 
#'   comparing 2 models the strength of the evidence against the model with the 
#'   higher BIC value can be summarized as follows according to Wikipedia
#'   (\url{https://en.wikipedia.org/wiki/Bayesian_information_criterion}):
#'   \tabular{ll}{
#'     ΔBIC \tab Evidence against higher BIC\cr
#'     0-2 \tab Not worth more than a bare mention\cr
#'     2-6 \tab Positive\cr
#'     6-10 \tab Strong\cr
#'     >10 \tab Very Strong\cr
#'   }
#' 
#' @examples
#' set.seed(1) 
#' x <- seq(-10,10,.1) 
#' b <- 3
#' 
#' w0 <- 0.2 ; w1 <- 3 ; w2 <- 10
#' 
#' y <- rnorm(201, mean = w0 + w1 * x + w2 *sin(x), sd = sqrt(1/b)) 
#' model <- blm(y ~ x + sin(x), prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' plot(model, xlim=c(-10,10)) 
#' bic(model) ## -1.582287
#' 
#' #another model removing the sinus term 
#' model2 <- blm(y ~ x, prior = NULL, beta = b, data = data.frame(x=x, y=y))
#' 
#' bic(model2) ## 785.2877
#'   
#' @export 
bic <- function(object, ...){
  k <- ncol(object$matrix)
  n <- nrow(object$frame)
  
  n * log( deviance(object) / n ) + k * log(n)
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
  print(x$posterior$means)
  cat('\nCovariance:\n')
  print(x$posterior$covar)
  cat('\n----\n')
  cat('Prior Coefficients:\n')
  cat('\nEstimate:\n')
  print(x$prior$means)
  cat('\nCovariance:\n')
  print(x$prior$covar)
}



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
  y_blm <- predict(object, data, report.var=T)
  
  #MAP fit lines
  y_blm_map <- y_blm$mean
  lines(y_blm_map~x_fit, col = blm_col, lty= blm_map_lty)
  
  #add quantiles of the fit
  if (show_blm_interval) {
    upper <- qnorm(blm_interval_level, mean = y_blm$mean, sd = sqrt(y_blm$var))
    lower <- qnorm(1-blm_interval_level, mean = y_blm$mean, sd = sqrt(y_blm$var))
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
  fits <- fitted(object, report.var = T)
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
  x <- fitted(object, report.var = FALSE)
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
  x <- fitted(object, report.var = FALSE)
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