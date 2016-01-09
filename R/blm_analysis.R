#### MODEL ANALYSIS AND COMPARISON ####
# Note: Experimental and not thoroughly tested!



#' Probability of Response 
#' 
#' Probability of observing the response given a (bayesian) linear regression
#' model.
#' 
#' @param object a \code{\link{blm}} object.
#' @param a_prior prior value for parameter a of inv. gamma distribution.
#' @param b_prior prior value for parameter b of inv. gamma distribution.
#' @param use_log compute the log-probability.
#'    
#' @details Computes as described in 
#'   \url{https://en.wikipedia.org/wiki/Bayesian_linear_regression}. By default
#'   the log-probability is computed as non log-transformed probabilites are not
#'   solvable.
#' 
#' @export
pblm <- function(object, a_prior = 1, b_prior = 1, use_log = TRUE){
  
  resp <- model.response(object$frame)
  
  mu_prior <- object$prior$means
  cov_prior <- object$prior$covar
  
  mu_post <- coef(object)
  cov_post <- covar.blm(object)
  
  n <- nrow(object$matrix)
  
  a_post <- a_prior + n/2
  b_post <- drop(b_prior + 0.5*( t(resp)%*%resp + 
                                   t(mu_prior) %*% cov_prior %*% mu_prior -
                                   t(mu_post) %*% cov_post %*% mu_post ) )
  
  if (use_log) {
    log(1/((2*pi)^(n/2))) + 
      log(sqrt( det(as.matrix(cov_prior))/det(as.matrix(cov_post)) )) + 
      a_prior*log(b_prior) - a_post*log(b_post) +
      log(gamma(a_post)) - log(gamma(a_prior))
  } else {
    1/((2*pi)^(n/2)) * 
      sqrt( det(as.matrix(cov_prior))/det(as.matrix(cov_post)) ) * 
      b_prior^a_prior / b_post^a_post *
      gamma(a_post) / gamma(a_prior)
  }
}



#' Bayes Factor
#' 
#' Computes the Bayes factor for comparison of 2 \code{\link{blm}} objects.
#' 
#' @param model1 first \code{blm} object
#' @param model2 second \code{blm} object
#' 
#' @details Returns the ratio of the probability of model1 divided by the
#'   probability of model2. A value of K > 1 means that M1 is more strongly
#'   supported by the data under consideration than M2
#'   (\url{https://en.wikipedia.org/wiki/Bayes_factor}). Wikipedia mentions
#'   additional scales for interpretation, ie by Harold Jeffreys:
#'    \itemize{
#'      \item K < 1: supports model2
#'      \item K > 10: strong support for model1
#'      \item K > 100: decisive support for model1
#'    }.
#'  An alternative table, widely cited, is provided by Kass and Raftery.
#'    \itemize{
#'      \item K 1-3: not worth more than a bare mention
#'      \item K 3-20: positive support for model1
#'      \item K 20-150: strong support for model1
#'      \item K >150: very strong support for model1
#'    }.
#' 
#' @examples
#' set.seed(1) 
#' x <- seq(-10,10,.1) 
#' b <- 0.3
#' 
#' w0 <- 20 ; w1 <- 3 ; w2 <- 10
#' 
#' y <- rnorm(201, mean = w0 + w1 * x + w2 *sin(x), sd = sqrt(1/b)) 
#' model <- blm(y ~ x + sin(x), prior = NULL, data = data.frame(x=x, y=y))
#' 
#' plot(model, xlim=c(-10,10)) 
#' pblm(model) ## -955.8591 this is the log-probability
#' 
#' #another model removing the sinus term 
#' model2 <- blm(y ~ x, prior = NULL, data = data.frame(x=x, y=y))
#' 
#' pblm(model2) ## -958.128 
#' bayes_factor(model, model2) #9.668859
#' 
#' model3 <- blm(y ~ x + 0, prior = NULL, data = data.frame(x=x, y=y))
#' pblm(model3) ## -960.2644
#' bayes_factor(model, model3) #137.592
#' 
#' @export
bayes_factor <- function(model1, model2) {
  exp( pblm(model1) - pblm(model2) )
}



#' Bayesian information criterion (BIC)
#' 
#' Bayesian information criterion (BIC) for the fit of a \code{blm} object.
#' 
#' @param object a \code{blm} object.
#' @param ... other arguments (currently ignored).
#' 
#' @details Bayesian information criterion (BIC) of the \code{blm} object. 
#'   \eqn{BIC = n * ln(RSS/n) + k * ln(n)}. Where n is the number of 
#'   observations, k is the number of free parameters to be estimated and RSS is
#'   the residuals sum of squares (ie the \code{\link{deviance}}). When 
#'   comparing 2 models the strength of the evidence against the model with the 
#'   higher BIC value can be summarized as follows according to Wikipedia
#'   (\url{https://en.wikipedia.org/wiki/Bayesian_information_criterion}):
#'   \tabular{ll}{
#'     \bold{\eqn{\Delta}BIC} \tab \bold{Evidence against higher BIC}\cr
#'     0-2 \tab Not worth more than a bare mention\cr
#'     2-6 \tab Positive\cr
#'     6-10 \tab Strong\cr
#'     >10 \tab Very Strong\cr
#'   }
#' 
#' @examples
#' set.seed(1) 
#' x <- seq(-10,10,.1) 
#' b <- 0.3
#' 
#' w0 <- 0.2 ; w1 <- 3 ; w2 <- 10
#' 
#' y <- rnorm(201, mean = w0 + w1 * x + w2 *sin(x), sd = sqrt(1/b)) 
#' model <- blm(y ~ x + sin(x), prior = NULL, data = data.frame(x=x, y=y))
#' 
#' plot(model, xlim=c(-10,10)) 
#' bic(model) ## 221.6777
#' 
#' #another model removing the sinus term 
#' model2 <- blm(y ~ x, prior = NULL, data = data.frame(x=x, y=y))
#' 
#' bic(model2) ## 805.4814
#'   
#' @export 
bic <- function(object, ...){
  k <- ncol(object$matrix)
  n <- nrow(object$frame)
  
  n * log( deviance(object) / n ) + k * log(n)
}