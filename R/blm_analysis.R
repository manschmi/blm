#### MODEL ANALYSIS AND COMPARISON ####
# Note: Experimental and not thoroughly tested!



#' Model Probability
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
#' w0 <- 0.2 ; w1 <- 3 ; w2 <- 10
#' 
#' y <- rnorm(201, mean = w0 + w1 * x + w2 *sin(x), sd = sqrt(1/b)) 
#' mod1 <- blm(y ~ x + sin(x))
#' 
#' pblm(mod1) ## -886.4801 this is the log-probability (log=TRUE by default)
#' 
#' #another model removing the sinus term 
#' mod2 <- blm(y ~ x)
#' 
#' pblm(mod2) ## -891.0581 
#' bayes_factor(mod1, mod2) #97.31074 #strong support for mod1 over mod2
#' 
#' mod3 <- blm(y ~ x + 0)
#' pblm(mod3) ## -891.842
#' bayes_factor(mod1, mod3) #213.128 #very strong support for mod1 over mod3
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
#' mod1 <- blm(y ~ x + sin(x))
#' 
#' \dontrun{
#'   plot(mod1, xlim=c(-10,10)) 
#' }
#' bic(mod1) ## 224.351
#' 
#' #another mod removing the sinus term, clearly less well fitting
#' mod2 <- blm(y ~ x)
#' 
#' bic(mod2) ## 805.5729 -> huge difference in bic
#' 
#' ##much less distinguished model
#' b <- 0.003
#' y <- rnorm(201, mean = w0 + w1 * x + w2 *sin(x), sd = sqrt(1/b)) 
#' mod1 <- blm(y ~ x + sin(x))
#' \dontrun{
#'   plot(mod1, xlim=c(-10,10)) 
#' }
#' mod2 <- blm(y ~ x)
#' bic(mod1) #1209.717
#' bic(mod2) #1215.66 ... still positive support but not strong
#'   
#' @export 
bic <- function(object, ...){
  k <- ncol(model.matrix(object))
  n <- nrow(model.frame(object))
  
  n * log( deviance(object) / n ) + k * log(n)
}



#' F-Test
#' 
#' Anova F-test for comparison of the fits of 2 \code{blm} objects.
#' 
#' @param model1 first, more simple, \code{blm} model, nested in model2.#
#' @param model2 second, more complex, \code{blm} object.
#' @param ... other arguments (currently ignored).
#' 
#' @details Extracts relevant features and compute the F value and the p-value
#'   for comparison 2 model fits.
#'   
#' @return Analysis of Variance Table as for \code{\link{anova.lm}}.
#' 
#' @examples
#' set.seed(1) 
#' x <- seq(-10,10,.1) 
#' b <- 0.3
#' 
#' w0 <- 0.2 ; w1 <- 3 ; w2 <- 10
#' 
#' y <- rnorm(201, mean = w0 + w1 * x + w2 *sin(x), sd = sqrt(1/b)) 
#' mod1 <- blm(y ~ x + sin(x))
#' 
#' \dontrun{
#'   plot(mod1, xlim=c(-10,10)) 
#' }
#' 
#' #another mod removing the sinus term, clearly less well fitting
#' mod2 <- blm(y ~ x)
#' 
#' anova(mod1, mod2)
#' 
#' ##much less distinguished model
#' b <- 0.003
#' y <- rnorm(201, mean = w0 + w1 * x + w2 *sin(x), sd = sqrt(1/b)) 
#' mod1 <- blm(y ~ x + sin(x))
#' \dontrun{
#'   plot(mod1, xlim=c(-10,10)) 
#' }
#' mod2 <- blm(y ~ x)
#' anova(mod1, mod2) #... still very low p-value !?
#'   
#' @export 
anova.blm <- function(model1, model2, ...){
  dev1 <- deviance(model1)
  dev2 <- deviance(model2)
  df1 <- df.residual(model1)
  df2 <- df.residual(model2)
  F <- ( (dev2 - dev1)/(df2 - df1) ) / ( dev2/(df2+1) )
  pF <- 1-pf(F, df1, df2, lower.tail=F)
  
  struct <- structure(data.frame(Res.Df = c(df1, df2),
                 RSS = c(dev1, dev2),
                 Df = c(NA, df2-df1),
                 'Sum of Sq' = c(NA, dev2-dev1),
                 F = c(NA, F),
                 pF = c(NA, pF)),
            class = c('anova.blm', 'anova', 'data.frame'))
  
  attr(struct, 'model1') <- mod1
  attr(struct, 'model2') <- mod2
  
  struct
}


#' Print 
#' 
#' Print a anova.blm object
#' 
#' @param x a \code{\link{anova.blm}} object.
#' @param ... other arguments (currently ignored).
#' 
#' @export
print.anova.blm <- function(x, ...){
  cat('Analysis of Variance Table:\n\n')
  cat('Model1: ')
  print(attr(x, 'model1')$call)
  cat('Model2: ')
  print(attr(x, 'model2')$call)
  cat('\n')
  print(t(x[,1:6]))
}