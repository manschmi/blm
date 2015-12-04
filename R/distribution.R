#' Multivariate distribution
#' 
#' Means and covariance matrix of a multivariate normal distribution .
#' 
#' @param means vector of means of the variables
#' @param covar covariance matrix of the variables
#' @param ... currently ignored
#'   
#' @return An object of class \code{distribution} containing: 
#' \describe{
#'    \item{means}{vector of means}
#'    \item{covar}{covariance matrix}
#' } 
#' 
#' @export
mv_dist <- function(means, covar){
  structure(list(means = drop(means), 
                 covar = drop(covar)),
            class = 'mv_dist')
}



#' Means
#' 
#' Means of multivariate normal distribution object of class \code{mv_dist}.
#' 
mv_dist.means <- function(d, ...){
  d$means
}



#' Covariance
#' 
#' Covariance matrix of multivariate normal distribution object of class 
#'  \code{mv_dist}.
#' 
mv_dist.covar <- function(d, ...){
  d$covar
}


#' Variable Names from a Multivariate Normal Distribution
#' 
#' Names of the variables of a multivariate normal distribution object of class 
#'  \code{mv_dist}.
#' 
mv_dist.var_names <- function(d, ...){
  names(d$means)
}

#' Variable Names from a Multivariate Normal Distribution
#' 
#' Names of the variables of a multivariate normal distribution object of class 
#'  \code{mv_dist}.
#' 
mv_dist.set_var_names <- function(d, var_names, ...){
  names(d$means) <- var_names
  colnames(d$covar) <- var_names
  rownames(d$covar) <- var_names
  d
}


#' Subset a Multivariate Normal Distribution
#' 
#' Return ubset of variables of a multivariate normal distribution object of 
#'  class \code{mv_dist}.
#' 
#' @param object 
subset.mv_dist <- function(object, subset, ...){
  object$means <- object$means[subset]
  object$covar <- object$covar[subset,subset]
  object
}

#' Distribution
#' 
#' Generic method returning the distribution of an object.
#' 
distribution <- function(x) UseMethod("distribution")
distribution.default <- function(x) x
distribution.blm <- function(x) x$posterior
