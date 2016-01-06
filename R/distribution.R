#' Multivariate distribution
#' 
#' Means and covariance matrix of a multivariate normal distribution .
#' 
#' @param means vector of means of the variables
#' @param covar covariance matrix of the variables
#' @param var_names optional names of the variables
#' @param ... other arguments (currently ignored).
#'   
#' @return An object of class \code{distribution} containing: 
#' \describe{
#'    \item{means}{vector of means}
#'    \item{covar}{covariance matrix}
#' } 
#' 
#' @examples
#' mv_dist(c(0,1),matrix(c(2,3,3,2),ncol=2))
#' 
#' @export
mv_dist <- function(means, covar, var_names, ...){
  
  d <- structure(list(means = drop(means), 
                 covar = drop(covar)),
            class = 'mv_dist')
  
  if ( !missing(var_names) ) {
    d <- mv_dist.set_var_names(d, var_names)
  }
  
  d
}



#' Means
#' 
#' Means of multivariate normal distribution object of class \code{mv_dist}.
#' 
#' @param d object of class \code{\link{mv_dist}}.
#' @param ... other arguments (currently ignored).
#' 
#' @export
mv_dist.means <- function(d, ...){
  d$means
}



#' Covariance
#' 
#' Covariance matrix of multivariate normal distribution object of class 
#'  \code{mv_dist}.
#' 
#' @param d object of class \code{\link{mv_dist}}.
#' @param ... other arguments (currently ignored).
#'
#' @export
mv_dist.covar <- function(d, ...){
  d$covar
}


#' Variable Names from a Multivariate Normal Distribution
#' 
#' Names of the variables of a multivariate normal distribution object of class 
#'  \code{mv_dist}.
#' 
#' @param d object of class \code{\link{mv_dist}}.
#' @param ... other arguments (currently ignored).
#'
#' @export
mv_dist.var_names <- function(d, ...){
  names(d$means)
}


#' Set variable Names from a Multivariate Normal Distribution
#' 
#' Set the names of the variables of a multivariate normal distribution object
#' of class \code{mv_dist}.
#' 
#' @param d object of class \code{\link{mv_dist}}.
#' @param var_names names of the variables of the distributions.
#' @param ... other arguments (currently ignored).
#'
#' @export
mv_dist.set_var_names <- function(d, var_names, ...){
  names(d$means) <- var_names
  
  if ( length(d$means) > 1 ) {
    colnames(d$covar) <- var_names
    rownames(d$covar) <- var_names
  } else {
    #univariate distribution case
    names(d$covar) <- var_names
  }
  
  d
}


#' Subset a Multivariate Normal Distribution
#' 
#' Return subset of variables of a multivariate normal distribution object of
#' class \code{mv_dist}.
#' 
#' @param object an object of class \code{mv_dist}.
#' @param subset names or indices of variables to extract from mv_dist.
#'   
#' @examples
#' d <- mv_dist(c(0,1),matrix(c(2,3,3,2),ncol=2), c('a', 'b'))
#' mv_dist.subset(d, 1)
#' mv_dist.subset(d, 'b')
#' 
#' @export
mv_dist.subset <- function(object, subset){
  object$means <- object$means[subset]
  object$covar <- object$covar[subset,subset]
  object
}


#' Distribution
#' 
#' Generic method returning the distribution of an object.
#' 
#' @param x object containing a distribution
distribution <- function(x) UseMethod("distribution")

#' Distribution
#' 
#' Default for generic method returning the distribution of an object.
#' 
#' @param x object containing a distribution
distribution.default <- function(x) x

#' Distribution
#' 
#' Distribution of a \code{blm} object.
#' 
#' @param x \code{blm} object.
distribution.blm <- function(x) x$posterior
