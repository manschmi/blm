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
#' mvnd(c(0,1),matrix(c(2,3,3,2),ncol=2))
#' 
#' @export
mvnd <- function(means, covar, var_names, ...){
  
  d <- structure(list(means = drop(means), 
                 covar = drop(covar)),
            class = 'mvnd')
  
  if ( !missing(var_names) ) {
    d <- mvnd.set_var_names(d, var_names)
  }
  
  d
}



#' Means
#' 
#' Means of multivariate normal distribution object of class \code{mvnd}.
#' 
#' @param d object of class \code{\link{mvnd}}.
#' @param ... other arguments (currently ignored).
#' 
#' @export
mvnd.means <- function(d, ...){
  d$means
}



#' Covariance
#' 
#' Covariance matrix of multivariate normal distribution object of class 
#'  \code{mvnd}.
#' 
#' @param d object of class \code{\link{mvnd}}.
#' @param ... other arguments (currently ignored).
#'
#' @export
mvnd.covar <- function(d, ...){
  d$covar
}


#' Variable Names from a Multivariate Normal Distribution
#' 
#' Names of the variables of a multivariate normal distribution object of class 
#'  \code{mvnd}.
#' 
#' @param d object of class \code{\link{mvnd}}.
#' @param ... other arguments (currently ignored).
#'
#' @export
mvnd.var_names <- function(d, ...){
  names(d$means)
}


#' Set variable Names from a Multivariate Normal Distribution
#' 
#' Set the names of the variables of a multivariate normal distribution object
#' of class \code{mvnd}.
#' 
#' @param d object of class \code{\link{mvnd}}.
#' @param var_names names of the variables of the distributions.
#' @param ... other arguments (currently ignored).
#'
#' @export
mvnd.set_var_names <- function(d, var_names, ...){
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
#' class \code{mvnd}.
#' 
#' @param object an object of class \code{mvnd}.
#' @param subset names or indices of variables to extract from mvnd.
#'   
#' @examples
#' d <- mvnd(c(0,1),matrix(c(2,3,3,2),ncol=2), c('a', 'b'))
#' mvnd.subset(d, 1)
#' mvnd.subset(d, 'b')
#' 
#' @export
mvnd.subset <- function(object, subset){
  object$means <- object$means[subset]
  object$covar <- object$covar[subset,subset]
  object
}



#' Mahalanobis Distance
#' 
#' Returns the Mahalanobis distance of a data point to a multivariate normal
#' distribution object of class \code{mvnd}.
#' 
#' @param object an object of class \code{mvnd}.
#' @param x a data point, provided as vector with same length as dimensions of
#'   the \code{mvnd} object.
#' 
#' @details Calls the base mahalanobis function \code{mahalanobis(x,
#'   object$means, object$covar)}
#'   
#' @examples
#' d <- mvnd(c(0,1),matrix(c(2,3,3,2),ncol=2), c('a', 'b'))
#' mahal(d, c(1,1))
#' mahal(d, c(0,2))
#' 
#' @export
mahal <- function(object, x){
  mahalanobis(x, object$means, object$covar)
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
