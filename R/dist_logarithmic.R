#' The Logarithmic distribution
#'
#' @description
#' `r lifecycle::badge('stable')`
#'
#' @inheritParams actuar::dlogarithmic
#'
#' @seealso [actuar::Logarithmic]
#'
#' @examples
#' dist <- dist_logarithmic(prob = c(0.33, 0.66, 0.99))
#' dist
#'
#' @examplesIf requireNamespace("actuar", quietly = TRUE)
#' mean(dist)
#' variance(dist)
#' support(dist)
#' generate(dist, 10)
#'
#' density(dist, 2)
#' density(dist, 2, log = TRUE)
#'
#' cdf(dist, 4)
#'
#' quantile(dist, 0.7)
#' @name dist_logarithmic
#' @export
dist_logarithmic <- function(prob){
  prob <- vec_cast(prob, double())
  if(any((prob < 0) | (prob > 1))){
    abort("The prob parameter of a Logarithmic distribution must be between 0 and 1.")
  }
  new_dist(p = prob, class = "dist_logarithmic")
}

#' @export
format.dist_logarithmic <- function(x, digits = 2, ...){
  sprintf(
    "Logarithmic(%s)",
    format(x[["p"]], digits = digits, ...)
  )
}

#' @export
density.dist_logarithmic <- function(x, at, ...){
  require_package("actuar")
  actuar::dlogarithmic(at, x[["p"]])
}

#' @export
log_density.dist_logarithmic <- function(x, at, ...){
  require_package("actuar")
  actuar::dlogarithmic(at, x[["p"]], log = TRUE)
}

#' @export
quantile.dist_logarithmic <- function(x, p, ...){
  require_package("actuar")
  actuar::qlogarithmic(p, x[["p"]])
}

#' @export
cdf.dist_logarithmic <- function(x, q, ...){
  require_package("actuar")
  actuar::plogarithmic(q, x[["p"]])
}

#' @export
generate.dist_logarithmic <- function(x, times, ...){
  require_package("actuar")
  actuar::rlogarithmic(times, x[["p"]])
}

#' @export
mean.dist_logarithmic <- function(x, ...){
  p <- x[["p"]]
  (-1/(log(1-p)))*(p/(1-p))
}

#' @export
covariance.dist_logarithmic <- function(x, ...){
  p <- x[["p"]]
  -(p^2 + p*log(1-p))/((1-p)*log(1-p))^2
}
