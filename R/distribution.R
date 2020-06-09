#' Create a new distribution
#'
#' @param ... Parameters of the distribution (named).
#' @param class The class of the distribution for S3 dispatch.
#' @param dimnames The names of the variables in the distribution (optional).
#'
#' @export
new_dist <- function(..., class = NULL, dimnames = NULL){
  args <- transpose(vctrs::vec_recycle_common(...))
  wrap_dist(
    lapply(args, structure, class = c(class, "dist_default")),
    dimnames = dimnames
  )
}

wrap_dist <- function(x, dimnames = NULL){
  vctrs::new_vctr(x, vars = dimnames, class = "distribution")
}

#' @export
vec_ptype_abbr.distribution <- function(x, ...){
  "dist"
}

#' @export
format.distribution <- function(x, ...){
  x <- vec_data(x)
  out <- vapply(x, format, character(1L), ...)
  out[vapply(x, is.null, logical(1L))] <- "?"
  out
}

#' @export
`dimnames<-.distribution` <- function(x, value){
  attr(x, "vars") <- value
  x
}
#' @export
dimnames.distribution <- function(x){
  attr(x, "vars")
}

#' The probability density/mass function
#'
#' \lifecycle{stable}
#'
#' Computes the probability density function for a continuous distribution, or
#' the probability mass function for a discrete distribution.
#'
#' @param x The distribution(s).
#' @param at The point at which to compute the density/mass.
#' @param ... Additional arguments passed to methods.
#'
#' @importFrom stats density
#' @export
density.distribution <- function(x, at, ...){
  vec_is(at, double(), 1L)
  dist_apply(x, density, at = at, ...)
}

#' Distribution Quantiles
#'
#' \lifecycle{stable}
#'
#' Computes the quantiles of a distribution.
#'
#' @param x The distribution(s).
#' @param p The probability of the quantile.
#' @param ... Additional arguments passed to methods.
#'
#' @importFrom stats quantile
#' @export
quantile.distribution <- function(x, p, ...){
  vec_is(p, double(), 1L)
  dist_apply(x, quantile, p = p, ...)
}

#' The cumulative distribution function
#'
#' \lifecycle{stable}
#'
#' @param x The distribution(s).
#' @param q The quantile at which the cdf is calculated.
#' @param ... Additional arguments used by methods.
#'
#' @name cdf
#' @export
cdf <- function (x, q, ...){
  ellipsis::check_dots_used()
  UseMethod("cdf")
}
#' @rdname cdf
#' @export
cdf.distribution <- function(x, q, ...){
  vec_is(q, double(), 1L)
  dist_apply(x, cdf, q = q, ...)
}

#' Randomly sample values from a distribution
#'
#' \lifecycle{stable}
#'
#' Generate random samples from probability distributions.
#'
#' @param x The distribution(s).
#' @param times The number of samples.
#' @param ... Additional arguments used by methods.
#'
#' @export
generate.distribution <- function(x, times, ...){
  times <- vec_cast(times, integer())
  times <- vec_recycle(times, size = length(x))
  mapply(generate, vec_data(x), times = times, ..., SIMPLIFY = FALSE)
  # dist_apply(x, generate, times = times, ...)
  # Needs work to structure MV appropriately.
}

#' Mean of a probability distribution
#'
#' \lifecycle{stable}
#'
#' Returns the empirical mean of the probability distribution. If the method
#' does not exist, the mean of a random sample will be returned.
#'
#' @param x The distribution(s).
#' @param ... Additional arguments used by methods.
#'
#' @export
mean.distribution <- function(x, ...){
  dist_apply(x, mean, ...)
}

#' Variance
#'
#' A generic function for computing the variance of an object. The default
#' method will use [`stats::var()`] to compute the variance.
#'
#' @param x An object.
#' @param ... Additional arguments used by methods.
#'
#' @seealso [`variance.distribution()`]
#'
#' @export
variance <- function(x, ...){
  UseMethod("variance")
}
#' @export
variance.default <- function(x, ...){
  stats::var(x, ...)
}

#' Variance of a probability distribution
#'
#' \lifecycle{stable}
#'
#' Returns the empirical mean of the probability distribution. If the method
#' does not exist, the mean of a random sample will be returned.
#'
#' @param x The distribution(s).
#' @param ... Additional arguments used by methods.
#'
#' @export
variance.distribution <- function(x, ...){
  dist_apply(x, variance, ...)
}

#' Median of a probability distribution
#'
#' \lifecycle{stable}
#'
#' Returns the median (50th percentile) of a probability distribution. This is
#' equivalent to `quantile(x, p=0.5)`.
#'
#' @inheritParams stats::median
#' @param x The distribution(s).
#' @param ... Additional arguments used by methods.
#'
#' @importFrom stats median
#' @export
median.distribution <- function(x, na.rm = FALSE, ...){
  quantile(x, p = 0.5, na.rm = na.rm, ...)
}

#' Probability intervals of a probability distribution
#'
#' \lifecycle{maturing}
#'
#' Returns a `hilo` central probability interval with probability coverage of
#' `size`. By default, the distribution's [`quantile()`] will be used to compute
#' the lower and upper bound for a centered interval
#'
#' @param x The distribution(s).
#' @param size The size of the interval (between 0 and 100).
#' @param ... Additional arguments used by methods.
#'
#' @seealso [`hdr.distribution()`]
#'
#' @importFrom stats median
#' @export
hilo.distribution <- function(x, size = 95, ...){
  dist_apply(x, hilo, size = size, ...)
}

#' Highest density regions of probability distributions
#'
#' \lifecycle{experimental}
#'
#' This function is highly experimental and will change in the future. In
#' particular, improved functionality for object classes and visualisation tools
#' will be added in a future release.
#'
#' Computes minimally sized probability intervals highest density regions.
#'
#' @param x The distribution(s).
#' @param size The size of the interval (between 0 and 100).
#' @param n The resolution used to estimate the distribution's density.
#' @param ... Additional arguments used by methods.
#'
#' @export
hdr.distribution <- function(x, size = 95, n = 512, ...){
  dist_x <- vapply(seq(0.5/n, 1 - 0.5/n, length.out = n), quantile, numeric(1L), x = x)
  dist_y <- vapply(dist_x, density, numeric(1L), x = x)
  alpha <- quantile(dist_y, probs = size/100)

  crossing_alpha <- function(alpha, x, y){
    it <- seq_len(length(y) - 1)
    dd <- y - alpha
    dd <- dd[it + 1] * dd[it]
    index <- it[dd <= 0]
    # unique() removes possible duplicates if sequential dd has same value.
    # More robust approach is required.
    unique(
      vapply(
        index,
        function(.x) stats::approx(y[.x + c(0,1)], x[.x + c(0,1)], xout = alpha)$y,
        numeric(1L)
      )
    )
  }

  # purrr::map(alpha, crossing_alpha, dist_x, dist_y)
  hdr <- crossing_alpha(alpha, dist_x, dist_y)
  lower_hdr <- seq_along(hdr)%%2==1
  hdr <- new_hilo(hdr[lower_hdr], hdr[!lower_hdr], size = size)
  new_hdr(list(hdr))
}

#' @export
sum.distribution <- function(x, ...){
  Reduce("+", x)
}

#' @method vec_arith distribution
#' @export
vec_arith.distribution <- function(op, x, y, ...){
  UseMethod("vec_arith.distribution", y)
}
#' @method vec_arith.distribution default
#' @export
vec_arith.distribution.default <- function(op, x, y, ...){
  if(is_empty(y)){
    out <- lapply(x, get(op))
  }
  else {
    x <- vec_recycle_common(x = x, y = y)
    y <- x[["y"]]
    x <- x[["x"]]
    out <- mapply(get(op), x = x, y = y, SIMPLIFY = FALSE)
  }
  vec_restore(out, x)
}

#' @method vec_arith.numeric distribution
#' @export
vec_arith.numeric.distribution <- function(op, x, y, ...){
  x <- vec_recycle_common(x = x, y = y)
  y <- x[["y"]]
  x <- x[["x"]]
  out <- mapply(get(op), x = x, y = y, SIMPLIFY = FALSE)
  vec_restore(out, y)
}

#' @method vec_math distribution
#' @export
vec_math.distribution <- function(.fn, .x, ...) {
  if(.fn %in% c("is.nan", "is.infinite")) return(rep_len(FALSE, length(.x)))
  if(.fn == "is.finite") return(rep_len(TRUE, length(.x)))
  out <- lapply(.x, get(.fn), ...)
  vec_restore(out, .x)
}

#' @export
vec_ptype2.distribution.distribution <- function(x, y, ...){
  if(!identical(dimnames(x), dimnames(y))){
    abort("Distributions must have the same `dimnames` to be combined.")
  }
  x
}
#' @export
vec_ptype2.double.distribution <- function(x, y, ...) new_dist()
#' @export
vec_ptype2.distribution.double <- function(x, y, ...) new_dist()
#' @export
vec_ptype2.integer.distribution <- function(x, y, ...) new_dist()
#' @export
vec_ptype2.distribution.integer <- function(x, y, ...) new_dist()

#' @export
vec_cast.distribution.distribution <- function(x, to, ...){
  dimnames(x) <- dimnames(to)
  x
}
#' @export
vec_cast.distribution.double <- function(x, to, ...){
  x <- dist_degenerate(x)
  dimnames(x) <- dimnames(to)
  x
}
#' @export
vec_cast.distribution.integer <- vec_cast.distribution.double