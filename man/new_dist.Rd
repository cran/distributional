% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{new_dist}
\alias{new_dist}
\title{Create a new distribution}
\usage{
new_dist(..., class = NULL, dimnames = NULL)
}
\arguments{
\item{...}{Parameters of the distribution (named).}

\item{class}{The class of the distribution for S3 dispatch.}

\item{dimnames}{The names of the variables in the distribution (optional).}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

Allows extension package developers to define a new distribution class
compatible with the distributional package.
}
