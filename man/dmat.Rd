% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bspline.R
\name{dmat}
\alias{dmat}
\title{Differentiation matrix of B-spline}
\usage{
dmat(nqw = NULL, xk = NULL, n = NULL, f = NULL, same_xk = FALSE, nderiv = 1L)
}
\arguments{
\item{nqw}{Integer scalar, row number of qw matrix (i.e. degree of freedom of a B-spline)}

\item{xk}{Numeric vector, knot positions}

\item{n}{Integer scalar, B-spline polynomial order}

\item{f}{Function from which previous parameters can be retrieved.
If both f and any of previous parameters are given then explicitly
set parameters take precedence over those retrieved from f.}

\item{same_xk}{Logical scalar, the same meaning as in \code{\link{dbsp}}}

\item{nderiv}{Integer scalar, order of differentiation (default 1)}
}
\value{
Numeric matrix of size \code{nqw-1 x nqw}
}
\description{
Calculate matrix for obtaining coefficients of first-derivative B-spline.
They can be calculated as \code{dqw=Md \%*\% qw}. Here, dqw are coefficients
of the first derivative,
Md is the matrix returned by this function, and qw are the coefficients
of differentiated B-spline.\cr
}
