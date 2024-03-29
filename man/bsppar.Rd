% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bspline.R
\name{bsppar}
\alias{bsppar}
\title{Retrieve parameters of B-splines}
\usage{
bsppar(f)
}
\arguments{
\item{f}{Function, B-splines such that returned by par3bsp(), smbsp(), ...}
}
\value{
List having components: n - polynomial order, qw - coefficients, xk -
 knots
}
\description{
Retrieve parameters of B-splines
}
