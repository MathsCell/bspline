% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bspline.R
\name{dbsp}
\alias{dbsp}
\title{Derivative of B-spline}
\usage{
dbsp(f, nderiv = 1L)
}
\arguments{
\item{f}{Function, B-spline such as returned by \code{smbsp()} or \code{par2bsp()}}

\item{nderiv}{Integer scalar >= 0, order of derivative to calculate (1 by default)}
}
\value{
Function calculating requested derivative
}
\description{
Derivative of B-spline
}
\examples{
  x=seq(0., 1., length.out=11L)
  y=sin(2*pi*x)
  f=smbsp(x, y, nki=2L)
  d_f=dbsp(f)
  xf=seq(0., 1., length.out=101) # fine grid for plotting
  plot(xf, d_f(xf)) # derivative estimated by B-splines
  lines(xf, 2.*pi*cos(2*pi*xf), col="blue") # true derivative
  xk=bsppar(d_f)$xk
  points(xk, d_f(xk), pch="x", col="red") # knot positions
}