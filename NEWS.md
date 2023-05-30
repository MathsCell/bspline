## v2.2.1 2023-??

 - fixed typo in rgl example of `bcurve()`
 - fixed NEWS mention of `parr()` for v2.2.

## v2.2 2023-05-26

 - added `bcurve()` for building nD curves from a set of control points
 - added `parr()` polynomial formulation of B-splines

## v2.1 2022-09-26

 - `smbsp()` can estimate covariance matrix of estimated coefficients
 - added `dmat()` for differentiation operator
 - `ipk()` is now exported and tested
 - fixed 0-value for abscissas on the right hand side of the knot interval

## v2.0.1 2022-04-19

 - `bsc()` can calculate Jacobian of basis vectors as function of knots.
    By consequence, the suggestion of `numDeriv` package is removed
 - added `jacw()` calculating Jacobian of B-spline with weights
 - added `ibsp()` for integration
 - now, `iknots()` can deal with repeated `x` values
 - added monotonicity and positivity optional constraints
 - added Copyright field to DESCRIPTION
 - fixed inequality generation in `fitsmbsp()`
 - fixed gcc-12 compile error

## v1.0.2 2022-03-17

 - fixed ipk() ASAN problem signaled by R CRAN team
 - fixed export of `par2bsp()`
 - fixed error in example of `fitsmbsp()`
 - fixed debug printing in `iknots()`
 - increased resolution from 5 to 10 intervals for tv fitting by linear
    B-splines in `iknots()`
 - forced monotonicity in `iknots()` despite possible round-off errors
 - added URL and BugRepport in DESCRIPTION

## v1.0.1 2022-03-11

 - First release on CRAN
