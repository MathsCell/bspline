# bspline

R package to build and use B-splines for interpolation and regression.
  In case of regression, equality constraints as well as monotonicity
  requirement can be imposed. Moreover, 
  knot positions (not only spline coefficients) can be part of 
  optimized parameters too. User is provided with 
  functions calculating spline values at arbitrary points. This 
  functions can be differentiated to obtain B-splines calculating 
  derivatives at any point. B-splines of this package can 
  simultaneously operate on a series of curves sharing the same set of 
  knots. 'bspline' is written with concern about computing 
  performance that's why the basis calculation is implemented in C++.
  The rest is implemented in R but without notable impact on computing speed.
