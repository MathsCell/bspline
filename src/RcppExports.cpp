// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/bspline.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ipk
umat ipk(const vec& x, const vec& xk);
RcppExport SEXP _bspline_ipk(SEXP xSEXP, SEXP xkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const vec& >::type xk(xkSEXP);
    rcpp_result_gen = Rcpp::wrap(ipk(x, xk));
    return rcpp_result_gen;
END_RCPP
}
// bsc
SEXP bsc(const vec& x, const vec& xk, const size_t n, const bool cjac);
RcppExport SEXP _bspline_bsc(SEXP xSEXP, SEXP xkSEXP, SEXP nSEXP, SEXP cjacSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const vec& >::type xk(xkSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< const bool >::type cjac(cjacSEXP);
    rcpp_result_gen = Rcpp::wrap(bsc(x, xk, n, cjac));
    return rcpp_result_gen;
END_RCPP
}
// jacw
cube jacw(const cube& jac, const RObject& qws);
RcppExport SEXP _bspline_jacw(SEXP jacSEXP, SEXP qwsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const cube& >::type jac(jacSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type qws(qwsSEXP);
    rcpp_result_gen = Rcpp::wrap(jacw(jac, qws));
    return rcpp_result_gen;
END_RCPP
}
// parr
cube parr(const vec& xk, const size_t n);
RcppExport SEXP _bspline_parr(SEXP xkSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type xk(xkSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(parr(xk, n));
    return rcpp_result_gen;
END_RCPP
}
// pbsc
mat pbsc(const vec& x, const vec& xk, const cube& coeffs);
RcppExport SEXP _bspline_pbsc(SEXP xSEXP, SEXP xkSEXP, SEXP coeffsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const vec& >::type xk(xkSEXP);
    Rcpp::traits::input_parameter< const cube& >::type coeffs(coeffsSEXP);
    rcpp_result_gen = Rcpp::wrap(pbsc(x, xk, coeffs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bspline_ipk", (DL_FUNC) &_bspline_ipk, 2},
    {"_bspline_bsc", (DL_FUNC) &_bspline_bsc, 4},
    {"_bspline_jacw", (DL_FUNC) &_bspline_jacw, 2},
    {"_bspline_parr", (DL_FUNC) &_bspline_parr, 2},
    {"_bspline_pbsc", (DL_FUNC) &_bspline_pbsc, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_bspline(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
