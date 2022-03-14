// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Intervals of points in knot intervals
//
// Find first and last+1 indexes iip s.t. x[iip] belongs to interval starting at xk[iik]
//
// @param x Numeric vector, abscissa points (must be non decreasing)
// @param xk Numeric vector, knots (must be non decreasing)
// @return Integer matrix of size \code{(2 x length(xk)-1)}. Indexes are 0-based
umat ipk(const vec& x, const vec& xk) {
    // find first and last+1 indexes iip s.t. x[iip] belongs to interval starting at xk[iik]
    // if not found, first==last+1, i.e. the length of interval is 0.
    auto nk=xk.n_elem, np=x.n_elem;
    umat ip(2, nk-1);
    double xkmax=xk[nk-1];
    //x.print("x");
    //xk.print("xk");
    //Rcout << "&x=" << &x[0] << std::endl;
    //Rcout << "np=" << np << std::endl;
    for (size_t iik=0, iip=0; iik < nk-1; iik++) {
        double xki=xk[iik];
        //Rcout << "beging: iik, iip=" << iik << ", " << iip << std::endl;
        while ((iip < np) && (xki == xkmax || (iip == np-1 && iip > 0) ? x[iip] <= xki : x[iip] < xki))
            iip++;
        //Rcout << "iik, iip=" << iik << ", " << iip << std::endl;
        ip(0, iik)=iip;
        if (iik > 0)
            ip(1, iik-1)=iip;
        if (iik == nk-2)
            ip(1, iik)=(iip == np-1 && iip > 0) ? iip : np;
        //ip.print("ip");
    }
    return ip;
}
//' Basis matrix for B-spline of order 0 (step function) and higher
//'
//' This function is analogous but not equivalent to \code{splines:bs()} and \code{splines2::bSpline()}.
//' It is also several times faster.
//'
//' @param x Numeric vector, abscissa points
//' @param xk Numeric vector, knots
//' @param n Integer scalar, polynomial order (3 by default)
//' @return Numeric matrix, each column correspond to a B-spline calculated on x
//' @details
//'   For n==0, step function is defined as constant on each interval
//'   \code{[xk[i]; xk[i+1][}, i.e. closed on the left and open on the right
//'   except for the last interval which is closed on the right too.
//' @seealso {splines::bs()}, {splines2::bSpline()}
//' @examples
//'   x=seq(0, 5, length.out=101)
//'   # cubic basis matrix
//'   n=3
//'   m=bsc(x, xk=c(rep(0, n+1), 1:4, rep(5, n+1)), n=n)
//'   matplot(x, m, t="l")
//'   stopifnot(all.equal.numeric(c(m), c(splines::bs(x, knots = 1:4, degree = n, intercept = TRUE))))
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
mat bsc(const vec& x, const vec &xk, const size_t n=3) {
    // matrix  of b-splines(n) for each running xk set of 2+n points.
    // first column is first spline, etc.
    // spline values for x outside of xk are set to 0.
    auto nk=xk.n_elem;
    if (nk < 2+n)
        stop("Knot number must be >= n+2=%i (got %i)", n+2, nk);
    auto np=x.n_elem;
    vec idxk;
    uvec::fixed<2> len;
    mat res;
    umat ip=ipk(x, xk);
    mat dd=zeros<mat>(np, 2);
    for (size_t nc=0; nc <= n; nc++) {
        // nc is the current n, B-spline order
        auto nw=nk-1-nc;
        if (nc == 0) {
            res=zeros<mat>(np, nk-1);
            for (size_t iw=0; iw < nw; iw++) {
                size_t len=ip(1, iw)-ip(0, iw);
                if (len)
                    res(ip(0, iw), iw, size(len, 1)).fill(1.);
            }
            continue;
        }
        idxk=1./(xk.tail(nw+1)-xk.head(nw+1));
        len[0]=ip(1, nc-1)-ip(0, 0);
        if (len[0] > 0)
            dd(0, 0, size(len[0], 1))=(x(span(ip(0, 0), ip(1, nc-1)-1))-xk[0])*idxk[0];
        for (size_t iw=0; iw < nw; iw++) {
            auto iw1=iw+1;
            int i0=iw%2, i1=iw1%2;
            if (len[i0])
                res(ip(0, iw), iw, size(len[i0], 1)) %= dd(0, i0, size(len[i0], 1));
            len[i1]=ip(1, iw1+nc-1)-ip(0, iw1);
            if (len[i1] > 0) {
                dd(0, i1, size(len[i1], 1))=(x(span(ip(0, iw1), ip(1, iw1+nc-1)-1))-xk[iw1])*idxk[iw1];
                res(ip(0, iw1), iw, size(len[i1], 1)) += (1.-dd(0, i1, size(len[i1], 1)))%res(ip(0, iw1), iw1, size(len[i1], 1));
            }
        }
    }
    res.reshape(res.n_rows, nk-1-n);
    return res;
}
