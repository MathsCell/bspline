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
//' Basis matrix and knot Jacobian for B-spline of order 0 (step function) and higher
//'
//' This function is analogous but not equivalent to \code{splines:bs()} and \code{splines2::bSpline()}.
//' It is also several times faster.
//'
//' @param x Numeric vector, abscissa points
//' @param xk Numeric vector, knots
//' @param n Integer scalar, polynomial order (3 by default)
//' @param cjac Logical scalar, if \code{TRUE} makes to calculate Jacobian of basis
//'   vectors as function of knot positions (FALSE by default)
//' @return Numeric matrix (for cjac=FALSE), each column correspond to a
//'   B-spline calculated on x; or List (for cjac=TRUE) with components \describe{
//'      \item{mat}{basis matrix of dimension \code{nx x nw}, where nx is the length
//'        of x and \code{nw=nk-n-1} is the number of basis vectors}
//'      \item{jac}{array of dimension \code{nx x (n+2) x nw} where n+2
//'        is the number of support knots for each basis vector
//'    }
//' }
//' @details
//'   For n==0, step function is defined as constant on each interval
//'   \code{[xk[i]; xk[i+1][}, i.e. closed on the left and open on the right
//'   except for the last interval which is closed on the right too. The
//'   Jacobian for step function is considered 0 in every x point even if
//'   in points where x=xk, the derivative is not defined.\cr
//'   For n==1, Jacobian is discontinuous in such points so for
//'   these points we take the derivative from the right.
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
SEXP bsc(const vec& x, const vec &xk, const size_t n=3, const bool cjac=false) {
    // matrix  of b-splines(n) for each running xk set of 2+n points.
    // first column is first spline, etc.
    // spline values for x outside of xk are set to 0.
    auto nk=xk.n_elem;
    if (nk < 2+n)
        stop("Knot number must be >= n+2=%i (got %i)", n+2, nk);
    auto np=x.n_elem;
    vec idxk, vtmp(np);
    size_t len[2];
    mat res;
    umat ip=ipk(x, xk);
    mat dd=zeros<mat>(np, 2);
    cube jac, jacprev;
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
            if (cjac && n == 0)
                jac=zeros(np, 2, nk-1);
            continue;
        } else if (cjac && nc > 1) {
            jacprev=jac;
            //jacprev.print("jacprev");
        }
        if (cjac)
            jac=zeros<cube>(np, nc+2, nw);
        idxk=1./(xk.tail(nw+1)-xk.head(nw+1));
        len[0]=ip(1, nc-1)-ip(0, 0);
        if (len[0] > 0)
            dd(0, 0, size(len[0], 1))=(x(span(ip(0, 0), ip(1, nc-1)-1))-xk[0])*idxk[0];
        for (size_t iw=0; iw < nw; iw++) {
            auto iw1=iw+1;
            int i0=iw%2, i1=1-i0;
            len[i1]=ip(1, iw1+nc-1)-ip(0, iw1);
            if (len[i1] > 0)
                dd(0, i1, size(len[i1], 1))=(x(span(ip(0, iw1), ip(1, iw1+nc-1)-1))-xk[iw1])*idxk[iw1];
            // first jacobian as it needs lower order B-splines
            if (cjac) {
                if (nc == 1) {
                    if (len[i0] > 0) {
                        jac(ip(0, iw), 1, iw, size(len[i0], 1, 1))=-dd(0, i0, size(len[i0], 1))*idxk[iw];
                        jac(ip(0, iw), 0, iw, size(len[i0], 1, 1))=-idxk[iw]-jac(ip(0, iw), 1, iw, size(len[i0], 1, 1));
                    }
                    if (len[i1] > 0) {
                        jac(ip(0, iw1), 2, iw, size(len[i1], 1, 1))=dd(0, i1, size(len[i1], 1))*idxk[iw1];
                        jac(ip(0, iw1), 1, iw, size(len[i1], 1, 1))=idxk[iw1]-jac(ip(0, iw1), 2, iw, size(len[i1], 1, 1));
                    }
                } else {
                    // recursive jacobian definition
                    //Rcout << "iw=" << iw << std::endl;
                    if (len[i0] > 0) {
                        vtmp.head(len[i0])=dd(0, i0, size(len[i0], 1))*idxk[iw];
                        jac(ip(0, iw), 0, iw, size(len[i0], 1, 1))=(-idxk[iw]+vtmp.head(len[i0]))%res(ip(0, iw), iw, size(len[i0], 1));
                        //jac.print("j1");
                        jac(ip(0, iw), 0, iw, size(len[i0], nc+1, 1)) += 
                            dd(0, i0, size(len[i0], 1))%jacprev.slice(iw)(ip(0, iw), 0, size(len[i0], nc+1)).each_col();
                        //jac.print("j2");
                        jac(ip(0, iw), nc, iw, size(len[i0], 1, 1)) -= vtmp.head(len[i0])%res(ip(0, iw), iw, size(len[i0], 1));
                        //jac.print("j3");
                    }
                    //Rcout << "iw1,nc=" << iw1 <<"," << nc << "; jac size=" << size(jac) << std::endl;
                    if (len[i1] > 0) {
                        vtmp.head(len[i1])=dd(0, i1, size(len[i1], 1))*idxk[iw1];
                        jac(ip(0, iw1), 1, iw, size(len[i1], 1, 1)) += (idxk[iw1]-vtmp.head(len[i1]))%res(ip(0, iw1), iw1, size(len[i1], 1));
                        //jac.print("j4");
                        //Rcout << "ip(0, iw1)=" << ip(0, iw1) << "; size(len[i1], nc+1)=" << size(len[i1], nc+1) << std::endl;
                        jac(ip(0, iw1), 1, iw, size(len[i1], nc+1, 1)) += 
                            (1.-dd(0, i1, size(len[i1], 1)))%jacprev.slice(iw1)(ip(0, iw1), 0, size(len[i1], nc+1)).each_col();
                        //jac.print("j5");
                        jac(ip(0, iw1), nc+1, iw, size(len[i1], 1, 1)) += vtmp.head(len[i1])%res(ip(0, iw1), iw1, size(len[i1], 1));
                        //jac.print("j6");
                    }
                }
            }
            // then B-splines them-selves
            if (len[i0])
                res(ip(0, iw), iw, size(len[i0], 1)) %= dd(0, i0, size(len[i0], 1));
            if (len[i1] > 0)
                res(ip(0, iw1), iw, size(len[i1], 1)) += (1.-dd(0, i1, size(len[i1], 1)))%res(ip(0, iw1), iw1, size(len[i1], 1));
        }
    }
    res.reshape(res.n_rows, nk-1-n);
    if (cjac)
        return wrap(List::create(_["mat"]=res, _["jac"]=jac));
    else
        return wrap(res);
}
//' Knot Jacobian of B-spline with weights
//'
//' @param jac Numeric array, such as returned by \code{bsc(..., cjac=TRUE)}
//' @param qws Numeric matrix, each column is a set of weights forming a
//'   B-spline. If qws is a vector, it is coerced to 1-column matrix.
//' @return Numeric array of size \code{nx x ncol(qw) x nk}, where \code{nx=dim(jac)[1]}
//'   and nk is the number of knots \code{dim(jac)[3]+n+1} (n being polynomial order).
//' @export
// [[Rcpp::export]]
cube jacw(const cube& jac, const RObject& qws) {
    mat qw;
    if (qws.hasAttribute("dim")) {
        //Rcout << "mat\n";
        qw=as<mat>(qws);
    } else {
        //Rcout << "vec\n";
        NumericVector qwv(qws);
        qw=mat(&(*qwv.begin()), qwv.size(), 1, false);
    }
    //qw.print("qw");
    //qw=mat(&qwm(0,0), qwm.nrow(), qwm.ncol(), false);
    auto si=size(jac);
    //Rcout << "si=" << si << std::endl;
    auto np1=si[1]-1, nw=si[2], nk=nw+np1;
    if (nw != qw.n_rows)
        stop("nrow(qw) (=%d) must be equal to %d to be consistent with dim(jac)=(%d,%d,%d).", qw.n_rows, nw, si[0], si[1], si[2]);
    cube res=zeros<cube>(si[0], qw.n_cols, nk);
    for (size_t iw=0; iw < nw; iw++) {
        for (size_t i=0; i <= np1; i++) {
            auto ik=iw+i;
            res.slice(ik) += jac.slice(iw).col(i)*qw.row(iw);
            //char buf[256];
            //sprintf(buf, "%d, %d, %d", iw, i, ik);
            //res.print(buf);
        }
    }
    return res;
}
