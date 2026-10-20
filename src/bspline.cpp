// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Intervals of points in knot intervals
//'
//' Find first and last+1 indexes iip s.t. x[iip] belongs to interval starting at xk[iik]
//'
//' @param x Numeric vector, abscissa points (must be non decreasing)
//' @param xk Numeric vector, knots (must be non decreasing)
//' @return Integer matrix of size \code{(2 x length(xk)-1)}. Indexes are 0-based
//' @export
// [[Rcpp::export]]
umat ipk(const vec& x, const vec& xk) {
    // find first and last+1 indexes iip s.t. x[iip] belongs to interval starting at xk[iik]
    // if not found, first==last+1, i.e. the length of interval is 0.
    auto nk=xk.n_elem, np=x.n_elem;
    umat ip(2, nk-1);
    double xkmax=xk[nk-1];
    auto ip_in=uvec(find(x <= xkmax, 1, "last"))[0];
    //x.print("x");
    //xk.print("xk");
    //Rcout << "&x=" << &x[0] << std::endl;
    //Rcout << "np=" << np << std::endl;
    for (size_t iik=0, iip=0; iik < nk-1; iik++) {
        double xki=xk[iik];
        //Rcout << "beging: iik, iip=" << iik << ", " << iip << std::endl;
        while ((iip < np && x[iip] <= xkmax) && (xki == xkmax || (iip == np-1 && iip > 0) ? x[iip] <= xki : x[iip] < xki))
            iip++; // find max iip s.t. x[iip]<[=] xk[iik] // "=" for the last iik
        //Rcout << "iik, iip=" << iik << ", " << iip << std::endl;
        ip(0, iik)=iip;
        if (iik > 0)
            ip(1, iik-1)=iip;
        if (iik == nk-2)
            ip(1, iik)=(iip == ip_in && iip > 0) ? iip : ip_in+1;
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
//' @seealso [splines::bs()], [splines2::bSpline()]
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
vec pmult(size_t n, vec v, double c0, double c1) {
    // v[0]*x^(n-1) + ... + v[n-1]
    vec res(n+1, fill::zeros);
    res.head(n)=c1*v.head(n);
    res.subvec(1, n)+=c0*v.head(n);
    return res;
}
//' Polynomial formulation of B-spline
//'
//' @param xk Numeric vector, knots
//' @param n Integer scalar, polynomial order (3 by default)
//' @param cjac Logical scalar, if \code{TRUE} makes to calculate Jacobian of polynomial
//'   coefficients as function of knot positions (FALSE by default)
//' @return Numeric 3D array (for cjac=FALSE), the first index runs through n+1 polynomial coefficients;
//'    the second -- through n+1 supporting intervals; and the last one through nk-n-1
//'    B-splines (here nk=length(xk)). Knot interval of length 0 will have corresponding
//'    coefficients set to 0; or List (for cjac=TRUE) with components \describe{
//'      \item{coeffs}{the 3D array described above}
//'      \item{jac}{4D array of dimension \code{(n+1) x (n+1) x nw x (n+2)} where
//'        \code{nw=nk-n-1} is the number of B-splines and n+2 is the number of support
//'        knots for each B-spline. The element \code{jac[i, j, k, l]} is the derivative of
//'        \code{coeffs[i, j, k]} with respect to the knot \code{xk[k+l-1]}
//'      }
//'    }
//' @details
//'  The polynomial coefficients on the j-th supporting interval of the k-th B-spline
//'  are those of a polynomial in the local variable \code{x-xk[k+j-1]}. The Jacobian is
//'  calculated for coefficients as such, i.e. keeping the local variable fixed. The
//'  dependence of the local variable on knot \code{xk[k+j-1]} itself is not included.
//'  Zero-length knot intervals as well as B-splines whose support knots are all identical
//'  have zero Jacobian (as they have zero coefficients).
//' @export
// [[Rcpp::export]]
SEXP parr(const vec &xk, const size_t n=3, const bool cjac=false) {
    auto nk=xk.n_elem;
    if (nk < 2+n)
        stop("Knot number must be >= n+2=%i (got %i)", n+2, nk);
    cube res(n+1, n+1, nk-1, fill::zeros); // will resize to nk-n-1 before return
    // jac[ik] is the derivative of res by xk[isp+ik], ik=0,...,n+1 (support knots of spline isp)
    std::vector<cube> jac;
    if (cjac)
        jac.assign(n+2, cube(n+1, n+1, nk-1, fill::zeros)); // will resize to nk-n-1 before return
    vec dk=diff(xk);
    vec dkn=dk;
    uvec inzk;
    double di0, di1, c00, c01;
    for (size_t nc=0; nc <= n; nc++) {
//Rcout << "nc=" << nc << std::endl;
//res.print("entry");
        // nc is the current n, B-spline order
        // We start with polynomials of degree 0 and apply Cox-de Boor recursive formula
        // by multiplying with polynomials of degree 1.
        inzk=find(dk);
        if (nc == 0) {
            for (auto k: inzk)
                res(0, 0, k)=1.;
            continue;
        }
        // multiply polynomials
        di0=dkn[0];
        if (di0 != 0.)
            di0=1./di0;
        for (size_t isp=0; isp < nk-nc-1; isp++) {
            // variable change on each interval: ksi in [0; dk[isp+isup]], ksi=x-xk[isp+isup]
//Rcout << "isp=" << isp << std::endl;
            di1=dkn[isp+1];
            if (di1 != 0.) {
                di1=1./di1;
                //c01=xk[isp+nc+1]*di1;
                c01=dkn[isp+1];
            }
            if (di0 != 0.)
                c00=0.;
                //c00=-xk[isp]*di0;
            for (size_t isup=0; isup <= nc; ++isup) {
//Rcout << "isup=" << isup << std::endl;
                if (di0 != 0. && isup < nc) {
                    if (isup > 0)
                        c00 += dk[isp+isup-1];
//Rcout << "c00=" << c00 << std::endl;
                    if (cjac) {
                        // d(M0*P)/dxk = M0*dP/dxk + dM0/dxk*P, where M0=c00*di0+di0*ksi
                        // depends on xk[isp] (ik=0), xk[isp+isup] (ik=isup) and xk[isp+nc] (ik=nc).
                        // It must be done before res is updated as old P is needed here
                        for (size_t ik=0; ik <= nc; ++ik)
                            jac[ik].slice(isp).col(isup).rows(0, nc)=pmult(nc, jac[ik].slice(isp).col(isup), c00*di0, di0);
                        jac[0].slice(isp).col(isup).rows(0, nc) += pmult(nc, res.slice(isp).col(isup), (c00*di0-1.)*di0, di0*di0);
                        jac[isup].slice(isp).col(isup).rows(0, nc) += pmult(nc, res.slice(isp).col(isup), di0, 0.);
                        jac[nc].slice(isp).col(isup).rows(0, nc) += pmult(nc, res.slice(isp).col(isup), -c00*di0*di0, -di0*di0);
                    }
                    res.slice(isp).col(isup).rows(0, nc)=pmult(nc, res.slice(isp).col(isup), c00*di0, di0);
//res.print("d0");
                }
                if (di1 != 0. && isup > 0) {
//Rcout << "c01=" << c01 << std::endl;
                    if (cjac) {
                        // d(M1*P)/dxk = M1*dP/dxk + dM1/dxk*P, where M1=c01*di1-di1*ksi
                        // depends on xk[isp+1] (ik=1), xk[isp+isup] (ik=isup) and xk[isp+nc+1] (ik=nc+1).
                        // P is taken from the spline isp+1, which is not yet updated.
                        for (size_t ik=1; ik <= nc+1; ++ik)
                            jac[ik].slice(isp).col(isup).rows(0, nc) += pmult(nc, jac[ik-1].slice(isp+1).col(isup-1), c01*di1, -di1);
                        jac[1].slice(isp).col(isup).rows(0, nc) += pmult(nc, res.slice(isp+1).col(isup-1), c01*di1*di1, -di1*di1);
                        jac[isup].slice(isp).col(isup).rows(0, nc) += pmult(nc, res.slice(isp+1).col(isup-1), -di1, 0.);
                        jac[nc+1].slice(isp).col(isup).rows(0, nc) += pmult(nc, res.slice(isp+1).col(isup-1), (1.-c01*di1)*di1, di1*di1);
                    }
                    res.slice(isp).col(isup).rows(0, nc) += pmult(nc, res.slice(isp+1).col(isup-1), c01*di1, -di1);
                    c01 -= dk[isp+isup];
//res.print("d1");
                }
            }
            di0=di1;
        }
        // increase interval size
        for (size_t ik=0; ik < nk-nc-1; ++ik)
            dkn[ik] += dk[ik+nc];
    }
    res.resize(n+1, n+1, nk-n-1);
    if (!cjac)
        return wrap(res);
    // pack jac into one 4D array (n+1, n+1, nk-n-1, n+2)
    size_t nw=nk-n-1;
    cube jacres(n+1, n+1, nw*(n+2));
    for (size_t ik=0; ik <= n+1; ++ik)
        jacres.slices(ik*nw, (ik+1)*nw-1)=jac[ik].slices(0, nw-1);
    NumericVector jr=wrap(jacres);
    jr.attr("dim")=IntegerVector::create(int(n+1), int(n+1), int(nw), int(n+2));
    return wrap(List::create(_["coeffs"]=res, _["jac"]=jr));
}

//' Polynomial B-spline Calculation of Basis Matrix
//'
//' @param x Numeric,vector, abscissa points
//' @param xk Numeric vector, knots
//' @param coeffs, Numeric 3D array, polynomial coefficients such as calculated by \code{\link{parr}}
//' @return Numeric matrix, basis vectors, one per column. Row number is \code{length(x)}.
//' @details
//'  Polynomials are calculated recursively by Cox-de Boor formula. However, it is not applied to
//'  final values but to polynomial coefficients. Multiplication by a linear functions gives
//'  a raise of polynomial degree by 1.\cr
//'  Polynomial coefficients stored in the first dimension of \code{coeffs} are used as in
//'  the following formula \code{p[1]*x^n + p[1]*x^(n-1) + ... + p[n+1]}. \cr
//'  Resulting matrix is the same as returned by \code{bsc(x, xk, n=dim(coeffs)[1]-1)}
//' @examples
//'   n=3
//'   x=seq(0, 5, length.out=101)
//'   xk=c(rep(0, n+1), 1:4, rep(5, n+1))
//'   # cubic polynomial coefficients
//'   coeffs=parr(xk)
//'   # basis matrix
//'   m=pbsc(x, xk, coeffs)
//'   matplot(x, m, t="l")
//'   stopifnot(all.equal.numeric(c(m), c(bsc(x, xk))))
//' @seealso \code{\link{bsc}}
//' @export
// [[Rcpp::export]]
mat pbsc(const vec& x, const vec& xk, const cube& coeffs) {
    size_t n=coeffs.n_rows-1;
    size_t nk=coeffs.n_slices+n+1;
    size_t nw=nk-n-1;
    size_t first, last;
    if (nk != xk.n_elem)
        stop("pbsc: the length(xk)=%d must be equal to dim(coeffs)[3]+dim(coeffs)[1]=%d+%d=%d", xk.n_elem, nw, n+1, nk);
    mat res(x.n_elem, nw, fill::zeros);
    umat ip=ipk(x, xk);
//ip.print("ip");
    for (size_t isp=0; isp < nw; ++isp) {
//Rcout << "isp=" << isp << std::endl;
        for (size_t isup=0; isup <= n; ++isup) {
//Rcout << "\tisup=" << isup << std::endl;
            first=ip(0, isp+isup);
            last=ip(1, isp+isup);
            if (last > first) {
                last -= 1;
                res.col(isp).subvec(first, last) += polyval(coeffs.slice(isp).col(isup), x.subvec(first, last)-xk[isp+isup]);
//res.print("res");
            }
        }
    }
    return res;
}

//' Knot Jacobian of polynomial B-spline with weights
//'
//' Analogous to \code{\link{jacw}} but based on polynomial coefficients of B-splines
//' (see \code{\link{parr}} and \code{\link{pbsc}}).
//'
//' @param x Numeric vector, abscissa points
//' @param xk Numeric vector, knots
//' @param coeffs Numeric 3D array, polynomial coefficients such as calculated by \code{\link{parr}}
//' @param jac Numeric 4D array, Jacobian of \code{coeffs} by knots such as calculated by
//'   \code{parr(xk, n, cjac=TRUE)$jac}
//' @param qws Numeric matrix, each column is a set of weights forming a
//'   B-spline. If qws is a vector, it is coerced to 1-column matrix.
//' @return Numeric array of size \code{nx x ncol(qws) x nk}, where \code{nx=length(x)}
//'   and \code{nk=length(xk)}. The element \code{[i, j, k]} is the derivative of the
//'   j-th weighted B-spline at \code{x[i]} with respect to \code{xk[k]}. It has the same layout
//'   as the result of \code{jacw(bsc(x, xk, n, cjac=TRUE)$jac, qws)}.
//' @details
//'  On the j-th supporting interval of the k-th B-spline, the polynomial is taken in the local
//'  variable \code{x-xk[k+j-1]}, so the derivative by \code{xk[k+j-1]} gets an additional term
//'  coming from the derivative of this polynomial with respect to its variable.
//'  For n=0 the Jacobian is 0 and for n=1 it is taken from the right in points where
//'  it is discontinuous (as in \code{\link{bsc}}).
//' @examples
//'   n=3
//'   x=seq(0, 5, length.out=101)
//'   xk=c(rep(0, n+1), 1:4, rep(5, n+1))
//'   qws=cbind(1:8, 8:1)
//'   pc=parr(xk, n, cjac=TRUE)
//'   jw=pjacw(x, xk, pc$coeffs, pc$jac, qws)
//'   stopifnot(all.equal(jw, jacw(bsc(x, xk, n, cjac=TRUE)$jac, qws)))
//' @seealso \code{\link{jacw}}, \code{\link{parr}}, \code{\link{pbsc}}
//' @export
// [[Rcpp::export]]
cube pjacw(const vec& x, const vec& xk, const cube& coeffs, const NumericVector& jac, const RObject& qws) {
    size_t n=coeffs.n_rows-1;
    size_t nw=coeffs.n_slices;
    size_t nk=nw+n+1;
    size_t first, last;
    if (nk != xk.n_elem)
        stop("pjacw: the length(xk)=%d must be equal to dim(coeffs)[3]+dim(coeffs)[1]=%d+%d=%d", xk.n_elem, nw, n+1, nk);
    if (size_t(jac.size()) != (n+1)*(n+1)*nw*(n+2))
        stop("pjacw: length(jac)=%d must be equal to (n+1)*(n+1)*nw*(n+2)=%d", jac.size(), (n+1)*(n+1)*nw*(n+2));
    mat qw;
    if (qws.hasAttribute("dim")) {
        qw=as<mat>(qws);
    } else {
        NumericVector qwv(qws);
        qw=mat(&(*qwv.begin()), qwv.size(), 1, false);
    }
    if (nw != qw.n_rows)
        stop("nrow(qw) (=%d) must be equal to %d to be consistent with dim(coeffs)=(%d,%d,%d).", qw.n_rows, nw, n+1, n+1, nw);
    // jac(i, j, isp+ik*nw) is the derivative of coeffs(i, j, isp) by xk[isp+ik]
    cube jc(const_cast<double*>(&jac[0]), n+1, n+1, nw*(n+2), false); // used as read-only
    cube res(x.n_elem, qw.n_cols, nk, fill::zeros);
    umat ip=ipk(x, xk);
    vec dc(n > 0 ? n : 1, fill::zeros);
    for (size_t isp=0; isp < nw; ++isp) {
        for (size_t isup=0; isup <= n; ++isup) {
            first=ip(0, isp+isup);
            last=ip(1, isp+isup);
            if (last > first) {
                last -= 1;
                vec ksi=x.subvec(first, last)-xk[isp+isup];
                // derivative of polynomial by its variable ksi, it is needed for
                // the dependence of ksi on xk[isp+isup]
                vec pd(ksi.n_elem, fill::zeros);
                if (n > 0) {
                    dc=coeffs.slice(isp).col(isup).head(n)%regspace(n, -1, 1);
                    pd=polyval(dc, ksi);
                }
                for (size_t ik=0; ik <= n+1; ++ik) {
                    vec v=polyval(jc.slice(isp+ik*nw).col(isup), ksi);
                    if (ik == isup)
                        v -= pd;
                    res.slice(isp+ik).rows(first, last) += v*qw.row(isp);
                }
            }
        }
    }
    return res;
}
//' Polynomial coefficients of B-splines of order n-1 from those of order n
//'
//' Calculate polynomial coefficients for B-splines of order n-1 on the knot grid
//' of the derivative of an order n B-spline (see \code{\link{dbsp}}), starting
//' from the polynomial coefficients of order n B-splines.
//'
//' @param xk Numeric vector, knots of order n B-splines
//' @param coeffs Numeric 3D array, polynomial coefficients of order n B-splines
//'   such as calculated by \code{\link{parr}(xk, n)}
//' @param same_xk Logical scalar, if \code{TRUE}, order n-1 B-splines are defined
//'   on the same knots \code{xk}, otherwise on \code{xk[-c(1, length(xk))]}
//'   (FALSE by default)
//' @return Numeric 3D array of dimension \code{n x n x nw1}, where \code{nw1} is
//'   \code{dim(coeffs)[3]+1} for \code{same_xk=TRUE} and \code{dim(coeffs)[3]-1} otherwise.
//'   It is the same as \code{parr(xk, n-1)} (\code{same_xk=TRUE}) or
//'   \code{parr(xk[-c(1, length(xk))], n-1)} (\code{same_xk=FALSE}) up to numerical precision.
//' @details
//'  Let B be a B-spline of order n on knots t[0], ..., t[n+1], and U, V be B-splines
//'  of order n-1 on knots t[0], ..., t[n] and t[1], ..., t[n+1] respectively.
//'  The Cox-de Boor formula and the derivative of B-spline
//'  give a linear system on U and V. Its solution is
//'  \code{V=(B-(x-t[0])*B'/n)*(t[n+1]-t[1])/(t[n+1]-t[0])} and
//'  \code{U=(B+(t[n+1]-x)*B'/n)*(t[n]-t[0])/(t[n+1]-t[0])}. These are polynomial
//'  operations applied on each supporting interval of B.
//' @seealso \code{\link{parr}}, \code{\link{dbsp}}
//' @export
// [[Rcpp::export]]
cube dparr(const vec &xk, const cube &coeffs, const bool same_xk=false) {
    size_t n=coeffs.n_rows-1;
    size_t nw=coeffs.n_slices;
    size_t nk=nw+n+1;
    if (n < 1)
        stop("dparr: polynomial order must be >= 1 (got %i)", n);
    if (nk != xk.n_elem)
        stop("dparr: the length(xk)=%d must be equal to dim(coeffs)[3]+dim(coeffs)[1]=%d+%d=%d", xk.n_elem, nw, n+1, nk);
    size_t off=same_xk ? 1 : 0;
    size_t nwn=same_xk ? nw+1 : nw-1; // number of order n-1 B-splines
    cube res(n, n, nwn, fill::zeros);
    vec p, dc, q;
    double den, fac, c;
    for (size_t isp=0; isp < nw; isp++) {
        // B-spline isp of order n is on knots xk[isp], ..., xk[isp+n+1]
        den=xk[isp+n+1]-xk[isp];
        if (den == 0.)
            continue; // zero polynomial
        // V=B-spline isp+1 of order n-1: its supporting interval isup is isup+1 of B-spline isp
        // (if same_xk is false, it is the B-spline isp as the first knot is removed and
        // the last one gets no V as it would need the removed last knot)
        fac=(xk[isp+n+1]-xk[isp+1])/den;
        for (size_t isup=0; isup < n && isp+off < nwn; isup++) {
            p=coeffs.slice(isp).col(isup+1);
            dc=p.head(n)%regspace(n, -1, 1)/n; // B'/n
            // local variable ksi=x-xk[isp+isup+1], so x-xk[isp]=ksi+c
            c=xk[isp+isup+1]-xk[isp];
            q=p-pmult(n, dc, c, 1.);
            // the leading coefficient is 0 up to numerical precision
            res.slice(isp+off).col(isup)=fac*q.tail(n);
        }
        if (same_xk && isp == 0) {
            // U=B-spline 0 of order n-1, its supporting interval isup is isup of B-spline 0
            fac=(xk[n]-xk[0])/den;
            for (size_t isup=0; isup < n; isup++) {
                p=coeffs.slice(0).col(isup);
                dc=p.head(n)%regspace(n, -1, 1)/n;
                // xk[n+1]-x=c-ksi, ksi=x-xk[isup]
                c=xk[n+1]-xk[isup];
                q=p+pmult(n, dc, c, -1.);
                res.slice(0).col(isup)=fac*q.tail(n);
            }
        }
    }
    return res;
}
