#' bspline: build and use B-splines for interpolation and regression.
#' 
#' Build and use B-splines for interpolation and regression.
#'   In case of regression, equality constraints as well as monotonicity
#'   requirement can be imposed. Moreover, 
#'   knot positions (not only spline coefficients) can be part of 
#'   optimized parameters too. User is provided with 
#'   functions calculating spline values at arbitrary points. This 
#'   functions can be differentiated to obtain B-splines calculating 
#'   derivatives at any point. B-splines of this package can 
#'   simultaneously operate on a series of curves sharing the same set of 
#'   knots. 'bspline' is written with concern about computing 
#'   performance that's why the basis calculation is implemented in C++.
#'   The rest is implemented in R but without notable impact on computing speed.
#'
#' @section bspline functions: \describe{
#'  \item{"\code{bsc}:"}{ basis matrix (implemented in C++)}
#'  \item{"\code{bsp}:"}{ values of B-spline from its coefficients}
#'  \item{"\code{dbsp}:"}{ derivative of B-spline}
#'  \item{"\code{par2bsp}:"}{ build B-spline function from parameters }
#'  \item{"\code{bsppar}:"}{ retrieve B-spline parameters from its function}
#'  \item{"\code{smbsp}:"}{ build smoothing B-spline}
#'  \item{"\code{fitsmbsp}:"}{ build smoothing B-spline with optimized knot positions}
#'  \item{"\code{diffn}:"}{ finite differences}
#' }
#'
#' @docType package
#' @name bspline
#' @useDynLib bspline, .registration=TRUE
NULL

#' Calculate B-spline values from their coefficients qw and knots xk
#'
#' @param x Numeric vector, abscissa points at which B-splines should be calculated.
#'   They are supposed to be non decreasing.
#' @param xk Numeric vector, knots of the B-splines. They are supposed to be non decreasing.
#' @param qw Numeric vector or matrix, coefficients of B-splines. \code{NROW(qw)}
#'   must be equal to \code{length(xk)-n-1} where \code{n} is the next parameter
#' @param n Integer scalar, polynomial order of B-splines, by default cubic splines
#'   are calculated.
#' @return Numeric matrix (column number depends on qw dimensions), B-spline values on x.
#' @details This function does nothing else than calculate a dot-product between
#'   a B-spline basis matrix calculated by \code{bsc()} and coefficients \code{qw}.
#'   If qw is a matrix, each
#'   column corresponds to a separate set of coefficients.
#'   For x values falling outside of xk range, the B-splines values are set to 0.
#'   To get a function calculating spline values at arbitrary points from \code{xk}
#'   and \code{qw}, cf. \code{par2bsp()}.
#' @seealso [bsc], [par2bsp]
#' @export
bsp=function(x, xk, qw, n=3L) {
    stopifnot(NROW(qw) == length(xk)-n-1L)
    bsc(x, xk, n=n)%*%qw
}

#' Smoothing B-spline of order n >= 0
#'
#' @param x Numeric vector, abscissa points
#' @param y Numeric vector or matrix or data.frame, ordinate values to be smoothed
#'   (one set per column in case of matrix or data.frame)
#' @param n Integer scalar, polynomial order of B-splines (3 by default)
#' @param xki Numeric vector, strictly internal B-spline knots, i.e. lying strictly
#'   inside of \code{x} bounds. If NULL (by default), they are
#'   estimated with the help of \code{iknots()}. This vector is used as initial approximation
#'   during optimization process. Must be non decreasing if not NULL.
#' @param nki Integer scalar, internal knot number (1 by default). When
#'   nki==0, it corresponds to polynomial regression. If \code{xki}
#'   is not NULL, this parameter is ignored.
#' @param lieq List, equality constraints to respect by the smoothing spline,
#'   one list item per y column. By default (NULL), no constraint is imposed.
#'   Constraints are given as a 2-column matrix \code{(xe, ye)} where
#'   for each xe, an ye value is imposed. If a list item is NULL, no constraint
#'   is imposed on corresponding y column.
#' @param monotone Numeric scalar or vector, if \code{monotone > 0}, resulting B-spline
#'   weights must be increasing;
#'   if \code{monotone < 0}, B-spline weights must be decreasing; if \code{monotone == 0} (default), no
#'   constraint on monotonicity is imposed. If 'monotone' is a vector it
#'   must be of length \code{ncol(y)}, in which case each component indicates
#'   the constraint for corresponding column of y.
#' @param positive Numeric scalar, if \code{positive > 0}, resulting B-spline weights
#'   must be >= 0;
#'   if \code{positive < 0}, B-spline weights must be decreasing; if \code{positive == 0} (default), no
#'   constraint on positivity is imposed. If 'positive' is a vector it
#'   must be of length \code{ncol(y)}, in which case each component indicates
#'   the constraint for corresponding column of y.
#' @param mat Numeric matrix of basis vectors, if NULL it is recalculated
#'   by \code{bsc()}. If provided, it is the responsibility of the user
#'   to ensure that this matrix be adequate to xki vector.
#' @param estSD Logical scalar, if TRUE, indicates to calculate: SD of each y column, covariance
#'   matrix and SD of spline coefficients. All these values can be retrieved
#'   with bsppar() call (FALSE by default). These estimations are made under assumption
#'   that all y points have uncorrelated noise. Optional constraints are not taken
#'   into account of SD.
#' @param tol Numerical scalar, relative tolerance for small singular values
#'   that should be considered as 0 if \code{s[i] <= tol*s[1]}. This parameter
#'   is ignored if estSD=FALSE (1.e-10 by default).
#' @return Function, smoothing B-splines
#'   respecting optional constraints (generated by \code{par2bsp()}).
#' @details
#'   If constraints are set, we use \code{nlsic::lsie_ln()} to solve a
#'   least squares
#'   problem with equality constraints in least norm sens for each y column.
#'   Otherwise, \code{nlsic::ls_ln_svd()} is used for the whole y matrix.
#'   The solution of least squares problem is a vector of B-splines coefficients \code{qw},
#'   one vector per \code{y} column. These vectors are used to define B-spline function
#'   which is returned as the result.\cr\cr
#'   NB. When \code{nki >= length(x)-n-1} (be it from direct setting or calculated
#'   from \code{length(xki)}), it corresponds
#'   to spline interpolation, i.e. the resulting spline will pass
#'   exactly by (x,y) points (well, up to numerical precision).
#' @seealso \code{bsppar} for retrieving parameters of B-spline functions; \code{par2bsp}
#'   for generating B-spline function; \code{iknots} for estimation of knot positions
#' @examples
#'   x=seq(0, 1, length.out=11)
#'   y=sin(pi*x)+rnorm(x, sd=0.1)
#'   # constraint B-spline to be 0 at the interval ends
#'   fsm=smbsp(x, y, nki=1, lieq=list(rbind(c(0, 0), c(1, 0))))
#'   # check parameters of found B-splines
#'   bsppar(fsm)
#'   plot(x, y) # original "measurements"
#'   # fine grained x
#'   xfine=seq(0, 1, length.out=101)
#'   lines(xfine, fsm(xfine)) # fitted B-splines
#'   lines(xfine, sin(pi*xfine), col="blue") # original function
#'   # visualize knot positions
#'   xk=bsppar(fsm)$xk
#'   points(xk, fsm(xk), pch="x", col="red")
#' @importFrom nlsic lsie_ln ls_ln_svd
#' @importFrom stats sd
#' @aliases fitsmbsp smbsp
#' @export
smbsp=function(x, y, n=3L, xki=NULL, nki=1L, lieq=NULL, monotone=0, positive=0, mat=NULL, estSD=FALSE, tol=1.e-10) {
    y=as.matrix(y)
    nx=length(x)
    nc=ncol(y)
    stopifnot(nx == nrow(y))

    if (is.null(xki)) {
        # build default xki
        #dxk=diff(range(x))/(nki+1)
        #xki=x[1L]+seq_len(nki)*dxk
        xki=iknots(x, y, nki=nki, n=n)
    } else {
        nki=length(xki)
    }
    #stopifnot(nki > 0)
    if (nki > 0L) {
        stopifnot("xki is not inside of x range"=all((range(xki)-range(x))*c(1,-1) > 0))
        xk=c(rep(x[1L], n+1), xki, rep(x[nx], n+1))
    } else {
        xk=c(rep(x[1L], n+1), rep(x[nx], n+1))
    }
    
    nk=length(xk)
    nw=nk-n-1L
    if (is.null(mat))
        mat=bsc(x, xk, n=n)
    # equality constraints
    if (!is.null(lieq)) {
        stopifnot(length(lieq) == ncol(y))
        # todo: bezn() only on unique x from lieq
        liece=lapply(lieq, function(meq) {
            if (is.null(meq) || length(meq) == 0L) {
                NULL
            } else {
                list(e=bsc(meq[,1L], n=n, xk), ce=meq[,2L])
            }
        })
    } else {
        liece=NULL
    }
    monotone[monotone > 0]=1.
    monotone[monotone < 0]=-1.
    positive[positive > 0]=1.
    positive[positive < 0]=-1.
    if (length(monotone) == 1L && nc > 1L)
        monotone=rep(monotone, nc)
    stopifnot(length(monotone) == ncol(y))
    if (length(positive) == 1L && nc > 1L)
        positive=rep(positive, nc)
    stopifnot(length(positive) == ncol(y))
    if (any(monotone != 0)) {
        # B-spline must be monotonously in- or de-creasing
        i=seq_len(nw-1L)
        um=cbind(diag(x=-1., nrow=nw-1L), double(nw-1L))
        um[cbind(i, i+1)]=1.
        com=double(nw-1)
    } else {
        um=NULL
        com=NULL
    }
    if (any(positive != 0)) {
        # B-spline must be positive or negative
        up=diag(x=1., nrow=nw)
        cop=double(nw)
    } else {
        up=NULL
        cop=NULL
    }
    # u, co for all present combinations of monotone and positive
    uco=list()
    for (m in names(table(monotone))) {
        if (is.null(uco[[m]]))
            uco[[m]]=list()
        if (m == "-1") {
            u1=-um
            co1=com
        } else if (m == "1") {
            u1=um
            co1=com
        } else {
            u1=NULL
            co1=NULL
        }
        for (p in names(table(positive))) {
            if (is.null(uco[[m]][[p]]))
                uco[[m]][[p]]=list()
            if (p == "-1") {
                u2=-up
                co1=cop
            } else if (p == "1") {
                u2=up
                co2=cop
            } else {
                u2=NULL
                co2=NULL
            }
            uco[[m]][[p]]$u=rbind(u1, u2)
            uco[[m]][[p]]$co=c(co1, co2)
        }
    }
    monotone=as.character(monotone)
    positive=as.character(positive)
    if (!is.null(lieq) || any(monotone != 0) || any(positive != 0)) {
        #browser()
        qra=qr(mat, LAPACK=TRUE)
        qw=structure(vapply(seq_len(nc), function(ic) {
            m=monotone[ic]
            p=positive[ic]
            nlsic::lsie_ln(qra, y[,ic], u=uco[[m]][[p]]$u, co=uco[[m]][[p]]$co, e=liece[[ic]]$e, ce=liece[[ic]]$ce)
        }, double(nw)), dim=c(nw, nc), dimnames=list(NULL, colnames(y)))
        #browser()
    } else {
        qw=nlsic::ls_ln_svd(mat, y)
    }
    f=par2bsp(n, qw, xk)
    if (estSD) {
        #browser()
        # estimate sd of the fit
        sdy=apply(f(x) - y, 2, sd, na.rm=TRUE)
        # estimate SD of qw
        s=svd(mat)
        s$d=s$d[s$d/s$d[1L] >= tol]
        s$d=1./s$d
        r=length(s$d)
        if (r < ncol(mat)) {
            ir=seq_len(r)
            s$u=s$u[,ir,drop=FALSE]
            s$v=s$v[,ir,drop=FALSE]
        }
        covqw=tcrossprod(arrApply::arrApply(s$v, 2, "multv", v=s$d*s$d), s$v)
        sdqw=sqrt(diag(covqw))%o%sdy
        e=environment(f)
        e$sdqw=sdqw
        e$sdy=sdy
        e$covqw=covqw
    }
    f
}

#' Smoothing B-spline with optimized knot positions
#'
#' Optimize smoothing B-spline coefficients (smbsp) and knot positions (fitsmbsp)
#' such that residual squared sum is minimized for all y columns.
#' @details
#' Border and external knots are fixed, only strictly internal knots can move
#' during optimization. The optimization process is constrained to respect a minimal
#' distance between knots as well as to bound them to x range.
#' This is done to avoid knots getting unsorted during iterations and/or going
#' outside of a meaningful range.
#' @param control List, passed through to \code{nlsic()} call
#' @rdname smbsp
#' @aliases fitsmbsp smbsp
#' @examples
#'  # fit broken line with linear B-splines
#'  x1=seq(0, 1, length.out=11)
#'  x2=seq(1, 3, length.out=21)
#'  x3=seq(3, 4, length.out=11)
#'  y1=x1+rnorm(x1, sd=0.1)
#'  y2=-2+3*x2+rnorm(x2, sd=0.1)
#'  y3=4+x3+rnorm(x3, sd=0.1)
#'  x=c(x1, x2, x3)
#'  y=c(y1, y2, y3)
#'  plot(x, y)
#'  f=fitsmbsp(x, y, n=1, nki=2)
#'  lines(x, f(x))
#' @importFrom nlsic nlsic lsi_ln
#' @export
fitsmbsp=function(x, y, n=3L, xki=NULL, nki=1L, lieq=NULL, monotone=0, positive=0, control=list(), estSD=FALSE, tol=1.e-10) {
    np=length(x)
    y=as.matrix(y)
    stopifnot(nrow(y) == np)
    if (is.null(xki)) {
        xki=iknots(x, y, n=n, nki=nki)
    }
    nki=length(xki)
    if (nki == 0L)
        return(smbsp(x, y, n=n, xki=NULL, nki=0, lieq, monotone, positive))
    # inequalities
    # p[1] >= tail(x_1, 1)+epsx; p[i]+epsx <= p[i+1]; p[nki]+epsx <= x_n
    dx=diff(x)
    dx=dx[dx != 0.]
    stopifnot(dx > 0.)
    epsx=min(dx)*0.1
    u=matrix(0., nrow=nki+1L, ncol=nki)
    co=double(nrow(u))
    u[1L,1L]=1.; co[1L]=x[1L]+epsx
    i=seq_len(nki-1L)
    u[cbind(i+1L, i)]=-1.
    u[cbind(i+1L, i+1L)]=1.; co[i+1L]=epsx
    u[nki+1L, nki]=-1.; co[nki+1]=-x[np]+epsx
    x1=rep(x[1L], n+1L)
    x2=rep(x[np], n+1L)
    fresid=function(xki, cjac) {
        if (cjac) {
            li=bsc(x, c(x1, xki, x2), n=n, cjac=TRUE)
            f=smbsp(x, y, n=n, xki=xki, nki=0, lieq=lieq, monotone=monotone, positive=positive, mat=li$mat)
            res=f(x)-y
            jac=jacw(li$jac, bsppar(f)$qw)[,,(n+1)+seq_len(nki), drop=FALSE]
            dim(jac)=c(np*ncol(y), nki)
            list(res=res, jacobian=jac)
        } else {
            list(res=smbsp(x, y, n=n, xki=xki, nki=0, lieq=lieq, monotone=monotone, positive=positive)(x)-y)
        }
    }
    if (is.null(control$errx))
        control$errx=epsx*0.1
    fit=nlsic::nlsic(xki, fresid,
        u=u, co=co, control=control, flsi=nlsic::lsi_ln)
    if (fit$error != 0)
        stop(fit$mes)
    smbsp(x, y, n=n, xki=fit$par, nki=0, lieq=lieq, monotone=monotone, positive=positive, estSD=estSD, tol=tol)
}

#' Derivative of B-spline
#'
#' @param f Function, B-spline such as returned by \code{smbsp()} or \code{par2bsp()}
#' @param nderiv Integer scalar >= 0, order of derivative to calculate (1 by default)
#' @param same_xk Logical scalar, if TRUE, indicates to calculate derivative
#'   on the same knot grid as original function. In this case, coefficient number
#'   will be incremented by 2. Otherwise, extreme knots are
#'   removed on each side of the grid and coefficient number is maintained (FALSE by default).
#' @return Function calculating requested derivative
#' @examples
#'   x=seq(0., 1., length.out=11L)
#'   y=sin(2*pi*x)
#'   f=smbsp(x, y, nki=2L)
#'   d_f=dbsp(f)
#'   xf=seq(0., 1., length.out=101) # fine grid for plotting
#'   plot(xf, d_f(xf)) # derivative estimated by B-splines
#'   lines(xf, 2.*pi*cos(2*pi*xf), col="blue") # true derivative
#'   xk=bsppar(d_f)$xk
#'   points(xk, d_f(xk), pch="x", col="red") # knot positions
#' @export
dbsp=function(f, nderiv=1L, same_xk=FALSE) {
    #base::.Deprecated("Dbsp")
    stopifnot(nderiv >= 0L)
    if (nderiv == 0L)
        return(f)
    e=environment(f)
    qw=e$qw
    xk=e$xk
    n=e$n
    stopifnot(n >= nderiv)
    dxi=n/diff(xk, lag=n)
    dxi[!is.finite(dxi)]=0.
    if (same_xk) {
        qwn=arrApply::arrApply(rbind(0., qw, 0.), 1L, "diff")*dxi
        colnames(qwn)=colnames(qw)
        #browser()
        res=par2bsp(n-1L, qwn, xk)
        if (nderiv == 1L) {
            res
        } else {
            dbsp(res, nderiv=nderiv-1L)
        }
    } else {
        dxi=dxi[c(-1L, -length(dxi))]
        qwn=arrApply::arrApply(qw, 1L, "diff")*dxi
        colnames(qwn)=colnames(qw)
        #browser()
        res=par2bsp(n-1L, qwn, xk[c(-1L, -length(xk))])
        if (nderiv == 1L) {
            res
        } else {
            dbsp(res, nderiv=nderiv-1L, same_xk)
        }
    }
}

#' Differentiation matrix
#'
#' Calculate matrix for obtaining coefficients of first-derivative B-spline.
#' They can be calculated as \code{dqw=Md \%*\% qw}. Here, dqw are coefficients
#' of the first derivative,
#' Md is the matrix returned by this function, and qw are the coefficients
#' of differentiated B-spline.\cr
#'
#' @param nqw Integer scalar, row number of qw matrix (i.e. degree of freedom of a B-spline)
#' @param xk Numeric vector, knot positions
#' @param n Integer scalar, B-spline polynomial order
#' @param f Function from which previous parameters can be retrieved.
#'   If both f and any of previous parameters are given then explicitly
#'   set parameters take precedence over those retrieved from f.
#' @param same_xk Logical scalar, the same meaning as in \code{\link{dbsp}}
#' @return Numeric matrix of size \code{nqw-1 x nqw}
#' @export
dmat=function(nqw=NULL, xk=NULL, n=NULL, f=NULL, same_xk=FALSE) {
    if (is.null(f) && is.function(nqw))
        f=nqw
    if (!is.null(f)) {
        e=environment(f)
        if (is.null(nqw))
            nqw=NROW(e$qw)
        if (is.null(xk))
            xk=e$xk
        if (is.null(n))
            n=e$n
    }
    if (is.null(nqw))
        stop("nqw must be deduced from f or given explicitly")
    if (is.null(xk))
        stop("xk must be deduced from f or given explicitly")
    if (is.null(n))
        stop("n must be deduced from f or given explicitly")
    stopifnot(nqw == length(xk)-n-1L)
    n1=nqw+1L
    i=seq_len(nqw)
    mat=diag(nrow=n1, ncol=nqw)
    mat[cbind(i+1L, i)]=-1.
    dxi=n/diff(xk, lag=n)
    dxi[!is.finite(dxi)]=0.
    mat=mat*dxi
    if (!same_xk)
        mat=mat[c(-1L,-nrow(mat)),,drop=FALSE]
    mat
}

#' Indefinite integral of B-spline
#'
#' @param f Function, B-spline such as returned by \code{smbsp()} or \code{par2bsp()}
#' @param const Numeric scalar or vector of length \code{ncol(qw)} where
#'   qw is weight matrix of f. Defines starting value of weights for indefinite
#'   integral (0 by default).
#' @param nint Integer scalar >= 0, defines how many times to take integral (1 by default)
#' @return Function calculating requested integral
#' @details
#' If f is B-spline, then following identity is held: Dbsp(ibsp(f)) is identical to f.
#' Generally, it does not work in the other sens: ibsp(Dbsp(f)) is not f
#' but not very far. If we can get an appropriate constant C=f(min(x)) then
#' we can assert that ibsp(Dbsp(f), const=C) is the same as f.
#' @export
ibsp=function(f, const=0, nint=1L) {
    stopifnot(nint >= 0L)
    if (nint == 0L)
        return(f)
    e=environment(f)
    qw=e$qw
    stopifnot(length(const) == 1L || length(const) == ncol(qw))
    xk=e$xk
    n=e$n
    n1=n+1L
    xkn=xk[c(1L, seq_along(xk), length(xk))]
    qwn=rbind(0, arrApply::arrApply(qw*diff(xk, lag=n1)/n1, 1L, "cumsum"))
    if (length(const) == 1L) {
        qwn=qwn+const
    } else {
        qwn=arrApply::arrApply(qwn, 2, "addv", v=const)
    }
    colnames(qwn)=colnames(qw)
    #browser()
    res=par2bsp(n+1L, qwn, xkn)
    if (nint == 1L) {
        res
    } else {
        ibsp(res, nint=nint-1L, const)
    }
}
#' Retrieve parameters of B-splines
#'
#' @param f Function, B-splines such that returned by par3bsp(), smbsp(), ...
#' @return List having components: n - polynomial order, qw - coefficients, xk -
#'  knots
#' @export
bsppar=function(f) {
    as.list(environment(f))
}
#' Convert parameters to B-spline function
#'
#' @param n Integer scalar, polynomial order of B-splines
#' @param qw Numeric vector or matrix, coefficients of B-splines, one set per
#'  column in case of matrix
#' @param xk Numeric vector, knots
#' @param covqw Numeric Matrix, covariance matrix of qw (can be estimated in \code{\link{smbsp}}).
#' @param sdy Numeric vector, SD of each y column (can be estimated in \code{\link{smbsp}}).
#' @param sdqw Numeric Matrix, SD of qw thus having the same dimension
#'   as qw (can be estimated in \code{\link{smbsp}}).
#' @return Function, calculating B-splines at arbitrary points and having
#'  interface \code{f(x, select)} where \code{x} is a vector of abscissa points.
#'  Parameter \code{select} is passed to
#'  \code{qw[, select, drop=FALSE]} and can be missing. This function will return
#'  a matrix of size \code{length(x) x ncol(qw)} if \code{select} is missing. Elsewhere,
#'  a number of column will depend on \code{select} parameter. Column names in
#'  the result matrix will be inherited from \code{qw}.
#' @export
par2bsp=function(n, qw, xk, covqw=NULL, sdy=NULL, sdqw=NULL)
    local({
        stopifnot(NROW(qw) == length(xk)-n-1)
        n=n
        qw=as.matrix(qw)
        xk=xk
        covqw=covqw
        sdy=sdy
        sdqw=sdqw
        function(x, select, fsd=0.) {
            if (fsd == 0.) {
                if (base::missing(select)) bspline::bsp(x, xk, qw, n=n) else
                bspline::bsp(x, xk, qw[, select, drop=FALSE], n=n)
            } else {
                if (is.null(sdqw))
                    stop("B-spline is asked for fsd != 0 but sdqw is NULL")
                if (base::missing(select)) bspline::bsp(x, xk, qw+fsd*sdqw, n=n) else
                bspline::bsp(x, xk, qw[, select, drop=FALSE]+fsd*sdqw[, select, drop=FALSE], n=n)
            }
        }
    })
#' Finite differences
#'
#' Calculate dy/dx where x,y are first and the rest of columns in the entry matrix 'm'
#' @param m 2- or more-column numeric matrix
#' @param ndiff Integer scalar, order of finite difference (1 by default)
#' @return Numeric matrix, first column is midpoints of x, the second
#'  and following are dy/dx
#' @importFrom arrApply arrApply
#' @export
diffn=function(m, ndiff=1L) {
    stopifnot(ndiff >= 0L)
    if (ndiff == 0L)
        return(m)
    stopifnot(nrow(m) > ndiff)
    nr=nrow(m)
    d=arrApply::arrApply(m, 1L, "diff")
    res=cbind(m[-nr,1L]+0.5*d[, 1L, drop=FALSE], d[, -1L, drop=FALSE]/d[,1L])
    if (ndiff > 1L) {
        res=diffn(res, ndiff-1L)
    }
    colnames(res)=colnames(m)
    res
}
#' Estimate internal knot positions equalizing jumps in n-th derivative
#'
#' Normalized total variation of n-th finite differences is calculated for each column in
#' \code{y} then averaged. These averaged values are fitted by a linear spline to
#' find knot positions that equalize the jumps of n-th derivative.\cr
#' NB. This function is used internally in \code{(fit)smbsp()} and a priori
#' has no interest to be called directly by user.
#' @param x Numeric vector
#' @param y Numeric vector or matrix
#' @param nki Integer scalar, number of internal knots to estimate (1 by default)
#' @param n Integer scalar, polynomial order of B-spline (3 by default)
#' @return Numeric vector, estimated knot positions
#' @export
iknots=function(x, y, nki=1L, n=3L) {
    stopifnot(nki >= 0L)
    if (nki == 0L)
        return(double(0L))
    y=as.matrix(y)
    xy=cbind(x, y)
    dubx=duplicated(x)
    if (any(dubx)) {
        ix=split(seq_along(x), x)
        xy=t(vapply(ix, function(ii) if (length(ii) > 1L) colMeans(xy[ii,,drop=FALSE]) else xy[ii,], xy[1L,]))
    }
    dxyn=diffn(xy, n)
    dxy=diffn(dxyn, 1L)
    dy=dxy[,-1L, drop=FALSE]
    xtv=dxyn[,1L]
    tv=rbind(0, arrApply::arrApply(abs(dy)*diff(xtv), 1, "cumsum"))
    tv=arrApply::arrApply(tv, 2L, "multv", v=1./tv[nrow(tv),])
    tv[!is.finite(tv)]=0.
    tv=rowSums(tv)
    tv=tv/tv[length(tv)]
    ra=range(xtv)
    ftv=smbsp(xtv, tv, xki=xtv[1]+diff(ra)*seq(0, 1, length.out=12)[2:11], n=1L, lieq=list(cbind(ra,0:1)), monotone=1)
    par=bsppar(ftv)
    if (FALSE) {
        print(c("knots qw=", par$qw))
        print(c("knots dqw=", diff(par$qw)))
    }
    etv=seq(0, 1, length.out=nki+2L)[c(-1L, -(nki+2L))] # equalized tv
    qw=par$qw
    dq=diff(qw)
    dq[dq < 0]=0 # to force monotonicity despite round off errors
    qw=c(0., cumsum(dq))
    qw=qw/qw[length(qw)] # to make qw end up in 1
    ik=findInterval(etv, qw, rightmost.closed=TRUE)
    i=seq_along(etv)
    k=ik[i]
    x1=par$xk[k+1L]
    x2=par$xk[k+2L]
    y1=qw[k]
    y2=qw[k+1L]
    x1+(x2-x1)*(etv-y1)/(y2-y1)
}

#' nD B-curve governed by (x,y,...) control points.
#'
#' @param xy Real matrix of (x,y,...) coordinates, one control point per row.
#' @param n Integer scalar, polynomial order of B-spline (3 by default)
#' @return Function of one argument calculating B-curve. The argument is supposed
#'   to be in [0, 1] interval.
#' @details
#'   The curve will pass by the first and the last points in 'xy'. The tangents at the
#'   first and last points will coincide with the first and last segments of
#'   control points. Example of signature is inspired from this \href{https://www.r-bloggers.com/2023/03/little-useless-useful-r-functions-using-xspline-to-create-wacky-signatures/}{blog}.
#' @examples
#'   # simulate doctor's signature ;)
#'   set.seed(71);
#'   xy=matrix(rnorm(16), ncol=2)
#'   tp=seq(0,1,len=301)
#'   doc_signtr=bcurve(xy)
#'   plot(doc_signtr(tp), t="l", xaxt='n',  yaxt='n', ann=FALSE, frame.plot=FALSE,
#'       xlim=range(xy[,1]), ylim=range(xy[,2]))
#'   # see where control points are
#'   text(xy, labels=seq(nrow(xy)), col=rgb(0, 0, 0, 0.25))
#'   # join them by segments
#'   lines(bcurve(xy, n=1)(tp), col=rgb(0, 0, 1, 0.25))
#'   
#'   # randomly curved wire in 3D space
#'\dontrun{
#'   if (requireNamespace("rgl", quietly=TRUE)) {
#'      xyz=matrix(rnorm(24),ncol=3)
#'      tp=seq(0,1,len=201)
#'      curv3d=bcurve(xyz)
#'      rgl::plot3d(curv3d(tp), t="l", decorate=FALSE)
#'   }
#'}
#' @export
bcurve=function(xy, n=3) {
    nq=nrow(xy)
    stopifnot(nq >= 2)
    bspline::par2bsp(n, xy, c(rep(0, n), seq(0, 1, length.out=nq-n+1), rep(1, n)))
}
