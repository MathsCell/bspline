library(RUnit)
library(bspline)
x010=c(seq(0, 1, by=0.2), seq(0.8, 0, by=-0.2))
x100=c(seq(1, 0, by=-0.2), 0, 0, 0, 0, 0)
x001=rev(x100)
nrm2=function(x) sum(x*x)
test.ipk=function() {
  x=seq(0, 2, len=11L)
  checkEqualsNumeric(ipk(x, 0:1), t(t(c(0, 6)))) # and not (0,11) as before
}
test.bsc0=function() {
  x=seq(0, 2, len=11L)
  v=c(rep(1, 5), rep(0, 6))
  checkEqualsNumeric(bsc(x, 0:1, n=0L), rev(1-v))
  checkEqualsNumeric(bsc(x, 0:2, n=0L), cbind(v, 1-v))
  checkEqualsNumeric(bsc(x, -1:3, n=0L), cbind(0,v, 1-v,0))
  checkEqualsNumeric(bsc(x, c(0, 0:2, 2), n=0L), cbind(0,v, 1-v,0))
  checkEqualsNumeric(bsc(x, c(0, 0:1, 1, 1:2, 2), n=0L), cbind(
    0,
    v,
    0.,
    0,
    1-v,
    0
  ))
}
test.bsc1=function() {
  x=seq(0, 2, len=11)
  checkEqualsNumeric(bsc(x, 0:2, n=1L), x010)
  checkEqualsNumeric(bsc(x, 0:2, n=1L, cjac=TRUE)$jac[,,1L], cbind(-x100, c(-x010[1:5], x010[6:11]), x001))
  checkEqualsNumeric(bsc(x, -1:3, n=1L), cbind(x100, x010, x001))
  checkEqualsNumeric(bsc(x, c(0, 0:2, 2), n=1L), cbind(x100, x010, x001))
  checkEqualsNumeric(bsc(x, c(0, 0:1, 1, 1:2, 2), n=1L), cbind(
    x100,
    c(seq(0, 0.8, by=0.2), rep(0, 6L)),
    0.,
    c(rep(0, 5L), seq(1, 0, by=-0.2)),
    x001
  ))
}
test.bsc2=function() {
  x=seq(0, 3, len=13L)
  x1=(x**2/2.)[1:5]
  x2=0.75-c(x[2]**2, 0., x[2]**2)
  checkEqualsNumeric(bsc(x, 0:3, n=2L), c(x1, x2, rev(x1)))
  j2=bsc(x, 0:3, n=2L, cjac=TRUE)$jac[,,1L]
  checkEqualsNumeric(j2[,1:2], -rev(j2[,3:4]))
  checkEqualsNumeric(bsc(x, c(1,1:3), n=2L), c(rep(0, 5), 0.40625, 0.625, 0.65625, rev(x1)))
}
test.smbsp3=function() {
  x=seq(0, 1, len=11L)
  y=1+2*x+3*x**2+4*x**3
  f=smbsp(x, y, nki=0L)
  checkEqualsNumeric(y, f(x))
  checkEqualsNumeric(1, f(0))
  checkEqualsNumeric(2, dbsp(f)(0))
  checkEqualsNumeric(6, dbsp(f, 2)(0))
  checkEqualsNumeric(24, dbsp(f, 3)(0))
}
test.jacw2=function() {
  x=seq(0, 4, len=13L)
  j2=bsc(x, 0:4, n=2L, cjac=TRUE)$jac
  checkEqualsNumeric(j2[1:10,1:2,1], -j2[10:1,4:3,1]) # must be anti-symmetric on appropriate interval
  checkException(jacw(j2, 1:4))
  jw=jacw(j2, c(1,1))
  checkEqualsNumeric(jw[1:13,1,1:5], -jw[13:1,1,5:1]) # must be anti-symmetric
}
test.ibsp=function() {
  x=seq(0, 1, len=11)
  f2=smbsp(x, x**2, nki=2)
  fi2=ibsp(f2)
  checkEqualsNumeric(fi2(x), x**3/3.)
  checkEqualsNumeric(bsppar(fi2)$qw, imat(f2)%*%bsppar(f2)$qw)
  f3=smbsp(x, x**3, nki=2)
  fi3=ibsp(f3)
  checkEqualsNumeric(fi3(x), x**4/4.)
  checkEqualsNumeric(bsppar(fi3)$qw, imat(f3)%*%bsppar(f3)$qw)
}
test.mnorm=function() {
  x=seq(0, 5, len=11)
  y=exp(-x)
  xpp=seq(0, 5, len=301);
  s=smbsp(x, y, n=3, nki=8, regular_grid=TRUE);
  checkEqualsNumeric(nrm2(s(xpp)-exp(-xpp)), 8.010032e-05, tolerance=1.e-7)
}
#----------------------------------
# doc examples
test.smbsp.ex=function() {
  x=seq(0, 1, length.out=11)
  y=sin(pi*x)+rnorm(x, sd=0.1)
  # constraint B-spline to be 0 at the interval ends
  fsm=smbsp(x, y, nki=1, lieq=list(rbind(c(0, 0), c(1, 0))))
  checkEqualsNumeric(c(0, 0), fsm(0:1))
}
test.fitsmbsp.ex=function() {
  # fit broken line with linear B-splines
  x1=seq(0, 1, length.out=11)
  x2=seq(1, 3, length.out=21)
  x3=seq(3, 4, length.out=11)
  set.seed(7)
  y1=x1+rnorm(x1, sd=0.1)
  y2=-2+3*x2+rnorm(x2, sd=0.1)
  y3=4+x3+rnorm(x3, sd=0.1)
  x=c(x1, x2, x3)
  y=c(y1, y2, y3)
  f=fitsmbsp(x, y, n=1, nki=2)
  r=f(x)-y
  #print(nrm2(r))
  if (FALSE) {
    pdf("~/tmp/tmp.pdf")
    plot(x, y)
    lines(x, f(x))
    dev.off()
  }
  #print(c("r2=", nrm2(r)))
  checkEqualsNumeric(0.3744, nrm2(r), tolerance=1.e-3)
}
test.dbsp.ex=function() {
  x=seq(0., 1., length.out=11L)
  y=sin(2*pi*x)
  f=smbsp(x, y, nki=2L)
  d_f=dbsp(f)
  r=d_f(x)-2*pi*cos(2*pi*x)
  checkEqualsNumeric(0.1992064, nrm2(r), tolerance=1.e-6)
}
test.bcurve.ex=function() {
  set.seed(71)
  xy=matrix(rnorm(16), ncol=2)
  #print(bcurve(xy)(0.5))
  checkEqualsNumeric(c(-0.03506473, -0.1508343), bcurve(xy)(0.5), tolerance=1.e-6)
}
test.pbsc.ex=function() {
  n=3
  x=seq(0, 5, length.out=101)
  xk=c(rep(0, n+1), 1:4, rep(5, n+1))
  # cubic polynomial coefficients
  coeffs=parr(xk, n)
  # basis matrix
  m=pbsc(x, xk, coeffs)
  checkEqualsNumeric(c(m), c(bsc(x, xk)))
}
test.parr.jac=function() {
  # compare analytical Jacobian of polynomial coefficients to central finite differences
  for (n in 1:4) {
    xk=c(0, 0.3, 1, 1.7, 2.5, 3.1, 4, 5.2, 6.1, 7, 8.3, 9)
    nk=length(xk)
    nw=nk-n-1
    r=parr(xk, n, cjac=TRUE)
    checkEqualsNumeric(parr(xk, n), r$coeffs)
    checkEquals(c(n+1, n+1, nw, n+2), dim(r$jac))
    h=1.e-6
    for (m in seq_len(nk)) {
      xp=xm=xk
      xp[m]=xk[m]+h
      xm[m]=xk[m]-h
      d=(parr(xp, n)-parr(xm, n))/(2*h)
      for (k in seq_len(nw)) {
        l=m-k+1 # xk[m]=xk[k+l-1]
        if (l < 1 || l > n+2)
          next
        checkEqualsNumeric(d[,,k], r$jac[,,k,l], tolerance=1.e-6)
      }
    }
  }
}
test.pjacw=function() {
  # must give the same result as jacw() applied to the Jacobian of bsc()
  x=seq(0, 5, length.out=101)
  for (n in 1:4) {
    xk=c(rep(0, n+1), 1:4, rep(5, n+1))
    nw=length(xk)-n-1
    pc=parr(xk, n, cjac=TRUE)
    set.seed(7)
    qws=matrix(rnorm(3*nw), ncol=3)
    ref=jacw(bsc(x, xk, n, cjac=TRUE)$jac, qws)
    checkEqualsNumeric(ref, pjacw(x, xk, pc$coeffs, pc$jac, qws), tolerance=1.e-8)
    # vector weights
    checkEqualsNumeric(jacw(bsc(x, xk, n, cjac=TRUE)$jac, qws[,1]),
      pjacw(x, xk, pc$coeffs, pc$jac, qws[,1]), tolerance=1.e-8)
  }
  # wrong number of weights
  pc=parr(xk, 3L, cjac=TRUE)
  checkException(pjacw(x, xk, pc$coeffs, pc$jac, 1:3))
}
test.dbsp.coeffs=function() {
  x=seq(0, 1, length.out=11)
  y=sin(2*pi*x)
  f=smbsp(x, y, nki=2L)
  e=environment(f)
  #e$coeffs=parr(e$xk, e$n)
  xf=seq(0.05, 0.95, length.out=41)
  for (same_xk in c(FALSE, TRUE)) {
    for (nd in 1:3) {
      d=dbsp(f, nd, same_xk)
      p=bsppar(d)
      checkEquals(e$n-nd, p$n)
      # coeffs deduced by dparr() must be the same as calculated by parr() for the order and knots of the derivative
      checkEquals(parr(p$xk, p$n), p$coeffs)
      # and give the same function values
      checkEqualsNumeric(d(xf), pbsc(xf, p$xk, p$coeffs) %*% p$qw)
    }
  }
  # coeffs of wrong dimension are refused by par2bsp
  checkException(par2bsp(3L, 1:5, c(0, 0, 0, 0, 1, 2, 2, 2, 2), coeffs=array(0, c(3, 3, 5))))
  # no coeffs in f => none in the derivative
  e$coeffs=NULL
  checkTrue(is.null(bsppar(dbsp(f))$coeffs))
}
test.dparr=function() {
  xk=c(0, 0, 0, 0, 1, 1, 2.5, 3, 5, 5, 5, 5)
  for (n in 1:4) {
    pc=parr(xk, n)
    checkEquals(parr(xk, n-1L), dparr(xk, pc, same_xk=TRUE))
    checkEquals(parr(xk[c(-1L, -length(xk))], n-1L), dparr(xk, pc, same_xk=FALSE))
  }
  checkException(dparr(xk, parr(xk, 3L)[,,-1L]))
  checkException(dparr(xk, parr(xk, 0L)))
}
