library(RUnit)
library(bspline)
x010=c(seq(0, 1, by=0.2), seq(0.8, 0, by=-0.2))
x100=c(seq(1, 0, by=-0.2), 0, 0, 0, 0, 0)
x001=rev(x100)
test.bsc0=function() {
  x=seq(0, 2, len=11L)
  v=c(rep(1, 5), rep(0, 6))
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
test.smbsp.ex=function() {
  x=seq(0, 1, length.out=11)
  y=sin(pi*x)+rnorm(x, sd=0.1)
  # constraint B-spline to be 0 at the interval ends
  fsm=smbsp(x, y, nki=1, lieq=list(rbind(c(0, 0), c(1, 0))))
  checkEqualsNumeric(c(0, 0), fsm(0:1))
}
test.fitsmbsp.ex=function() {
  if (requireNamespace("numDeriv", quietly=TRUE)) {
    # fit broken line with linear B-splines
    x1=seq(0, 1, length.out=11)
    x2=seq(1, 3, length.out=21)
    x3=seq(3, 4, len=11)
    set.seed(7)
    y1=x1+rnorm(x1, sd=0.1)
    y2=-2+3*x2+rnorm(x2, sd=0.1)
    y3=4+x3+rnorm(x3, sd=0.1)
    x=c(x1, x2[-1], x3[-1])
    y=c(y1, y2[-1], y3[-1])
    f=fitsmbsp(x, y, n=1, nki=2)
    r=f(x)-y
    #print(sum(r*r))
    if (FALSE) {
      pdf("tmp.pdf")
      plot(x, y)
      lines(x, f(x))
      dev.off()
    }
    checkEqualsNumeric(0.4565935, sum(r*r))
  }
}
test.dbsp.ex=function() {
  x=seq(0., 1., length.out=11L)
  y=sin(2*pi*x)
  f=smbsp(x, y, nki=2L)
  d_f=dbsp(f)
  r=d_f(x)-2*pi*cos(2*pi*x)
  checkEqualsNumeric(1.494275547, sum(r*r))
}
