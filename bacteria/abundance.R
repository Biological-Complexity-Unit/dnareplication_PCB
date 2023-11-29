library(statmod)
source("parmlapply.R")

# ************************************************************************
# **** Model v = v(x) is a harmonic oscillation **************************
# ************************************************************************

# Velocity is strictly a function of position here. This is similar, but
# not identical to the model of Bhat et al. where velocity is actually a
# function of distance traveled. See below for the v = v(d) model for an
# implementation of the model of Bhat et al.
#
# See Mathematica notebook elife_oscillation.nb and the TeX file
# harmonic_vx.tex, the code is a one-to-one translation of those.

vx.pp <- function(x, p)
  1 + p$delta * cos(x)

vx.p <- function(x, p)
  p$v0 * vx.pp(x * p$omega + p$phi, p)

vx <- function(x, p)
  vx.p(p$L/2 - abs(p$L/2 - x), p)

vx.Vpp1 <- function(x, p)
  2 * atan((1 - p$delta) * tan(x/2) / sqrt(1 - p$delta^2)) /
  sqrt(1 - p$delta^2)

vx.Vpp1OnePeriod <- function(p)
  2 * pi / sqrt(1 - p$delta^2)

vx.Ipp1 <- function(x, p)
  vx.Vpp1((x + pi) %% (2*pi) - pi, p) + vx.Vpp1OnePeriod(p) * ((x + pi) %/% (2*pi))

vx.Ip1 <- function(x, p) {
  if (p$omega > 0)
    vx.Ipp1(x*p$omega + p$phi, p) / (p$v0 * p$omega) -
    vx.Ipp1(p$phi, p) / (p$v0 * p$omega)
  else
    x / vx.p(0, p)
}

vx.I1 <- function(x, p)
  ifelse(x <= p$L/2, vx.Ip1(x, p), 2*vx.Ip1(p$L/2, p) - vx.Ip1(p$L-x, p))

vx.mu1 <- function(x, p)
  vx.I1(x, p)

vx.mu1.numeric <- Vectorize(function(x, p, dx=100)
  sum(dx/vx(seq(from=0, to=x, by=dx), p))
  , vectorize.args="x")

vx.mu2 <- function(x, p)
  vx.I1(p$L, p) - vx.I1(x, p)

vx.mu2.numeric <- Vectorize(function(x, p, dx=100)
  sum(dx/vx(seq(from=x, to=p$L, by=dx), p))
  , vectorize.args="x")

vx.Vpp3 <- function(x, p)
  (2 + p$delta^2) * atan((1 - p$delta) * tan(x/2) / sqrt(1 - p$delta^2)) /
  (1 - p$delta^2)^(5/2) +
  p$delta * (p$delta^2 - 3 * p$delta * cos(x) - 4) * sin(x) /
  (2 * (p$delta^2 - 1)^2 * (1 + p$delta * cos(x))^2)

vx.Vpp3OnePeriod <- function(p)
  2 * pi * (2 + p$delta^2) / (2 * (1 - p$delta^2)^(5/2))

vx.Ipp3 <- function(x, p)
  vx.Vpp3((x + pi) %% (2*pi) - pi, p) + vx.Vpp3OnePeriod(p) * ((x + pi) %/% (2*pi))

vx.Ip3 <- function(x, p) {
  if (p$omega > 0)
    vx.Ipp3(x*p$omega + p$phi, p) / (p$v0^3 * p$omega) -
    vx.Ipp3(p$phi, p) / (p$v0^3 * p$omega)
  else
    x / vx.p(0, p)^3
}

vx.I3 <- function(x, p)
  ifelse(x <= p$L/2, vx.Ip3(x, p), 2*vx.Ip3(p$L/2, p) - vx.Ip3(p$L-x, p))

vx.w1 <- function(x, p)
  vx.I3(x, p)

vx.w1.numeric <- Vectorize(function(x, p, dx=100)
  sum(dx/vx(seq(from=0, to=x, by=dx), p)^3)
  , vectorize.args="x")

vx.w2 <- function(x, p)
  vx.I3(p$L, p) - vx.I3(x, p)

vx.w2.numeric <- Vectorize(function(x, p, dx=100)
  sum(dx/vx(seq(from=x, to=p$L, by=dx), p)^3)
  , vectorize.args="x")

vx.f <- function(t, x, p)
  1 - pnorm(t, mean=vx.mu1(x, p), sd=sqrt(2*p$D*vx.w1(x, p)), lower.tail=FALSE) *
      pnorm(t, mean=vx.mu2(x, p), sd=sqrt(2*p$D*vx.w2(x, p)), lower.tail=FALSE)

vx.A <- Vectorize(function(x, p, force.numeric=FALSE) {
  # We have to compute the integral
  #    Int_0^inf lambda * exp(-lambda*t) * f(t, x) dt
  # where f is as defined above, i.e.
  #    f(t, x) = 1 - g_1(x,t) * g_2(x,t)
  # with
  #    g_i(x, t) = 1 - f_i(x, t) = S[(t - mu_i(x)) / 2*D*w_i(x)]
  # and S the survival function of the standard normal distribution.
  #
  R <- 5
  m <- c(vx.mu1(x, p), vx.mu2(x, p))
  sd <- sqrt(2 * p$D * c(vx.w1(x, p), vx.w2(x, p)))
  if ((m[1] + R*sd[1] < m[2] - R*sd[2]) && !force.numeric) {
    # mu_1 < mu_2, and the distributions don't significantly overlap
    # We treat this as a single replisome problem and use the analytic solution
    return(exp(-p$lambda * m[1] + p$lambda^2 * sd[1] ^ 2 / 2))
  } else if ((m[2] + R*sd[2] < m[1] - R*sd[1]) && !force.numeric) {
    # mu_2 < mu_1, and the distributions don't significantly overlap
    # We treat this as a single replisome problem and use the analytic solution
    return(exp(-p$lambda * m[2] + p$lambda^2 * sd[2] ^ 2 / 2))
  } else {
    # For small enough x we have that g_1(x, t) =~ 1, g_2(x, t) =~ 1 and
    # hence f(t, x) =~ 0, and for sufficiently large x we have that at least
    # one of g_1(x, t) =~ 0 or g_2(x, t) =~ 0 and hence f(t, x) =~ 1.
    #
    # We exploit that and split the integration domain into three regions,
    #    [0, T1], [T1, T2], [T2, inf]
    # such that
    #    f(t, x) =~ 0 on [0, T1] and
    #    f(t, x) =~ 1 on [T2, inf].
    # To obtain T1 and T2, we assume that S(-5) =~ 1 and S(5) =~ 0.
    T1 <- max(min(m - R*sd), 0)
    T2 <- max(min(m + R*sd), 0)
    # We then have that Int_D lambda * exp(-lambda*t) * f(t. x) dt is
    #    0 for D = [0, T1],
    #    exp(-lambda*T2) for D = [T2, infinity]
    # and only have to numerically integrate over [T2, T2]
    q2 <- integrate(function(t) p$lambda * exp(-p$lambda*t) * vx.f(t, x, p),
                    lower=T1, upper=T2)$value
    q3 <- exp(-p$lambda * T2)
    return (q2 + q3)
  }
}, vectorize.args=c("x"))

vx.harmonic.model <- list(
  name="v(x) harmonic",
  v1=function(x, t, p) vx(x, p),
  v2=function(x, t, p) -vx(x, p),
  D1=function(x, t, p) p$D,
  D2=function(x, t, p) p$D,
  reflect=FALSE,
  position="max(x(t))",
  f=vx.f,
  A=vx.A,
  scales=c(v0=500, delta=0.5, omega=1e-6, phi=1, D=1e6),
  upper.bounds=c(v0=Inf, delta=1, omega=Inf, phi=Inf, D=Inf),
  initial.value.grid=function(L, lambda) list(
    v0=10^seq(from=-2, to=1, by=0.33) * L * lambda,
    D=c(0, 10^seq(from=0, to=10, by=0.5)),
    delta=c(0, 0.5),
    omega=2^seq(from=-2, to=4, by=0.25) * 2 * pi / L,
    phi=seq(from=0, to=1-0.125, by=0.125) * 2 * pi
  )
)

vx.harmonic.reflect.backtrack.model <- list(
  name="v(x) harmonic w/ reflection & backtracking (simulation only)",
  v1=function(x, t, p) vx(x, p),
  v2=function(x, t, p) -vx(x, p),
  D1=function(x, t, p) p$D,
  D2=function(x, t, p) p$D,
  reflect=TRUE,
  position="x(t)"
)

# ************************************************************************
# **** Model v = v(d) is a harmonic oscillation **************************
# ************************************************************************
#
# The difference between this and v = v(x) is that d is the distance from
# the origin, *measured as experienced by each replisome*. Once a relisome
# crosses the terminus, it's d keeps increasing. The velocity here is thus
# not strictly speaking a function of genomic position, but rather of
# distance traveled.

vd <- function(d, p)
  p$v0 * (1 + p$delta * cos(d * p$omega + p$phi))

vd.hatG1 <- function(u, p)
  2 * atan((1 - p$delta) * tan(u/2) / sqrt(1 - p$delta^2)) /
  sqrt(1 - p$delta^2)

vd.barG1 <- function(p)
  2 * pi / sqrt(1 - p$delta^2)

vd.G1 <- function(u, p)
  vd.hatG1(u, p) + vd.barG1(p) * round(u / (2*pi))

vd.V1 <- function(d, p) {
  if (p$omega > 0)
    (1 / (p$v0 * p$omega)) * (vd.G1(d * p$omega + p$phi, p) - vd.G1(p$phi, p))
  else
    d / vd(0, p)
}

vd.V1.numeric <- Vectorize(function(d, p, dx=100)
  sum(dx/vd(seq(from=0, to=d, by=dx), p))
  , vectorize.args="d")

vd.hatG3 <- function(u, p)
  (2 + p$delta^2) * atan((1 - p$delta) * tan(u/2) / sqrt(1 - p$delta^2)) /
  (1 - p$delta^2)^(5/2) +
  p$delta * (p$delta^2 - 3 * p$delta * cos(u) - 4) * sin(u) /
  (2 * (p$delta^2 - 1)^2 * (1 + p$delta * cos(u))^2)

vd.barG3 <- function(p)
  2 * pi * (2 + p$delta^2) / (2 * (1 - p$delta^2)^(5/2))

vd.G3 <- function(u, p)
  vd.hatG3(u, p) + vd.barG3(p) * round(u / (2*pi))

vd.V3 <- function(d, p) {
  if (p$omega > 0)
    (1 / (p$v0^3 * p$omega)) * (vd.G3(d * p$omega + p$phi, p) - vd.G3(p$phi, p))
  else
    d / vd(0, p)^3
}

vd.V3.numeric <- Vectorize(function(d, p, dx=100)
  sum(dx/(vd(seq(from=0, to=d, by=dx), p) ^ 3))
  , vectorize.args="d")

vd.f <- function(t, x, p)
  1 - pnorm(t, mean=vd.V1(x, p), sd=sqrt(2*p$D*vd.V3(x, p)), lower.tail=FALSE) *
  pnorm(t, mean=vd.V1(p$L-x, p), sd=sqrt(2*p$D*vd.V3(p$L-x, p)), lower.tail=FALSE)

vd.A <- Vectorize(function(x, p, force.numeric=FALSE) {
  # We have to compute the integral
  #    Int_0^inf lambda * exp(-lambda*t) * f(t, x) dt
  # where f is as defined above, i.e.
  #    f(t, x) = 1 - g_1(x,t) * g_2(x,t)
  # with
  #    g_i(x, t) = 1 - f_i(x, t) = S[(t - mu_i(x)) / 2*D*w_i(x)]
  # and S the survival function of the standard normal distribution.
  #
  R <- 5
  m <- c(vd.V1(x, p), vd.V1(p$L-x, p))
  sd <- sqrt(2 * p$D * c(vd.V3(x, p), vd.V3(p$L-x, p)))
  if ((m[1] + R*sd[1] < m[2] - R*sd[2]) && !force.numeric) {
    # mu_1 < mu_2, and the distributions don't significantly overlap
    # We treat this as a single replisome problem and use the analytic solution
    return(exp(-p$lambda * m[1] + p$lambda^2 * sd[1] ^ 2 / 2))
  } else if ((m[2] + R*sd[2] < m[1] - R*sd[1]) && !force.numeric) {
    # mu_2 < mu_1, and the distributions don't significantly overlap
    # We treat this as a single replisome problem and use the analytic solution
    return(exp(-p$lambda * m[2] + p$lambda^2 * sd[2] ^ 2 / 2))
  } else {
    # For small enough x we have that g_1(x, t) =~ 1, g_2(x, t) =~ 1 and
    # hence f(t, x) =~ 0, and for sufficiently large x we have that at least
    # one of g_1(x, t) =~ 0 or g_2(x, t) =~ 0 and hence f(t, x) =~ 1.
    #
    # We exploit that and split the integration domain into three regions,
    #    [0, T1], [T1, T2], [T2, inf]
    # such that
    #    f(t, x) =~ 0 on [0, T1] and
    #    f(t, x) =~ 1 on [T2, inf].
    # To obtain T1 and T2, we assume that S(-5) =~ 1 and S(5) =~ 0.
    T1 <- max(min(m - R*sd), 0)
    T2 <- max(min(m + R*sd), 0)
    # We then have that Int_D lambda * exp(-lambda*t) * f(t. x) dt is
    #    0 for D = [0, T1],
    #    exp(-lambda*T2) for D = [T2, infinity]
    # and only have to numerically integrate over [T2, T2]
    q2 <- integrate(function(t) p$lambda * exp(-p$lambda*t) * vd.f(t, x, p),
                    lower=T1, upper=T2)$value
    q3 <- exp(-p$lambda * T2)
    return (q2 + q3)
  }
}, vectorize.args=c("x"))

vd.harmonic.model <- list(
  name="v(d) harmonic",
  v1=function(x, t, p) vd(x, p),
  v2=function(x, t, p) -vd(p$L - x, p),
  D1=function(x, t, p) p$D,
  D2=function(x, t, p) p$D,
  reflect=FALSE,
  position="max(x(t))",
  f=vd.f,
  A=vd.A,
  scales=c(v0=500, delta=0.5, omega=1e-6, phi=1, D=1e6),
  upper.bounds=c(v0=Inf, delta=1, omega=Inf, phi=Inf, D=Inf),
  initial.value.grid=function(L, lambda) list(
    v0=10^seq(from=-2, to=1, by=0.33) * L * lambda,
    D=c(0, 10^seq(from=0, to=10, by=0.5)),
    delta=c(0, 0.5),
    omega=2^seq(from=-2, to=4, by=0.25) * 2 * pi / L,
    phi=seq(from=0, to=1-0.125, by=0.125) * 2 * pi
  )
)

vd.harmonic.reflect.backtrack.model <- list(
  name="v(d) harmonic w/ reflection & backtracking (simulation only)",
  v1=function(x, t, p) vd(x, p),
  v2=function(x, t, p) -vd(p$L-x, p),
  D1=function(x, t, p) p$D,
  D2=function(x, t, p) p$D,
  reflect=TRUE,
  position="x(t)"
)

# ************************************************************************
# **** Model v(x,t) = v(t) a harmonic oscillation ************************
# ************************************************************************

# Mean and variance of the normal approximation for the time
# a replisome with harmonic time-dependent speed takes to reach
# position x.

vt.h <- function(t, p)
  1 + p$delta * cos(t * p$omega + p$phi)

vt.H <- function(t, p)
  if (p$omega != 0) {
    t + p$delta * (sin(t * p$omega + p$phi) - sin(p$phi)) / p$omega
  } else t

vt <- function(t, p)
  p$v0 * vt.h(t, p)

vt.D <- function(t, p)
  p$D * vt.h(t, p)

vt.Hmean <- function(x, p)
  x / p$v0

vt.Hvar <- function(x, p)
  2 * p$D * x / p$v0 ^ 3

vt.Hshape <- function(x, p)
  vt.Hmean(x, p)^3 / vt.Hvar(x, p)

vt.f <- function(t, x, p) {
  stopifnot(p$D >= 0)
  if (p$D > 0)
    1 - pinvgauss(vt.H(t, p), mean=vt.Hmean(x, p), shape=vt.Hshape(x, p), lower.tail=FALSE) *
        pinvgauss(vt.H(t, p), mean=vt.Hmean(p$L - x, p), shape=vt.Hshape(p$L - x, p), lower.tail=FALSE)
  else
    ifelse((x <= p$v0 * vt.H(t, p)) | (x >= p$L - p$v0 * vt.H(t, p)), 1, 0)
}

vt.A <- Vectorize(function(x, p, tol=1e-3) {
  # We want to numerically compute
  #   A(x) = int_0^inf Lambda * exp(-Lambda * t) f(t, x) dt
  # We do that by finding t0, t1 such that the integral over [0, t1] is
  # negligible and the integral over [t1, inf] is approximately exp(-Lambda * t1),
  # and then compute A(x) as
  #  A(x) = exp(-Lambda * t1) + int_t0^t1 Lambda * exp(-Lambda * t) f(t, x) dt

  # (1) Find max. t0 such that
  #   I0 = int_0^t0 Lambda * exp(-Lambda * t) f(t, x) dt <= tol/3.
  # We use that I0 <= B0(t0) = f(t0, x) * (1 - exp(-Lambda * t0)) and thus
  # search for the downward zero crossing of tol/2 - B0(t0).
  g0 <- function(t) tol/3 - vt.f(t, x, p) * (1 - exp(-p$lambda * t))
  r0 <- uniroot(g0, lower=0, upper=1, extendInt="downX")
  t0 <- r0$root
  
  # (2) Find min. t1 such that the error of replacing f <= 1 with constant 1 in
  #   I1 = int_t1^inf Lambda * exp(-Lambda * t) f(t, x) dt does not exceed tol/3,
  # or in other words that
  #   I1 = int_t1^inf Lambda * exp(-Lambda * t) (1 - f(t, x)) dt <= tol/3.
  # We use that I1 <= B1(t1) = (1 - f(t, x)) * exp(-Lambda * t1) and thus
  # search for the upward zero crossing of tol/2 - B1(t1).
  g1 <- function(t) tol/3 - (1 - vt.f(t, x, p)) * exp(-p$lambda * t)
  t1 <- if (g1(t0) < 0) {
    r1 <- uniroot(g1, lower=t0, upper=2*t0, extendInt="upX")
    r1$root
  } else {
    r1 <- NULL
    t0
  }
  #print(cbind(x=x, t0=t0, t1=t1, d0=r0$f.root, d1=r1$f.root))

  # (3) Integrate int_t0^t1 Lambda * exp(-Lambda * t) f(t, x) dt
  i12 <- if (t0 < t1)
    integrate(function(t) p$lambda * exp(-p$lambda*t) * vt.f(t, x, p),
              lower=t0, upper=t1, abs.tol=tol/3)$value
  else 0
  
  # (4) Add contribution of [t1, inf] and return
  return(i12 + exp(-p$lambda * t1))
}, vectorize.args=c("x"))

vt.harmonic.model <- list(
  name="v(t) harmonic",
  v1=function(x, t, p) vt(t, p),
  v2=function(x, t, p) -vt(t, p),
  D1=function(x, t, p) vt.D(t, p),
  D2=function(x, t, p) vt.D(t, p),
  reflect=FALSE,
  position="max(x(t))",
  f=vt.f,
  A=vt.A,
  A.terminus=function(eps, p) exp(- p$lambda * p$L
                                  / (2 * p$v0)
                                  + p$lambda * sqrt(p$D * p$L)
                                  / sqrt(pi * p$v0^3)
                                  + p$lambda ^ 2 * p$D * p$L * (pi - 1)
                                  / (2 * pi * p$v0^3)),
  A.terminus.curv=function(p) 2 * p$lambda / (sqrt(pi * p$D * p$L * p$v0)),
  Z.uncertainty=function(p) p$D * p$L / (2 * p$v0),
  scales=c(v0=500, delta=0.5, omega=2*pi/1800, phi=1, D=1e6),
  upper.bounds=c(v0=Inf, delta=1, omega=Inf, phi=Inf, D=Inf),
  initial.value.grid=function(L, lambda) list(
    v0=10^seq(from=-2, to=1, by=0.33) * L * lambda,
    D=c(0, 10^seq(from=0, to=10, by=0.5)),
    delta=c(0, 0.5),
    omega=2^seq(from=-2, to=4, by=0.25) * 2 * pi * lambda,
    phi=seq(from=0, to=1-0.125, by=0.125)*2*pi
  )
)

vt.harmonic.reflect.backtrack.model <- list(
  name="v(t) harmonic w/ reflection & backtracking (simulation only)",
  v1=function(x, t, p) vt(t, p),
  v2=function(x, t, p) -vt(t, p),
  D1=function(x, t, p) vt.D(t, p),
  D2=function(x, t, p) vt.D(t, p),
  reflect=TRUE,
  position="x(t)"
)

# ************************************************************************
# **** Computing f(t, x) and A(x) ****************************************
# ************************************************************************

f <- function(t, x, p, ...)
  p$model$f(t, x, p, ...)

A <- function(x, p, ...)
  p$model$A(x, p, ...)
  
# ************************************************************************
# **** Parameter Estimation **********************************************
# ************************************************************************

read.counts.deepak <- function(file, genome.length, genome.origin) {
  data <- data.table(read.table(file, header=FALSE, sep="\t", col.names=c("bp", "count", "sqrt")))
  data[, x := ifelse(bp >= genome.origin, bp - genome.origin, genome.length + bp - genome.origin) ]
  return(data)
}

normalize.counts <- function(x, sample, reference, dx=NA) {
  # Find window size dx, assuming a fixed window size except possibly for the last window
  if (is.na(dx)) {
    dx <- median(diff(x))
    cat(paste0("Detected window size: ", dx, "bp\n"))
  }
  # Compute normalized abundances
  y <- sample / reference
  y <- y / (dx * sum(y))
  stopifnot(is.finite(y))
  # Compute standard deviation using error propagation and Poissonian variances
  sd <- y * sqrt(1/sample + 1/reference)
  stopifnot(is.finite(sd) & (sd > 0))
  # Compute constant part of log-likelihood. For model predictions m, the full
  # log-likelihood assuming normality is then: logl.const + w * (a - m)^2
  logl.const <- -sum(log(sd)) -0.5*log(2*pi)*length(x)
  return(list(x=x, y=y, sd=sd, w=1/(2*sd^2), logl.const=logl.const))
}

evaluate.model <- function(p, data_or_x, y=NULL, w=NULL, logl.const=NA) {
  data <- NULL
  if (is.list(data_or_x) && is.null(y) && is.null(w)) {
    data <- data_or_x
    x <- data$x
    y <- data$y
    w <- data$w
    logl.const <- data$logl.const
  } else {
    x <- data_or_x
    if (is.null(y))
      stop("observed values y missing")
    if (is.null(w))
      w <- rep(1, length(y))
  }
  # Evaluate abundances
  a <- A(x, p)
  # Normalize to mach data
  ap <- sum(y) * a / sum(a)
  # Compute weighted squared distances
  sq <- sum(w * (y - ap)^2)
  # Compute log-likelihood
  logl <- logl.const - sq
  return(list(logl=logl, residual=sq, a=a, ap=ap, x=x))
}

fit.abundance <- function(x, y, w=NULL, vars, p0, scales.override=numeric(), rscale=NA,
                          method="Nelder-Mead", control=list(), ...) {
  # Determine variable scales.
  # If no explicit scales was specified, use the initial value if its finite and
  # non-zero, otherwise use the default scale for that parameters.
  scales.override <- scales.override[intersect(names(scales.override), names(vars))]
  scales <- p0$model$scales[names(vars)]
  vf <- is.finite(vars) & !is.na(vars) & (vars != 0)
  scales[names(vars[vf])] <- vars[vf]
  scales[names(scales.override)] <- scales.override
  stopifnot(is.finite(scales))
  stopifnot(length(scales) == length(vars))
  ub <- p0$model$upper.bounds[names(vars)] / scales

  # Generate goal function
  stopifnot(length(x) == length(y))
  s <- sum(y)
  if (is.null(w))
    w <- rep(1, length(x))
  g <- function(theta) {
    stopifnot(length(theta) == length(scales))
    if (any(theta < 0) || any(is.na(theta)) || any(theta >= ub))
      return(NA)

    # Merge constant parameters p0 with parameter vector theta
    p <- p0
    p[names(scales)] <- theta * scales
    # Evaluate abundances
    a <- A(x, p)
    # Normalize to mach data
    ap <- s * a / sum(a)
    # Compute weighted squared distances
    sq <- sum(w * (y - ap)^2)
    return(sq / rscale)
  }

  # Determine initial values
  # Use specified initial values or if NA use the parameter scale
  theta0 <- ifelse(is.na(vars), 1, vars / scales)

  # Determine residual scale by evaluating at the initial parameters
  if (is.na(rscale)) {
    rscale <- 1
    rscale <- g(theta0)
  }
  
  # Optimize
  r <- optim(par=theta0, fn=g, method=method, control=control, ...)
  stopifnot(names(r$par) == names(vars))
  
  # Return full parameter vector
  p <- p0
  p[names(scales)] <- r$par * scales
  attr(p, "optim") <- r
  attr(p, "residual") <- r$value * rscale
  return(p)
}

fit.grid <- function(x, y, w=NULL, grid, p0, mc.cores=1, ...) {
  mylapply <- if (mc.cores > 1) {
    function(...) mclapply(..., mc.cores=mc.cores)
  } else {
    lapply
  }
  
  rs <- mylapply(split(grid, 1:nrow(grid)), function(p) {
    fit.abundance(x, y, w, unlist(p), p0, ...)
  })
  r <- rs[[which.min(lapply(rs, function(p) attr(p, "residual")))]]
  if (attr(r, "optim")$convergence != 0)
    warning("Optimization did not converge!")
  return(r)
}


fit.abundance.stepwise <- function(x, y, w=NULL, L, lambda, model, mc.cores=1, reltol=1e-3) {
  report <- function(p, vars) {
    message("  Found ", paste0(vars, " = ", signif(unlist(p[vars]), 3), collapse=", "),
            " (residual ", signif(attr(p, "residual"), 6), ")")
  }

  # Grid of initial values (only D, omega, phi are used!)
  grid <- model$initial.value.grid(L, lambda)

  message("Fitting for genome length L = ", L, " and lambda = ", lambda)
  p0 <- list(L=L, lambda=lambda, model=model, D=0, v0=0, delta=0, omega=0, phi=0)

  message("Fitting v0 for D=0")
  p.vc.Dz <- fit.abundance(x, y, w, c(v0=NA), p0, method="BFGS")
  r1 <- attr(p.vc.Dz, "residual")
  if (attr(p.vc.Dz, "optim")$convergence != 0)
    warning("Optimization did not converge!")
  report(p.vc.Dz, "v0")

  message("Fitting v0, delta, omega and phi on a grid for D=0")
  p.vosc.Dz <- fit.grid(x, y, w, expand.grid(v0=p.vc.Dz$v0,
                                             delta=grid$delta,
                                             omega=grid$omega,
                                             phi=grid$phi),
                        p0, mc.cores=mc.cores, rscale=r1, control=list(reltol=reltol))
  report(p.vosc.Dz, c("v0", "delta", "omega", "phi"))

  message("Fitting v0 and D")
  p.vc.Dnz <- fit.grid(x, y, w, expand.grid(v0=p.vc.Dz$v0,
                                            D=grid$D),
                       p0, mc.cores=mc.cores, rscale=r1, control=list(reltol=reltol))
  report(p.vc.Dnz, c("v0", "D"))

  message("Fitting v0, delta, omega and phi on a grid for constant D")
  p.vosc.Dnz <- fit.grid(x, y, w, expand.grid(v0=p.vc.Dnz$v0,
                                              delta=grid$delta,
                                              omega=grid$omega,
                                              phi=grid$phi),
                         p.vc.Dnz, mc.cores=mc.cores, rscale=r1, control=list(reltol=reltol))
  report(p.vosc.Dnz, c("v0", "delta", "omega", "phi"))

  message("Re-fitting v0, delta, omega, phi and D")
  p.vosc.Dnz <- fit.grid(x, y, w, rbind(as.data.frame(p.vc.Dz[c("v0", "delta", "omega", "phi", "D")]),
                                        as.data.frame(p.vc.Dnz[c("v0", "delta", "omega", "phi", "D")]),
                                        as.data.frame(p.vosc.Dz[c("v0", "delta", "omega", "phi", "D")]),
                                        as.data.frame(p.vosc.Dnz[c("v0", "delta", "omega", "phi", "D")])),
                         p0, mc.cores=mc.cores, rscale=r1, control=list(maxit=1000, reltol=reltol))
  report(p.vosc.Dnz, c("v0", "delta", "omega", "phi", "D"))
  
  return(list(p.vc.Dz=p.vc.Dz, p.vosc.Dz=p.vosc.Dz, p.vc.Dnz=p.vc.Dnz, p.vosc.Dnz=p.vosc.Dnz))
}

