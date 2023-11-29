# Simulation of a growing population using its age distribution

# Samples replisome positions for given genome ages tau
sample.positions <- function(p, tau, dt, output="x(t)") {
  # Check whether to output raw x(tau) or the record max(x(t) | 0 <= t <= tau)
  output <- match.arg(output, c("x(t)", "max(x(t))"))
  # Number of parallel simulations
  n <- length(tau)
  # Initial positions, exponentially distributed positions
  X1 <- rep(0, n)
  X2 <- rep(p$L, n)
  M1 <- rep(0, n)
  M2 <- rep(p$L, n)
  # State of simulation, running or stopped
  R <- rep(TRUE, n)
  t <- 0
  tmax <- max(tau)
  while ((t < tmax) && any(R)) {
    t <- t + dt
    k <- sum(R)
    # Do time step dx = v(x) * dt + sqrt(2*D*dt) * eps
    X1[R] <- (X1[R] + p$model$v1(X1[R], rep(t, k), p)*dt +
                sqrt(2*p$model$D1(X1[R], rep(t, k), p)*dt) * rnorm(k))
    X2[R] <- (X2[R] + p$model$v2(X2[R], rep(t, k), p)*dt +
                sqrt(2*p$model$D2(X2[R], rep(t, k), p)*dt) * rnorm(k))
    # Reflecting boundary at 0 respectively L
    if (p$model$reflect) {
      X1[R] <- abs(X1[R])
      X2[R] <- p$L - abs(p$L - X2[R])
    }
    # Update maximal position
    M1[R] <- pmax(M1[R], X1[R])
    M2[R] <- pmin(M2[R], X2[R])
    # Stop if trajectories cross
    if (output == "x(t)") {
      s <- R & (X1 >= X2)
      if (any(s)) {
        m <- X1[s]*0.5 + X2[s]*0.5
        X1[s] <- m
        X2[s] <- m
      }
    } else if (output == "max(x(t))") {
      s <- R & (M1 >= M2)
      if (any(s)) {
        m <- M1[s]*0.5 + M2[s]*0.5
        M1[s] <- m
        M2[s] <- m
      }
    } else stop("unknown output variable")
    R[s] <- FALSE
    # Stop if time has reached tau
    s <- R & (t >= tau)
    R[s] <- FALSE
    #print(cbind(t=t, tau=tau[1:5], x1=X1[1:5], x2=X2[1:5]))
  }
  if (output == "x(t)")
    return(cbind(X1, X2))
  else if (output == "max(x(t))")
    return(cbind(M1, M2))
  else stop("unknown output variable")
}

# Determines the position-wise genome abundance in a population with ages tau
simulate.fxt <- function(p, tau, dx, dt, normalize=TRUE, cluster=CLUSTER, batchsize=1000) {
  mylapply <- if(is.null(cluster)) lapply else function(...) parLapply(cl=cluster, ...)
  if (missing(dx))
    dx <- 10e3
  if (missing(dt))
    dt <- 1e-3 / p$lambda
  # Number of bins and indices
  k <- ceiling(p$L / dx)
  i <- 1:k
  # Bin lower- and upper bounds
  j1 <- dx * (0:(k-1))
  j2 <- pmin(j1 + dx, p$L)
  # Split tau into batches
  tau <- as.numeric(tau)
  batches <- splitIndices(length(tau), ceiling(length(tau) / batchsize))
  # Sample replisome positions and add coverages
  a <- Reduce(`+`, mylapply(batches, function(b) {
    # Get batch
    t <- tau[b]
    # Sample replisome positions for batch
    x <- sample.positions(p, t, dt, output=p$model$position)
    # Add coverages for batch
    Reduce(`+`, lapply(1:nrow(x), function(i) {
      # Sampled replisome positions
      x1 <- x[i,1]
      x2 <- x[i,2]
      # Compute weights (0..1) for how much of each bin the replicated ranges
      # [0, x1] and [x2, L] cover
      w1 <- pmax(0, pmin(x1 - j1, dx)) / dx
      w2 <- pmax(0, pmin(j2 - x2, dx)) / dx
      # Return total abundance contribution of the replisome pair
      w1 + w2
    }))
  }))
  # Normalize abundance
  if (normalize)
    a <- a / (dx * sum(a))
  return(list(x=j1*0.5+j2*0.5, a=a, dx=dx))
}

simulate.uncertainty <- function(p, N, dt, cluster=CLUSTER) {
  mylapply <- if(is.null(cluster)) lapply else function(...) parLapply(cl=cluster, ...)
  if (missing(dt))
    dt <- 1e-3 / p$lambda
  batches <- splitIndices(N, length(cluster))
  do.call(mean, mylapply(batches, function(b) {
    # Get batch size
    n <- length(b)
    # Simulate
    x <- sample.positions(p, rep(Inf, n), dt, output=p$model$position)
    # Compute uncertainty
    stopifnot(all(x[, 1] == x[, 2]))
    mean((x[, 1] - p$L/2)^2)
  }))
}

# Determines the position-wise genome abundance in a growing populations
simulate.abundance <- function(p, dx, N, dt, normalize=TRUE, cluster=CLUSTER) {
  if (missing(N))
    N <- 1e4
  # Sample genome ages
  tau <- rexp(N, p$lambda)
  # Simulate
  simulate.fxt(p, tau, dx, dt, normalize, cluster)
}
