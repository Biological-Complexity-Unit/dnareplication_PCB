f <- function(xgrid, taugrid, origins, rates, v) {
  stopifnot(!is.unsorted(xgrid))
  stopifnot(!is.unsorted(taugrid))
  tau <- matrix(taugrid, nrow=length(taugrid), ncol=length(origins), byrow=FALSE)
  rates <- matrix(rates, nrow=length(taugrid), ncol=length(origins), byrow=TRUE)
  1 - exp(-sapply(xgrid, function(x) {
    dt <- matrix(abs(origins - x) / v, nrow=length(taugrid), ncol=length(origins), byrow=TRUE)
    rowSums(pmax((tau - dt), 0) * rates)
  }))
}

A <- function(x, origins, rates, lambda, v) {
  # Validate parameters
  stopifnot(is.numeric(x), is.numeric(origins) && is.numeric(rates) && is.numeric(lambda) && is.numeric(v))
  stopifnot(length(origins) == length(rates))
  stopifnot(all(is.finite(origins)) && all(!is.na(rates)) && (all(!is.na(lambda))) && (all(is.finite(v))) && (all(rates > 0)) && (all(v > 0)))
  stopifnot((length(lambda) == 1) && (length(v) == 1))

  # Compute weights
  N <- length(origins)
  W <- 1 + cumsum(rates) / lambda
  
  # Evaluate abundance at every position
  sapply(x, function(x) {
    # Compute absolute distances between current position x and origins
    d <- abs(origins - x)
    # Make sure origins are in order of increasing absolute distance
    if (is.unsorted(d)) {
      # Sort in order of increasing distance
      o <- order(d)
      d <- d[o]
      origins <<- origins[o]
      rates <<- rates[o]
      # Update weights
      W <<- 1 + cumsum(rates) / lambda
    }
    # Compute travel times from origins to x
    tau <- d / v
    # Compute effective replication times
    #  Tau_n = tau_n + sum[j=1..n] (tau_n - tau_j) * I_j / lambda
    #        = tau_n * W_n - sum[j=1..n] tau_j * I_j / lambda
    Tau <- tau * W - cumsum(tau * rates) / lambda
    # Compute abundance
    sum((1 / c(1, W[1:(N-1)]) - 1/W) * exp(-lambda*Tau))
  })
}
