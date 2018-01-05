# Takes `y`, a stereographic projection of S^N into the plane in R^(N+1) at
# `x[1] = 0`, and inverts the projection. If `y` has length `N`, then the
# returned vector has length `N+1` and sum-of-squares 1.
#
# https://en.wikipedia.org/wiki/Stereographic_projection
unproject_stereographic <- function(y) {
  d <- sum(y ^ 2) + 1
  x <- c( (d - 2) / d, 2 * y / d)
  return(x)
}

# Given a vector of integers, `n`, returns a vector `e` of length `length(n)`.
# If `H` is the physicists' Hermite polynomial of degree `n[i]`, then `e[i]` is
# the log of the integral of (H(x))^2 exp(-x^2) over the real line.
log_normalizer <- function(n) 0.5 * log(pi) + n * log(2) + lfactorial(n)

# Returns `x` evaluated at each Hermite polynomial of degree 0 to degree
# `degree`. Returns a `length(x)` by `degree+1` matrix.
evaluate_H_at <- function(x, degree) {
  stopifnot(degree >= 0)
  stopifnot(round(degree) == degree)
  H <- matrix(nrow = length(x), ncol = (1 + degree))
  H[, 1] <- 1.0
  if (degree == 0) {
    return(H)
  }
  H[, 2] <- 2.0 * x
  if (degree == 1) {
    return(H)
  }
  for (d in 2:degree) {
    H[, d + 1] <- 2 * x * H[, d] - 2 * (d - 1) * H[, d - 1]
  }
  return(H)
}

# Returns `x` evaluated in the normalized Hermite basis of the specified
# degree. Returns a `length(x)` by `degree+1` matrix.
#
# Distinct from `evaluate_H_at` for ease of testing.
evaluate_in_basis <- function(x, degree) {
  norm <- exp(0.5 * log_normalizer(0:degree))
  M <- evaluate_H_at(x, degree) %*% diag(1 / norm)
  return(M)
}

assert_sum_sq_is_one <- function(w) {
  if (abs(sum(w ^ 2) - 1.0) > 0.0001) {
    stop("The sum of squares of 'w' must be one; is is ", sum(w ^ 2))
  }
}

# exported
wavefunction_fit <- function(x, degree) {
  M <- evaluate_in_basis(x, degree)

  # Returns the negative log likelihood of `x` for the transformed coefficients
  # `gamma`.
  neg_log_lik_fn <- function(gamma) {
    w <- unproject_stereographic(gamma)  # untransformed coefficients
    amp <- rowSums(M %*% w)              # amp[i] is the amplitude for x[i]
    stopifnot(length(amp) == length(x))
    loglik <- 2.0 * log(abs(amp)) - x ^ 2
    return(-sum(loglik))
  }

  # Fit the coefficients, `w` by minimizing the negative log likelihood of the
  # stereographic projection of the coefficients.
  gamma0 <- rep(0, degree)
  opt <- stats::optim(gamma0, neg_log_lik_fn, method = "L-BFGS")
  w <- unproject_stereographic(opt$par)

  # The density implied by `w` and `-w` is the same. Force the first
  # coefficient to be non-negative for uniqueness.
  if (w[1] < 0) {
    w <- w * -1
  }

  assert_sum_sq_is_one(w)
  return(w)
}

# exported
dwavefunction <- function(x, w, log = FALSE, amplitude = FALSE) {
  stopifnot(is.numeric(w))
  assert_sum_sq_is_one(w)
  degree <- length(w) - 1
  M <- evaluate_in_basis(x, degree)
  amp <- as.numeric(M %*% w)
  e <- ifelse(amplitude, 1.0, 2.0)  # exponent of `amp` to return
  if (log) {
    return(e * base::log(abs(amp)) - (e / 2) * x ^ 2)
  } else {
    return(amp ^ e * exp(-1 * (e / 2) * x ^ 2))
  }
}
