\name{dwavefunction}
\alias{dwavefunction}
\title{Wave Function Density}
\description{Evaluate the density of a wave function model}
\usage{
dwavefunction(x, w, log = FALSE, amplitude = FALSE)
}
\arguments{
  \item{x}{a numeric vector}
  \item{w}{a vector of coefficients from \code{\link{wavefunction_fit}}}
  \item{log}{if \code{TRUE}, returns the log density instead of the density}
  \item{amplitude}{if \code{TRUE}, returns the amplitude (or the log of the
    absolute value of the amplitude) instead of the density. The density is the
    squared amplitude, but the amplitude may be positive or negative.}
}
\value{a numeric vector of the same length as \code{x}}
\details{
  The elements of the returned vector \eqn{p} are (when \code{log} and
  \code{amplitude} are \code{FALSE}):
  \deqn{
    p_i = \left(\sum_{k=0}^K \frac{w_{k+1}}{(\sqrt{\pi}2^k k!)^{1/2}} H_k(x_i) \right)^2 e^{-x_i^2}
  }{%
    p[i] = (w[1] H[0](x) / e[1] + ... + w[K+1] H[k](x) / e[K+1])^2 * exp(-x^2)
    where e[k] = sqrt(sqrt(pi) * 2^k * k!)
  }
  Here, \eqn{K} is the maximum degree, equal to \code{length(w)-1}, and
  \eqn{H_k}{H[k]} is the Hermite polynomial of degree \eqn{k}. Note that
  \code{w}, being an R vector, is one-indexed, so \eqn{w_k}{w[k]} is associated
  with the Hermite polynomial of degree \eqn{k-1}.
}
\examples{
  x <- rnorm(100)
  w <- wavefunction_fit(x, degree = 6)
  p <- dwavefunction(x, w)
}
\seealso{
  Madeleine B. Thompson, \dQuote{Wave function representation of probability
  distributions,} 2017, \url{https://arxiv.org/abs/1712.07764}.
}
