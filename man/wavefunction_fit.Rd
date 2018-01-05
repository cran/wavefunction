\name{wavefunction_fit}
\alias{wavefunction_fit}
\title{Fit Wave Function}
\description{Fit wave function coefficients from a sample}
\usage{
wavefunction_fit(x, degree)
}
\arguments{
  \item{x}{a sample from a distribution on the reals}
  \item{degree}{the Hermite polynomial degree to fit}
}
\value{a numeric vector of coefficients of length \code{degree+1}}
\details{
  Fits a Hermite wave function density of degree \code{degree}. The values will
  maximize the likelihood under the density specified under
  \code{\link{dwavefunction}}. A more accurate representation is obtained for a
  low \code{degree} if the sample is standardized to have mean zero and
  variance one-half. There are diminishing returns to \code{degree} greater
  than 20 or so due to floating point limitations.
}
\examples{
  x <- rt(100, df = 5)
  w <- wavefunction_fit(x, degree = 6)
}
\seealso{
  Madeleine B. Thompson, \dQuote{Wave function representation of probability
  distributions,} 2017, \url{https://arxiv.org/abs/1712.07764}.
}
