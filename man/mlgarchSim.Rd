\name{mlgarchSim}
\alias{mlgarchSim}

\title{ Simulate from a multivariate log-GARCH(1,1) model }

\description{
Simulate the y series (typically a collection of financial returns or regression errors) from a log-GARCH model. Optionally, the conditional standard deviation and the standardised error, together with their logarithmic transformations, are also returned.
}
\usage{
mlgarchSim(n, constant = c(0,0), arch = diag(c(0.1, 0.05)),
  garch = diag(c(0.7, 0.8)), xreg = NULL,
  backcast.values = list(lnsigma2 = NULL, lnz2 = NULL, xreg = NULL),
  innovations = NULL, innovations.vcov = diag(rep(1,
  length(constant))), check.stability = TRUE, verbose = FALSE)
}

\arguments{
  \item{n}{integer, i.e. number of observations}
  \item{constant}{vector with the values of the intercepts in the log-volatility specification}
  \item{arch}{matrix with the arch coefficients}
  \item{garch}{matrix with the garch coefficients}
  \item{xreg}{a vector (of length n) or matrix (with rows n) with the values of the conditioning variables. The first column enters the first equation, the second enters the second equation, and so on}
  \item{backcast.values}{backcast values for the recursion (chosen automatically if NULL)}
  \item{check.stability}{logical. If TRUE (default), then the system is checked for stability}
  \item{innovations}{Either NULL (default) or a vector or matrix of length n with the standardised errors. If NULL, then the innovations are multivariate N(0,1) with correlations equal to zero}
  \item{innovations.vcov}{numeric matrix, the variance-covariance matrix of the standardised multivariate normal innovations. Only applicable if innovations = NULL}
  \item{verbose}{logical. If FALSE (default), then only the matrix with the y series is returned. If TRUE, then also additional information is returned}
}
\details{Empty}
\value{
A \code{\link{zoo}} matrix with n rows.
}
\references{
Sucarrat, Gronneberg and Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/} \cr
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{ \code{\link{lgarchSim}}, \code{\link{mlgarch}} and \code{\link{zoo}} }
\examples{
##simulate 1000 observations from a multivariate
##ccc-log-garch(1,1) w/default parameter values:
set.seed(123)
y <- mlgarchSim(1000)

##simulate the same series, but with more output:
set.seed(123)
y <- mlgarchSim(1000, verbose=TRUE)
head(y)

##plot the simulated values:
plot(y)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
