\name{lgarchSim}
\alias{lgarchSim}

\title{ Simulate from a log-GARCH model }

\description{
Simulate the y series (typically interpreted as a financial return or the error in a regression) from a log-GARCH model. Optionally, the conditional standard deviation and the standardised error (z), together with their logarithmic transformations, are also returned.
}
\usage{
lgarchSim(n, constant = 0, arch = 0.05, garch = 0.9, xreg = NULL,
  backcast.values = list(lnsigma2 = NULL, lnz2 = NULL, xreg = NULL),
  check.stability = TRUE, innovations = NULL, verbose = FALSE)
}

\arguments{
  \item{n}{integer, length of y (i.e. no of observations)}
  \item{constant}{the value of the intercept in the log-volatility specification}
  \item{arch}{a vector with the arch coefficients}
  \item{garch}{a vector with the garch coefficients}
  \item{xreg}{a vector or a matrix with the conditioning variables}
  \item{backcast.values}{backcast values for the recursion (chosen automatically if NULL)}
  \item{check.stability}{logical. If TRUE (default) then the roots of arch+garch are checked for stability}
  \item{innovations}{Etiher NULL (default) or a vector of length n with the standardised errors (i.e. z). If NULL, then the innovations are N(0,1)}
  \item{verbose}{logical. If FALSE (default), then only the vector y is returned. If TRUE, then a matrix with all the output is returned}
}
\details{Empty}
\value{
A \code{\link{zoo}} vector of length n if verbose = FALSE, or a \code{\link{zoo}} matrix with n rows if verbose = TRUE.
}
\references{
Sucarrat, Gronneberg and Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/} \cr
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{ \code{\link{mlgarchSim}}, \code{\link{lgarch}}, \code{\link{zoo}} }
\examples{
##simulate 500 observations w/default parameter values:
set.seed(123)
y <- lgarchSim(500)

##simulate the same series, but with more output:
set.seed(123)
y <- lgarchSim(500, verbose=TRUE)
head(y)

##plot the simulated values:
plot(y)

##simulate from an integrated log-GARCH:
y <- lgarchSim(500, arch=0.1, garch=0.9)

##simulate w/conditioning variable:
x <- rnorm(500)
y <- lgarchSim(500, xreg=0.05*x)

##simulate from a log-GARCH with a simple form of leverage:
z <- rnorm(500)
zneg <- as.numeric(z < 0)
zneglagged <- glag(zneg, pad=TRUE, pad.value=0)
y <- lgarchSim(500, xreg=0.05*zneglagged, innovations=z)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
