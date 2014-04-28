\name{mlgarchSim}
\alias{mlgarchSim}

\title{ Simulate from a multivariate log-GARCH model }
\description{ Simulate from a multivariate log-GARCH(1,1) model (higher orders will be available in future versions) }
\usage{
mlgarchSim(n, constant = c(0, 0), arch = diag(c(0.1, 0.05)), garch = diag(c(0.7, 0.8)),
  xreg = NULL, backcast.values = list(lnsigma2 = NULL, lnz2 = NULL, xreg = NULL),
  innovations = NULL, innovations.vcov = diag(c(1, 1)), verbose = FALSE)
}
\arguments{
  \item{n}{integer, i.e. the number of observations}
  \item{constant}{vector with the log-volatility intercepts}
  \item{arch}{matrix with arch terms}
  \item{garch}{matrix with garch terms}
  \item{xreg}{matrix with the X-variables. Column 1 should contain the x-regressor of equation 1, column 2 that of equation 2, and so on}
  \item{backcast.values}{backcast values for the recursion (chosen automatically if NULL)}
  \item{innovations}{Etiher NULL (default) or a matrix (rows n and columns equal to model-dimension) containing the standardised innovations. If NULL, then the innovations are N(0,1)}
  \item{innovations.vcov}{This argument is available only if innovations=NULL. In that case, the innovations.vcov argument contains the variance-covariance matrix of the normal innovations}
  \item{verbose}{logical. If FALSE (default), then only the matrix y is returned. If TRUE, then all the output is returned}
}
\details{Empty}

\value{
A \code{\link{zoo}} matrix with rows n. If verbose = TRUE, then all the output is returned.
}
\references{
Sucarrat, Gronneberg and Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/} \cr
}
\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{ \code{\link{lgarchSim}}, \code{\link{lgarch}}, \code{\link{zoo}} }

\examples{
##simulate 500 observations w/default parameter values from
##a 2-dimensional log-garch(1,1):
set.seed(123)
y <- mlgarchSim(500)

##simulate the same series, but with more output:
set.seed(123)
y <- mlgarchSim(500, verbose=TRUE)
head(y)

##plot the simulated values:
plot(y)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
