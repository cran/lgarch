\name{coef.mlgarch}
\alias{coef.mlgarch}
\alias{fitted.mlgarch}
\alias{logLik.mlgarch}
\alias{print.mlgarch}
\alias{residuals.mlgarch}
\alias{vcov.mlgarch}

\title{ Extraction methods for 'mlgarch' objects }
\description{
Extraction methods for objects of class 'mlgarch' (i.e. the result of estimating a multivariate CCC-log-GARCH model)
}

\usage{
\method{coef}{mlgarch}(object, varma = FALSE, ...)
\method{fitted}{mlgarch}(object, varma = FALSE, verbose = FALSE, ...)
\method{logLik}{mlgarch}(object, varma = FALSE, ...)
\method{print}{mlgarch}(x, varma = FALSE, ...)
\method{residuals}{mlgarch}(object, varma = FALSE, ...)
\method{vcov}{mlgarch}(object, varma = FALSE, ...)
}

\arguments{
  \item{object}{an object of class 'mlgarch'}
  \item{x}{an object of class 'mlgarch'}
  \item{verbose}{logical. If FALSE (default) then only basic information is returned}
  \item{varma}{logical. If FALSE (default), then information relating to the multivariate CCC-log-GARCH model is returned. If TRUE, then information relating to the VARMA representation is returned}
  \item{...}{additional arguments}
}

\details{Empty}

\value{
  \item{coef:}{A numeric vector containing the parameter estimates}
  \item{fitted:}{A \code{\link{zoo}} object (a matrix). If verbose = FALSE (default), then the zoo object contains the fitted conditional standard deviations of each equation. If verbose = TRUE, then the zoo object also contains additional information}
  \item{logLik:}{The value of the log-likelihood (contributions at zeros excluded) at the maximum}
  \item{print:}{Prints the most important parts of the estimation results}
  \item{residuals:}{A \code{\link{zoo}} object (a matrix) with the residuals. If varma = FALSE (default), then the standardised residuals are returned. If varma = TRUE, then the residuals of the VARMA representation is returned}
  \item{vcov:}{The variance-covariance matrix}
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}
\seealso{ \code{\link{mlgarch}} }

\examples{
##simulate 1000 observations from 2-dimensional model w/default parameter values:
set.seed(123)
y <- mlgarchSim(1000)

##estimate a 2-dimensional ccc-log-garch(1,1):
mymod <- mlgarch(y)

##print results:
print(mymod)

##extract ccc-log-garch coefficients:
coef(mymod)

##extract Gaussian log-likelihood (zeros excluded) of the ccc-log-garch model:
logLik(mymod)

##extract Gaussian log-likelihood (zeros excluded) of the varma representation:
logLik(mymod, varma=TRUE)

##extract variance-covariance matrix:
vcov(mymod)

##extract and plot the fitted conditional standard deviations:
sdhat <- fitted(mymod)
plot(sdhat)

##extract and plot standardised residuals:
zhat <- residuals(mymod)
plot(zhat)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
