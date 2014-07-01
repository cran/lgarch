\name{coef.lgarch}
\alias{coef.lgarch}
\alias{fitted.lgarch}
\alias{logLik.lgarch}
\alias{print.lgarch}
\alias{rss}
\alias{residuals.lgarch}
\alias{vcov.lgarch}

\title{ Extraction methods for 'lgarch' objects }
\description{
Extraction methods for objects of class 'lgarch' (i.e. the result of estimating a log-GARCH model)
}

\usage{
\method{coef}{lgarch}(object, arma = FALSE, ...)
\method{fitted}{lgarch}(object, verbose = FALSE, ...)
\method{logLik}{lgarch}(object, arma = FALSE, ...)
\method{print}{lgarch}(x, arma = FALSE, ...)
## informal method for class 'lgarch'
rss(object, ...)
\method{residuals}{lgarch}(object, arma = FALSE, ...)
\method{vcov}{lgarch}(object, arma = FALSE, ...)
}

\arguments{
  \item{object}{an object of class 'lgarch'}
  \item{x}{an object of class 'lgarch'}
  \item{verbose}{logical. If FALSE (default), then only basic information is returned}
  \item{arma}{logical. If FALSE (default), then information relating to the log-garch model is returned. If TRUE, then information relating to the ARMA representation is returned}
  \item{...}{additional arguments}
}

\details{Note: The rss function is not a formal S3 method.}

\value{
  \item{coef:}{A numeric vector containing the parameter estimates}
  \item{fitted:}{A \code{\link{zoo}} object. If verbose = FALSE (default), then the zoo object is a vector containing the fitted conditional standard deviations. If verbose = TRUE, then the zoo object is a matrix containing the conditional standard deviations and additional information}
  \item{logLik:}{The value of the log-likelihood (contributions at zeros excluded) at the maximum}
  \item{print:}{Prints the most important parts of the estimation results}
  \item{residuals:}{A \code{\link{zoo}} object with the residuals. If arma = FALSE (default), then the standardised residuals are returned. If arma = TRUE, then the residuals of the ARMA representation is returned}
  \item{rss:}{A numeric; the Residual Sum of Squares of the ARMA representation}
  \item{vcov:}{The variance-covariance matrix}
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}
\seealso{ \code{\link{lgarch}} }

\examples{
##simulate 500 observations w/default parameter values:
set.seed(123)
y <- lgarchSim(500)

##estimate a log-garch(1,1):
mymod <- lgarch(y)

##print results:
print(mymod)

##extract coefficients:
coef(mymod)

##extract Gaussian log-likelihood (zeros excluded) of the log-garch model:
logLik(mymod)

##extract the Residual Sum of Squares of the ARMA representation:
rss(mymod)

##extract log-likelihood (zeros excluded) of the ARMA representation:
logLik(mymod, arma=TRUE)

##extract variance-covariance matrix:
vcov(mymod)

##extract and plot the fitted conditional standard deviation:
sdhat <- fitted(mymod)
plot(sdhat)

##extract and plot standardised residuals:
zhat <- residuals(mymod)
plot(zhat)

##extract and plot all the fitted series:
myhat <- fitted(mymod, verbose=TRUE)
plot(myhat)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
