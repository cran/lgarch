\name{coef.lgarch}
\alias{coef.lgarch}
\alias{fitted.lgarch}
\alias{logLik.lgarch}
\alias{print.lgarch}
\alias{residuals.lgarch}
\alias{rss.lgarch}
\alias{vcov.lgarch}

\title{ Extraction methods for 'lgarch' objects }
\description{
Extraction methods for objects of class 'lgarch' (i.e. the result of estimating a log-GARCH model)
}

\usage{
\method{coef}{lgarch}(object, arma = FALSE, ...)
\method{fitted}{lgarch}(object, verbose = FALSE, ...)
\method{logLik}{lgarch}(object, arma = FALSE, ...)
\method{print}{lgarch}(x, arma.version = FALSE, ...)
\method{residuals}{lgarch}(object, arma = FALSE, ...)
\method{rss}{lgarch}(object, ...)
\method{vcov}{lgarch}(object, arma.version = FALSE, full.matrix = FALSE, ...)
}

\arguments{
  \item{object}{an object of class 'lgarch'}
  \item{x}{an object of class 'lgarch'}
  \item{verbose}{logical. If FALSE (default) then only basic information is returned}
  \item{arma}{logical. If TRUE, then information relating to the arma representation is returned. If FALSE (default), then information relating to the log-garch model is returned}
  \item{arma.version}{logical. If TRUE, then information relating to the arma representation is returned. If FALSE (default), then information relating to the log-garch model is returned}
  \item{full.matrix}{logical. If TRUE, then the full variance-covariance matrix is returned}
  \item{...}{additional arguments}
}

\details{Empty}

\value{
  \item{coef:}{A numeric vector containing the parameter estimates}
  \item{fitted:}{A \code{\link{zoo}} object. If verbose=FALSE (default), then the zoo object is a vector containing the fitted conditional standard deviations. If verbose = TRUE, then the zoo object is a matrix containing the return series y, fitted conditional standard deviation (sigma)and more}
  \item{logLik:}{The value of the log-likelihood (contributions at zeros excluded) at the maximum}
  \item{print:}{Prints the most important parts of the estimation results}
  \item{residuals:}{A \code{\link{zoo}} object}
  \item{rss:}{Residual Sum of Squares of the estimated ARMA representation (zeros excluded)}
  \item{vcov:}{the variance-covariance matrix}
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

##extract Gaussian log-likelihood (zeros excluded) of the arma representation:
logLik(mymod, arma=TRUE)

##extract variance-covariance matrix:
vcov(mymod)

##extract and plot the fitted conditional standard deviation:
sigmahat <- fitted(mymod)
plot(sigmahat)

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
