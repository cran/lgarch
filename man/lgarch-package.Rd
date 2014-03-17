\name{lgarch-package}
\alias{lgarch-package}
\docType{package}
\title{Simulation and estimation of log-GARCH models}

\description{
This package provides facilities for the simulation and estimation of log-GARCH models, see Sucarrat, Gronneberg and Escribano (2013).\cr

Let y[t] denote a financial return or the error of a regression at time t such that

y[t] = sigma[t]*z[t],

where sigma[t] > 0 is the volatility or conditional standard deviation, and where z[t] is IID with mean zero and unit variance. The log-volatility specifiction of the log-GARCH-X model is given by

ln sigma[t]^2 = constant + Sum_i alpha_i * ln y[t-i]^2 + Sum_j beta_j *ln sigma[t-1]^2 + Sum_k lambda_k * x[t]_k,\cr

where the x-variables may be contemporaneous or lagged. The lgarch package estimates this model via its ARMA-X representation, see Sucarrat, Gronneberg and Escribano (2013), and treats zeros as missing values, see Sucarrat and Escribano (2013).

}

\details{
\tabular{ll}{
Package: \tab lgarch\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-03-17\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The two main functions of the package are:

\code{\link{lgarchSim}}
\code{\link{lgarch}}

The first function simulates a log-GARCH model whereas the second estimates one. The second function returns an object (a list) of class 'lgarch', and a collection of methods can be applied to this class: \code{\link{coef.lgarch}}, \code{\link{fitted.lgarch}}, \code{\link{logLik.lgarch}}, \code{\link{print.lgarch}}, \code{\link{residuals.lgarch}} and \code{\link{vcov.lgarch}}. The output produced by the \code{\link{lgarchSim}} function and the \code{\link{fitted.lgarch}} and \code{\link{residuals.lgarch}} methods are of the Z's ordered observations (\code{\link{zoo}}) class, which means a range of time-series methods are available for these objects.
}
\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}
\references{
Sucarrat and Escribano (2013), 'Unbiased QML Estimation of Log-GARCH Models in the Presence of Zero Returns', MPRA Paper 50699: \url{http://mpra.ub.uni-muenchen.de/50699/} \cr

Sucarrat, Gronneberg and Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/} \cr
}
\seealso{
\code{\link{lgarchSim}}, \code{\link{lgarch}}, \code{\link{coef.lgarch}}, \code{\link{fitted.lgarch}}, \code{\link{logLik.lgarch}}, \code{\link{print.lgarch}}, \code{\link{residuals.lgarch}} and \code{\link{vcov.lgarch}}
}
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
