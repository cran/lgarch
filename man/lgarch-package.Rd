\name{lgarch-package}
\alias{lgarch-package}
\docType{package}
\title{Simulation and estimation of log-GARCH models}

\description{
This package provides facilities for the simulation and estimation of univariate log-GARCH models, and for the multivariate CCC-log-GARCH(1,1) model, see Sucarrat, Gronneberg and Escribano (2013), Sucarrat and Escribano (2013), and Francq and Sucarrat (2013).\cr

Let y[t] denote a financial return or the error of a regression at time t such that

y[t] = sigma[t]*z[t],

where sigma[t] > 0 is the conditional standard deviation or volatility at t, and where z[t] is an IID innovation with mean zero and unit variance. The log-volatility specifiction of the log-GARCH-X model is given by

ln sigma[t]^2 = intercept + Sum_i alpha_i * ln y[t-i]^2 + Sum_j beta_j *ln sigma[t-1]^2 + Sum_k lambda_k * x[t]_k,\cr

where the conditioning x-variables can be contemporaneous and/or lagged. The lgarch package estimates this model via its ARMA-X representation, see Sucarrat, Gronneberg and Escribano (2013), and treats zeros on y as missing values, see Sucarrat and Escribano (2013).

}

\details{
\tabular{ll}{
Package: \tab lgarch\cr
Type: \tab Package\cr
Version: \tab 0.6-2\cr
Date: \tab 2015-09-14\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The main functions of the package are: \code{\link{lgarchSim}}, \code{\link{mlgarchSim}}, \code{\link{lgarch}} and \code{\link{mlgarch}}. The first two functions simulate from a univariate and a multivariate log-GARCH model, respectively, whereas the latter two estimate a univariate and a multivariate log-GARCH model, respectively.

The lgarch and mlgarch functions return an object (a list) of class 'lgarch' and 'mlgarch', respectively. In both cases a collection of methods can be applied to each of them: coef, fitted, logLik, print, residuals, summary and vcov. In addition, the function \code{\link{rss}} can be used to extract the Residual Sum of Squares of the estimated ARMA representation from an lgarch object.

The output produced by the \code{\link{lgarchSim}} and \code{\link{mlgarchSim}} functions, and by the fitted and residuals methods, are of the Z's ordered observations (\code{\link{zoo}}) class, see Zeileis and Grothendieck (2005), and Zeileis, Grothendieck and Ryan (2014). This means a range of time-series and plotting methods are available for these objects.
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\references{
Francq, C. and G. Sucarrat (2013), 'An Exponential Chi-Squared QMLE for Log-GARCH Models via the ARMA Representation', MPRA Paper 51783: \url{http://mpra.ub.uni-muenchen.de/51783/} \cr

Sucarrat, G. and A. Escribano (2013), 'Unbiased QML Estimation of Log-GARCH Models in the Presence of Zero Returns', MPRA Paper 50699: \url{http://mpra.ub.uni-muenchen.de/50699/} \cr

Sucarrat, G., S. Gronneberg and A. Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/} \cr

Zeileis, A. and G. Grothendieck (2005), 'zoo: S3 Infrastructure for Regular and Irregular Time Series', Journal of Statistical Software 14, pp. 1-27 \cr

Zeileis, A., G. Grothendieck, J.A. Ryan and F. Andrews(2014), 'zoo: S3 Infrastructure for Regular and Irregular Time Series (Z's ordered observations)', R package version 1.7-11, http://CRAN.R-project.org/package=zoo/
}
\seealso{
\code{\link{lgarchSim}}, \code{\link{mlgarchSim}}, \code{\link{lgarch}}, \code{\link{mlgarch}}, \code{\link{coef.lgarch}}, \code{\link{coef.mlgarch}}, \code{\link{fitted.lgarch}}, \cr
\code{\link{fitted.mlgarch}}, \code{\link{logLik.lgarch}}, \code{\link{logLik.mlgarch}}, \code{\link{print.lgarch}}, \code{\link{print.mlgarch}}, \cr
\code{\link{residuals.lgarch}}, \code{\link{residuals.mlgarch}}, \code{\link{rss}}, \code{\link{summary.mlgarch}}, \code{\link{summary.mlgarch}}, \code{\link{vcov.lgarch}}, \cr
\code{\link{vcov.mlgarch}} and \code{\link{zoo}}
}

\examples{
##simulate 500 observations w/default parameter values from
##a univariate log-garch(1,1):
set.seed(123)
y <- lgarchSim(500)

##estimate a log-garch(1,1):
mymod <- lgarch(y)

##print results:
print(mymod)

##extract coefficients:
coef(mymod)

##extract Gaussian log-likelihood (zeros excluded, if any) of the log-garch model:
logLik(mymod)

##extract Gaussian log-likelihood (zeros excluded, if any) of the arma representation:
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

##simulate 1000 observations from a two-dimensional
##ccc-log-garch(1,1) w/default parameter values:
set.seed(123)
yy <- mlgarchSim(1000)

##estimate a 2-dimensional ccc-log-garch(1,1):
myymod <- mlgarch(yy)

##print results:
print(myymod)

}
