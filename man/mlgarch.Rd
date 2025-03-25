\name{mlgarch}
\alias{mlgarch}

\title{ Estimate a multivariate CCC-log-GARCH(1,1) model }

\description{
Fit a multivariate Constant Conditional Correlation (CCC) log-GARCH(1,1) model with multivariate Gaussian Quasi Maximum Likelihood (QML) via the VARMA representation, see Sucarrat, Gronneberg and Escribano (2013). Zero-values on y are treated as missing values, as suggested in Sucarrat and Escribano (2013). Estimation is via the \code{\link[stats]{nlminb}} function, whereas a numerical estimate of the Hessian is obtained with \code{\link[stats]{optimHess}} for the computation of the variance-covariance matrix
}

\usage{
mlgarch(y, arch = 1, garch = 1, xreg = NULL, initial.values = NULL,
  lower = NULL, upper = NULL, nlminb.control = list(), vcov = TRUE,
  objective.penalty = NULL, solve.tol = .Machine$double.eps, c.code = TRUE)
}

\arguments{
  \item{y}{a numeric matrix, typically financial returns or regression errors}
  \item{arch}{the arch order (i.e. an integer equal to or greater than 0). The default is 1. NOTE: in the current version the order cannot be greater than 1}
  \item{garch}{the garch order (i.e. an integer equal to or greater than 0). The default is 1. NOTE: in the current version the order cannot be greater than 1}
  \item{xreg}{a vector or a matrix with the conditioning variables. The x-variables enter in each of the equations}
  \item{initial.values}{NULL (default) or a vector with the initial values of the VARMA representation}
  \item{lower}{NULL (default) or a vector with the lower bounds of the parameter space (of the VARMA representation). If NULL, then the values are automatically chosen}
  \item{upper}{NULL (default) or a vector with the upper bounds of the parameter space (of the VARMA representation). If NULL, then the values are automatically chosen}
  \item{nlminb.control}{list of control options passed on to the \code{\link[stats]{nlminb}} optimiser}
  \item{vcov}{logical. If TRUE (default), then the variance-covariance matrix is computed. The FALSE options makes estimation faster, but the variance-covariance matrix cannot be extracted subsequently}
  \item{objective.penalty}{NULL (default) or a numeric value. If NULL, then the log-likelihood value associated with the initial values is used. Sometimes estimation can result in NA and/or +/-Inf values (this can be fatal for simulations). The value objective.penalty is the value returned by the log-likelihood function \code{\link[lgarch]{lgarchObjective}} in the presence of NA or +/-Inf values}
  \item{solve.tol}{The function \code{\link[base]{solve}} is used for the inversion of the Hessian in computing the variance-covariance matrix. The value solve.tol is passed on to \code{\link[base]{solve}}, and is the tolerance for detecting linear dependencies in the columns}
  \item{c.code}{logical. TRUE (default) is (much) faster, since it makes use of compiled C-code}
}

\value{ A list of class 'mlgarch' }

\references{
G. Sucarrat and A. Escribano (2018), 'Estimation of Log-GARCH Models in the Presence of Zero Returns', European Journal of Finance 24, pp. 809-827, \doi{https://doi.org/10.1080/1351847X.2017.1336452}

G. Sucarrat, S. Gronneberg and A. Escribano (2016), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', Computational Statistics and Data Analysis 100, pp. 582-594, \doi{http://dx.doi.org/10.1016/j.csda.2015.12.005}
}

\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\note{Empty}

\seealso{
\code{\link[lgarch]{lgarchSim}}, \code{\link[lgarch]{coef.lgarch}}, \code{\link[lgarch]{fitted.lgarch}}, \code{\link[lgarch]{logLik.lgarch}}, \code{\link[lgarch]{print.lgarch}}, \code{\link[lgarch]{residuals.lgarch}} and \code{\link[lgarch]{vcov.lgarch}}
}

\examples{

##simulate from a 2-dimensional ccc-log-garch(1,1) w/defaults:
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
