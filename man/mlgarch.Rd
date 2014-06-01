\name{mlgarch}
\alias{mlgarch}

\title{ Estimate a multivariate CCC-log-GARCH(1,1) model }

\description{
Fit a multivariate Constant Conditional Correlation (CCC) log-GARCH(1,1) model with multivariate Gaussian Quasi Maximum Likelihood (QML) via the VARMA representation, see Sucarrat, Gronneberg and Escribano (2013). Zero-values on y are treated as missing values, as suggested in Sucarrat and Escribano (2013). Estimation is via the \code{\link{nlminb}} function, whereas a numerical estimate of the Hessian is obtained with \code{\link{optimHess}} for the computation of the variance-covariance matrix
}

\usage{
mlgarch(y, arch = 1, garch = 1, xreg = NULL, initial.values = NULL,
  backcast.values = list(lny2 = NULL, lnz2 = NULL), lower = NULL,
  upper = NULL, nlminb.control = list(), vcov = TRUE,
  objective.penalty = NULL, solve.tol = .Machine$double.eps, c.code = TRUE)
}

\arguments{
  \item{y}{a numeric matrix, typically consisting of financial returns or errors from several regressions}
  \item{arch}{the arch order (i.e. an integer equal to or greater than 0). The default is 1. NOTE: in the current version the order canno be greater than 1}
  \item{garch}{the garch order (i.e. an integer equal to or greater than 0). The default is 1. NOTE: in the current version the order canno be greater than 1}
  \item{xreg}{a vector or a matrix with the conditioning variables. It should be noted that the x-variables all enter in each of the equations}
  \item{initial.values}{NULL (default) or a vector with the initial values of the VARMA representation}
  \item{backcast.values}{initial values used for the recursions during optimisation. If the values are NULL, then the values are automatically chosen}
  \item{lower}{NULL (default) or a vector with the lower bounds of the parameter space (of the VARMA representation). If NULL, then the values are automatically chosen}
  \item{upper}{NULL (default) or a vector with the upper bounds of the parameter space (of the VARMA representation). If NULL, then the values are automatically chosen}
  \item{nlminb.control}{list of control options passed on to the \code{\link{nlminb}} optimiser}
  \item{vcov}{logical. If TRUE (default), then the variance-covariance matrix is computed. The FALSE options makes estimation faster, but the variance-covariance matrix cannot be extracted subsequently}
  \item{objective.penalty}{NULL (default) or a numeric value. If NULL, then the log-likelihood value associated with the initial values is used. Sometimes estimation can result in NA and/or +/-Inf values (this can be fatal for simulations). The value objective.penalty is the value returned by the log-likelihood function \code{\link{lgarchObjective}} in the presence of NA or +/-Inf values}
  \item{solve.tol}{The function \code{\link{solve}} is used for the inversion of the Hessian in computing the variance-covariance matrix. The value solve.tol is passed on to \code{\link{solve}}, and is the tolerance for detecting linear dependencies in the columns}
  \item{c.code}{logical. TRUE (default) is (much) faster, since it makes use of compiled C-code}
}

\value{ A list of class 'mlgarch' }

\references{
Sucarrat and Escribano (2013), 'Unbiased QML Estimation of Log-GARCH Models in the Presence of Zero Returns', MPRA Paper 50699: \url{http://mpra.ub.uni-muenchen.de/50699/} \cr

Sucarrat, Gronneberg and Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/} \cr
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\note{Empty}

\seealso{
\code{\link{lgarchSim}}, \code{\link{coef.lgarch}}, \code{\link{fitted.lgarch}}, \code{\link{logLik.lgarch}}, \code{\link{print.lgarch}}, \code{\link{residuals.lgarch}} and \code{\link{vcov.lgarch}}
}

\examples{
##simulate 1000 observations from a 2-dimensional
##ccc-log-garch(1,1) w/default parameter values:
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
