\name{lgarch}
\alias{lgarch}

\title{ Estimate a log-GARCH model }

\description{
Fit a log-GARCH model by either (nonlinear) Least Squares (LS) or Quasi Maximum Likelihood (QML) via the ARMA representation. For QML either the Gaussian or centred exponential chi-squared distribution can be used as instrumental density, see Sucarrat, Gronneberg and Escribano (2013), and Francq and Sucarrat (2013). Zero-values on the dependent variable y are treated as missing values, as suggested in Sucarrat and Escribano (2013). Estimation is via the \code{\link[stats]{nlminb}} function, whereas a numerical estimate of the Hessian is obtained with \code{\link[stats]{optimHess}} for the computation of the variance-covariance matrix
}

\usage{
lgarch(y, arch = 1, garch = 1, xreg = NULL, initial.values = NULL,
  lower = NULL, upper = NULL, nlminb.control = list(), vcov = TRUE,
  method=c("ls","ml","cex2"), mean.correction=FALSE,
  objective.penalty = NULL, solve.tol = .Machine$double.eps,
  c.code = TRUE)
}

\arguments{
  \item{y}{numeric vector, typically a financial return series or the error of a regression}
  \item{arch}{the arch order (i.e. an integer equal to or greater than 0). The default is 1. NOTE: in the current version the order canno be greater than 1}
  \item{garch}{the garch order (i.e. an integer equal to or greater than 0). The default is 1. NOTE: in the current version the order canno be greater than 1}
  \item{xreg}{vector or matrix with conditioning variables}
  \item{initial.values}{NULL (default) or a vector with the initial values of the ARMA-representation}
  \item{lower}{NULL (default) or a vector with the lower bounds of the parameter space (of the ARMA-representation). If NULL, then the values are automatically chosen}
  \item{upper}{NULL (default) or a vector with the upper bounds of the parameter space (of the ARMA-representation). If NULL, then the values are automatically chosen}
  \item{nlminb.control}{list of control options passed on to the \code{\link[stats]{nlminb}} optimiser}
  \item{vcov}{logical. If TRUE (default), then the variance-covariance matrix is computed. The FALSE options makes estimation faster, but the variance-covariance matrix cannot be extracted subsequently}
  \item{method}{Estimation method to use. Either "ls", i.e. Nonlinear Least Squares (default), "ml", i.e. Gaussian QML or "cex2", i.e. Centred exponential Chi-squared QML, see Francq and Sucarrat (2013). Note: For the cex2 method mean-correction = FALSE is not available}
  \item{mean.correction}{Whether to mean-correct the ARMA representation. Mean-correction is usually faster, but not always recommended if covariates are added (i.e. if xreg is not NULL)}
  \item{objective.penalty}{NULL (default) or a numeric value. If NULL, then the log-likelihood value associated with the initial values is used. Sometimes estimation can result in NA and/or +/-Inf values (this can be fatal for simulations). The value objective.penalty is the value returned by the objective function \code{\link[lgarch]{lgarchObjective}} in the presence of NA or +/-Inf values}
  \item{solve.tol}{The function \code{\link[base]{solve}} is used for the inversion of the negative of the Hessian in computing the variance-covariance matrix. The value solve.tol is passed on to \code{\link[base]{solve}}, and is the tolerance for detecting linear dependencies in the columns}
  \item{c.code}{logical. TRUE (default) is (much) faster, since it makes use of compiled C-code in the recursions}
}

\value{ A list of class 'lgarch' }

\references{
C. Francq and G. Sucarrat (2018), 'An Exponential Chi-Squared QMLE for Log-GARCH Models Via the ARMA Representation', Journal of Financial Econometrics 16, pp. 129-154 \doi{https://doi.org/10.1093/jjfinec/nbx032}

G. Sucarrat and A. Escribano (2018), 'Estimation of Log-GARCH Models in the Presence of Zero Returns', European Journal of Finance 24, pp. 809-827, \doi{https://doi.org/10.1080/1351847X.2017.1336452}

G. Sucarrat, S. Gronneberg and A. Escribano (2016), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', Computational Statistics and Data Analysis 100, pp. 582-594, \doi{http://dx.doi.org/10.1016/j.csda.2015.12.005}

}

\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\note{Empty}

\seealso{
\code{\link[lgarch]{lgarchSim}}, \code{\link[lgarch]{coef.lgarch}}, \code{\link[lgarch]{fitted.lgarch}}, \code{\link[lgarch]{logLik.lgarch}}, \code{\link[lgarch]{print.lgarch}}, \code{\link[lgarch]{residuals.lgarch}} and \code{\link[lgarch]{vcov.lgarch}}
}

\examples{

##simulate 500 observations w/default parameter values:
set.seed(123)
y <- lgarchSim(500)

##estimate a log-garch(1,1) w/least squares:
mymod <- lgarch(y)

##estimate the same model, but w/cex2 method:
mymod2 <- lgarch(y, method="cex2")

##print results:
print(mymod); print(mymod2)

##extract coefficients:
coef(mymod)

##extract Gaussian log-likelihood (zeros excluded) of the log-garch model:
logLik(mymod)

##extract Gaussian log-likelihood (zeros excluded) of the arma representation:
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
