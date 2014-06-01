\name{lgarchObjective}
\alias{lgarchObjective}
\alias{lgarchRecursion1}

\title{Auxiliary functions}

\description{
lgarchObjective and lgarchRecursion1 are auxiliary functions called by \code{\link{lgarch}}. The functions are not intended for the average user.
}

\usage{
lgarchObjective(pars, aux)
lgarchRecursion1(pars, aux)
}
\arguments{
  \item{pars}{numeric vector of parameters of the ARMA representation}
  \item{aux}{auxiliary list}
}
\details{To understand the structure and content of pars and aux, see the source code of the \code{\link{lgarch}} function}

\value{
lgarchObjective returns the value of the objective function (either the log-likelihood or the residual sum of squares) used in estimating the ARMA representation. lgarchRecursion1 returns the residuals of the ARMA representation associated with the ARMA parameters pars
}

\references{
Sucarrat, Gronneberg and Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/}
}
\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{ \code{\link{lgarch}} }

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}