\name{lgarchLogl}
\alias{lgarchLogl}
\alias{lgarchRecursion1}
\title{Auxiliary functions}
\description{
lgarchLogl and lgarchRecursion1 are auxiliary functions called by \code{\link{lgarch}}, and are not intended to be used by the average user.
}
\usage{
lgarchLogl(pars, aux)
lgarchRecursion1(pars, aux)
}
\arguments{
  \item{pars}{numeric vector of parameters}
  \item{aux}{auxiliary list}
}
\details{To understand the structure and content of pars and aux, see the source code of the \code{\link{lgarch}} function}
\value{
lgarchLogl returns the log-likelihood of the ARMA representation, whereas lgarchRecursion1 returns the residuals of the ARMA representation associated with parameters pars
}
\references{
Sucarrat, Gronneberg and Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/}
}
\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{ \code{\link{lgarch}} }

\keyword{Statistical Models}