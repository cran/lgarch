\name{mlgarchObjective}
\alias{mlgarchObjective}
\alias{mlgarchRecursion1}
\title{Auxiliary functions}
\description{
mlgarchObjective and mlgarchRecursion1 are auxiliary functions called by \code{\link{mlgarch}}. The functions are not intended for the average user.
}
\usage{
mlgarchObjective(pars, aux)
mlgarchRecursion1(pars, aux)
}
\arguments{
  \item{pars}{numeric vector of VARMA parameters}
  \item{aux}{auxiliary list}
}
\details{To understand the structure and content of pars and aux, see the source code of the \code{\link{mlgarch}} function}
\value{
mlgarchObjective returns the log-likelihood of the VARMA representation, whereas mlgarchRecursion1 returns the residuals of the VARMA representation associated with the VARMA parameters pars
}
\references{
Sucarrat, Gronneberg and Escribano (2013), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', MPRA Paper 49344: \url{http://mpra.ub.uni-muenchen.de/49344/}
}
\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{ \code{\link{mlgarch}} and \code{\link{lgarchObjective}} }

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
