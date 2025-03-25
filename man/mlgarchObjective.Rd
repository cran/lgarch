\name{mlgarchObjective}
\alias{mlgarchObjective}
\alias{mlgarchRecursion1}
\title{Auxiliary functions}
\description{
mlgarchObjective and mlgarchRecursion1 are auxiliary functions called by \code{\link[lgarch]{mlgarch}}. The functions are not intended for the average user.
}
\usage{
mlgarchObjective(pars, aux)
mlgarchRecursion1(pars, aux)
}
\arguments{
  \item{pars}{numeric vector of VARMA parameters}
  \item{aux}{auxiliary list}
}
\details{To understand the structure and content of pars and aux, see the source code of the \code{\link[lgarch]{mlgarch}} function}
\value{
mlgarchObjective returns the log-likelihood of the VARMA representation, whereas mlgarchRecursion1 returns the residuals of the VARMA representation associated with the VARMA parameters pars
}
\references{
G. Sucarrat, S. Gronneberg and A. Escribano (2016), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', Computational Statistics and Data Analysis 100, pp. 582-594, \doi{http://dx.doi.org/10.1016/j.csda.2015.12.005}
}
\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\seealso{ \code{\link[lgarch]{mlgarch}} and \code{\link[lgarch]{lgarchObjective}} }

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
