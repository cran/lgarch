\name{lgarchObjective}
\alias{lgarchObjective}
\alias{lgarchRecursion1}

\title{Auxiliary functions}

\description{
lgarchObjective and lgarchRecursion1 are auxiliary functions called by \code{\link[lgarch]{lgarch}}. The functions are not intended for the average user.
}

\usage{
lgarchObjective(pars, aux)
lgarchRecursion1(pars, aux)
}
\arguments{
  \item{pars}{numeric vector with the parameters of the ARMA representation}
  \item{aux}{auxiliary list}
}
\details{To understand the structure and content of pars and aux, see the source code of the \code{\link[lgarch]{lgarch}} function}

\value{
\code{lgarchObjectivek()} returns the value of the objective function (either the log-likelihood or the residual sum of squares) used in estimating the ARMA representation. lgarchRecursion1 returns the residuals of the ARMA representation associated with the ARMA parameters pars
}

\references{
C. Francq and G. Sucarrat (2018), 'An Exponential Chi-Squared QMLE for Log-GARCH Models Via the ARMA Representation', Journal of Financial Econometrics 16, pp. 129-154 \doi{https://doi.org/10.1093/jjfinec/nbx032}

G. Sucarrat, S. Gronneberg and A. Escribano (2016), 'Estimation and Inference in Univariate and Multivariate Log-GARCH-X Models When the Conditional Density is Unknown', Computational Statistics and Data Analysis 100, pp. 582-594, \doi{http://dx.doi.org/10.1016/j.csda.2015.12.005}
}
\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\seealso{ \code{\link[lgarch]{lgarch}} }

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}