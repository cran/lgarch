\name{rmnorm}
\alias{rmnorm}

\title{ Multivariate normal random number generation }

\description{Similar to the \code{\link{rmnorm}} function from the mnormt package, but slightly faster.}

\usage{
rmnorm(n, mean = NULL, vcov = 1)
}

\arguments{
  \item{n}{number of observations}
  \item{mean}{vector of means}
  \item{vcov}{variance-covariance matrix}
}
\details{Empty}
\value{A matrix with the random numbers}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{mlgarchSim}}, \code{\link{rnorm}}, \code{\link{chol}}
}
\examples{
##generate from 2-dimensional multivariate normal
##(standardised and independent variables):
y <- rmnorm(100, vcov=diag(c(1,1)))

##generate from 2-dimensional multivariate normal
##(standardised and dependent variables):
y <- rmnorm(100, vcov=cbind(c(1,0.5),c(0.5,1)))

##generate from 3-dimensional multivariate normal
##(independent variables with variances 1, 2 and 3):
y <- rmnorm(100, vcov=diag(c(1,2,3)))

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
